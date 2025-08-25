library(maps)
library(ggplot2)
library(sf)
library(terra)
library(tidyterra) # raster plotting
library(tidyr)
library(scales)
library(dplyr)
library(INLA)
library(inlabru)
library(fmesher)
library(akima)

# in the fullness of time replace with cifti commands
Coords <- read.table('/Users/robleech/Dropbox/Solaris/data/fslr_vCoordsMid.txt', header = FALSE, sep=" ")/100
Tris <- read.table('/Users/robleech/Dropbox/Solaris/data/fslr_vTrisMid.txt', header = FALSE, sep=" ")

max_dist <- apply(Coords, 2, function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) #MAX distance within the mesh (Euclidean, a bad proxy for geodesic)
range2 <- max(max_dist)/2 #this is our upper limit for the spatial correlation range
tri_dist <- apply(
  matrix(1:nrow(Tris), ncol=1), 1, #for each face/triangle #MIN distance within the mesh
  function(x) { mean(dist(Coords[unlist(Tris[1,]),])) } #compute the avg distance among vertices in the triangle
)
min_dist <- min(tri_dist) #distance between vertices in the smallest triangle
range1 <- min_dist*2
range0 <- range1*5
range_vec <- c(range1, range2, range0)

logkappa_vec <- log(sqrt(8*nu)/range_vec)
logkappa0 <- log(sqrt(8*nu)/range0)
var0 <- 0.1 #we usually expect most of the activation amplitudes to be between (-2 and 2) --> SD ~= 0.33, Var ~= 0.1
var_vec <- c(0.01, 1, var0) #reasonable range for variance
var2logtau <- function(var, d, kappa){
  nu <- 2 - d/2 #nu = alpha - d/2, alpha = 2
  tausq <- gamma(nu)/(var * (4*pi)^(d/2) * kappa^(2*nu))
  log(sqrt(tausq))
}
logtau_vec <- var2logtau(var_vec, d, exp(logkappa0))
logtau0 <- var2logtau(var0, d, exp(logkappa0)) #we will use a good starting value even with the default prior

Elog.tau <- logtau_vec[3]
prior_sd <- abs(diff(logtau_vec[1:2]))/4 #so that 95% of the prior density is within the range
Qlog.tau <- 1/(prior_sd^2)

Elog.kappa <- logkappa_vec[3]
prior_sd <- abs(diff(logkappa_vec[1:2]))/4
Qlog.kappa <- 1/(prior_sd^2)
  
# log_kappa_tau<-list(
#   logkappa_vec=logkappa_vec,
#   logkappa0=logkappa0,
#   logtau0=logtau0,
#   logtau_vec=logtau_vec
# )

names(Coords) <- c("l1", "l2", "l3")
names(Tris) <- c("i1", "i2", "i3")

Sig1 <- read.csv('/Users/robleech/Dropbox/Solaris/data/FC1Vals.txt', header = FALSE)
Sig2 <- read.csv('/Users/robleech/Dropbox/Solaris/data/MyelinVals.txt', header = FALSE)

names(Sig1) <- "Sig1"
names(Sig2) <- "Sig2"

df <- cbind(Coords, Sig1, Sig2)

# Create mesh with original 3D coordinates
mesh <- INLA::inla.mesh.create(
  loc = as.matrix(Coords),
  tv = as.matrix(Tris)
)


spde <- INLA::inla.spde2.matern(
  mesh,
  prior.tau = exp(Elog.tau),
  prior.kappa = exp(Elog.kappa),
  theta.prior.prec = diag(c(Qlog.tau, Qlog.kappa))
)
# Set up the SPDE

#spde <- inla.spde2.pcmatern(
#  mesh = mesh,
#  prior.range = c((max(Coords)-min(Coords))/5, 0.5),
#  prior.sigma = c(sd(train_data$Sig1),0.5)
#)

#sigma0<-sd(Sig1$Sig1)
#size <- max(Coords)-min(Coords)
#range0 <- size/5
#kappa0 <- sqrt(8)/range0
#tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)

#spde <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1), B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0), theta.prior.prec = c(0.1, 1))

# Alternative way to set up SPDE

#Elog.kappa <- Elog.tau <- 0
#Qlog.kappa <- Qlog.tau <- 0.1


#spde <- INLA::inla.spde2.matern(
#  mesh,
#  prior.tau = exp(Elog.tau),
#  prior.kappa = exp(Elog.kappa),
#  theta.prior.prec = diag(c(Qlog.tau, Qlog.kappa))
#)

# Set seed for reproducibility
set.seed(23)

# Determine split ratio (e.g., 80% train, 20% test)
# We'll do 10% for training to keep things a little faster
train_ratio <- 0.9
n_vertices <- nrow(df)
train_size <- round(train_ratio * n_vertices)
#train_size=20 # choose number of vertices by hand
# Create index for random sampling

train_idx <- sample(1:n_vertices, train_size)
test_idx <- setdiff(1:n_vertices, train_idx)

# Split the data

train_idx=train_idx[df$Sig1[train_idx]!= 0]
train_data <- df[train_idx, ]

test_idx=test_idx[df$Sig1[test_idx]!= 0]
test_data <- df[test_idx, ]

field_index <- inla.spde.make.index(name = "field", n.spde = spde$n.spde)

# Create projection matrix for training data 
A_train <- Matrix::sparseMatrix(
  i = 1:length(train_idx),
  j = train_idx,
  x = rep(1, length(train_idx)),
  dims = c(length(train_idx), spde$n.spde)
)

################################################### 
st.est_base <- inla.stack(data = list(y = train_data$Sig1), A = list(A_train), effects = list(c(field_index, 
  list(intercept = 1))), tag = "est")
  
formula <- y ~ -1 + intercept + f(field, model = spde)

# Fit model on combined stack
model_1 <- INLA::inla(
  formula,
  data = INLA::inla.stack.data(st.est_base),
  control.predictor = list(A = INLA::inla.stack.A(st.est_base), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  family = "gaussian"
)
summary(model_1)

####

demeanSig2=(train_data$Sig2-mean(train_data$Sig2))/sd(train_data$Sig2)

st.est_cov <- inla.stack(data = list(y = train_data$Sig1), A = list(A_train, 1), effects = list(c(field_index, 
  list(intercept = 1)), list(cov = demeanSig2)), tag = "est")

formula <- y ~ -1 + intercept + cov + f(field, model = spde)

model_2 <- INLA::inla(
  formula,
  data = INLA::inla.stack.data(st.est_cov),
  control.predictor = list(A = INLA::inla.stack.A(st.est_cov), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  family = "gaussian"
)
summary(model_2)

####


svc_index <- inla.spde.make.index(name = "svc", n.spde = spde$n.spde)

demeanSig2=(train_data$Sig2-mean(train_data$Sig2))/sd(train_data$Sig2)

A_train_svc <- Matrix::sparseMatrix(
  i = 1:length(train_idx),
  j = train_idx,
  x = demeanSig2,
  dims = c(length(train_idx), spde$n.spde)
)

st.est_svc <- inla.stack(data = list(y = train_data$Sig1), A = list(A_train,A_train_svc), effects = list(c(field_index, 
  list(intercept = 1)),svc_index), tag = "est")

formula <- y ~ -1 + intercept + f(svc, model=spde) + f(field, model = spde)

model_3 <- INLA::inla(
  formula,
  data = INLA::inla.stack.data(st.est_svc),
  control.predictor = list(A = INLA::inla.stack.A(st.est_svc), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  family = "gaussian"
)
summary(model_3)


