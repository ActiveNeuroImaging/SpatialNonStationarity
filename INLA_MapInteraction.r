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
library(gstat)
library(ciftiTools)

ciftiTools.setOption('wb_path', '/Users/robleech/workbench/macosxub_apps/wb_command.app/Contents/MacOS//')

all_myelin_xifti <- read_xifti("/Users/robleech/Dropbox/HCP_S900_GroupAvg_v1/S900.All.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii") 
all_myelin_xifti<-move_from_mwall(all_myelin_xifti)
all_thick_xifti <- read_xifti("/Users/robleech/Dropbox/HCP_S900_GroupAvg_v1/S900.All.corrThickness_MSMAll.32k_fs_LR.dscalar.nii") 
all_thick_xifti<-move_from_mwall(all_thick_xifti)
all_curve_xifti <- read_xifti("/Users/robleech/Dropbox/HCP_S900_GroupAvg_v1/S900.All.curvature_MSMAll.32k_fs_LR.dscalar.nii")
all_curve_xifti<-move_from_mwall(all_curve_xifti)

mwall <- all_curve_xifti$meta$cortex$medial_wall_mask$left

# in the fullness of time replace with cifti commands
Coords <- read.table('/Users/robleech/Dropbox/Solaris/data/fslr_vCoordsMid.txt', header = FALSE, sep=" ")/100
Tris <- read.table('/Users/robleech/Dropbox/Solaris/data/fslr_vTrisMid.txt', header = FALSE, sep=" ")

d<-2
nu <- 2 - d/2 #nu = alpha - d/2, alpha = 2

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

names(Coords) <- c("l1", "l2", "l3")
names(Tris) <- c("i1", "i2", "i3")

Sig1 <- all_myelin_xifti$data$cortex_left[,1]
Sig2 <- all_thick_xifti$data$cortex_left[,1]
Sig3 <- all_curve_xifti$data$cortex_left[,1]

names(Sig1) <- "Sig1"
names(Sig2) <- "Sig2"
names(Sig3) <- "Sig3"

df <- cbind(Coords, Sig1, Sig2,Sig3)

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


# Set seed for reproducibility
set.seed(23)

# Determine split ratio (e.g., 80% train, 20% test)
# We'll do 10% for training to keep things a little faster
train_ratio <- 0.1
n_vertices <- nrow(df)
train_size <- round(train_ratio * n_vertices)
#train_size=20 # choose number of vertices by hand
# Create index for random sampling

train_idx <- sample(1:n_vertices, train_size)
test_idx <- setdiff(1:n_vertices, train_idx)

# Split the data

#train_idx=train_idx[df$Sig1[train_idx]!= NaN]
train_idx=train_idx[!is.na(df$Sig1[train_idx])]
train_data <- df[train_idx, ]

#test_idx=test_idx[df$Sig1[test_idx]!= NaN]
test_idx=test_idx[!is.na(df$Sig1[test_idx])]
test_data <- df[test_idx, ]

field_index <- inla.spde.make.index(name = "field", n.spde = spde$n.spde)

# Create projection matrix for training data 
A_train <- Matrix::sparseMatrix(
  i = 1:length(train_idx),
  j = train_idx,
  x = rep(1, length(train_idx)),
  dims = c(length(train_idx), spde$n.spde)
)

A_test <- Matrix::sparseMatrix(
  i = 1:length(test_idx),
  j = test_idx,
  x = rep(1, length(test_idx)),
  dims = c(length(test_idx), spde$n.spde)
)

###################

# Baseline model covariates but no spatial field

st.est_cov <- inla.stack(data = list(y = train_data$Sig1), A = list(A_train, 1,1), effects = list(c(list(intercept = 1)), list(cov1 = train_data$Sig2),list(cov2 = train_data$Sig3)), tag = "est")

formula <- y ~ -1 + intercept + cov1 +  cov2

model_0 <- INLA::inla(
  formula,
  data = INLA::inla.stack.data(st.est_cov),
  control.predictor = list(A = INLA::inla.stack.A(st.est_cov), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  family = "gaussian",
  inla.mode = "compact"
)
summary(model_0)


################################################### 

# Baseline model with no covariates but spatial field

st.est_base <- inla.stack(data = list(y = train_data$Sig1), A = list(A_train), effects = list(c(field_index, 
  list(intercept = 1))), tag = "est")

st.pred_base <- inla.stack(data = list(y = NA), A = list(A_test), effects = list(c(field_index, 
  list(intercept = 1))), tag = "pred")


formula <- y ~ -1 + intercept + f(field, model = spde)


model_1 <- INLA::inla(
  formula,
  data = INLA::inla.stack.data(st.est_base),
  control.predictor = list(A = INLA::inla.stack.A(st.est_base), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  family = "gaussian",
  inla.mode = "compact"
)


summary(model_1)



# Model with two covariates for thickness and curvature 

st.est_cov <- inla.stack(data = list(y = train_data$Sig1), A = list(A_train, 1,1), effects = list(c(field_index, 
                                                                                                  list(intercept = 1)), list(cov1 = train_data$Sig2),list(cov2 = train_data$Sig3)), tag = "est")

formula <- y ~ -1 + intercept + cov1 +  cov2 + f(field, model = spde)

model_2 <- INLA::inla(
  formula,
  data = INLA::inla.stack.data(st.est_cov),
  control.predictor = list(A = INLA::inla.stack.A(st.est_cov), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  family = "gaussian",
  inla.mode = "compact"
)
summary(model_2)


# Have two spatially varying coefficients. This is slower
svc1_index <- inla.spde.make.index(name = "svc1", n.spde = spde$n.spde)
svc2_index <- inla.spde.make.index(name = "svc2", n.spde = spde$n.spde)

demeanSig2=(train_data$Sig2-mean(train_data$Sig2))/sd(train_data$Sig2)
demeanSig3=(train_data$Sig3-mean(train_data$Sig3))/sd(train_data$Sig3)


A_train_svc1 <- Matrix::sparseMatrix(
  i = 1:length(train_idx),
  j = train_idx,
  x = demeanSig2,
  dims = c(length(train_idx), spde$n.spde)
)
A_train_svc2 <- Matrix::sparseMatrix(
  i = 1:length(train_idx),
  j = train_idx,
  x = demeanSig3,
  dims = c(length(train_idx), spde$n.spde)
)


st.est_svc <- inla.stack(data = list(y = train_data$Sig1), A = list(A_train,A_train_svc1,A_train_svc2), effects = list(c(field_index, 
                                                                                                           list(intercept = 1)),svc1_index,svc2_index), tag = "est")

formula <- y ~ -1 + intercept + f(svc1, model=spde)+ f(svc2, model=spde) + f(field, model = spde)

model_3 <- INLA::inla(
  formula,
  data = INLA::inla.stack.data(st.est_svc),
  control.predictor = list(A = INLA::inla.stack.A(st.est_svc), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  family = "gaussian",
  inla.mode = "compact"
)

summary(model_3)

# Have an iteraction term between the (non spatially varying) covariates 

demeanSig2=(train_data$Sig2-mean(train_data$Sig2))
demeanSig3=(train_data$Sig3-mean(train_data$Sig3))
interactSig=demeanSig2*demeanSig3

st.est_cov <- inla.stack(data = list(y = train_data$Sig1), A = list(A_train, 1,1,1), effects = list(c(field_index, 
                                                                                                    list(intercept = 1)), list(cov1 = train_data$Sig2),list(cov2 = train_data$Sig3),list(inter = interactSig)), tag = "est")

formula <- y ~ -1 + intercept + cov1 +  cov2 + inter+ f(field, model = spde)

model_4 <- INLA::inla(
  formula,
  data = INLA::inla.stack.data(st.est_cov),
  control.predictor = list(A = INLA::inla.stack.A(st.est_cov), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  family = "gaussian",
  inla.mode = "compact"
)
summary(model_4) # doesn't work very well



# Model with one spatially varying coefficient for the interaction

svcInteraction_index <- inla.spde.make.index(name = "svc_int", n.spde = spde$n.spde)
demeanSig2=(train_data$Sig2-mean(train_data$Sig2))
demeanSig3=(train_data$Sig3-mean(train_data$Sig3))
interactSig=demeanSig2*demeanSig3


A_train_svcInt <- Matrix::sparseMatrix(
  i = 1:length(train_idx),
  j = train_idx,
  x = interactSig,
  dims = c(length(train_idx), spde$n.spde)
)


st.est_svc <- inla.stack(data = list(y = train_data$Sig1), A = list(A_train,A_train_svcInt,1,1), effects = list(c(field_index, list(intercept = 1)),svcInteraction_index,list(cov1 = train_data$Sig2),list(cov2 = train_data$Sig3)), tag = "est")

formula <- y ~ -1 + intercept + cov1 + cov2 +f(svc_int, model=spde) + f(field, model = spde)

model_5 <- INLA::inla(
  formula,
  data = INLA::inla.stack.data(st.est_svc),
  control.predictor = list(A = INLA::inla.stack.A(st.est_svc), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  family = "gaussian",
  inla.mode = "compact"
)

summary(model_5)
