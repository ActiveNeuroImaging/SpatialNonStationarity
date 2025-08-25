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
#  prior.sigma = c(sd(train_data$Sig1,0.5)
#)


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
train_ratio <- 0.001
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

A_test <- Matrix::sparseMatrix(
  i = 1:length(test_idx),
  j = test_idx,
  x = rep(1, length(test_idx)),
  dims = c(length(test_idx), spde$n.spde)
)

################################################### 
st.est_base <- inla.stack(data = list(y = train_data$Sig1), A = list(A_train), effects = list(c(field_index, 
  list(intercept = 1))), tag = "est")

st.pred_base <- inla.stack(data = list(y = NA), A = list(A_test), effects = list(c(field_index, 
  list(intercept = 1))), tag = "pred")

stack_full <- INLA::inla.stack(st.est_base, st.pred_base)

formula <- y ~ -1 + intercept + f(field, model = spde)

# Fit model on combined stack
model_1 <- INLA::inla(
  formula,
  data = INLA::inla.stack.data(stack_full),
  control.predictor = list(A = INLA::inla.stack.A(stack_full), compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
  family = "gaussian"
)
summary(model_1)

idx_train <- INLA::inla.stack.index(stack_full, "train")$data
idx_pred <- INLA::inla.stack.index(stack_full, "pred")$data

predictions <- model_1$summary.fitted.values$mean[idx_pred]

# Calculate mean squared error on test set
mse <- mean((predictions - test_data$Sig1)^2)
rmse <- sqrt(mse)
cat("Root Mean Squared Error on test set:", rmse, "\n")

# Calculate correlation between predicted and actual values
corr <- cor(predictions, test_data$Sig1)
cat("Correlation between predicted and actual values:", corr, "\n")

cat("WAIC:", model_1$waic$waic, "\n") 





results_all <- data.frame(
  x = Coords$l1,
  y = Coords$l2,
  z = Coords$l3,
  actual = NA,
  predicted = NA,
  training =NA
)

# Fill in actual values
results_all$training[train_idx] <- train_data$Sig1
results_all$actual[train_idx] <- train_data$Sig1
results_all$actual[test_idx] <- test_data$Sig1

# Fill in predicted values for test set
results_all$predicted[test_idx] <- predictions
# For training set, use fitted values

# Visualize predictions on the brain surface
library(rgl)

# Plot actual values
open3d()
bg3d("white")
material3d(specular = "black")
points3d(results_all$x, results_all$y, results_all$z, 
         color = heat.colors(100)[cut(results_all$training, 100)],
         size = 3)
title3d("Training Values")

# Plot actual values
open3d()
bg3d("white")
material3d(specular = "black")
points3d(results_all$x, results_all$y, results_all$z, 
         color = heat.colors(100)[cut(results_all$actual, 100)],
         size = 3)
title3d("Actual Values")

# Plot test points only
open3d()
bg3d("white")
material3d(specular = "black")
points3d(test_data$l1, test_data$l2, test_data$l3, 
         color = heat.colors(100)[cut(predictions, 100)],
         size = 5)
title3d("Predictions on Test Vertices")

# if you're interested int the spatial range etc after fitting the model (I think). This data is very very smooth.

spde.result <- inla.spde2.result(model_1, "field", spde)
range.marginal <- spde.result$marginals.range.nominal[[1]]

# Get posterior mean of range
fit_range.mean <- inla.emarginal(function(x) x, range.marginal)

# Get credible intervals
fit_range.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), range.marginal)

# Practical range (where correlation drops close to uncorrelated)

fit_range.practical<- range.mean*sqrt(8)
