
rm(list=ls())

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

if(Sys.info()['sysname'] == "Linux"){INLA:::inla.dynload.workaround()}
library(lattice)
library(fields)
library(GEOmap)
library(ggplot2)
library(gridExtra)
library(geoR)
library(plotly)


bru_safe_inla(multicore = FALSE, quietly = FALSE, minimum_version = "23.1.31")
output_path <- "InlaGeneratedNS"

generate_random_points_on_sphere <- function(num_points) {
  theta <- runif(num_points, 0, 2 * pi)  # Azimuthal angle (longitude)
  phi <- acos(runif(num_points, -1, 1))  # Polar angle (latitude)
  x <- sin(phi) * cos(theta)
  y <- sin(phi) * sin(theta)
  z <- cos(phi)
  return(data.frame(x = x, y = y, z = z))
}


nseg=30


mesh_fine <- inla.mesh.create(globe = nseg)


variable <- cbind(mesh_fine$loc[,1])
binary_variable <- ifelse(variable >= 0, 1, 0)-0.5




nu <- 2
alpha <- nu + 1 / 2

# log(kappa)
logkappa0 <- log(8 * nu) / 2

# log(tau); in two lines to keep code width within range
logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 
logtau0 <- logtau0 - logkappa0

theta <- c(-1, 0, 1)
# SPDE model
matern_fine <- inla.spde2.matern(mesh_fine, 
	
	
	B.tau = cbind(logtau0, -1, nu, nu * (binary_variable*1)), 
	B.kappa = cbind(logkappa0, 0, -1, -1 * (binary_variable*2)),
	
	theta.prior.mean = rep(0, 3), 
	theta.prior.prec = rep(0.01, 3))

 
  Q <- inla.spde.precision(matern_fine, theta = theta)
  

true_Q <- inla.spde.precision(matern_fine, theta = theta)
true_field <- inla.qsample(1, true_Q)[, 1]

true_field_alt <- inla.qsample(1, true_Q)[, 1]







plot_ly(x = mesh_fine$loc[,1], y = mesh_fine$loc[,2], z = mesh_fine$loc[,3], color = true_field,sizes=0.1)

num_points <- 1000 # Number of random points
random_points <- generate_random_points_on_sphere(num_points) # Generate random points on the unit sphere
mydata <- st_as_sf(random_points, coords = c("x", "y", "z"))


file_name <- file.path(output_path,'GroundTruth_trainCoords.txt')  
write.table(random_points, file = file_name, sep = " ",row.names = FALSE,col.names = FALSE) 


mydata$observed <-
  fm_evaluate(
    mesh_fine,
    loc = mydata,
    field = true_field
  ) 


mydata$observedAlt <-
  fm_evaluate(
    mesh_fine,
    loc = mydata,
    field = true_field_alt
  )



nseg_sample=20

mesh <- inla.mesh.create(globe = nseg_sample)


matern <-
  inla.spde2.pcmatern(mesh,
    prior.sigma = c(1, 0.01),
    prior.range = c(1, 0.01)
  )


nu <- 2
alpha <- nu + 1 / 2

# log(kappa)
logkappa0 <- log(8 * nu) / 2

# log(tau); in two lines to keep code width within range
logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 
logtau0 <- logtau0 - logkappa0

theta <- c(-1, 0, -1)
# SPDE model

variable <- cbind(mesh$loc[,1])
binary_variable <- ifelse(variable >= 0, 1, 0)-0.5

maternNS <- inla.spde2.matern(mesh, 
	
	
	B.tau = cbind(logtau0, -1, nu, nu * (binary_variable*1)), 
	B.kappa = cbind(logkappa0, 0, -1, -1 * (binary_variable*2)),
	
	theta.prior.mean = rep(0, 3), 
	theta.prior.prec = rep(1, 3))


cmp1 <- observed ~ -1 + field(geometry, model = matern)
fit1 <- bru(cmp1, mydata, family = "gaussian")


cmp2 <- observed ~ -1 + observedAlt + field(geometry, model = matern)

fit2 <- bru(cmp2, mydata, family = "gaussian")

cmp3 <- observed ~  -1 + field(geometry, model = maternNS)

fit3 <- bru(cmp3, mydata, family = "gaussian")

#summary(fit)

deltaIC(fit1,fit2,fit3,criterion="DIC")

num_evaluate_points <- 5000 # Number of random points
random_points_evaluate <- generate_random_points_on_sphere(num_evaluate_points) # Generate random points on the unit sphere
mydata_evaluate <- st_as_sf(random_points_evaluate, coords = c("x", "y", "z"))

file_name <- file.path(output_path,'GroundTruth_EvaluateCoords.txt')  
write.table(random_points_evaluate, file = file_name, sep = " ",row.names = FALSE,col.names = FALSE) 

mydata_evaluate$observedAlt <-
  fm_evaluate(
    mesh_fine,
    loc = mydata_evaluate,
    field = true_field_alt
  )

mydata_evaluate$observed <-
  fm_evaluate(
    mesh_fine,
    loc = mydata_evaluate,
    field = true_field
  )

pred1 <- predict(
  fit1, mydata_evaluate,
  ~ field
)

pred2 <- predict(
  fit2, mydata_evaluate,
  ~ field + observedAlt
)

pred3 <- predict(
  fit3, mydata_evaluate,
  ~ field
)

plot_ly(x = random_points_evaluate$x, y = random_points_evaluate$y, z = random_points_evaluate$z, color = pred1$mean,sizes=0.1)

plot_ly(x = random_points_evaluate$x, y = random_points_evaluate$y, z = random_points_evaluate$z, color = pred2$mean,sizes=0.1)

plot_ly(x = random_points_evaluate$x, y = random_points_evaluate$y, z = random_points_evaluate$z, color = pred3$mean,sizes=0.1)



file_name <- file.path(output_path,'SphereCoords_FineMesh.txt')  
write.table(mesh_fine$loc, file = file_name, row.names = FALSE, col.names = FALSE)
file_name <- file.path(output_path,'SphereTrisFineMesh.txt')  
write.table(mesh_fine$graph$tv, file = file_name, row.names = FALSE, col.names = FALSE)


file_name <- file.path(output_path,'SphereCoords_GroundTruthEval.txt')  
write.table(mesh$loc, file = file_name, row.names = FALSE, col.names = FALSE)
file_name <- file.path(output_path,'SphereTris_GroundTruthEval.txt')  
write.table(mesh$graph$tv, file = file_name, row.names = FALSE, col.names = FALSE)


file_name <- file.path(output_path,'GroundTruth_evaluate.txt')  
write.table(mydata_evaluate, file = file_name, sep = " ",row.names = TRUE,col.names = TRUE) 

file_name <- file.path(output_path,'GroundTruthPredCorrect.txt')  
write.table(pred3$median, file = file_name, sep = " ",row.names = TRUE,col.names = TRUE) 


file_name <- file.path(output_path,'GroundTruth_Train.txt')  
write.table(mydata, file = file_name, sep = " ",row.names = TRUE,col.names = TRUE) 




