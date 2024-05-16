
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


generate_random_points_on_sphere <- function(num_points) {
  theta <- runif(num_points, 0, 2 * pi)  # Azimuthal angle (longitude)
  phi <- acos(runif(num_points, -1, 1))  # Polar angle (latitude)
  x <- sin(phi) * cos(theta)
  y <- sin(phi) * sin(theta)
  z <- cos(phi)
  return(data.frame(x = x, y = y, z = z))
}


nseg=30
#output_path <- "/Users/robleech/Dropbox/BrainSmash/FlatMaps/InlaGeneratedNS"

mesh_fine <- inla.mesh.create(globe = nseg)


dic <- data.frame()



#for(i in 1:1){

#variable <- cbind(mesh_fine$loc[,1])
#binary_variable <- ifelse(variable >= 0, 1, 0)-0.5

#nu <- 2
#alpha <- nu + 1 / 2


#logkappa0 <- log(8 * nu) / 2

# log(tau); in two lines to keep code width within range
#logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 
#logtau0 <- logtau0 - logkappa0



matern_fine <-
  inla.spde2.pcmatern(mesh_fine,
    prior.sigma = c(1, 0.01),
    prior.range = c(1, 0.01)
  )
  
true_range <- 0.5
true_sigma <- 1
true_Q <- inla.spde.precision(matern_fine, theta = log(c(true_range, true_sigma)))
true_field <- inla.qsample(1, true_Q)[, 1]

true_field_alt <- inla.qsample(1, true_Q)[, 1]




#plot_ly(x = mesh$loc[,1], y = mesh$loc[,2], z = mesh$loc[,3], color = sample1,sizes=0.1)


plot_ly(x = mesh_fine$loc[,1], y = mesh_fine$loc[,2], z = mesh_fine$loc[,3], color = true_field,sizes=0.1)

num_points <- 500 # Number of random points
random_points <- generate_random_points_on_sphere(num_points) # Generate random points on the unit sphere
mydata <- st_as_sf(random_points, coords = c("x", "y", "z"))

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

#mesh<-fm_rcdt_2d(mydata$geometry)

nseg_sample=20
#output_path <- "/Users/robleech/Dropbox/BrainSmash/FlatMaps/InlaGeneratedNS"
mesh <- inla.mesh.create(globe = nseg_sample)


matern <-
  inla.spde2.pcmatern(mesh,
    prior.sigma = c(10, 0.01),
    prior.range = c(1, 0.01)
  )

cmp1 <- observed ~ field(geometry, model = matern) + Intercept(1)
fit1 <- bru(cmp1, mydata, family = "gaussian")


cmp2 <- observed ~  observedAlt + field(geometry, model = matern)  + Intercept(1)

fit2 <- bru(cmp2, mydata, family = "gaussian")

#summary(fit)

deltaIC(fit1,fit2,criterion="DIC")

num_evaluate_points <- 5000 # Number of random points
random_points_evaluate <- generate_random_points_on_sphere(num_evaluate_points) # Generate random points on the unit sphere
mydata_evaluate <- st_as_sf(random_points_evaluate, coords = c("x", "y", "z"))


pred <- predict(
  fit1, mydata_evaluate,
  ~ field + Intercept
)



plot_ly(x = random_points_evaluate$x, y = random_points_evaluate$y, z = random_points_evaluate$z, color = pred$mean,sizes=0.1)


#df <- cbind(mesh$loc, sample1,sample2)
#colnames(df) <- c("V1", "V2", "V3", "Sig1", "Sig2")

#df <- as.data.frame(df) 

#mydata <- sf::st_as_sf(
#  df,
#  coords = c("V1", "V2","V3")
#)





