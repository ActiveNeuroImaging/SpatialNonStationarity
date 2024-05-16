
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



nseg=30
output_path <- "/Users/robleech/Dropbox/BrainSmash/FlatMaps/InlaGeneratedNS"

mesh <- inla.mesh.create(globe = nseg)

file_name <- file.path(output_path,'SphereCoords.txt')  
write.table(mesh$loc, file = file_name, row.names = FALSE, col.names = FALSE)
file_name <- file.path(output_path,'SphereTris.txt')  

write.table(mesh$graph$tv, file = file_name, row.names = FALSE, col.names = FALSE)

for(j in 1:5){
for(i in 1:1000){

variable <- cbind(mesh$loc[,1])
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
spde <- inla.spde2.matern(mesh, 
	
	B.tau = cbind(logtau0, -1, nu, nu*binary_variable),  # Change for non-stationary
	B.kappa = cbind(logkappa0, 0, -1, -1 * binary_variable*j/2), # Change for non-stationary
	#B.tau = cbind(logtau0, -1, nu, nu * (binary_variable*0+1)), 
	#B.kappa = cbind(logkappa0, 0, -1, -1 * (binary_variable*0+1)),
	
	theta.prior.mean = rep(0, 3), 
	theta.prior.prec = rep(1, 3))
  Q <- inla.spde.precision(spde, theta = theta)


#spde <- 
  
  sample1 <- as.vector(inla.qsample(1, Q))  
  sample2 <- as.vector(inla.qsample(1, Q))
 samples<-cbind(sample1,sample2)

 #file_name <- file.path(output_path, paste0("NS_Maps_", i, ".txt"))
 file_name <- file.path(output_path, paste0("NS_Maps_",j ,"_", i, ".txt"))
  
  # Save the corresponding vector element to the file
  write.table(samples, file = file_name, row.names = FALSE, col.names = FALSE)
  
  
  #mvar <- diag(inla.qinv(Q))
  
  
  #image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(mvar)), 
  #           col = topo.colors(20),
  #           xlab = "", ylab="",
  #           main = paste("edge length = ", 1/nseg))
             


}

}
####### JUST FOR TESTING



variable <- cbind(mesh$loc[,1])
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
spde <- inla.spde2.matern(mesh, 
	#B.tau = cbind(logtau0, -1, nu, nu * mesh$loc[,1]), 
	#B.kappa = cbind(logkappa0, 0, -1, -1 * mesh$loc[,1]),
	B.tau = cbind(logtau0, -1, nu, nu * (binary_variable*(j-1)/4)), 
	B.kappa = cbind(logkappa0, 0, -1, -1 * (binary_variable*(j-1)/2)),
	
	theta.prior.mean = rep(0, 3), 
	theta.prior.prec = rep(1, 3))
  Q <- inla.spde.precision(spde, theta = theta)
  
  
  sample1 <- as.vector(inla.qsample(1, Q))  
  sample2 <- as.vector(inla.qsample(1, Q))
 samples<-cbind(sample1,sample2)

plot_ly(mtcars, x = mesh$loc[,1], y = mesh$loc[,2], z = mesh$loc[,3], color = sample1,sizes=0.1)
plot_ly(mtcars, x = mesh$loc[,1], y = mesh$loc[,2], z = mesh$loc[,3], color = sample2,sizes=0.1)

