
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

nseg=10
output_path <- "/Users/robleech/Dropbox/BrainSmash/FlatMaps/InlaGeneratedNS"

mesh <- inla.mesh.create(globe = nseg)

#file_name <- file.path(output_path,'SphereCoords.txt')  
#write.table(mesh$loc, file = file_name, row.names = FALSE, col.names = FALSE)

dic <- data.frame()



for(i in 1:1){

variable <- cbind(mesh$loc[,1])
binary_variable <- ifelse(variable >= 0, 1, 0)-0.5


#nu <- 2
#alpha <- nu + 1 / 2

# log(kappa)
#logkappa0 <- log(8 * nu) / 2

# log(tau); in two lines to keep code width within range
#logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 
#logtau0 <- logtau0 - logkappa0

#theta <- c(-1, 0, 1)
# SPDE model
#spde <- inla.spde2.matern(mesh, 
	#B.tau = cbind(logtau0, -1, nu, nu * mesh$loc[,1]*1), 
	#B.kappa = cbind(logkappa0, 0, -1, -1 * mesh$loc[,1]*1),
#	B.tau = cbind(logtau0, -1, 0, 0 * (binary_variable*0)), 
#	B.kappa = cbind(logkappa0, 0, -1, -1 * (binary_variable*0)),
	
#	theta.prior.mean = rep(0, 3), 
#	theta.prior.prec = rep(1, 3))


 
  #Q <- inla.spde.precision(spde, theta = theta)
  
  #Q <- inla.spde.precision(spde, theta = c(log(3), 0))
  
  samples <- read.table(paste0("/Users/robleech/Dropbox/BrainSmash/FlatMaps/InlaGeneratedNS/NS_V3_Maps_", i, ".txt"),sep=" ",header=F)
  print(samples[1,1])
  sample1 <- samples[,1]
  sample2 <-  samples[,2] #as.vector(inla.qsample(1, Q))
  #samples<-cbind(sample1,sample2)

print(sample1[1])
# file_name <- file.path(output_path, paste0("NS_Maps_", i, ".txt"))
  
  # Save the corresponding vector element to the file
#  write.table(samples, file = file_name, row.names = FALSE, col.names = FALSE)
  
  
  #mvar <- diag(inla.qinv(Q))
  
  
  #image.plot(proj$x, proj$y, inla.mesh.project(proj, as.vector(mvar)), 
  #           col = topo.colors(20),
  #           xlab = "", ylab="",
  #           main = paste("edge length = ", 1/nseg))
             


#sigma0<-sd(sample1)
#size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
#range0 <- size/5
#kappa0 <- sqrt(8)/range0
#tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
#matern <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1), B.kappa = cbind(log(kappa0), 
#  0, -1), theta.prior.mean = c(0, 0), theta.prior.prec = c(0.1, 1))

spde <-inla.spde2.pcmatern(mesh,prior.sigma = c(0.5, 0.1),prior.range = c(0.5, 0.1))
#spde <- inla.spde2.matern(mesh)

#matern <-inla.spde2.pcmatern(mesh,prior.sigma = c(0.1, 0.01),prior.range = c(0.1, 0.1))
#mydata <- sf::st_as_sf(
#  sampled_df2,
#  coords = c("V1", "V2","V3")
#)



#names(CoordsNew)<-list("V1","V2","V3")

df <- cbind(mesh$loc, sample1,sample2)
colnames(df) <- c("V1", "V2", "V3", "Sig1", "Sig2")

df <- as.data.frame(df) 

mydata <- sf::st_as_sf(
  df,
  coords = c("V1", "V2","V3")
)



#ggplot() + gg(mesh_fine) +   gg(mydata, aes(col = Sig1))

#met<-(mydata$Sig1-mean(mydata$Sig1))/sd(mydata$Sig1)

#cmp <- Sig3 ~ -1 +  Met(geometry, weights = Sig1, model = matern) + field(geometry, model = matern)

#cmp1 <- Sig1 ~  -1 + Met(geometry, weights = Sig3, model = matern) + field(geometry, model = matern)
#fit1 <- bru(cmp1, mydata, family = "gaussian")
###,options = list(control.compute = list(cpo = TRUE,config = TRUE)))

cmp2 <- Sig1 ~ Intercept(1) + Sig2 + field(geometry, model = spde)

fit2 <- bru(cmp2, mydata, family = "gaussian")

cmp3 <- Sig1 ~ Intercept(1) + field(geometry, model = spde)

fit3 <- bru(cmp3, mydata, family = "gaussian")

cmp4 <- Sig1 ~  Sig2 + Intercept(1)

fit4 <- bru(cmp4, mydata, family = "gaussian")


nu <- 2
alpha <- nu + 1 / 2

# log(kappa)
logkappa0 <- log(8 * nu) / 2

# log(tau); in two lines to keep code width within range
logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 
logtau0 <- logtau0 - logkappa0
print(logtau0)

theta <- c(-1, 0, -1)

# SPDE model
maternNS <- inla.spde2.matern(mesh, 
	B.tau = cbind(logtau0, -1, nu, nu * binary_variable), 
	B.kappa = cbind(logkappa0, 0, -1, -1 * binary_variable*2),
	theta.prior.mean = theta, 
	theta.prior.prec = rep(0.1, 3))

#maternNS <- inla.spde2.matern(mesh, 
#	B.tau = cbind(logtau0, -1, nu, nu * mesh$loc[,1]), 
#	B.kappa = cbind(logkappa0, 0, -1, -1 * mesh$loc[,1]),
#	theta.prior.mean = rep(0, 3), 
#	theta.prior.prec = rep(1, 3))


cmp5 <- Sig1 ~   Intercept(1) + Sig2 + field(geometry, model = maternNS)
fit5 <- bru(cmp5, mydata, family = "gaussian")


cmp6 <- Sig1 ~   Intercept(1) + field(geometry, model = maternNS)
fit6 <- bru(cmp6, mydata, family = "gaussian")

deltaIC(fit2,fit3,fit4,fit5, fit6)


fit2dic <- fit2$dic$dic
fit3dic <- fit3$dic$dic
fit4dic <- fit4$dic$dic
fit5dic <- fit5$dic$dic
fit6dic <- fit6$dic$dic

dic_temp<-cbind(fit2dic,fit3dic,fit4dic,fit5dic,fit6dic)
print(dic_temp)

}


#plot_ly(mtcars, x = mesh$loc[,1], y = mesh$loc[,2], z = mesh$loc[,3], color = sample1,sizes=0.1)

