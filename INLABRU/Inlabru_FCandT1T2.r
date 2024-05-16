

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

colsc <- function(...) {
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    limits = range(..., na.rm = TRUE)
  )
}


Coords=read.table('FlatCoords.txt',header = FALSE,sep=" ")/100


Sig1=read.csv('FlatFCGradient_0.txt',header = FALSE)

Sig2=read.csv('Myelin.txt',header = FALSE)


names(Coords)<-list("l1","l2","l3")


Coords <- subset(Coords, select = -c(l3))
names(Sig1)<-list("Sig1")

names(Sig2)<-list("Sig2")





df <- cbind(Coords, Sig1,Sig2)





# Remove rows with any 0 values
df_no_zeros <- df %>% 
  filter_all(all_vars(. != 0))
  
set.seed(23)

df_no_zeros$id <- 1:nrow(df_no_zeros)


yourData<-df_no_zeros[sample(nrow(df_no_zeros)),]



#Create 10 equally size folds
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)



rmse <- data.frame()
dic <- data.frame()

#Perform 10 fold cross validation
for(i in 1:10){
 
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- yourData[testIndexes, ]
    trainData <- yourData[-testIndexes, ]
	



train <- trainData %>% dplyr::sample_frac(0.3) 
train_mesh <- train %>% dplyr::sample_frac(0.1) 
train_fit <- train %>% anti_join(train_mesh,train, by='id')



# Take a random subset of rows
sampled_df <- train_mesh

sampled_df2 <- train_fit

sampled_df3 <- testData

coords=cbind(sampled_df$l1,sampled_df$l2)

boundary=fm_extensions(coords, c(50/100, 100/100))

mesh <- fm_mesh_2d(boundary = boundary,max.edge = c(10/100, 30/100))
plot(mesh)

#matern <-inla.spde2.matern(mesh)

#matern <-inla.spde2.matern(mesh)
#matern <-inla.spde2.pcmatern(mesh,prior.range=c(i,0.1), prior.sigma=c(100,0.01))


sigma0<-sd(Sig1$Sig1)
size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
range0 <- size/5
kappa0 <- sqrt(8)/range0
tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
matern <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1), B.kappa = cbind(log(kappa0), 
  0, -1), theta.prior.mean = c(0, 0), theta.prior.prec = c(0.1, 1))

mydata <- sf::st_as_sf(
  sampled_df2,
  coords = c("l1", "l2")
)

mydataPred <- sf::st_as_sf(
  sampled_df3,
  coords = c("l1", "l2")
)

print("new") # delete me
cmp1 <- Sig1 ~  -1 + Met(geometry, weights = Sig2, model = matern) + field(geometry, model = matern)
fit1 <- bru(cmp1, mydata, family = "gaussian")
#print(fit1$dic$dic) # delete me


cmp2 <- Sig1 ~ -1 + Sig2 + field(geometry, model = matern)
fit2 <- bru(cmp2, mydata, family = "gaussian")


cmp3 <- Sig1 ~ -1 + field(geometry, model = matern)

fit3 <- bru(cmp3, mydata, family = "gaussian")

cmp4 <- Sig1 ~  Sig2 + Intercept(1)

fit4 <- bru(cmp4, mydata, family = "gaussian")

cmp6 <- Sig1 ~  -1 + Sig2 + Met(geometry, weights = Sig2, model = matern) + field(geometry, model = matern)
fit6 <- bru(cmp6, mydata, family = "gaussian")

cmp7 <- Sig1 ~  1 + Sig2 + Met(geometry, weights = Sig2, model = matern)
fit7 <- bru(cmp7, mydata, family = "gaussian")

cmp8 <- Sig1 ~  1 + Met(geometry, weights = Sig2, model = matern)
fit8 <- bru(cmp8, mydata, family = "gaussian")

# Interpolate the covariate so it is at every point on the mesh
MeshInterp<-akima::interpp(x=sampled_df2$l1,y=sampled_df2$l2,z=sampled_df2$Sig2,xo=mesh$loc[,1],yo=mesh$loc[,2],extrap=TRUE,linear=FALSE)
Sig2_MeshInterp<-(MeshInterp$z-mean(MeshInterp$z))/sd(MeshInterp$z)

nu <- 1 
alpha <- nu + 2 / 2
# log(kappa)
logkappa0 <- log(8 * nu) / 2
# log(tau); in two lines to keep code width within range
logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 
logtau0 <- logtau0 - logkappa0
# SPDE model
maternNS <- inla.spde2.matern(mesh, 
  B.tau = cbind(logtau0, -1, nu, nu * Sig2_MeshInterp), 
  B.kappa = cbind(logkappa0, 0, -1, -1 * Sig2_MeshInterp),
  theta.prior.mean = rep(0, 3), 
  theta.prior.prec = rep(1, 3))


cmp5 <- Sig1 ~  -1 + Sig2 + field(geometry, model = maternNS)
fit5 <- bru(cmp5, mydata, family = "gaussian")

sigma0 <- 1
size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
range0 <- size/5
kappa0 <- sqrt(8)/range0
tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
spde <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1),B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0),theta.prior.prec = c(0.1, 1))

nu <- 1 
alpha <- nu + 2 / 2
# log(kappa)
logkappa0 <- log(8 * nu) / 2
# log(tau); in two lines to keep code width within range
logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2 
logtau0 <- logtau0 - logkappa0
# SPDE model
maternNS <- inla.spde2.matern(mesh, 

  B.tau = cbind(logtau0, -1, nu, nu * mesh$loc[,1]), 
  B.kappa = cbind(logkappa0, 0, -1, -1 * mesh$loc[,1]),

  theta.prior.mean = rep(0, 3), 
  theta.prior.prec = rep(1, 3))

cmp9 <- Sig1 ~  -1 + field(geometry, model = maternNS)
fit9 <- bru(cmp9, mydata, family = "gaussian")

pred1 <- predict(
  fit1, mydataPred,~(field+Met)
)
pred2 <- predict(
  fit2, mydataPred, ~(field+Sig2))

pred3 <- predict(
  fit3, mydataPred, ~field)

pred4 <- predict(
  fit4, mydataPred, ~Sig2)
  
pred5 <- predict(
  fit5, mydataPred, ~(field+Sig3))

pred7 <- predict(
  fit7, mydataPred, ~(Sig2 + Met))	

pred8 <- predict(
  fit8, mydataPred, ~(Met))
  
pred9 <- predict(
  fit9, mydataPred, ~(field))  



pred1rmse<- sqrt(mean((mydataPred$Sig1 - pred1$mean)^2))
pred2rmse<- sqrt(mean((mydataPred$Sig1 - pred2$mean)^2))
pred3rmse<- sqrt(mean((mydataPred$Sig1 - pred3$mean)^2))
pred4rmse<- sqrt(mean((mydataPred$Sig1 - pred4$mean)^2))
pred5rmse<- sqrt(mean((mydataPred$Sig1 - pred5$mean)^2))
pred7rmse<- sqrt(mean((mydataPred$Sig1 - pred7$mean)^2))
pred8rmse<- sqrt(mean((mydataPred$Sig1 - pred8$mean)^2))
pred9rmse<- sqrt(mean((mydataPred$Sig1 - pred9$mean)^2))

rmse_temp<-cbind(pred1rmse,pred2rmse,pred3rmse,pred4rmse,pred5rmse,pred7rmse,pred8rmse,pred9rmse)

deltadf<-deltaIC(fit1,fit2,fit3,fit4,fit5,fit7,fit8,fit9, criterion="DIC")

fit1dic <- fit1$dic$dic
fit2dic <- fit2$dic$dic
fit3dic <- fit3$dic$dic
fit4dic <- fit4$dic$dic
fit5dic <- fit5$dic$dic
fit7dic <- fit7$dic$dic
fit8dic <- fit8$dic$dic
fit9dic <- fit9$dic$dic




dic_temp<-cbind(fit1dic,fit2dic,fit3dic,fit4dic,fit5dic,fit7dic,fit8dic,fit9dic)

dic <- rbind(dic, dic_temp)
rmse <- rbind(rmse,rmse_temp)
}

write.table(dic, file = "dic.txt", sep = " ",row.names = TRUE,col.names = TRUE)
write.table(rmse, file = "rmse.txt", sep = " ",row.names = TRUE,col.names = TRUE)


