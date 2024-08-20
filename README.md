# R Code for paper on the impact of heterogeneous spatial autocorrelation on comparisons of brain maps

![alt text](https://github.com/ActiveNeuroImaging/SpatialNonStationarity/blob/main/Spinning.gif "Spinning")

### INLABRU documentation https://cran.r-project.org/web/packages/inlabru/inlabru.pdf

### This is heavily adapted from the tutorial at https://inlabru-org.github.io/inlabru/articles/svc.html

### Below is a simple example of using INLABRU with two brain maps to assess spatially varying coefficient.
```
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

set.seed(23) # using 23 rather than 42 because of https://en.wikipedia.org/wiki/23_enigma
```





Prepare the data: in this case the xyz coordinates of each vertex were taken from the fslr 32k map downloaded from https://www.nature.com/articles/s41592-022-01625-w 
The vertex values were from the first principal gradient (Margulies et al, 2016) and T1/T2 ratio map from the Human Connectome Project (Glasser et al 2013) also downloaded from Markello et al, (2022) and converted to a text file.
```
Coords=read.table('FlatCoords.txt',header = FALSE,sep=" ")/100


Sig1=read.csv('FC_0.txt',header = FALSE) 


Sig2=read.csv('Myelin.txt',header = FALSE)


names(Coords)<-list("l1","l2","l3") # just keep xy coordinates because it's a flat map projection of cortical surface

Coords <- subset(Coords, select = -c(l3))

names(Sig1)<-list("Sig1")

names(Sig2)<-list("Sig2")

df <- cbind(Coords, Sig1,Sig2)
```

Basic preproceessing of the data (i.e., removing zeros, splitting into folds for cross-validation to assess out-of-sample predictive performance).
```
df_no_zeros <- df %>% 
  filter_all(all_vars(. != 0))

df_no_zeros$id <- 1:nrow(df_no_zeros)

Data<-df_no_zeros[sample(nrow(df_no_zeros)),]

folds <- cut(seq(1,nrow(Data)),breaks=10,labels=FALSE) # 10 folds, although in this example we will only test one iteration testing on fold 1 and training on the remaining data. 
 
testIndexes <- which(folds==1,arr.ind=TRUE)
testData <- Data[testIndexes, ]
trainData <- Data[-testIndexes, ]
	
train <- trainData %>% dplyr::sample_frac(1) # proportion of the data to be selected at random for training (here we use all of it, but can reduce to reduce computational burden)
train_mesh <- train %>% dplyr::sample_frac(0.08) # proportion of the data to be used to make the mesh (here we only use a small proportion of the data to reduce memory usage).
train_fit <- train %>% anti_join(train_mesh,train, by='id') # use the data not used to create the mesh for training.

sampled_df <- train_mesh
sampled_df2 <- train_fit
sampled_df3 <- testData


coords=cbind(sampled_df$l1,sampled_df$l2)
```

Make a mesh to use for fitting the model (see https://github.com/inlabru-org/fmesher)
```
boundary=fm_extensions(coords, c(50/100, 100/100))

mesh <- fm_mesh_2d(boundary = boundary,max.edge = c(10/100, 30/100))
plot(mesh)

```

Set up the Matern spde (see https://inlabru-org.github.io/inlabru/articles/svc.html)
```
sigma0<-sd(Sig1$Sig1)
size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
range0 <- size/5
kappa0 <- sqrt(8)/range0
tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
matern <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1), B.kappa = cbind(log(kappa0), 
  0, -1), theta.prior.mean = c(0, 0), theta.prior.prec = c(0.1, 1))
```

Convert the fitting and testing data into spatial dataframes
```
mydata <- sf::st_as_sf(
  sampled_df2,
  coords = c("l1", "l2")
)

mydataPred <- sf::st_as_sf(
  sampled_df3,
  coords = c("l1", "l2")
)
```

Set up and fit the model, in this case predicting the Principal gradient with a spatially varying coefficient (of the T1/T2 map) and a spatial field instead of a constant intercept
```
cmp1 <- Sig1 ~  -1 + Svc(geometry, weights = Sig2, model = matern) + field(geometry, model = matern)
fit1 <- bru(cmp1, mydata, family = "gaussian")
```
Evaluate model performance
```
print(fit1$dic$dic) # Deviance information criteria for the model fit

pred1 <- predict(
  fit1, mydataPred,~(field+Svc)
)

print(sqrt(mean((mydataPred$Sig1 - pred1$mean)^2))) # RMSE for out of sample prediction

```


