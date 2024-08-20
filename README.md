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
```





Prepare the data: in this case the xyz coordinates were downloaded from the fslr map from https://www.nature.com/articles/s41592-022-01625-w 

```
Coords=read.table('FlatCoords.txt',header = FALSE,sep=" ")/100

Sig1=read.csv('FlatFCGradient_0.txt',header = FALSE) 

Sig2=read.csv('Myelin.txt',header = FALSE)


names(Coords)<-list("l1","l2","l3")


Coords <- subset(Coords, select = -c(l3))
names(Sig1)<-list("Sig1")

names(Sig2)<-list("Sig2")

df <- cbind(Coords, Sig1,Sig2)
```




\# Remove rows with any 0 values
df_no_zeros <- df %>% 
  filter_all(all_vars(. != 0))
  
set.seed(23)

df_no_zeros$id <- 1:nrow(df_no_zeros)


yourData<-df_no_zeros[sample(nrow(df_no_zeros)),]



\# Create 10 equally size folds
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)



rmse <- data.frame()
dic <- data.frame()



 
testIndexes <- which(folds==1,arr.ind=TRUE)
testData <- yourData[testIndexes, ]
trainData <- yourData[-testIndexes, ]
	



train <- trainData %>% dplyr::sample_frac(1) 
train_mesh <- train %>% dplyr::sample_frac(0.08) 
train_fit <- train %>% anti_join(train_mesh,train, by='id')



\# Take a random subset of rows
sampled_df <- train_mesh

sampled_df2 <- train_fit

sampled_df3 <- testData

coords=cbind(sampled_df$l1,sampled_df$l2)

boundary=fm_extensions(coords, c(50/100, 100/100))

mesh <- fm_mesh_2d(boundary = boundary,max.edge = c(10/100, 30/100))
plot(mesh)


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
cmp1 <- Sig1 ~  -1 + Svc(geometry, weights = Sig2, model = matern) + field(geometry, model = matern)
fit1 <- bru(cmp1, mydata, family = "gaussian")


pred1 <- predict(
  fit1, mydataPred,~(field+Svc)
)


