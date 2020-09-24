##################
# Load packages #
#################
library(sp)
library(spatstat)
library(dplyr)
library(maptools)
library(glmnet)
library(stringr)
library(here)
library(doParallel)
######################
lassoBTCv <- readRDS(here::here("RegularizedRegResults", "lassoBT.rds"))
lassoPoisCv <- readRDS(here::here("RegularizedRegResults", "lassoPois.rds"))
lassoLogitCv <- readRDS(here::here("RegularizedRegResults", "lassoLogit.rds"))
######################
X44forBT <- readRDS(here::here("DerivedData", "covMatForBTModel.rds"))
yResponseforBT <-  readRDS(here::here("DerivedData", "responseForBTModel.rds"))
wtsforBT<- readRDS(here::here("DerivedData", "weightsForBTModel.rds"))
X44forDiscretizedModel <- readRDS(here::here("DerivedData", "covMatForDiscretizedModels.rds"))
yResponseforLogit <-  readRDS(here::here("DerivedData", "indicatorResponseForDiscretizedModels.rds"))
wtsforLogit<- readRDS(here::here("DerivedData", "centredlogLengthForDiscretizedModels.rds"))
yResponseforPois <-  readRDS(here::here("DerivedData", "countResponseForDiscretizedModels.rds"))
wtsforPois<- readRDS(here::here("DerivedData", "logLengthForDiscretizedModels.rds"))
######################
lassoBTFit <- glmnet::glmnet(x = X44forBT, 
                             y = yResponseforBT,
                             weights = wtsforBT,
                             family = "poisson", 
                             alpha = 1,
                             standardize = FALSE)
saveRDS(lassoBTFit, here::here("RegularizedRegResults", "lassoBTFit.rds"))
######################
lassoPoisFit <- glmnet::glmnet(x = X44forDiscretizedModel, 
                               y = yResponseforPois,
                               family = "poisson", 
                               alpha = 1,
                               offset = wtsforPois,
                               standardize = FALSE)
saveRDS(lassoPoisFit, here::here("RegularizedRegResults", "lassoPoisFit.rds"))
########################
lassoLogitFit <- glmnet::glmnet(x = X44forDiscretizedModel, 
                                y = yResponseforLogit,
                                family = "binomial", 
                                alpha = 1,
                                offset = wtsforLogit,
                                standardize = FALSE)
saveRDS(lassoLogitFit, here::here("RegularizedRegResults", "lassoLogitFit.rds"))
########################
options(scipen = 999)
par(mfrow = c(3,2))
plot(lassoBTCv, xlab=expression(log(gamma)), ylab="Cross-validation MSE", cex.lab=1.4, cex.axis=1.2,
     ylim = c(0.0000335, 0.0000365))
plot(lassoBTFit, xvar = "dev", cex.lab=1.4, cex.axis=1.2)
plot(lassoLogitCv, xlab=expression(log(gamma)), ylab="Cross-validation MSE", cex.lab=1.4, cex.axis=1.2)
plot(lassoLogitFit, xvar = "dev", cex.lab=1.4, cex.axis=1.2)
plot(lassoPoisCv, xlab=expression(log(gamma)), ylab="Cross-validation MSE", cex.lab=1.4, cex.axis=1.2)
plot(lassoPoisFit, xvar = "dev", cex.lab=1.4, cex.axis=1.2)
# Figure saved 10x12 inches in pdf
########################





