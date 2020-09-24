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
######################
lassoMarkedBTCv <- readRDS(here::here("RegularizedRegResults", "lassoMarkedBT.rds"))
ridgeMarkedBTCv <- readRDS(here::here("RegularizedRegResults", "ridgeMarkedBT.rds"))
enetMarkedBTCv <- readRDS(here::here("RegularizedRegResults", "enetMarkedBT.rds"))
######################
X89forMarkedBT <- readRDS(here::here("DerivedData", "covMatForMarkedBTModel.rds"))
yResponseforMarkedBT <-  readRDS(here::here("DerivedData", "responseForMarkedBTModel.rds"))
wtsforMarkedBT <- readRDS(here::here("DerivedData", "weightsForMarkedBTModel.rds"))
######################
lassoMarkedBTFit <- glmnet::glmnet(x = X89forMarkedBT, 
                                   y = yResponseforMarkedBT, 
                                   weights = wtsforMarkedBT,
                                   family = "poisson", 
                                   alpha = 1,
                                   standardize = FALSE)
saveRDS(lassoMarkedBTFit, here::here("RegularizedRegResults", "lassoMarkedBTFit.rds"))

ridgeMarkedBTFit <- glmnet::glmnet(x = X89forMarkedBT, 
                                   y = yResponseforMarkedBT, 
                                   weights = wtsforMarkedBT,
                                   family = "poisson",
                                   alpha = 0,
                                   standardize = FALSE)
saveRDS(ridgeMarkedBTFit, here::here("RegularizedRegResults", "ridgeMarkedBTFit.rds"))

enetMarkedBTFit <- glmnet::glmnet(x = X89forMarkedBT, 
                                  y = yResponseforMarkedBT, 
                                  weights = wtsforMarkedBT,
                                  family = "poisson",
                                  alpha = 0.5,
                                  standardize = FALSE)
saveRDS(enetMarkedBTFit, here::here("RegularizedRegResults", "enetMarkedBTFit.rds"))
#######################
lassoMarkedBTFit <- readRDS(here::here("RegularizedRegResults", "lassoMarkedBTFit.rds"))
ridgeMarkedBTFit <- readRDS(here::here("RegularizedRegResults", "ridgeMarkedBTFit.rds"))
enetMarkedBTFit <- readRDS(here::here("RegularizedRegResults", "enetMarkedBTFit.rds"))
options(scipen = 999)
par(mfrow = c(3,2))
plot(lassoMarkedBTCv, xlab=expression(log(gamma)), ylab="Cross-validation MSE", cex.lab=1.4, cex.axis=1.2)
plot(lassoMarkedBTFit, xvar = "dev", cex.lab=1.4, cex.axis=1.2)
plot(ridgeMarkedBTCv, xlab=expression(log(gamma)), ylab="Cross-validation MSE", cex.lab=1.4, cex.axis=1.2)
plot(ridgeMarkedBTFit, xvar = "dev", cex.lab=1.4, cex.axis=1.2)
plot(enetMarkedBTCv, xlab=expression(log(gamma)), ylab="Cross-validation MSE", cex.lab=1.4, cex.axis=1.2)
plot(enetMarkedBTFit, xvar = "dev", cex.lab=1.4, cex.axis=1.2)
# Figure saved 10x12 inches in pdf




