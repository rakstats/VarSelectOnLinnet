###########################################################################################
#                                                                                         #
# This R-code is for fitting the three regularized models (lasso, ridge and enet)         #
# using the B-T approximation.                                                            #
#                                                                                         #
###########################################################################################
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
############################
# Load the useful functions#
############################
source(here::here("RCode", "VarSelectUtils.R"))
#############################################################
#                                                           #
# Fit the regularized models using B-T apporximation        #
#                                                           #
#############################################################
X89forMarkedBT <- readRDS(here::here("DerivedData", "covMatForMarkedBTModel.rds"))
yResponseforMarkedBT <-  readRDS(here::here("DerivedData", "responseForMarkedBTModel.rds"))
wtsforMarkedBT<- readRDS(here::here("DerivedData", "weightsForMarkedBTModel.rds"))
#############################################################
# Fitting Lasso                                             #
#############################################################
registerDoParallel(detectCores())
lassoMarkedBT <- glmnet::cv.glmnet(x = X89forMarkedBT, 
                             y = yResponseforMarkedBT, 
                             weights = wtsforMarkedBT,
                             family = "poisson", 
                             type.measure = "mse",
                             alpha = 1,
                             standardize = FALSE,
                             parallel = TRUE)
saveRDS(lassoMarkedBT, here::here("RegularizedRegResults", "lassoMarkedBT.rds"))
plot(lassoMarkedBT)
#############################################################
# Fitting Ridge                                             #
#############################################################
registerDoParallel(detectCores())
ridgeMarkedBT <- glmnet::cv.glmnet(x = X89forMarkedBT, 
                             y = yResponseforMarkedBT, 
                             weights = wtsforMarkedBT,
                             family = "poisson", 
                             type.measure = "mse",
                             alpha = 0,
                             standardize = FALSE,
                             parallel = TRUE)
saveRDS(ridgeMarkedBT, here::here("RegularizedRegResults", "ridgeMarkedBT.rds"))
rm(ridgeMarkedBT)
#############################################################
# Fitting Enet                                              #
#############################################################
registerDoParallel(detectCores())
enetMarkedBT <- glmnet::cv.glmnet(x = X89forMarkedBT, 
                            y = yResponseforMarkedBT, 
                            weights = wtsforMarkedBT,
                            family = "poisson", 
                            type.measure = "mse",
                            alpha = 0.5,
                            standardize = FALSE,
                            parallel = TRUE)
saveRDS(enetMarkedBT, here::here("RegularizedRegResults", "enetMarkedBT.rds"))
plot(enetMarkedBT)
#############################################################