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
X44forBT <- readRDS(here::here("DerivedData", "covMatForBTModel.rds"))
yResponseforBT <-  readRDS(here::here("DerivedData", "responseForBTModel.rds"))
wtsforBT<- readRDS(here::here("DerivedData", "weightsForBTModel.rds"))
#############################################################
# Fitting Lasso                                             #
#############################################################
registerDoParallel(detectCores())
lassoBT <- glmnet::cv.glmnet(x = X44forBT, 
                             y = yResponseforBT, 
                             weights = wtsforBT,
                             family = "poisson", 
                             type.measure = "mse",
                             alpha = 1,
                             standardize = FALSE,
                             parallel = TRUE)
saveRDS(lassoBT, here::here("RegularizedRegResults", "lassoBT.rds"))
plot(lassoBT)
#############################################################
# Fitting Ridge                                             #
#############################################################
registerDoParallel(detectCores())
ridgeBT <- glmnet::cv.glmnet(x = X44forBT, 
                             y = yResponseforBT, 
                             weights = wtsforBT,
                             family = "poisson", 
                             type.measure = "mse",
                             alpha = 0,
                             standardize = FALSE,
                             parallel = TRUE)
saveRDS(ridgeBT, here::here("RegularizedRegResults", "ridgeBT.rds"))
plot(ridgeBT)
#############################################################
# Fitting Enet                                              #
#############################################################
registerDoParallel(detectCores())
enetBT <- glmnet::cv.glmnet(x = X44forBT, 
                            y = yResponseforBT, 
                            weights = wtsforBT,
                            family = "poisson", 
                            type.measure = "mse",
                            alpha = 0.5,
                            standardize = FALSE,
                            parallel = TRUE)
saveRDS(enetBT, here::here("RegularizedRegResults", "enetBT.rds"))
plot(enetBT)
#############################################################