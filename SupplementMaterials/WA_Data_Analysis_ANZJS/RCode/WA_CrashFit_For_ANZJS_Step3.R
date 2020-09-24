###########################################################################################
#                                                                                         #
# This R-code is for fitting the three regularized models (lasso, ridge and enet)         #
# using the logistic regression.                                                          #
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
##################
############################
# Load the useful functions#
############################
source(here::here("RCode", "VarSelectUtils.R"))
#############################################################
#                                                           #
# Fit the regularized models using logistic regression      #
#                                                           #
#############################################################
X44forLogit <- readRDS(here::here("DerivedData", "covMatForDiscretizedModels.rds"))
yResponseforLogit <-  readRDS(here::here("DerivedData", "indicatorResponseForDiscretizedModels.rds"))
wtsforLogit<- readRDS(here::here("DerivedData", "centredlogLengthForDiscretizedModels.rds"))
#############################################################
# Fitting Lasso                                             #
#############################################################
registerDoParallel(detectCores())
lassoLogit <- glmnet::cv.glmnet(x = X44forLogit, 
                                y = yResponseforLogit, 
                                offset = wtsforLogit,
                                family = "binomial", 
                                type.measure = "mse",
                                alpha = 1,
                                standardize = FALSE,
                                parallel = TRUE)
saveRDS(lassoLogit, here::here("RegularizedRegResults", "lassoLogit.rds"))
plot(lassoLogit)
#############################################################
# Fitting Ridge                                             #
#############################################################
registerDoParallel(detectCores())
ridgeLogit <- glmnet::cv.glmnet(x = X44forLogit, 
                                y = yResponseforLogit, 
                                offset = wtsforLogit,
                                family = "binomial", 
                                type.measure = "mse",
                                alpha = 0,
                                standardize = FALSE,
                                parallel = TRUE)
saveRDS(ridgeLogit, here::here("RegularizedRegResults", "ridgeLogit.rds"))
plot(ridgeLogit)
#############################################################
# Fitting Enet                                              #
#############################################################
registerDoParallel(detectCores())
enetLogit <- glmnet::cv.glmnet(x = X44forLogit, 
                               y = yResponseforLogit, 
                               offset = wtsforLogit,
                               family = "binomial", 
                               type.measure = "mse",
                               alpha = 0.5,
                               standardize = FALSE,
                               parallel = TRUE)
saveRDS(enetLogit, here::here("RegularizedRegResults", "enetLogit.rds"))
plot(enetLogit)
#############################################################