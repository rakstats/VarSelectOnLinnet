###########################################################################################
#                                                                                         #
# This R-code is for fitting the three regularized models (lasso, ridge and enet)         #
# using the Poisson Count regression.                                                     #
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
X44forPois <- readRDS(here::here("DerivedData", "covMatForDiscretizedModels.rds"))
yResponseforPois <-  readRDS(here::here("DerivedData", "countResponseForDiscretizedModels.rds"))
wtsforPois<- readRDS(here::here("DerivedData", "logLengthForDiscretizedModels.rds"))
#############################################################
# Fitting Lasso                                             #
#############################################################
registerDoParallel(detectCores())
lassoPois <- glmnet::cv.glmnet(x = X44forPois, 
                                y = yResponseforPois, 
                                offset = wtsforPois,
                                family = "poisson", 
                                type.measure = "mse",
                                alpha = 1,
                                standardize = FALSE,
                                parallel = TRUE)
saveRDS(lassoPois, here::here("RegularizedRegResults", "lassoPois.rds"))
plot(lassoPois)
#############################################################
# Fitting Ridge                                             #
#############################################################
registerDoParallel(detectCores())
ridgePois <- glmnet::cv.glmnet(x = X44forPois, 
                                y = yResponseforPois, 
                                offset = wtsforPois,
                                family = "poisson", 
                                type.measure = "mse",
                                alpha = 0,
                                standardize = FALSE,
                                parallel = TRUE)
saveRDS(ridgePois, here::here("RegularizedRegResults", "ridgePois.rds"))
plot(ridgePois)
#############################################################
# Fitting Enet                                              #
#############################################################
registerDoParallel(detectCores())
enetPois <- glmnet::cv.glmnet(x = X44forPois, 
                               y = yResponseforPois, 
                               offset = wtsforPois,
                               family = "poisson", 
                               type.measure = "mse",
                               alpha = 0.5,
                               standardize = FALSE,
                               parallel = TRUE)
saveRDS(enetPois, here::here("RegularizedRegResults", "enetPois.rds"))
plot(enetPois)
#############################################################