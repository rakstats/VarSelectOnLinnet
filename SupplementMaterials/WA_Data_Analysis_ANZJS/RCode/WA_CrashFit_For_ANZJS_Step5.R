#############################################################################################
#                                                                                           #
# This R-code is for preparing the design matrix (X) and and response vector (y) for        #
# fitting the three regularized models (lasso, ridge and enet) to the marked point pattern  #
# using the B-T approximation, logistic regression and Poisson count regression.            #
#                                                                                           #
#############################################################################################
#################
# Load packages #
#################
library(sp)
library(spatstat)
library(dplyr)
library(maptools)
library(glmnet)
library(stringr)
library(here)
############################
# Load the useful functions#
############################
source(here::here("RCode", "VarSelectUtils.R"))
############################
# Load necessary datasets  #
############################
waCrash <- readRDS(here::here("InitialData", "Marked_WA_Crash.rds"))
Qmarked <- linequad(waCrash, nd=10000)
linfunListScaled <- readRDS(here::here("DerivedData", "linfun_list_scaled.rds"))
vnames <- names(linfunListScaled)
form_final <- as.formula(paste("~ " , paste0(vnames, collapse = " + ")))
prepMarkedWAData <- mpl.engine(Qmarked, trend = form_final, interaction = Poisson(), 
                               covariates = linfunListScaled, preponly = TRUE)
saveRDS(prepMarkedWAData, here::here("DerivedData", "prepMarkedWAData.rds"))

#####################################################################
glmDataMarkedWa <- prepMarkedWAData$glmdata
yResponseforMarkedBT <- glmDataMarkedWa$.mpl.Y
wtsforMarkedBT <- glmDataMarkedWa$.mpl.W
covarForMarkedBT <- glmDataMarkedWa[3:20]
all89CovarsForMarkedBTMethod <- get_interaction_terms_for_marked_pattern(covarForMarkedBT)
X89forMarkedBT <- data.matrix(all89CovarsForMarkedBTMethod)


saveRDS(X89forMarkedBT, here::here("DerivedData", "covMatForMarkedBTModel.rds"))
saveRDS(yResponseforMarkedBT, here::here("DerivedData", "responseForMarkedBTModel.rds"))
saveRDS(wtsforMarkedBT, here::here("DerivedData", "weightsForMarkedBTModel.rds"))












