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
library(xtable)
############################
# Load the useful functions#
############################
source(here::here("RCode", "VarSelectUtils.R"))
#############################
#####################################################
# Load Cross-validation and regularized fit results #
#####################################################
#################################################
#1. Load lasso Cross-validation and lasso fit   #
##############################################################################
# 1a. Load Lasso CV and fit results for Berman-Turner (BT) approximation
lassoMarkedBTCV <- readRDS(here::here("RegularizedRegResults", "lassoMarkedBT.rds"))
lassoMarkedBTFit <- readRDS(here::here("RegularizedRegResults", "lassoMarkedBTFit.rds"))
print(lassoMarkedBTCV)
print(lassoMarkedBTFit)
coef(lassoMarkedBTCV, s="lambda.min") 
coef(lassoMarkedBTCV, s="lambda.1se")
lambdaMin <- lassoMarkedBTCV$lambda.min
lambda1Se <- lassoMarkedBTCV$lambda.1se
lambdaAvg <- mean(c(lambdaMin,lambda1Se))
lambdaGmean <- sqrt(lambdaMin*lambda1Se)
log(lambdaGmean)

number_of_nonzero_coef(coef(lassoMarkedBTCV, s=lambdaMin))
number_of_nonzero_coef(coef(lassoMarkedBTCV, s=lambda1Se))
number_of_nonzero_coef(coef(lassoMarkedBTCV, s=lambdaAvg))
number_of_nonzero_coef(coef(lassoMarkedBTCV, s=lambdaGmean))

#coef(lassoMarkedBTCV, s=lambdaAvg)

lassoBTCoef00 <- coef(lassoMarkedBTCV, s=lambdaGmean) # Selects 16 variables

lassoBTCoef0 <- as.data.frame(as.matrix(lassoBTCoef00))
names(lassoBTCoef0) <- c("B-T (lasso)")
lassoBTCoef0$Variable <- rownames(lassoBTCoef0)
lassoBTCoef <- lassoBTCoef0[c(2,1)]
Col1 <- lassoBTCoef[1:45,]
Col4 <- lassoBTCoef[46:90, ]
###############################################################################
ridgeMarkedBTCV <- readRDS(here::here("RegularizedRegResults", "ridgeMarkedBT.rds"))
plot(ridgeMarkedBTCV)
lambdaRidgeMin <- ridgeMarkedBTCV$lambda.min
lambdaRidge1Se <- ridgeMarkedBTCV$lambda.1se
lambdaRidgeAvg <- mean(c(lambdaRidgeMin, lambdaRidge1Se))
lambdaRidgeGmean <- sqrt(lambdaRidgeMin*lambdaRidge1Se)

get_max_deviance(ridgeMarkedBTCV)
options(scipen=999)
sapply(c(lambdaRidgeMin, lambdaRidge1Se, lambdaRidgeAvg, lambdaRidgeGmean),
       FUN = function(z){get_fraction_deviance_explained(X=ridgeMarkedBTCV, gamma=z)})

ridgeBTCoef00 <- coef(ridgeMarkedBTCV, s=lambdaRidgeMin)
ridgeBTCoef0 <- as.data.frame(as.matrix(ridgeBTCoef00))
names(ridgeBTCoef0) <- c("B-T (ridge)")
ridgeBTCoef0$Variable <- rownames(ridgeBTCoef0)
ridgeBTCoef <- ridgeBTCoef0[c(2,1)]
Col2 <- ridgeBTCoef[1:45,]
Col5 <- ridgeBTCoef[46:90, ]
###################################################################################
enetMarkedBTCV <- readRDS(here::here("RegularizedRegResults", "enetMarkedBT.rds"))
plot(enetMarkedBTCV)
lambdaenetMin <- enetMarkedBTCV$lambda.min
lambdaenet1Se <- enetMarkedBTCV$lambda.1se
lambdaenetAvg <- mean(c(lambdaenetMin,lambdaenet1Se))
lambdaenetGmean <- sqrt(lambdaenetMin*lambdaenet1Se)
log(lambdaenetGmean)

number_of_nonzero_coef(coef(enetMarkedBTCV, s=lambdaenetMin))
number_of_nonzero_coef(coef(enetMarkedBTCV, s=lambdaenet1Se))
number_of_nonzero_coef(coef(enetMarkedBTCV, s=lambdaenetAvg))
number_of_nonzero_coef(coef(enetMarkedBTCV, s=lambdaenetGmean))

enetBTCoef0 <- as.data.frame(as.matrix(coef(enetMarkedBTCV, s=lambdaenetGmean)))
names(enetBTCoef0) <- c("B-T (e-net)")
enetBTCoef0$Variable <- rownames(enetBTCoef0)
enetBTCoef <- enetBTCoef0[c(2,1)]
Col3 <- enetBTCoef[1:45,]
Col6 <- enetBTCoef[46:90, ]

###################################################################################
Table8a <- Col1 %>% left_join(Col2, by = "Variable") %>% left_join(Col3, by="Variable")
Table8b <- Col4 %>% left_join(Col5, by = "Variable") %>% left_join(Col6, by="Variable")
###################################################################################
Table8a[Table8a$Variable == "SPLI_SPEED", "Variable"] <- "SPD_LIM"
Table8a[Table8a$Variable == "HOAL_CURVE", "Variable"] <- "H_CURVE"
Table8a[Table8a$Variable == "TOTAL_PAVE", "Variable"] <- "TOT_P"
Table8a[Table8a$Variable == "TOTAL_SEAL", "Variable"] <- "TOT_S"
Table8a[Table8a$Variable == "TRAFFICABL", "Variable"] <- "TRFABL"
Table8a[Table8a$Variable == "NO_OF_LANE", "Variable"] <- "N_LANE"

Table8a[Table8a$Variable == "SPLI_SPEED2", "Variable"] <- "SPD_LIM2"
Table8a[Table8a$Variable == "HOAL_CURVE2", "Variable"] <- "H_CURVE2"
Table8a[Table8a$Variable == "TOTAL_PAVE2", "Variable"] <- "TOT_P2"
Table8a[Table8a$Variable == "TOTAL_SEAL2", "Variable"] <- "TOT_S2"
Table8a[Table8a$Variable == "TRAFFICABL2", "Variable"] <- "TRFABL2"
Table8a[Table8a$Variable == "NO_OF_LANE2", "Variable"] <- "N_LANE2"

Table8a[Table8a$Variable == "SHOULDER_S", "Variable"] <- "SHLDR"
Table8a[Table8a$Variable == "FLOODWAY", "Variable"] <- "FLDWY"
Table8a[Table8a$Variable == "BRIDGE", "Variable"] <- "BRDG"


Table8a[Table8a$Variable == "SPLI_SPEEDxKERB_L", "Variable"] <- "SPD_LIMxKERB_L"
Table8a[Table8a$Variable == "HOAL_CURVExKERB_L", "Variable"] <- "H_CURVExKERB_L"
Table8a[Table8a$Variable == "TOTAL_PAVExKERB_L", "Variable"] <- "TOT_PxKERB_L"
Table8a[Table8a$Variable == "TOTAL_SEALxKERB_L", "Variable"] <- "TOT_SxKERB_L"
Table8a[Table8a$Variable == "TRAFFICABLxKERB_L", "Variable"]<- "TRFABLxKERB_L"
Table8a[Table8a$Variable == "NO_OF_LANExKERB_L", "Variable"]<- "N_LANExKERB_L"
Table8a[Table8a$Variable == "SHOULDER_SxKERB_L", "Variable"]<- "SHLDRxKERB_L"
Table8a[Table8a$Variable == "FLOODWAYxKERB_L", "Variable"]<- "FLDWYxKERB_L"
Table8a[Table8a$Variable == "BRIDGExKERB_L", "Variable"]<- "BRDGxKERB_L"

Table8a[Table8a$Variable == "SPLI_SPEEDxSHOULDER_S", "Variable"] <- "SPD_LIMxSHLDR"
Table8a[Table8a$Variable == "HOAL_CURVExSHOULDER_S", "Variable"] <- "H_CURVExSHLDR"
Table8a[Table8a$Variable == "TOTAL_PAVExSHOULDER_S", "Variable"] <- "TOT_PxSHLDR"
Table8a[Table8a$Variable == "TOTAL_SEALxSHOULDER_S", "Variable"] <- "TOT_SxSHLDR"
Table8a[Table8a$Variable == "TRAFFICABLxSHOULDER_S", "Variable"] <- "TRFABLxSHLDR"
Table8a[Table8a$Variable == "NO_OF_LANExSHOULDER_S", "Variable"] <- "N_LANExSHLDR"
Table8a[Table8a$Variable == "KERB_RxSHOULDER_S", "Variable"] <- "KERB_RxSHLDR"
Table8a[Table8a$Variable == "FLOODWAYxSHOULDER_S", "Variable"] <- "FLDWYxSHLDR"
Table8a[Table8a$Variable == "BRIDGExSHOULDER_S", "Variable"] <- "BRDGxSHLDR"

Table8a[Table8a$Variable == "SPLI_SPEEDxKERB_R", "Variable"] <- "SPD_LIMxKERB_R"
Table8a[Table8a$Variable == "HOAL_CURVExKERB_R", "Variable"] <- "H_CURVExKERB_R"
Table8a[Table8a$Variable == "TOTAL_PAVExKERB_R", "Variable"] <- "TOT_PxKERB_R"
Table8a[Table8a$Variable == "TOTAL_SEALxKERB_R", "Variable"] <- "TOT_SxKERB_R"
Table8a[Table8a$Variable == "TRAFFICABLxKERB_R", "Variable"] <- "TRFABLxKERB_R"
Table8a[Table8a$Variable == "NO_OF_LANExKERB_R", "Variable"] <- "N_LANExKERB_R"
Table8a[Table8a$Variable == "FLOODWAYxKERB_R", "Variable"] <- "FLDWYxKERB_R"
Table8a[Table8a$Variable == "BRIDGExKERB_R", "Variable"] <- "BRDGxKERB_R"
##############################################################################
# Table8b
##############################################################################
Table8b[Table8b$Variable == "MARKSxSPLI_SPEED", "Variable"] <- "MarkxSPD_LIM"
Table8b[Table8b$Variable == "MARKSxHOAL_CURVE", "Variable"] <- "MarkxH_CURVE"
Table8b[Table8b$Variable == "MARKSxTOTAL_PAVE", "Variable"] <- "MarkxTOT_P"
Table8b[Table8b$Variable == "MARKSxTOTAL_SEAL", "Variable"] <- "MarkxTOT_S"
Table8b[Table8b$Variable == "MARKSxTRAFFICABL", "Variable"] <- "MarkxTRFABL"
Table8b[Table8b$Variable == "MARKSxNO_OF_LANE", "Variable"] <- "MarkxN_LANE"

Table8b[Table8b$Variable == "MARKSxSPLI_SPEED2", "Variable"] <- "MarkxSPD_LIM2"
Table8b[Table8b$Variable == "MARKSxHOAL_CURVE2", "Variable"] <- "MarkxH_CURVE2"
Table8b[Table8b$Variable == "MARKSxTOTAL_PAVE2", "Variable"] <- "MarkxTOT_P2"
Table8b[Table8b$Variable == "MARKSxTOTAL_SEAL2", "Variable"] <- "MarkxTOT_S2"
Table8b[Table8b$Variable == "MARKSxTRAFFICABL2", "Variable"] <- "MarkxTRFABL2"
Table8b[Table8b$Variable == "MARKSxNO_OF_LANE2", "Variable"] <- "MarkxN_LANE2"

Table8b[Table8b$Variable == "MARKSxSHOULDER_S", "Variable"] <- "MarkxSHLDR"
Table8b[Table8b$Variable == "MARKSxFLOODWAY", "Variable"] <- "MarkxFLDWY"
Table8b[Table8b$Variable == "MARKSxBRIDGE", "Variable"] <- "MarkxBRDG"


Table8b[Table8b$Variable == "MARKSxSPLI_SPEEDxKERB_L", "Variable"] <- "MarkxSPD_LIMxKERB_L"
Table8b[Table8b$Variable == "MARKSxHOAL_CURVExKERB_L", "Variable"] <- "MarkxH_CURVExKERB_L"
Table8b[Table8b$Variable == "MARKSxTOTAL_PAVExKERB_L", "Variable"] <- "MarkxTOT_PxKERB_L"
Table8b[Table8b$Variable == "MARKSxTOTAL_SEALxKERB_L", "Variable"] <- "MarkxTOT_SxKERB_L"
Table8b[Table8b$Variable == "MARKSxTRAFFICABLxKERB_L", "Variable"]<- "MarkxTRFABLxKERB_L"
Table8b[Table8b$Variable == "MARKSxNO_OF_LANExKERB_L", "Variable"]<- "MarkxN_LANExKERB_L"
Table8b[Table8b$Variable == "MARKSxSHOULDER_SxKERB_L", "Variable"]<- "MarkxSHLDRxKERB_L"
Table8b[Table8b$Variable == "MARKSxFLOODWAYxKERB_L", "Variable"]<- "MarkxFLDWYxKERB_L"
Table8b[Table8b$Variable == "MARKSxBRIDGExKERB_L", "Variable"]<- "MarkxBRDGxKERB_L"

Table8b[Table8b$Variable == "MARKSxSPLI_SPEEDxSHOULDER_S", "Variable"] <- "MarkxSPD_LIMxSHLDR"
Table8b[Table8b$Variable == "MARKSxHOAL_CURVExSHOULDER_S", "Variable"] <- "MarkxH_CURVExSHLDR"
Table8b[Table8b$Variable == "MARKSxTOTAL_PAVExSHOULDER_S", "Variable"] <- "MarkxTOT_PxSHLDR"
Table8b[Table8b$Variable == "MARKSxTOTAL_SEALxSHOULDER_S", "Variable"] <- "MarkxTOT_SxSHLDR"
Table8b[Table8b$Variable == "MARKSxTRAFFICABLxSHOULDER_S", "Variable"] <- "MarkxTRFABLxSHLDR"
Table8b[Table8b$Variable == "MARKSxNO_OF_LANExSHOULDER_S", "Variable"] <- "MarkxN_LANExSHLDR"
Table8b[Table8b$Variable == "MARKSxKERB_RxSHOULDER_S", "Variable"] <- "MarkxKERB_RxSHLDR"
Table8b[Table8b$Variable == "MARKSxFLOODWAYxSHOULDER_S", "Variable"] <- "MarkxFLDWYxSHLDR"
Table8b[Table8b$Variable == "MARKSxBRIDGExSHOULDER_S", "Variable"] <- "MarkxBRDGxSHLDR"

Table8b[Table8b$Variable == "MARKSxSPLI_SPEEDxKERB_R", "Variable"] <- "MarkxSPD_LIMxKERB_R"
Table8b[Table8b$Variable == "MARKSxHOAL_CURVExKERB_R", "Variable"] <- "MarkxH_CURVExKERB_R"
Table8b[Table8b$Variable == "MARKSxTOTAL_PAVExKERB_R", "Variable"] <- "MarkxTOT_PxKERB_R"
Table8b[Table8b$Variable == "MARKSxTOTAL_SEALxKERB_R", "Variable"] <- "MarkxTOT_SxKERB_R"
Table8b[Table8b$Variable == "MARKSxTRAFFICABLxKERB_R", "Variable"] <- "MarkxTRFABLxKERB_R"
Table8b[Table8b$Variable == "MARKSxNO_OF_LANExKERB_R", "Variable"] <- "MarkxN_LANExKERB_R"
Table8b[Table8b$Variable == "MARKSxFLOODWAYxKERB_R", "Variable"] <- "MarkxFLDWYxKERB_R"
Table8b[Table8b$Variable == "MARKSxBRIDGExKERB_R", "Variable"] <- "MarkxBRDGxKERB_R"
#######################################################################################
Table8 <- cbind(Table8a, Table8b)
#######################################################################################
print(xtable(Table8, digits = 3), include.rownames = FALSE)




