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
#####################################################
# Load Cross-validation and regularized fit results #
#####################################################
#################################################
#1. Load lasso Cross-validation and lasso fit   #
##############################################################################
# 1a. Load Lasso CV and fit results for Berman-Turner (BT) approximation
lassoBTCV <- readRDS(here::here("RegularizedRegResults", "lassoBT.rds"))
lassoBTFit <- readRDS(here::here("RegularizedRegResults", "lassoBTFit.rds"))
print(lassoBTFit)
print(lassoBTCV)
coef(lassoBTCV, s="lambda.min") # Selected 36 variables
coef(lassoBTCV, s="lambda.1se") # Selected 2 variables
lambdaMin <- lassoBTCV$lambda.min
lambda1Se <- lassoBTCV$lambda.1se
lambdaAvg <- 0.5*(lambdaMin+lambda1Se)
lambdaGmean <- sqrt(lambdaMin*lambda1Se)

coef(lassoBTCV, s=lam) # Selects 18 variables
coef(lassoBTCV, s=lambdaAvg) # Selects 3 variables
coef(lassoBTFit, s=lam) # Selects 18 variables
coef(lassoBTFit, s=lambdaAvg) # Selects 3 variables

lassoBTCoef0 <- as.data.frame(as.matrix(coef(lassoBTCV, s=lambdaGmean)))
names(lassoBTCoef0) <- c("B-T (lasso)")
lassoBTCoef0$Variable <- rownames(lassoBTCoef0)
lassoBTCoef <- lassoBTCoef0[c(2,1)]
###############################################################################
# 1b. Load Lasso CV and fit results for logistic regeression approximation
lassoLogitCV <- readRDS(here::here("RegularizedRegResults", "lassoLogit.rds"))
lassoLogitFit <- readRDS(here::here("RegularizedRegResults", "lassoLogitFit.rds"))
lambdaMinLogit <- lassoLogitCV$lambda.min
lambda1SeLogit <- lassoLogitCV$lambda.1se
lambdaAvgLogit <- 0.5*(lambdaMinLogit+lambda1SeLogit)
coefLamMinLassoLogit <- coef(lassoLogitCV, s="lambda.min") # Selected 43 Variables
coefLam1SeLassoLogit <- coef(lassoLogitCV, s="lambda.1se") # Selected 15 Variables
coefLamAvgLassoLogit <- coef(lassoLogitCV, s=lambdaAvgLogit) # Selected 22 Variables
number_of_nonzero_coef(coefLamMinLassoLogit)
number_of_nonzero_coef(coefLam1SeLassoLogit)
number_of_nonzero_coef(coefLamAvgLassoLogit)
print(lassoLogitCV)

lassoLogitCoef0 <- as.data.frame(as.matrix(coef(lassoLogitCV, s=lambdaAvgLogit)))
names(lassoLogitCoef0) <- c("Logistic (lasso)")
lassoLogitCoef0$Variable <- rownames(lassoLogitCoef0)
lassoLogitCoef <- lassoLogitCoef0[c(2,1)]
# 1c. Load Lasso CV and fit results for poisson count regeression approximation
lassoPoisCV <- readRDS(here::here("RegularizedRegResults", "lassoPois.rds"))
lassoPoisFit <- readRDS(here::here("RegularizedRegResults", "lassoPoisFit.rds"))

lambdaMinPois <- lassoPoisCV$lambda.min
lambda1SePois <- lassoPoisCV$lambda.1se
lambdaAvgPois <- 0.5*(lambdaMinPois+lambda1SePois)
coefLamMinLassoPois <- coef(lassoPoisCV, s="lambda.min") # Selected 33 Variables
coefLam1SeLassoPois <- coef(lassoPoisCV, s="lambda.1se") # Selected 24 Variables
coefLamAvgLassoPois <- coef(lassoPoisCV, s=lambdaAvgPois) # Selected 27 Variables
number_of_nonzero_coef(coefLamMinLassoPois)
number_of_nonzero_coef(coefLam1SeLassoPois)
number_of_nonzero_coef(coefLamAvgLassoPois)
print(lassoPoisCV)
lassoPoisCoef0 <- as.data.frame(as.matrix(coef(lassoPoisCV, s=lambdaAvgPois)))
names(lassoPoisCoef0) <- c("Poisson (lasso)")
lassoPoisCoef0$Variable <- rownames(lassoPoisCoef0)
lassoPoisCoef <- lassoPoisCoef0[c(2,1)]
#################################################
#################################################
#2. Load ridge Cross-validation and ridge fit   #
##############################################################################
# 2a. Load ridge CV and fit results for Berman-Turner (BT) approximation
ridgeBTCV <- readRDS(here::here("RegularizedRegResults", "ridgeBT.rds"))
plot(ridgeBTCV)
coef(ridgeBTCV, s="lambda.min") # Selected 36 variables
coef(ridgeBTCV, s="lambda.1se") # Selected 2 variables
lambdaMinRidge <- ridgeBTCV$lambda.min
lambda1SeRidge <- ridgeBTCV$lambda.1se


coefRidgeBT <- coef(ridgeBTCV, s=lambdaMinRidge)
number_of_nonzero_coef(coefRidgeBT)

ridgeBTCoef0 <- as.data.frame(as.matrix(coefRidgeBT))
names(ridgeBTCoef0) <- c("B-T (ridge)")
ridgeBTCoef0$Variable <- rownames(ridgeBTCoef0)
ridgeBTCoef <- ridgeBTCoef0[c(2,1)]
###############################################################################
###############################################################################
# 2b. Load Ridge CV and fit results for logistic regeression approximation
ridgeLogitCV <- readRDS(here::here("RegularizedRegResults", "ridgeLogit.rds"))
lambdaMinRidgeLogit <- ridgeLogitCV$lambda.min
lambda1SeRidgeLogit <- ridgeLogitCV$lambda.1se
lambdaAvgRidgeLogit <- 0.5*(lambdaMinRidgeLogit+lambda1SeRidgeLogit)

coefRidgeLogit <- coef(ridgeLogitCV, s=lambdaAvgRidgeLogit)
number_of_nonzero_coef(coefRidgeLogit)

ridgeLogitCoef0 <- as.data.frame(as.matrix(coefRidgeLogit))
names(ridgeLogitCoef0) <- c("Logistic (ridge)")
ridgeLogitCoef0$Variable <- rownames(ridgeLogitCoef0)
ridgeLogitCoef <- ridgeLogitCoef0[c(2,1)]
# 2c. Load Ridge CV and fit results for poisson count regeression approximation
ridgePoisCV <- readRDS(here::here("RegularizedRegResults", "ridgePois.rds"))

lambdaMinRidgePois <- ridgePoisCV$lambda.min
lambda1SeRidgePois <- ridgePoisCV$lambda.1se
lambdaAvgRidgePois <- 0.5*(lambdaMinRidgePois+lambda1SeRidgePois)

coefRidgePois <- coef(ridgePoisCV, s=lambdaAvgRidgePois)
number_of_nonzero_coef(coefRidgePois)

ridgePoisCoef0 <- as.data.frame(as.matrix(coefRidgePois))
names(ridgePoisCoef0) <- c("Poisson (ridge)")
ridgePoisCoef0$Variable <- rownames(ridgePoisCoef0)
ridgePoisCoef <- ridgePoisCoef0[c(2,1)]
############################################################################
############################################################################
#################################################
#3. Load e-net Cross-validation and enet fit   #
##############################################################################
# 3a. Load e-net CV and fit results for Berman-Turner (BT) approximation
enetBTCV <- readRDS(here::here("RegularizedRegResults", "enetBT.rds"))
plot(enetBTCV)

lambdaMinEnetBT <- enetBTCV$lambda.min
lambda1SeEnetBT <- enetBTCV$lambda.1se
lambdaGmeanEnetBT <- sqrt(lambdaMinEnetBT*lambda1SeEnetBT)

coefEnetBT <- coef(enetBTCV, s=lambdaGmeanEnetBT)
number_of_nonzero_coef(coefEnetBT)

enetBTCoef0 <- as.data.frame(as.matrix(coefEnetBT))
names(enetBTCoef0) <- c("B-T (e-net)")
enetBTCoef0$Variable <- rownames(enetBTCoef0)
enetBTCoef <- enetBTCoef0[c(2,1)]
###############################################################################
###############################################################################
# 3b. Load e-net CV and fit results for logistic regeression approximation
enetLogitCV <- readRDS(here::here("RegularizedRegResults", "enetLogit.rds"))
lambdaMinEnetLogit <- enetLogitCV$lambda.min
lambda1SeEnetLogit <- enetLogitCV$lambda.1se
lambdaAvgEnetLogit <- 0.5*(lambdaMinEnetLogit+lambda1SeEnetLogit)

coefEnetLogit <- coef(enetLogitCV, s=lambdaAvgEnetLogit)
number_of_nonzero_coef(coefEnetLogit)

enetLogitCoef0 <- as.data.frame(as.matrix(coefEnetLogit))
names(enetLogitCoef0) <- c("Logistic (e-net)")
enetLogitCoef0$Variable <- rownames(enetLogitCoef0)
enetLogitCoef <- enetLogitCoef0[c(2,1)]
# 3c. Load enet CV and fit results for poisson count regeression approximation
enetPoisCV <- readRDS(here::here("RegularizedRegResults", "enetPois.rds"))

lambdaMinEnetPois <- enetPoisCV$lambda.min
lambda1SeEnetPois <- enetPoisCV$lambda.1se
lambdaAvgEnetPois <- 0.5*(lambdaMinEnetPois+lambda1SeEnetPois)

coefEnetPois <- coef(enetPoisCV, s=lambdaAvgEnetPois)
number_of_nonzero_coef(coefEnetPois)

enetPoisCoef0 <- as.data.frame(as.matrix(coefEnetPois))
names(enetPoisCoef0) <- c("Poisson (e-net)")
enetPoisCoef0$Variable <- rownames(enetPoisCoef0)
enetPoisCoef <- enetPoisCoef0[c(2,1)]
############################################################################
logLength <- readRDS(here::here("DerivedData", "logLengthForDiscretizedModels.rds"))
meanLogLength <- mean(logLength)
lassoLogitCoef[1,2] <- lassoLogitCoef[1,2] - meanLogLength
ridgeLogitCoef[1,2] <- ridgeLogitCoef[1,2] - meanLogLength
enetLogitCoef[1,2] <- enetLogitCoef[1,2] - meanLogLength

Col1 <- lassoBTCoef
Col2 <- lassoLogitCoef
Col3 <- lassoPoisCoef
Col4 <- ridgeBTCoef
Col5 <- ridgeLogitCoef
Col6 <- ridgePoisCoef
Col7 <- enetBTCoef
Col8 <- enetLogitCoef
Col9 <- enetPoisCoef
############################################################################
Table6 <- Col1 %>% left_join(Col2, by = "Variable") %>%
          left_join(Col3, by = "Variable") %>%
          left_join(Col4, by = "Variable") %>%
          left_join(Col5, by = "Variable") %>%
          left_join(Col6, by = "Variable") %>%
          left_join(Col7, by = "Variable") %>%
          left_join(Col8, by = "Variable") %>%
          left_join(Col9, by = "Variable")
Table6[Table6$Variable == "SPLI_SPEED", "Variable"] <- "SPD_LIM"
Table6[Table6$Variable == "HOAL_CURVE", "Variable"] <- "H_CURVE"
Table6[Table6$Variable == "TOTAL_PAVE", "Variable"] <- "TOT_P"
Table6[Table6$Variable == "TOTAL_SEAL", "Variable"] <- "TOT_S"
Table6[Table6$Variable == "TRAFFICABL", "Variable"] <- "TRFABL"
Table6[Table6$Variable == "NO_OF_LANE", "Variable"] <- "N_LANE"

Table6[Table6$Variable == "SPLI_SPEED2", "Variable"] <- "SPD_LIM2"
Table6[Table6$Variable == "HOAL_CURVE2", "Variable"] <- "H_CURVE2"
Table6[Table6$Variable == "TOTAL_PAVE2", "Variable"] <- "TOT_P2"
Table6[Table6$Variable == "TOTAL_SEAL2", "Variable"] <- "TOT_S2"
Table6[Table6$Variable == "TRAFFICABL2", "Variable"] <- "TRFABL2"
Table6[Table6$Variable == "NO_OF_LANE2", "Variable"] <- "N_LANE2"

Table6[Table6$Variable == "SHOULDER_S", "Variable"] <- "SHLDR"
Table6[Table6$Variable == "FLOODWAY", "Variable"] <- "FLDWY"
Table6[Table6$Variable == "BRIDGE", "Variable"] <- "BRDG"


Table6[Table6$Variable == "SPLI_SPEEDxKERB_L", "Variable"] <- "SPD_LIMxKERB_L"
Table6[Table6$Variable == "HOAL_CURVExKERB_L", "Variable"] <- "H_CURVExKERB_L"
Table6[Table6$Variable == "TOTAL_PAVExKERB_L", "Variable"] <- "TOT_PxKERB_L"
Table6[Table6$Variable == "TOTAL_SEALxKERB_L", "Variable"] <- "TOT_SxKERB_L"
Table6[Table6$Variable == "TRAFFICABLxKERB_L", "Variable"]<- "TRFABLxKERB_L"
Table6[Table6$Variable == "NO_OF_LANExKERB_L", "Variable"]<- "N_LANExKERB_L"
Table6[Table6$Variable == "SHOULDER_SxKERB_L", "Variable"]<- "SHLDRxKERB_L"
Table6[Table6$Variable == "FLOODWAYxKERB_L", "Variable"]<- "FLDWYxKERB_L"
Table6[Table6$Variable == "BRIDGExKERB_L", "Variable"]<- "BRDGxKERB_L"

Table6[Table6$Variable == "SPLI_SPEEDxSHOULDER_S", "Variable"] <- "SPD_LIMxSHLDR"
Table6[Table6$Variable == "HOAL_CURVExSHOULDER_S", "Variable"] <- "H_CURVExSHLDR"
Table6[Table6$Variable == "TOTAL_PAVExSHOULDER_S", "Variable"] <- "TOT_PxSHLDR"
Table6[Table6$Variable == "TOTAL_SEALxSHOULDER_S", "Variable"] <- "TOT_SxSHLDR"
Table6[Table6$Variable == "TRAFFICABLxSHOULDER_S", "Variable"] <- "TRFABLxSHLDR"
Table6[Table6$Variable == "NO_OF_LANExSHOULDER_S", "Variable"] <- "N_LANExSHLDR"
Table6[Table6$Variable == "KERB_RxSHOULDER_S", "Variable"] <- "KERB_RxSHLDR"
Table6[Table6$Variable == "FLOODWAYxSHOULDER_S", "Variable"] <- "FLDWYxSHLDR"
Table6[Table6$Variable == "BRIDGExSHOULDER_S", "Variable"] <- "BRDGxSHLDR"

Table6[Table6$Variable == "SPLI_SPEEDxKERB_R", "Variable"] <- "SPD_LIMxKERB_R"
Table6[Table6$Variable == "HOAL_CURVExKERB_R", "Variable"] <- "H_CURVExKERB_R"
Table6[Table6$Variable == "TOTAL_PAVExKERB_R", "Variable"] <- "TOT_PxKERB_R"
Table6[Table6$Variable == "TOTAL_SEALxKERB_R", "Variable"] <- "TOT_SxKERB_R"
Table6[Table6$Variable == "TRAFFICABLxKERB_R", "Variable"] <- "TRFABLxKERB_R"
Table6[Table6$Variable == "NO_OF_LANExKERB_R", "Variable"] <- "N_LANExKERB_R"
Table6[Table6$Variable == "FLOODWAYxKERB_R", "Variable"] <- "FLDWYxKERB_R"
Table6[Table6$Variable == "BRIDGExKERB_R", "Variable"] <- "BRDGxKERB_R"



############################################################################
print(xtable(Table6, digits = 3), include.rownames = FALSE)
##############################################################################
##############################################
# glmnet test for offset effect on intercept #
##############################################
set.seed(100)
res <- sample(x=c(0,1), size=100,  replace=TRUE)
res <- factor(sort(res))
cov <- c(rnorm(50, mean=2), rnorm(50, mean=5))
cov2 <- matrix(cbind(rep(1,100), cov), ncol=2)
offst <- log(abs(rnorm(100)))
meanOffst <- mean(offst)
centeredOffst <- offst - meanOffst
glmFit <- glmnet(x=cov2, y=res, offset = offst, family = "binomial")
coef(glmFit, s=0.0007790)
glmFit2 <- glmnet(x=cov2, y=res, offset = centeredOffst, family = "binomial")
coef(glmFit2, s=0.0007790)
-10.269886 # Intercept without centering
-10.819482 # Intercept with centering
-10.819482 - meanOffst # Reltionship: Intercept with and without centering
#########################################################
# Let's compute the fraction deviance
# explained by different regularization parameters
########################################################
cvVec <- list(lassoBTCV, lassoLogitCV, lassoPoisCV,
              ridgeBTCV, ridgeLogitCV, ridgePoisCV,
              enetBTCV, enetLogitCV, enetPoisCV)
methodVec <- rep(c("lasso", "ridge", "e-net"), each=3)
approxVec <- rep(c("B-T", "Logistic", "Poisson"), 3)
########
modelfitList <- mapply(FUN = get_model_fit_info, cvVec, methodVec, approxVec, SIMPLIFY = FALSE)
########
modelfitDf <- do.call("rbind", modelfitList)
########
print(xtable(modelfitDf, digits = 2), include.rownames = FALSE)





