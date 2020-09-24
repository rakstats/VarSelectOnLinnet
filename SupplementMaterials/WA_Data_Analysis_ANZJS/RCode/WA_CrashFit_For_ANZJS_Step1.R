#############################################################################################
#                                                                                           #
# This R-code is for preparing the design matrix (X) and and response vector (y) for        #
# fitting the three regularized models (lasso, ridge and enet) using the B-T approximation, #
# logistic regression and Poisson count regression.                                         #
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
####################################################
# Load the 2011 crash data on the WA road network  #
####################################################
netLpp <- readRDS(here::here("InitialData", "net_lpp.rds"))
linnetWithCovars <- readRDS(here::here("InitialData","WA_Linnet_Covar_Data.rds"))
#########################
# Extract the variables #
#########################
roadVars <- c("SPLI_SPEED", "HOAL_CURVE", "TOTAL_PAVE", 
              "TOTAL_SEAL","TRAFFICABL", "NO_OF_LANE",
              "SHOULDER_S", "KERB_L", "KERB_R", "FLOODWAY", "BRIDGE")
###########################
lmarks <- marks(linnetWithCovars$lines)
covarsOnSeg <- lmarks[roadVars]
###########################
##############################################################################################################
# Variables need scaling
#############################################################################################################
# Compute the max and min of the continuous variables. These values will be used form scaling               #
#############################################################################################################
contVars <- covarsOnSeg[c("SPLI_SPEED", "HOAL_CURVE", "TOTAL_PAVE", 
                          "TOTAL_SEAL","TRAFFICABL", "NO_OF_LANE")]
contVars2Names <- paste0(names(contVars), "2")
contVars2 <- as.data.frame(sapply(contVars, function(x){x^2}))
names(contVars2) <- contVars2Names
contVarsWithQuadTerms <- cbind(contVars, contVars2)
minContVars <- sapply(contVarsWithQuadTerms, "min")
maxContVars <- sapply(contVarsWithQuadTerms, "max")
scaleContVars <- maxContVars - minContVars
scaledContCovars <- as.data.frame(mapply(function(x, xmin, xscale){(x-xmin)/xscale},
                                         contVarsWithQuadTerms, minContVars, scaleContVars))
summary(scaledContCovars)
# Dummy coding of the factors
facVars <- covarsOnSeg[c("SHOULDER_S", "KERB_L", "KERB_R", "FLOODWAY", "BRIDGE")]
dummyCodedVars <- as.data.frame(sapply(facVars, function(x){as.integer(x)-1}))
# Combine continuous and factor variables
waCovarsOnSeg <- cbind(scaledContCovars, dummyCodedVars)
#####################################################################################################################
# Create the lpp object using Linnet_with_covars
#####################################################################################################################
pts <- as.ppp(netLpp)
scaledLinnetWithCovars <- linnetWithCovars
marks(scaledLinnetWithCovars$lines) <- waCovarsOnSeg
newNetLpp <- lpp(X=pts, L=scaledLinnetWithCovars)
nvertices(scaledLinnetWithCovars)
nsegments(newNetLpp)
sum(lengths.psp(as.psp(scaledLinnetWithCovars)))/1000
###################################################################################################################
funcListScaled <- lapply(waCovarsOnSeg, function(z){function(x,y,seg,tp){z[seg]}})
linfunList <- lapply(funcListScaled, function(z, net){linfun(z, net)}, net=scaledLinnetWithCovars)
saveRDS(linfunList, here::here("DerivedData", "linfun_list_scaled.rds"))
###################################################################################################################
# Get the data using the B-T quadrature scheme                                                                    # 
###################################################################################################################
vnames <- names(linfunList)
form1 <- as.formula(paste("~ " , paste0(vnames, collapse = " + ")))
nd <- 10000 
Q <- linequad(newNetLpp, nd=nd)
prepAllWAMassive <- mpl.engine(Q, trend = form1, 
                               interaction = Poisson(), 
                               covariates = linfunList, 
                               preponly = TRUE)
saveRDS(prepAllWAMassive, here::here("DerivedData", "prepAllWAMassive.rds"))
# Covariates, responses and weights for the B-T approximation
glmDataAllWa <- prepAllWAMassive$glmdata
yResponseforBT <- glmDataAllWa$.mpl.Y
wtsforBT <- glmDataAllWa$.mpl.W
covarForBTMethod <- glmDataAllWa[3:19]
all44CovarsForBTMethod <- get_interaction_terms(covarForBTMethod)
X44forBT <- data.matrix(all44CovarsForBTMethod)
saveRDS(X44forBT, here::here("DerivedData", "covMatForBTModel.rds"))
saveRDS(yResponseforBT, here::here("DerivedData", "responseForBTModel.rds"))
saveRDS(wtsforBT, here::here("DerivedData", "weightsForBTModel.rds"))
##################################################################################################################
# Get the data for fitting discretized model                                                                     #
##################################################################################################################
waAll44CovarsOnSeg <- get_interaction_terms(waCovarsOnSeg)
countsAllWa <- get_counts_on_segments(newNetLpp)
countsVecWA <- countsAllWa[["ncount"]]
lengthVecWA <- countsAllWa[["Length"]]
indiVecWa <- countsAllWa[["Indicator"]]
logLengthAllWa <- log(lengthVecWA)
logLengthCenteredAllWa <- logLengthAllWa - mean(logLengthAllWa)

X44 <- data.matrix(waAll44CovarsOnSeg)

saveRDS(X44, here::here("DerivedData", "covMatForDiscretizedModels.rds"))
saveRDS(countsVecWA, here::here("DerivedData", "countResponseForDiscretizedModels.rds"))
saveRDS(indiVecWa, here::here("DerivedData", "indicatorResponseForDiscretizedModels.rds"))
saveRDS(lengthVecWA, here::here("DerivedData", "segLengthForDiscretizedModels.rds"))
saveRDS(logLengthAllWa, here::here("DerivedData", "logLengthForDiscretizedModels.rds"))
saveRDS(logLengthCenteredAllWa, here::here("DerivedData", "centredLogLengthForDiscretizedModels.rds"))
#################################################################################################################





