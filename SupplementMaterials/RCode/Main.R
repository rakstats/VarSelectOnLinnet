#########################################################################
# Load the necessary R-packages                                         #
#########################################################################
library(spatstat)
library(spatstat.utils)
library(sp)
library(glmnet)
library(dplyr)
library(stringr)
library(gstat)
#########################################################################
# Source R-script with utility functions                                #
#########################################################################
source("Utils.R")
#########################################################################
#                                                                       #
#                R-Code to create Figure-1 in the paper                 #
#                                                                       #
#########################################################################
## define a function to produce plot files in EPS and PDF
Plot <- function(fname, expr, w=6, h=6, ..., figdir="./figures") {
  fname <- paste0(figdir, .Platform$file.sep, fname)
  pdf(paste0(fname, ".pdf"), width=w, height=h)
  eval(expr)
  dev.off()
  postscript(paste0(fname, ".eps"), width=w, height=h, horiz=FALSE)
  eval(expr)
  dev.off()
  return(paste0(fname, c(".pdf", ".eps")))
}
#####################################
# Load the WA nework and crash data #
#####################################
wa_crash <- readRDS("../Data/Marked_WA_Crash_NetRepaired.rds")

PerthBox <- owin(c(358428.5, 426701), c(6389677, 6523267))
unitname(wa_crash) <- unitname(PerthBox) <- c("metre", "metres")
perth_crash <- wa_crash[PerthBox]

wa_crash_unmarked <- unmark(wa_crash)
perth_crash_unmarked <- unmark(perth_crash)
################################
# Create Figure 1 in the paper #
################################
netargs <- list(col="grey", cols="red", cex=0.4, pch=16)
yard1000 <- yardstick(160000, 8245000, 1160000, 8245000, txt="1000 km")
yard10 <- yardstick(338000, 6500000, 348000, 6500000, txt="10 km")
fig1stuff <- solist(layered(wa_crash_unmarked, PerthBox, yard1000,
                            plotargs=list(netargs,
                                          list(border="green"),
                                          list(angle=90, frac=1/16))),
                    layered(PerthBox, perth_crash_unmarked, yard10,
                            plotargs=list(list(lty=3), netargs, list(angle=90))))
Plot("fig1", {
  par(mar=rep(0.1, 4))
  plot(fig1stuff, main="", main.panel="")
}, 7, 5)
#########################################################################
#                                                                       #
#                R-Code to create Figure-2 in the paper                 #
#                                                                       #
#########################################################################
#####################################
# Load the Chicago nework           #
#####################################
data("chicago")
plot(domain(chicago), main="")
#########################################################################
#                                                                       #
#                R-Code to create Figure-3 in the paper                 #
#                                                                       #
#########################################################################
########################################################
# Create a continuous covariate on the chicago network #
########################################################
sample_linim <- covariate.as.linim(X=domain(chicago), Psill=2.4, Model='Exp',
                                   Range=100, Dummy=TRUE, Beta=1, Nmax=20)
################################
# Plot the simulated covariate #
################################
plot(sample_linim, main = "Simulated covariate")
#########################################################################
# Load the scaled covariates used in the paper for the simulation study #
#########################################################################
cov_linim_scaled <- readRDS("../Data/CovariatesForTheSimulation.rds")
#############################################################
# Load the intensity function for generating point patterns #
#############################################################
lamb <- readRDS("../Data/Lambda.rds")
###############################################################
# Generate a realization of the inhomogeneous Poisson process #
###############################################################
pp_on_chicago <- rpoislpp(lambda = lamb, L = domain(chicago))
################################
# Create Figure-3 in the paper #
################################
op <- par(no.readonly = TRUE)
linim_brk <- get_unified_ribbon_breaks(cov_linim_scaled[1:5])
linim_test <- plot(cov_linim_scaled[[2]], main="")
linim_stuff <- attr(linim_test, "stuff")
linim_cols <- linim_stuff$outputs
col_map <- colourmap(col=linim_cols, breaks = linim_brk)

grid <- rbind(c(1,1,1,1, 2,2,2,2, 3,3,3,3,  8),
              c(1,1,1,1, 2,2,2,2, 3,3,3,3,  4),
              c(1,1,1,1, 2,2,2,2, 3,3,3,3,  4),
              c(10,10,10,10, 11,11,11,11,12,12,12,12, 4), 
              c(5,5,5,5, 6,6,6,6, 7,7,7,7,  4),
              c(5,5,5,5, 6,6,6,6, 7,7,7,7,  4),
              c(5,5,5,5, 6,6,6,6, 7,7,7,7,  4),
              c(13,13,13,13, 14,14,14,14,15,15,15,15,9))
layout(grid)
par(mar = c(0, 1,0,1))
for(pid in 1:3){
  linImage <- cov_linim_scaled[[pid]]
  # index <- as.character(pid)
  # lab <- expression(Z[1])
  plot(linImage, main="", ribbon=FALSE)
}
par(mar = c(0,0,0,2))

plot(col_map, main="", vertical = TRUE, cex.axis = 1.5)

par(mar = c(0, 1,0,1))
for(pid in 4:5){
  linImage <- cov_linim_scaled[[pid]]
  plot(linImage, main="", ribbon=FALSE)
}

plot(pp_on_chicago, main="", ribbon=FALSE, cex = 0.5, pch = 19, alpha = 0.5, col = "gray")
frame()
frame()
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

text(x = 0.5, y = 0.5, labels=expression(italic("Z")["1"]), cex = 1.5, col = "black")
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, labels=expression(italic("Z")["2"]), cex = 1.5, col = "black")
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, labels=expression(italic("Z")["3"]), cex = 1.5, col = "black")
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, labels=expression(italic("Z")["4"]), cex = 1.5, col = "black")
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.5, y = 0.5, labels=expression(italic("Z")["5"]), cex = 1.5, col = "black")
par(op)
#######################################################################################
#########################################################################
#                                                                       #
#                R-Code to create Table-1 and Table-2 in the paper      #
#                                                                       #
#####################################################################################
# We have saved the simulation results as RDS files; these files can be loaded to   #
# reproduce the results of Tables 1 and 2.                                          #
# However, the results can be reproduced by running the following block of codes.   #
############################ Takes long time to run ##################################################
# lamb <- readRDS("../Data/Lambda.rds")                                                              #
# covLinImages_trans <- readRDS("../Data/CovariatesForTheSimulation.rds")                            #
# seed_vec <- seq(100, by=25, length=1000)                                                           #
# lasso_est_df <- list()                                                                             #
# ridge_est_df <- list()                                                                             #
# enet_est_df <- list()                                                                              #
# for(s in seq_along(seed_vec)){                                                                     #
#  sid <- seed_vec[s]                                                                                #
#  set.seed(sid)                                                                                     #
#  ntlpp <- rpoislpp(lambda = lamb)                                                                  #
#  res <- varselect_single_simulate(net_lpp = ntlpp, covLinImages_trans = covLinImages_trans)        #
#  lasso_est_df[[s]] <- res$LASSO_EST                                                                #
#  ridge_est_df[[s]] <- res$RIDGE_EST                                                                #
#  enet_est_df[[s]] <- res$ENET_EST                                                                  #
# }                                                                                                  #
# reg_coef_lasso_all <- do.call("rbind", lasso_est_df)                                               #
# reg_coef_ridge_all <- do.call("rbind", ridge_est_df)                                               #
# reg_coef_enet_all <- do.call("rbind", enet_est_df)                                                 #
######################################################################################################
# Load the LASSO results in Table-1  and Table-2 #
##################################################
reg_coef_lasso_all <- readRDS("../Data/lasso_est_for_paper.rds")
reg_coef_lasso <- reg_coef_lasso_all %>% filter(TuningCoef == "Min")
reg_coef_lasso2 <- get_bias_se_utils(reg_coef_lasso)
######################################
# Print the LASSO results in Table-1 #
######################################
lasso_b1tob5 <- reg_coef_lasso2$bias_se_b1_to_b5
print(lasso_b1tob5)
######################################
# Print the LASSO results in Table-2 #
######################################
lasso_b6tob10 <- reg_coef_lasso2$bias_se_b6_to_b10 %>% 
               rowwise() %>% 
               mutate(Model = paste0(Model, ifelse(is.na(lixel),"",lixel))) %>%
               select(-one_of(c("lixel")))
print(lasso_b6tob10)
##################################################
# Load the RIDGE results in Table-1  and Table-2 #
##################################################
reg_coef_ridge_all <- readRDS("../Data/ridge_est_for_paper.rds")
reg_coef_ridge <- reg_coef_ridge_all %>% filter(TuningCoef == "Min")
reg_coef_ridge2 <- get_bias_se_utils(reg_coef_ridge)
######################################
# Print the RIDGE results in Table-1 #
######################################
ridge_b1tob5 <- reg_coef_ridge2$bias_se_b1_to_b5
print(ridge_b1tob5)
######################################
# Print the RIDGE results in Table-2 #
######################################
ridge_b6tob10 <- reg_coef_ridge2$bias_se_b6_to_b10 %>% 
  rowwise() %>% 
  mutate(Model = paste0(Model, ifelse(is.na(lixel),"",lixel))) %>%
  select(-one_of(c("lixel")))
print(ridge_b6tob10)
##################################################
# Load the ENET results in Table-1  and Table-2  #
##################################################
reg_coef_enet_all <- readRDS("../Data/enet_est_for_paper.rds")
reg_coef_enet <- reg_coef_enet_all %>% filter(TuningCoef == "Min")
reg_coef_enet2 <- get_bias_se_utils(reg_coef_enet)
######################################
# Print the ENET results in Table-1  #
######################################
enet_b1tob5 <- reg_coef_enet2$bias_se_b1_to_b5
print(enet_b1tob5)
######################################
# Print the ENET results in Table-2  #
######################################
enet_b6tob10 <- reg_coef_enet2$bias_se_b6_to_b10 %>% 
  rowwise() %>% 
  mutate(Model = paste0(Model, ifelse(is.na(lixel),"",lixel))) %>%
  select(-one_of(c("lixel")))
print(enet_b6tob10)
#########################################################################
#                                                                       #
#                R-Code to create Table-3 in the paper                  #
#                                                                       #
#########################################################################  
##################################################
# Calculate the proportions for the LASSO method #
##################################################
reg_coef_lasso_prop0 <- get_coef_inclusion_prop(reg_coef_lasso_all)
#########################################################
# Print the proportions in Table-3 for the LASSO method #
#########################################################
print(reg_coef_lasso_prop0)
##################################################
# Calculate the proportions for the RIDGE method #
##################################################
reg_coef_ridge_prop0 <- get_coef_inclusion_prop(reg_coef_ridge_all)
#########################################################
# Print the proportions in Table-3 for the RIDGE method #
#########################################################
print(reg_coef_ridge_prop0)
##################################################
# Calculate the proportions for the ENET method  #
##################################################
reg_coef_enet_prop0 <- get_coef_inclusion_prop(reg_coef_enet_all)
#########################################################
# Print the proportions in Table-3 for the ENET method  #
#########################################################
print(reg_coef_enet_prop0)
########### END of SIMULATION #############################################
########### START APPLICATION #############################################
#########################################################################
#                                                                       #
#                R-Code to create Figure-4 in the paper                 #
#                                                                       #
######################################################################### 
wa_roads <- readRDS("../Data/WA_Linnet_Covar_Data.rds")
L <- wa_roads$lines
M <- marks(L)
road_vars <- c("SPLI_SPEED", "KERB_L", "SHOULDER_S")
V0 <- M[road_vars]
V <- V0
V$KERB_L <- factor(V0$KERB_L)
V$SHOULDER_S <- factor(V0$SHOULDER_S)
func_list <- lapply(V, function(z){function(x,y,seg,tp){z[seg]}})
linfun_list <- lapply(func_list, function(z, net){linfun(z, net)}, net=wa_roads)
#########################################
# Create 3 covariates plots in Figure-4 #
#########################################
plot(linfun_list[["SPLI_SPEED"]], main="", leg.side="top", window=FALSE, leg.wid = 0.025, leg.sep=0.001, zlim=c(40,110))
plot(linfun_list[["KERB_L"]], main="", leg.side="top", window=FALSE, leg.wid=0.025, leg.sep=0.001)
plot(linfun_list[["SHOULDER_S"]], main="", leg.side="top", window=FALSE, leg.wid=0.025, leg.sep=0.001)
######################################################################################
#                                                                                    #
#                R-Code to create Figure-5 and Table-6 in the paper                  #
#                                                                                    #
########################################################################################
# Results in Table-6 are saved in RDS files and to simply reproduce the results        #
# do not run the following block of code; the code to reproduce the results in the     #
# table are provided after the next block. However, to repeat the computation for the  #
# WA dataset, one can use the following block of codes.                                #
########################################################################################
##########################################################################################
#                              R-code block                                              #
##########################################################################################
# road_vars <- c("SPLI_SPEED", "HOAL_CURVE", "TOTAL_PAVE", "TOTAL_SEAL",                 #
#                "TRAFFICABL", "SHOULDER_S", "KERB_L", "KERB_R",                         #
#                "NO_OF_LANE", "FLOODWAY", "BRIDGE")                                     #
# wa_crash <- readRDS("../Data/Marked_WA_Crash_NetRepaired.rds")                         #
# Linnet_with_covars <- readRDS("../Data/WA_Linnet_Covar_Data.rds")                      #
# Lmarks <- marks(Linnet_with_covars$lines)                                              #
# covars.on.seg <- Lmarks[road_vars]                                                     #
# covars_scaled <- covars.on.seg                                                         #
# covars_scaled$SHOULDER_S <- factor(as.integer(covars.on.seg$SHOULDER_S) -1)            #
# covars_scaled$KERB_L <- factor(as.integer(covars.on.seg$KERB_L) -1)                    #
# covars_scaled$KERB_R <- factor(as.integer(covars.on.seg$KERB_R) -1)                    #
# covars_scaled$FLOODWAY <- factor(as.integer(covars.on.seg$FLOODWAY) -1)                #
# covars_scaled$BRIDGE <- factor(as.integer(covars.on.seg$BRIDGE) -1)                    #
# covars_scaled$SPLI_SPEED <- covars.on.seg$SPLI_SPEED/max(covars.on.seg$SPLI_SPEED)     #
# covars_scaled$HOAL_CURVE <- covars.on.seg$HOAL_CURVE/max(covars.on.seg$HOAL_CURVE)     #
# covars_scaled$TOTAL_PAVE <- covars.on.seg$TOTAL_PAVE/max(covars.on.seg$TOTAL_PAVE)     #
# covars_scaled$TOTAL_SEAL <- covars.on.seg$TOTAL_SEAL/max(covars.on.seg$TOTAL_SEAL)     #
# covars_scaled$TRAFFICABL <- covars.on.seg$TRAFFICABL/max(covars.on.seg$TRAFFICABL)     #
# func_list_scaled <- lapply(covars_scaled, function(z){function(x,y,seg,tp){z[seg]}})   #
# linfun_list_scaled <- lapply(func_list_scaled,                                         #
#                              function(z, net){linfun(z, net)},                         #
#                              net=Linnet_with_covars)                                   #  
# form1 <- as.formula(paste("~ " , paste0(names(linfun_list_scaled), collapse = " + "))) #
# nd <- 10000                                                                            #
# Q <- linequad(net.lpp.new, nd=nd)                                                      #
# X <- Q$data                                                                            #
# P <- union.quad(Q)                                                                     #
# prep_all_wa <- mpl.engine(Q, X, P, trend = form1,                                      #
#                           interaction = Poisson(),                                     #
#                           covariates = linfun_list_scaled,                             #  
#                           preponly = TRUE)                                             #
# glmData_all_wa <- prep_all_wa$glmdata                                                  #
# Y_response_all_wa <- glmData_all_wa$.mpl.Y                                             #
# Wts_all_wa <- glmData_all_wa$.mpl.W                                                    #
# covars_wa_massive_df <- glmData_all_wa[3:13]                                           #
# covars_wa_massive_df <- get_derived_covars(covars_wa_massive_df)                       #
# X_cov_all_wa <- as.matrix(covars_wa_massive_df)                                        #
#                                                                                        #
# ipp_lasso_cv_for_paper <- glmnet::cv.glmnet(x = X_cov_all_wa,                          #
#                                             y = Y_response_all_wa,                     #
#                                             weights = Wts_all_wa,                      #
#                                             family = "poisson",                        #
#                                             type.measure = "mse",                      #
#                                             alpha = 1)                                 #
#                                                                                        #
# ipp_ridge_cv_for_paper <- glmnet::cv.glmnet(x = X_cov_all_wa,                          #
#                                             y = Y_response_all_wa,                     #
#                                             weights = Wts_all_wa,                      #
#                                             family = "poisson",                        #
#                                             type.measure = "mse",                      #
#                                             alpha = 0)                                 #
#                                                                                        #
# ipp_enet_cv_for_paper <- glmnet::cv.glmnet(x = X_cov_all_wa,                           #
#                                             y = Y_response_all_wa,                     #
#                                             weights = Wts_all_wa,                      #
#                                             family = "poisson",                        #
#                                             type.measure = "mse",                      #
#                                             alpha = 0.5)                               #
#                                                                                        #
#                                                                                        #
# pts <- as.ppp(wa_crash)                                                                #
# wa_crash2 <- lpp(X=pts, L=Linnet_with_covars)                                          #
# counts_data_all_wa <- get_counts_on_segments(wa_crash2)                                #
# counts_vec_all_wa <- counts_data_all_wa[["ncount"]]                                    #
# length_vec_all_wa <- counts_data_all_wa[["Length"]]                                    #
# indi_vec_all_wa <- counts_data_all_wa[["Indicator"]]                                   #
# log_length_all_wa <- log(length_vec_all_wa)                                            #  
# log_length_centered_all_wa <- log_length_all_wa - mean(log_length_all_wa)              #
# cov_mat_all_Wa <- as.matrix(get_derived_covars(covars_scaled))                         #
#                                                                                        #
# pois_lasso_cv_for_paper <- glmnet::cv.glmnet(x = cov_mat_all_Wa,                       #        
#                                              y = counts_vec_all_wa,                    #
#                                              family = "poisson",                       #
#                                              type.measure = "mse",                     #
#                                              alpha = 1,                                #
#                                              offset = log_length_all_wa)               #
#                                                                                        #
# pois_ridge_cv_for_paper <- glmnet::cv.glmnet(x = cov_mat_all_Wa,                       #        
#                                              y = counts_vec_all_wa,                    #
#                                              family = "poisson",                       #
#                                              type.measure = "mse",                     #
#                                              alpha = 0,                                #
#                                              offset = log_length_all_wa)               #
#                                                                                        #
# pois_enet_cv_for_paper <- glmnet::cv.glmnet(x = cov_mat_all_Wa,                        #        
#                                              y = counts_vec_all_wa,                    #
#                                              family = "poisson",                       #
#                                              type.measure = "mse",                     #
#                                              alpha = 0.5,                              #
#                                              offset = log_length_all_wa)               #
#                                                                                        #
# logit_lasso_cv_for_paper <- glmnet::cv.glmnet(x = cov_mat_all_Wa,                      #
#                                               y = indi_vec_all_wa,                     #
#                                               family = "binomial",                     #
#                                               type.measure = "mse",                    #               
#                                               alpha = 1,                               #
#                                               offset = log_length_centered_all_wa)     #
#                                                                                        #
# logit_ridge_cv_for_paper <- glmnet::cv.glmnet(x = cov_mat_all_Wa,                      #
#                                               y = indi_vec_all_wa,                     #
#                                               family = "binomial",                     #
#                                               type.measure = "mse",                    #               
#                                               alpha = 0,                               #
#                                               offset = log_length_centered_all_wa)     #
#                                                                                        #
# logit_enet_cv_for_paper <- glmnet::cv.glmnet(x = cov_mat_all_Wa,                       #
#                                               y = indi_vec_all_wa,                     #
#                                               family = "binomial",                     #
#                                               type.measure = "mse",                    #               
#                                               alpha = 0.5,                             #
#                                               offset = log_length_centered_all_wa)     #
##########################################################################################
########################################################################
options(scipen = 999)
########################################################################
# B-T method: load the LASSO (cross-validation) result for the WA data #
########################################################################
ipp_lasso_cv_for_paper <- readRDS("../Data/ipp_lasso_cv_for_paper.rds")
############################################################
# Create the top-left plot of the Figure-5                 #
############################################################
plot(ipp_lasso_cv_for_paper, xlab=expression(log(gamma)), ylab="Cross-validation MSE", cex.lab=1.0, cex.axis=1.0)
###################################################
# B-T method: load the LASSO fit for the WA data  #
###################################################
ipp_lasso_fit_for_paper <- readRDS("../Data/ipp_lasso_fit_for_paper.rds")
############################################################
# Create the top-right plot of the Figure-5                #
############################################################
plot(ipp_lasso_fit_for_paper, xvar = "dev", cex.lab=1.0, cex.axis=1.0)
################################################################################
# The column numbers refer to the columns with the coefficient estimates.      #
# Note that the 1st column in Table-6 provides the names of the variables, so  #
# this column is considered 0-th column of the table.                          #
################################################################################
################################################################################
# Print 1st column of Table-6: B-T coefficients computed using LASSO           #
################################################################################
lam1 <- exp(-11)
coef(ipp_lasso_cv_for_paper, s=lam1)
##########################################################################
# LOGIT method: load the LASSO (cross-validation) result for the WA data #
##########################################################################
logit_lasso_cv_for_paper <- readRDS("../Data/logit_lasso_cv_for_paper.rds")
############################################################
# Create the middle-left plot of the Figure-5              #
############################################################
plot(logit_lasso_cv_for_paper, xlab=expression(log(gamma)), ylab="Cross-validation MSE", cex.lab=1.0, cex.axis=1.0)
#####################################################
# LOGIT method: load the LASSO fit for the WA data  #
#####################################################
logit_lasso_fit_for_paper <- readRDS("../Data/logit_lasso_fit_for_paper.rds")
############################################################
# Create the middle-right plot of the Figure-5             #
############################################################
plot(logit_lasso_fit_for_paper, xvar = "dev", cex.lab=1.0, cex.axis=1.0)
##################################################################################
# Print 2nd column of Table-6: LOGIT coefficients computed using LASSO           #
##################################################################################
coef(logit_lasso_cv_for_paper, s="lambda.1se")
#############################################################################
# POISSON method: load the LASSO (cross-validation) result for the WA data  #
#############################################################################
pois_lasso_cv_for_paper <- readRDS("../Data/pois_lasso_cv_for_paper.rds")
############################################################
# Create the bottom-left plot of the Figure-5              #
############################################################
plot(pois_lasso_cv_for_paper, xlab=expression(log(gamma)), ylab="Cross-validation MSE", cex.lab=1.0, cex.axis=1.0)
#######################################################
# POISSON method: load the LASSO fit for the WA data  #
#######################################################
pois_lasso_fit_for_paper <- readRDS("../Data/pois_lasso_fit_for_paper.rds")
############################################################
# Create the bottom-right plot of the Figure-5             #
############################################################
plot(pois_lasso_fit_for_paper, xvar = "dev", cex.lab=1.0, cex.axis=1.0)
##################################################################################
# Print 3rd column of Table-6 : Poisson coefficients computed using LASSO        #
##################################################################################
coef(pois_lasso_cv_for_paper, s="lambda.1se")
##################################################################################
# Print 4th column of Table-6 : B-T coefficients computed using RIDGE            #
##################################################################################
ipp_ridge_cv_for_paper <- readRDS("../Data/ipp_ridge_cv_for_paper.rds")
ipp_ridge_coef_df <- as.data.frame(as.matrix(coef(ipp_ridge_cv_for_paper, s="lambda.min")))
ipp_ridge_coef_df1 <- data.frame(varname = rownames(ipp_ridge_coef_df), ipp_ridge = ipp_ridge_coef_df$`1`)
print(ipp_ridge_coef_df1)
##################################################################################
# Print 5th column of Table-6 : LOGIT coefficients computed using RIDGE          #
##################################################################################
logit_ridge_cv_for_paper <- readRDS("../Data/logit_ridge_cv_for_paper.rds")
logit_ridge_coef_df <- as.data.frame(as.matrix(coef(logit_ridge_cv_for_paper, s="lambda.min")))
logit_ridge_coef_df1 <- data.frame(varname = rownames(logit_ridge_coef_df), logit_ridge = logit_ridge_coef_df$`1`)
print(logit_ridge_coef_df1)
##################################################################################
# Print 6th column of Table-6 : POISSON coefficients computed using RIDGE        #
##################################################################################
pois_ridge_cv_for_paper <- readRDS("../Data/pois_ridge_cv_for_paper.rds")
pois_ridge_coef_df <- as.data.frame(as.matrix(coef(pois_ridge_cv_for_paper, s="lambda.min")))
pois_ridge_coef_df1 <- data.frame(varname = rownames(pois_ridge_coef_df), pois_ridge = pois_ridge_coef_df$`1`)
print(pois_ridge_coef_df1)
##################################################################################
# Print 7th column of Table-6 : B-T coefficients computed using ENET (0.5)       #
##################################################################################
ipp_enet_cv_for_paper <- readRDS("../Data/ipp_enet_cv_for_paper.rds")
ipp_enet_coef_df <- as.data.frame(as.matrix(coef(ipp_enet_cv_for_paper, s=lam1)))
ipp_enet_coef_df1 <- data.frame(varname = rownames(ipp_enet_coef_df), ipp_enet = ipp_enet_coef_df$`1`)
print(ipp_enet_coef_df1)
##################################################################################
# Print 8th column of Table-6 : LOGIT coefficients computed using ENET (0.5)     #
##################################################################################
logit_enet_cv_for_paper <- readRDS("../Data/logit_enet_cv_for_paper.rds")
logit_enet_coef_df <- as.data.frame(as.matrix(coef(logit_enet_cv_for_paper, s="lambda.1se")))
logit_enet_coef_df1 <- data.frame(varname = rownames(logit_enet_coef_df), logit_enet = logit_enet_coef_df$`1`)
print(logit_enet_coef_df1)
###########################################################################################
# Print 9th (last) column of Table-6 : POISSON coefficients computed using ENET (0.5)     #
###########################################################################################
pois_enet_cv_for_paper <- readRDS("../Data/pois_enet_cv_for_paper.rds")
pois_enet_coef_df <- as.data.frame(as.matrix(coef(pois_enet_cv_for_paper, s="lambda.1se")))
pois_enet_coef_df1 <- data.frame(varname = rownames(pois_enet_coef_df), pois_enet = pois_enet_coef_df$`1`)
print(pois_enet_coef_df1)
######################################################################################
#                                                                                    #
#                R-Code to create Figure-6 in the paper                              #
#                                                                                    #
######################################################################################
chic_data <- readRDS("../Data/chicago_new.rds")
plot(chic_data, main="", col="grey", pch=c(1, 24))
######################################################################################
#                                                                                    #
#                R-Code to create Figure-7                                           #
#                                                                                    #
######################################################################################
wa_crash_spilt <- split(wa_crash)
wa_crash_hi <- wa_crash_spilt$High
wa_cash_lo <- wa_crash_spilt$Low
plot(wa_crash_hi, main="High", col="gray", cols="red", pch="+", cex=0.6)
plot(wa_cash_lo, main="Low", col="gray", cols="blue", pch="+", cex=0.6)

wa_split <- split(wa_crash)
pa <- function(i) {
  list(list(pch=16, cols=2),
       list(pch=3, cols=3))[[i]]
}

Plot("fig7", {
  par(mar=rep(0.1, 4))
  plot(wa_split, main="", main.panel="", panel.args=pa, col="grey", cex=0.4)
}, 8, 5)

BB <- owin(c(377728.8, 402213.1),c(6448259, 6468609))
wa_split_BB <- split(wa_crash[BB])
Plot("fig7zoom", {
  par(mar=rep(0.1, 4))
  plot(wa_split_BB, main="", main.panel="", panel.args=pa, col="grey", cex=0.7)
}, 8, 4)
######################################################################################
#                                                                                    #
#                R-Code to create the Table with regression estimates                #
#                for the marked analysis of the WA data (Table-8)                    #
#                                                                                    #
################################################################################################
# We have saved the results that are presented in Table-8 as RDS files. Therefore,             #
# these results can be quickly reproduced by running the R-codes provided after the following  #
# block of R-codes. The next block of R-code shows how to generate these results and might     #
# take few hours to run.                                                                       #
################################################################################################
##########################################################################################
#                              R-code block                                              #
################################################################################################
# Qmarked <- linequad(wa_crash, nd=nd)                                                         #
# vnames <- names(linfun_list_scaled)                                                          #
# form_final <- as.formula(paste("~ " , paste0(vnames, collapse = " + ")))                     #
# ipp.marked.fit.data <- mpl.engine(Qmarked, trend = form_final, interaction = Poisson(),      # 
#                                   covariates = linfun_list_scaled, preponly = TRUE)          #
# glmData_marked_for_ipp <- ipp.marked.fit.data$glmdata                                        #
# Y_response_marked_ipp <- glmData_marked_for_ipp$.mpl.Y                                       #
# Wts_marked_ipp <- glmData_marked_for_ipp$.mpl.W                                              #
# initial_covars_marked_ipp <- glmData_marked_for_ipp[3:14]                                    #   
# glmData_marked_for_ipp <- ipp.marked.fit.data$glmdata                                        #
# Y_response_marked_ipp <- glmData_marked_for_ipp$.mpl.Y                                       #
# Wts_marked_ipp <- glmData_marked_for_ipp$.mpl.W                                              #
# initial_covars_marked_ipp <- glmData_marked_for_ipp[3:14]                                    #
# all_covars_marked_ipp <- get_derived_covars_for_marked_pattern(initial_covars_marked_ipp)    #
# X_marked_ipp <- data.matrix(all_covars_marked_ipp)                                           #
#                                                                                              #
# ipp_marked_lasso_89vars_cv_for_paper <- glmnet::cv.glmnet(x = X_marked_ipp,                  #
#                                                           y = Y_response_marked_ipp,         #
#                                                           weights = Wts_marked_ipp,          #
#                                                           family = "poisson",                #
#                                                           type.measure = "mse",              #  
#                                                           alpha = 1)                         #
#                                                                                              #
# ipp_marked_lasso_89vars_fit_for_paper <-  glmnet::glmnet(x = X_marked_ipp,                   #
#                                                          y = Y_response_marked_ipp,          #
#                                                          weights = Wts_marked_ipp,           #
#                                                          family = "poisson",                 #
#                                                          alpha = 1)                          #
#                                                                                              #
# ipp_marked_ridge_89vars_cv_for_paper <- glmnet::cv.glmnet(x = X_marked_ipp,                  #
#                                                           y = Y_response_marked_ipp,         #
#                                                           weights = Wts_marked_ipp,          #
#                                                           family = "poisson",                #
#                                                           type.measure = "mse",              #
#                                                           alpha = 0)                         #
#                                                                                              #
# ipp_marked_ridge_89vars_fit_for_paper <-  glmnet::glmnet(x = X_marked_ipp,                   #
#                                                          y = Y_response_marked_ipp,          #
#                                                          weights = Wts_marked_ipp,           #  
#                                                          family = "poisson",                 #
#                                                          alpha = 0)                          #
#                                                                                              #
# ipp_marked_enet_89vars_cv_for_paper <- glmnet::cv.glmnet(x = X_marked_ipp,                   #
#                                                           y = Y_response_marked_ipp,         #
#                                                           weights = Wts_marked_ipp,          #
#                                                           family = "poisson",                #
#                                                           type.measure = "mse",              #
#                                                           alpha = 0.5)                       #
#                                                                                              #
# ipp_marked_enet_89vars_fit_for_paper <-  glmnet::glmnet(x = X_marked_ipp,                    # 
#                                                          y = Y_response_marked_ipp,          #
#                                                          weights = Wts_marked_ipp,           #
#                                                          family = "poisson",                 #
#                                                          alpha = 0.5)                        #  
################################################################################################
###########################################################################################
# Print B-T method LASSO results for the marked analysis of the WA data                   #
###########################################################################################
ipp_marked_lasso_cv_for_paper <- readRDS("../Data/ipp_marked_lasso_89vars_cv_for_paper.rds")
coef(ipp_marked_lasso_cv_for_paper, s=0.0000086)
###########################################################################################
# Print B-T method RIDGE results for the marked analysis of the WA data                   #
###########################################################################################
ipp_marked_ridge_cv_for_paper <- readRDS("../Data/ipp_marked_ridge_89vars_cv_for_paper.rds")
coef(ipp_marked_ridge_cv_for_paper, s="lambda.min")
###########################################################################################
# Print B-T method ENET results for the marked analysis of the WA data                    #
###########################################################################################
ipp_marked_enet_cv_for_paper <- readRDS("../Data/ipp_marked_enet_89vars_cv_for_paper.rds")
coef(ipp_marked_enet_cv_for_paper, s=0.0000086)
######################################## END ##############################################














