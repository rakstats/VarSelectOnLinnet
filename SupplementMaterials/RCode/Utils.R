########################################################################################
# Codes for the variable selection methods for point patterns on a linear network      #
# Author: Suman Rakshit                                                                #                                                                    #
# Curtin University                                                                    #
########################################################################################
# varselect_simulate_all
# no_of_cols = (Method, lixel, Model, MeasureType, TuningCoef, b1:b10)
# no_of_rows will be for lixels (discretised): 2 discretised model * 3 Eps value *
#   3 Tuning coef = 18
# no_of_rows for the IPP model is 3 IPP model * 3 Tuning coef = 9
# 18 + 9 = 27 rows in single simulation for 1 method
# 27 * 3 = 81 rows in single simulation for 3 methods of var selection
# 81,000 rows will be in total for the entire simulation of size 1000
########################################################################################
varselect_single_simulate <- function(net_lpp,
                                      covLinImages_trans,
                                      no_of_rows = 27,
                                      no_of_dismod_rows = 6,
                                      retcolnames = c("Method", "lixel", "Model",
                                                      "TuningCoef",
                                                      "b1", "b2", "b3", "b4", "b5",
                                                      "b6", "b7", "b8", "b9", "b10"),
                                      EPS = c(20, 10, 5),
                                      ND = c(1000, 2500, 5000))
{
  no_of_cols <- length(retcolnames)
  retdf_lasso <- as.data.frame(matrix(NA, nrow = no_of_rows, ncol = no_of_cols))
  retdf_ridge <- as.data.frame(matrix(NA, nrow = no_of_rows, ncol = no_of_cols))
  retdf_enet <- as.data.frame(matrix(NA, nrow = no_of_rows, ncol = no_of_cols))
  
  names(retdf_lasso) <- retcolnames
  names(retdf_ridge) <- retcolnames
  names(retdf_enet) <- retcolnames
  
  retdf_lasso[["Method"]] <- "Lasso"
  retdf_ridge[["Method"]] <- "Ridge"
  retdf_enet[["Method"]] <- "E-net"
  
  for(e in seq_along(EPS)){
    ep <- EPS[e]
    net.lpp.eps <- lixellate(net_lpp, eps = ep)
    nsegs <- nsegments(as.linnet(net.lpp.eps))
    
    istart <- (e-1) * no_of_dismod_rows + 1
    iend <- e * no_of_dismod_rows
    
    retdf_lasso[istart:iend, c("lixel")] <- nsegs
    retdf_ridge[istart:iend, c("lixel")] <- nsegs
    retdf_enet[istart:iend, c("lixel")] <- nsegs
    
    dom <- domain(net.lpp.eps)
    midpts <- get_midpt_of_segments(dom)
    cov_values <- lapply(covLinImages_trans, FUN = function(z){z[midpts]})
    cov_mat <- do.call("cbind", cov_values)
    
    counts_data <- get_counts_on_segments(net.lpp.eps)
    counts_vec <- counts_data[["ncount"]]
    length_vec <- counts_data[["Length"]]
    indi_vec <- counts_data[["Indicator"]]
    log_length <- log(length_vec)
    log_length_centered <- log_length - mean(log_length)
    ################# Start Lasso #####################################
    # Fit Poisson Count Model
    glm_dev_pois_vselect <- glmnet::cv.glmnet(x = cov_mat, 
                                              y = counts_vec,
                                              family = "poisson", 
                                              type.measure = "mse",
                                              alpha = 1,
                                              offset = log_length)
    
    retdf_lasso[(istart:(istart+2)), c("Model")] <- "poisson"
    get_coef <- get_coef_indicator(glm_dev_pois_vselect)
    retdf_lasso[istart, c("TuningCoef")] <- "Min"
    retdf_lasso[istart, 5:no_of_cols] <- get_coef[["coef_min"]]
    retdf_lasso[(istart+1), c("TuningCoef")] <- "1se"
    retdf_lasso[(istart+1), 5:no_of_cols] <- get_coef[["coef_1se"]]
    retdf_lasso[(istart+2), c("TuningCoef")] <- "Avg"
    retdf_lasso[istart+2, 5:no_of_cols] <- get_coef[["coef_avg"]]
    # Fit Logistic Regression Model
    glm_dev_bino_vselect <- glmnet::cv.glmnet(x = cov_mat, 
                                              y = as.factor(indi_vec),
                                              family = "binomial", 
                                              type.measure = "mse",
                                              offset = log_length_centered,
                                              alpha = 1)
    istart <- (istart+3)
    retdf_lasso[(istart:iend), c("Model")] <- "logit"
    get_coef <- get_coef_indicator(glm_dev_bino_vselect)
    retdf_lasso[istart, c("TuningCoef")] <- "Min"
    retdf_lasso[istart, 5:no_of_cols] <- get_coef[["coef_min"]]
    retdf_lasso[(istart+1), c("TuningCoef")] <- "1se"
    retdf_lasso[(istart+1), 5:no_of_cols] <- get_coef[["coef_1se"]]
    retdf_lasso[(istart+2), c("TuningCoef")] <- "Avg"
    retdf_lasso[istart+2, 5:no_of_cols] <- get_coef[["coef_avg"]]
    ######### LASSO Finished ####################################################
    ######### Start Ridge #######################################################
    # Fit Poisson Count Model
    glm_dev_pois_vselect <- glmnet::cv.glmnet(x = cov_mat, 
                                              y = counts_vec,
                                              family = "poisson", 
                                              type.measure = "mse",
                                              alpha = 0,
                                              offset = log_length)
    istart <- (e-1) * no_of_dismod_rows + 1
    
    retdf_ridge[(istart:(istart+2)), c("Model")] <- "poisson"
    get_coef <- get_coef_indicator(glm_dev_pois_vselect)
    retdf_ridge[istart, c("TuningCoef")] <- "Min"
    retdf_ridge[istart, 5:no_of_cols] <- get_coef[["coef_min"]]
    retdf_ridge[(istart+1), c("TuningCoef")] <- "1se"
    retdf_ridge[(istart+1), 5:no_of_cols] <- get_coef[["coef_1se"]]
    retdf_ridge[(istart+2), c("TuningCoef")] <- "Avg"
    retdf_ridge[istart+2, 5:no_of_cols] <- get_coef[["coef_avg"]]
    # Fit Logistic Regression Model
    glm_dev_bino_vselect <- glmnet::cv.glmnet(x = cov_mat, 
                                              y = as.factor(indi_vec),
                                              family = "binomial", 
                                              type.measure = "mse",
                                              offset = log_length_centered,
                                              alpha = 0)
    istart <- (istart+3)
    retdf_ridge[(istart:iend), c("Model")] <- "logit"
    get_coef <- get_coef_indicator(glm_dev_bino_vselect)
    retdf_ridge[istart, c("TuningCoef")] <- "Min"
    retdf_ridge[istart, 5:no_of_cols] <- get_coef[["coef_min"]]
    retdf_ridge[(istart+1), c("TuningCoef")] <- "1se"
    retdf_ridge[(istart+1), 5:no_of_cols] <- get_coef[["coef_1se"]]
    retdf_ridge[(istart+2), c("TuningCoef")] <- "Avg"
    retdf_ridge[istart+2, 5:no_of_cols] <- get_coef[["coef_avg"]]
    ######### Ridge Finished ####################################################
    ######### Start Enet #######################################################
    # Fit Poisson Count Model
    glm_dev_pois_vselect <- glmnet::cv.glmnet(x = cov_mat, 
                                              y = counts_vec,
                                              family = "poisson", 
                                              type.measure = "mse",
                                              alpha = 0.5,
                                              offset = log_length)
    istart <- (e-1) * no_of_dismod_rows + 1
    
    retdf_enet[(istart:(istart+2)), c("Model")] <- "poisson"
    get_coef <- get_coef_indicator(glm_dev_pois_vselect)
    retdf_enet[istart, c("TuningCoef")] <- "Min"
    retdf_enet[istart, 5:no_of_cols] <- get_coef[["coef_min"]]
    retdf_enet[(istart+1), c("TuningCoef")] <- "1se"
    retdf_enet[(istart+1), 5:no_of_cols] <- get_coef[["coef_1se"]]
    retdf_enet[(istart+2), c("TuningCoef")] <- "Avg"
    retdf_enet[istart+2, 5:no_of_cols] <- get_coef[["coef_avg"]]
    # Fit Logistic Regression Model
    glm_dev_bino_vselect <- glmnet::cv.glmnet(x = cov_mat, 
                                              y = as.factor(indi_vec),
                                              family = "binomial", 
                                              type.measure = "mse",
                                              offset = log_length_centered,
                                              alpha = 0.5)
    istart <- (istart+3)
    retdf_enet[(istart:iend), c("Model")] <- "logit"
    get_coef <- get_coef_indicator(glm_dev_bino_vselect)
    retdf_enet[istart, c("TuningCoef")] <- "Min"
    retdf_enet[istart, 5:no_of_cols] <- get_coef[["coef_min"]]
    retdf_enet[(istart+1), c("TuningCoef")] <- "1se"
    retdf_enet[(istart+1), 5:no_of_cols] <- get_coef[["coef_1se"]]
    retdf_enet[(istart+2), c("TuningCoef")] <- "Avg"
    retdf_enet[istart+2, 5:no_of_cols] <- get_coef[["coef_avg"]]
  }
  
  ## Start fitting lppm models
  # Fit inhomogeneous Poisson process model
  for(k in seq_along(ND)){
    nd <- ND[k]
    Q <- linequad(net_lpp, nd = nd)
    X <- Q$data
    P <- union.quad(Q)
    form1 <- as.formula(paste("~", paste0("v", 1:10, collapse = "+")))
    prep <- mpl.engine(Q, X, P, trend = form1, interaction = Poisson(),
                       covariates = covLinImages_trans,
                       preponly = TRUE)
    glmData <- prep$glmdata
    Y_response <- glmData$.mpl.Y
    Wts <- glmData$.mpl.W
    X_cov <- data.matrix(glmData[3:(10+2)])
    idk <- (length(EPS) * no_of_dismod_rows) + ((k-1)*3 + 1)
    ### LASSO #################################################
    retdf_lasso[(idk:(idk+2)), c("Model")] <- paste0("lppm",nd)
    glm_dev_ihom <- glmnet::cv.glmnet(x = X_cov, 
                                      y = Y_response, 
                                      weights = Wts,
                                      family = "poisson", 
                                      type.measure = "mse",
                                      alpha = 1)
    get_coef <- get_coef_indicator(glm_dev_ihom)
    retdf_lasso[idk, c("TuningCoef")] <- "Min"
    retdf_lasso[idk, 5:no_of_cols] <- get_coef[["coef_min"]]
    retdf_lasso[(idk+1), c("TuningCoef")] <- "1se"
    retdf_lasso[(idk+1), 5:no_of_cols] <- get_coef[["coef_1se"]]
    retdf_lasso[(idk+2), c("TuningCoef")] <- "Avg"
    retdf_lasso[(idk+2), 5:no_of_cols] <- get_coef[["coef_avg"]]
    ### RIDGE ####################################################
    retdf_ridge[(idk:(idk+2)), c("Model")] <- paste0("lppm",nd)
    glm_dev_ihom <- glmnet::cv.glmnet(x = X_cov, 
                                      y = Y_response, 
                                      weights = Wts,
                                      family = "poisson", 
                                      type.measure = "mse",
                                      alpha = 0)
    get_coef <- get_coef_indicator(glm_dev_ihom)
    retdf_ridge[idk, c("TuningCoef")] <- "Min"
    retdf_ridge[idk, 5:no_of_cols] <- get_coef[["coef_min"]]
    retdf_ridge[(idk+1), c("TuningCoef")] <- "1se"
    retdf_ridge[(idk+1), 5:no_of_cols] <- get_coef[["coef_1se"]]
    retdf_ridge[(idk+2), c("TuningCoef")] <- "Avg"
    retdf_ridge[(idk+2), 5:no_of_cols] <- get_coef[["coef_avg"]]
    ### ENET ####################################################
    retdf_enet[(idk:(idk+2)), c("Model")] <- paste0("lppm",nd)
    glm_dev_ihom <- glmnet::cv.glmnet(x = X_cov, 
                                      y = Y_response, 
                                      weights = Wts,
                                      family = "poisson", 
                                      type.measure = "mse",
                                      alpha = 0.5)
    get_coef <- get_coef_indicator(glm_dev_ihom)
    retdf_enet[idk, c("TuningCoef")] <- "Min"
    retdf_enet[idk, 5:no_of_cols] <- get_coef[["coef_min"]]
    retdf_enet[(idk+1), c("TuningCoef")] <- "1se"
    retdf_enet[(idk+1), 5:no_of_cols] <- get_coef[["coef_1se"]]
    retdf_enet[(idk+2), c("TuningCoef")] <- "Avg"
    retdf_enet[(idk+2), 5:no_of_cols] <- get_coef[["coef_avg"]]
  }
  ret_list <- list(LASSO_EST = retdf_lasso,
                   RIDGE_EST = retdf_ridge,
                   ENET_EST = retdf_enet)
  return(ret_list)
}
###################### covariate.as.linim ##############################################
varselect_simulate <- function(nsimulate, 
                               lamb, 
                               no_of_rows, 
                               no_of_cols,
                               meth,
                               alph,
                               type_measure,
                               coef_what){
  reg_coef_estimate_dlist <- list()
  for(i in 1:nsimulate){
    net.lpp <- rpoislpp(lambda = lamb)
    retdf <- as.data.frame(matrix(NA, nrow = no_of_rows, ncol = no_of_cols))
    names(retdf) <- retcolnames
    retdf[["Method"]] <- meth
    
    for(e in seq_along(EPS)){
      ep <- EPS[e]
      net.lpp.eps <- lixellate(net.lpp, eps = ep)
      nsegs <- nsegments(as.linnet(net.lpp.eps))
      
      istart <- (e-1) * no_of_dismod_rows + 1
      iend <- e * no_of_dismod_rows
      retdf[istart:iend, c("lixel")] <- nsegs
      
      dom <- domain(net.lpp.eps)
      midpts <- get_midpt_of_segments(dom)
      cov_values <- lapply(covLinImages_trans, FUN = function(z){z[midpts]})
      cov_mat <- do.call("cbind", cov_values)
      
      counts_data <- get_counts_on_segments(net.lpp.eps)
      counts_vec <- counts_data[["ncount"]]
      length_vec <- counts_data[["Length"]]
      indi_vec <- counts_data[["Indicator"]]
      log_length <- log(length_vec)
      log_length_centered <- log_length - mean(log_length)
      # Fit Poisson Count Model
      glm_dev_pois_vselect <- glmnet::cv.glmnet(x = cov_mat, 
                                                y = counts_vec,
                                                family = "poisson", 
                                                type.measure = type_measure,
                                                alpha = alph,
                                                offset = log_length)
      retdf[istart, c("Model")] <- dismod[1]
      retdf[istart, 4:no_of_cols] <- get_coef_indicator(glm_dev_pois_vselect)[[coef_what]]
      # Fit Logistic Regression Model
      glm_dev_bino_vselect <- glmnet::cv.glmnet(x = cov_mat, 
                                                y = as.factor(indi_vec),
                                                family = "binomial", 
                                                type.measure = type_measure,
                                                offset = log_length_centered,
                                                alpha = alph)
      retdf[(istart+1), c("Model")] <- dismod[2]
      retdf[(istart+1), 4:no_of_cols] <- get_coef_indicator(glm_dev_bino_vselect)[[coef_what]]
    }
    
    #ihomidstart <- length(EPS) * no_of_dismod_rows + 1
    # Fit inhomogeneous Poisson process model
    for(k in seq_along(ND)){
      nd <- ND[k]
      Q <- linequad(net.lpp, nd = nd)
      X <- Q$data
      P <- union.quad(Q)
      form1 <- as.formula(paste("~", paste0("v", 1:10, collapse = "+")))
      prep <- mpl.engine(Q, X, P, trend = form1, interaction = Poisson(),
                         covariates = covLinImages_trans,
                         preponly = TRUE)
      glmData <- prep$glmdata
      Y_response <- glmData$.mpl.Y
      Wts <- glmData$.mpl.W
      X_cov <- data.matrix(glmData[3:(10+2)])
      idk <- (length(EPS) * no_of_dismod_rows) + k
      retdf[idk, c("Model")] <- ihomod[k]
      glm_dev_ihom <- glmnet::cv.glmnet(x = X_cov, 
                                        y = Y_response, 
                                        weights = Wts,
                                        family = "poisson", 
                                        type.measure = type_measure,
                                        alpha = alph)
      retdf[idk, 4:no_of_cols] <- get_coef_indicator(glm_dev_ihom)[[coef_what]]
    }
    reg_coef_estimate_dlist[[i]] <- retdf
    print(i)
  }
  
  reg_coef_est <- do.call("rbind", reg_coef_estimate_dlist)
  return(reg_coef_est)
}
##############################
covariate.as.linim <- function(X,
                               Psill=NA,
                               Model='Sph',
                               Range=NA,
                               Dummy=TRUE,
                               Beta=1,
                               Nmax,
                               ...){
  
  stopifnot(!missing(X) && inherits(X, "linnet"))
  stopifnot(isNamespaceLoaded(name="gstat"))
  argh <- list(...)
  # Get the window corrensponding to the network.
  Win <- X$window
  # Create a mask of Win.
  Mask <- do.call.matched(as.mask, resolve.defaults(list(w=Win), argh))
  # Create a grid structure of Win (data.frame with x and y columns).
  xyGrid <- expand.grid(Mask$xcol,Mask$yrow)
  names(xyGrid) <- c('x','y')
  # Defining the Gaussian Random field model.
  rfModel <- do.call.matched(gstat::vgm, 
                             resolve.defaults(argh, list(psill=Psill,model=Model,range=Range)))
  
  g.dummy <- do.call.matched(gstat::gstat, 
                             resolve.defaults(list(formula=z~1, locations=~x+y, model=rfModel, ...), 
                                              list(dummy=Dummy, beta=Beta,  nmax=Nmax)))
  
  predxy <- predict(g.dummy, newdata=xyGrid, nsim=1)
  # Create a pixel image from the data.frame storing the values of the Gaussian Random Field
  grfImage <- as.im(predxy)
  # Create a linim from the grfImage
  Mpsp <- as.mask.psp(as.psp(X), W=as.owin(grfImage))
  grfLinImage <- linim(X, grfImage[Mpsp, drop=FALSE])
  
  return(grfLinImage)
}

#######################################################################################################
covariate_lpp <- function(expr.object=NULL,nsim=1,env=parent.frame()){
  stopifnot(!is.null(expr.object))
  stopifnot(inherits(expr.object,"expression"))
  
  linImages <- list()
  for(i in 1:nsim){
    linImages[[i]] <- eval(expr.object, env)
  }
  
  imnames <- sapply(seq(1,nsim), function(num){paste0("v", num)})
  names(linImages) <-  imnames
  return(linImages)
}
########################################################################################################
# Compute the proportion of LASSO solutions with nonzero regression coefficients for three regularization parameters
fraction_nonzero_coefs <- function(nsimulate = 1000, net, nvar=10, regcoeff=coefs0, 
                                   lambtimes=1, type_measure="deviance", ND=1000){
  
  frac_min <- frac_1se <- frac_avg <- matrix(0, nrow=nsimulate, ncol=nvar)
  npoints_avg <- rep(0, nsimulate)
  
  for(i in 1:nsimulate){
    lasso_fit <- simulate_lasso(net=net, nvar=nvar, regcoeff=regcoeff, lambtimes=lambtimes, 
                                type_measure=type_measure, ND=ND)
    
    coefs <- get_coef_indicator(lasso_fit$glm_dev)
    frac_min[i,] <- coefs$indi_min
    frac_1se[i,] <- coefs$indi_1se
    frac_avg[i,] <- coefs$indi_avg
    npoints_avg[i] <- lasso_fit$npts
    print(i)
  }
  
  retval <- list()
  
  retval$frac_min <- colMeans(frac_min)
  retval$frac_1se <- colMeans(frac_1se)
  retval$frac_avg <- colMeans(frac_avg)
  retval$npts_avg <- round(mean(npoints_avg))
  
  return(retval)
}


simulate_lasso <- function(net, nvar, regcoeff, lambtimes, type_measure, ND){
  newenv1 <- new.env()
  callObj <- expression(covariate.as.linim(X=net, Psill=0.25, Model='Exp', Range=100, Dummy=TRUE, Beta=1, Nmax=20))
  covLinImages <- covariate_lpp(callObj, nsim=nvar)
  
  newenv1$covLinImages_trans <- lapply(covLinImages, FUN = function(z){z/10})
  newenv1$beta0 <- regcoeff
  newenv1$lambtimes <- lambtimes
  
  expr <-  expression(lambtimes*0.00005 *exp(beta0[1] * v1 + beta0[2] * v2 + beta0[3] * v3 + beta0[4] * v4 + beta0[5] * v5))
  
  Lamb <- eval(expr, newenv1$covLinImages_trans, newenv1)
  
  net.lpp <- rpoislpp(lambda=Lamb, drop=T)
  
  Q <- linequad(net.lpp, nd=ND)
  X <- Q$data
  P <- union.quad(Q)
  form1 <- as.formula(paste("~", paste0("v", 1:nvar, collapse="+")))
  
  prep <- mpl.engine(Q, X, P, trend = form1, interaction=Poisson(), covariates = newenv1$covLinImages_trans, preponly = TRUE)
  
  glmData <- prep$glmdata
  y_response <- glmData$.mpl.Y
  wts <- glmData$.mpl.W
  x_cov <- data.matrix(glmData[3L:(nvar + 2L)])
  
  ## Deviance loss for Cross-Validation
  glm_dev <- cv.glmnet(x_cov, y_response, family = "poisson", weights=wts, type.measure = type_measure, alpha=1)
  
  retval <- list()
  retval$npts <- npoints(net.lpp)
  retval$glm_dev <- glm_dev
  
  
  return(retval)
}


get_coef_indicator <- function(glm_fit){
  stopifnot(inherits(glm_fit, "cv.glmnet"))
  
  nvar <- glm_fit$glmnet.fit$dim[1]
  lam_min_coefs <- lam_min_indi <- coef(glm_fit, s="lambda.min")[2L : (nvar+1), ]
  lam_min_indi[lam_min_indi != 0] <- 1
  
  lam_1se_coefs <- lam_1se_indi <- coef(glm_fit, s="lambda.1se")[2L : (nvar+1), ]
  lam_1se_indi[lam_1se_indi != 0] <- 1
  
  lam_avg <- 0.5 *(glm_fit$lambda.min + glm_fit$lambda.1se)
  lam_avg_coefs <- lam_avg_indi <- coef(glm_fit, s=lam_avg)[2L : (nvar+1), ]
  lam_avg_indi[lam_avg_indi != 0] <- 1
  
  retval <- list()
  
  retval$coef_min <- lam_min_coefs
  retval$indi_min <- lam_min_indi
  retval$coef_1se <- lam_1se_coefs
  retval$indi_1se <- lam_1se_indi
  retval$coef_avg <- lam_avg_coefs
  retval$indi_avg <- lam_avg_indi
  
  return(retval)
}

############################################################################################
get_counts_on_segments <- function(X, getIndicator=TRUE, getLengths=TRUE){
  stopifnot(inherits(X, "lpp"))
  stopifnot(isNamespaceLoaded("spatstat"))
  stopifnot(isNamespaceLoaded("dplyr"))
  
  nseg <- nsegments(X)
  seg_list <- data.frame(seg=1:nseg)
  events_df <- data.frame(seg=coords(X)$seg, freq=rep(1,npoints(X)))
  
  seg_with_counts <- summarize(group_by(events_df, seg), ncount=sum(freq))
  seg_with_counts <- left_join(seg_list, seg_with_counts, by="seg")
  seg_with_counts$ncount[is.na(seg_with_counts$ncount)] <- 0
  
  if(getIndicator){
    indi <- as.numeric((seg_with_counts$ncount > 0))
    seg_with_counts <- data.frame(seg_with_counts, Indicator=indi)
  }
  if(getLengths){
    segLengths <- lengths.psp(as.psp(X))
    seg_with_counts <- data.frame(seg_with_counts, Length=segLengths)
  }
  
  return(seg_with_counts)
}
###########################################################################################
compute_bias_se_lpp_models <- function(net, nsimulate=100, nvar=5, regcoeff=seq(10,6), 
                                       lambtimes=1L, ND=c(1000,2500,5000), EPS){
  
  ## All the checks
  stopifnot(inherits(net, "linnet"))
  stopifnot((length(nsimulate) == 1) && (nsimulate > 4))
  stopifnot((length(nvar) == 1) && (nvar > 0))
  stopifnot(is.vector(regcoeff) && (length(regcoeff) == nvar))
  stopifnot((length(lambtimes) == 1) && (lambtimes > 0))
  stopifnot(is.vector(ND))
  stopifnot((!missing(EPS)) && (is.vector(EPS)))
  
  n_d_model <- 3L # Number of discretised models (poisson, cloglog, and logit) considered. 
  # There are 3 * length(eps) discretised models and length(ND) lppm models
  neps <- length(EPS) # For each EPS value, all 3 discretised models will be fitted
  
  n_models <- ((neps + 1) * n_d_model) + length(ND) # Total number of models
  
  coef_estimates <- vector("list", length=nsimulate)
  npts <- rep(0, length=nsimulate)
    
  for(i in 1:nsimulate){
    newenv1 <- new.env()
    callObj <- expression(covariate.as.linim(X=net, Psill=0.25, Model='Exp', Range=100, Dummy=TRUE, Beta=1, Nmax=20))
    covLinImages <- covariate_lpp(callObj, nsim=nvar)
  
    newenv1$covLinImages_trans <- lapply(covLinImages, FUN = function(z){z/10})
    newenv1$beta0 <- regcoeff
    newenv1$lambtimes <- lambtimes
  
    expr <-  expression(lambtimes*0.00005 *exp(beta0[1] * v1 + beta0[2] * v2 + beta0[3] * v3 + beta0[4] * v4 + beta0[5] * v5))
    Lamb <- eval(expr, newenv1$covLinImages_trans, newenv1)
    net.lpp <- rpoislpp(lambda=Lamb, drop=T)
    npts[i] <- npoints(net.lpp)
    
    all_model_coefs <- matrix(0, nrow=nvar, ncol=n_models)
    
    ## Coefficient estimates of all 3 discretised models using the original network, 
    ## rendering number of lixels equal to number of segments.
    all_model_coefs[, 1:n_d_model] <- as.matrix(get_coef_of_discretised_models(net.lpp, newenv1$covLinImages_trans))
    
    # Fit all the discretised models using lixellated network.
    for(j in seq_along(EPS)){
      net.lpp.eps <- lixellate(net.lpp, eps=EPS[j])
      model_coefs_eps <- get_coef_of_discretised_models(net.lpp.eps, newenv1$covLinImages_trans)
      all_model_coefs[ ,(n_d_model + (j-1)*n_d_model + 1) : (n_d_model + j*n_d_model)] <- as.matrix(model_coefs_eps)
    }
    # 3+neps
    ## Coefficient estimates of the lppm model
    model_coefs2 <- get_coef_of_lppm(net.lpp, newenv1$covLinImages_trans, ND)
    all_model_coefs[ ,(n_d_model*(neps+1)+1):(n_models)] <- as.matrix(model_coefs2)
    
    coef_estimates[[i]] <- all_model_coefs
  }
  
  reg_coeff <- matrix(rep(regcoeff,n_models), nrow=nvar, ncol=n_models)
  # Estimate the bias
  coef_mean <- Reduce('+', coef_estimates)/nsimulate
  
  coef_bias <- lapply(coef_estimates, function(z, z0){z-z0}, z0=reg_coeff)
  
  coef_bias <- Reduce('+', coef_bias)/nsimulate
  # Estimate se
  coef_centerd_sq <- lapply(coef_estimates, function(z, z0){(z-z0)^2}, z0=coef_mean)
  coef_se <- sqrt(Reduce('+', coef_centerd_sq)/nsimulate)
  npts <- round(mean(npts))
  return(list(bias=coef_bias, se=coef_se, npts=npts))
} 

#################################################################################
get_coef_of_lppm <- function(net_lpp, covar, ND){
  ## IPP model fit
  nND <- length(ND)
  nvar <- length(covar)
  lppm_coef <- matrix(0, ncol=nND, nrow=nvar)
  for(j in seq_along(ND)){
    lppm_fit <- lppm(net_lpp ~ ., data=covar, nd=ND[j])
    lppm_coef[,j] <- coef(lppm_fit)[-1]
  }
 return(as.data.frame(lppm_coef))
}

get_coef_of_discretised_models <- function(X, cov_images){
  stopifnot(inherits(X, "lpp"))
  
  net <- domain(X)
  mid_pts <- get_midpt_of_segments(net)
  cov_values <- lapply(cov_images, FUN = function(z){z[mid_pts]})
  
  cov_mat <- do.call("cbind", cov_values)
  counts_data <- get_counts_on_segments(X)
  counts_vec <- counts_data[["ncount"]]
  length_vec <- counts_data[["Length"]]
  indi_vec <- counts_data[["Indicator"]]
  
  cov_frame <- as.data.frame(cov_mat)
  log_length <- log(length_vec)
  ## Fit Poisson model to count data
  pois_data <- data.frame(cov_frame, y = counts_vec)
  pois_fit <- glm(y ~ ., data=pois_data, family="poisson", offset=log_length)
  pois_coef <- coef(pois_fit)[-1]
  ## Fit Complementary log-log model
  model_data <- data.frame(cov_frame, y=as.factor(indi_vec))
  cloglog_fit <- glm(y ~., data=model_data, family=binomial(link = "cloglog"), offset=log_length)
  cloglog_coef <- coef(cloglog_fit)[-1]
  ## Fit a logistic-regression model
  logit_fit <- glm(y ~., data=model_data, family=binomial(link = "logit"), offset=log_length)
  logit_coef <- coef(logit_fit)[-1]
  
  return(data.frame(pois_coef = pois_coef, cloglog_coef = cloglog_coef, logit_coef = logit_coef))
}

####################################################################################
get_midpt_of_segments <- function(X){
  
  stopifnot(inherits(X, "linnet"))
  net_ends <- as.psp(X)$ends
  mid_pts_x <- 0.5 * (net_ends$x0 + net_ends$x1)
  mid_pts_y <- 0.5 * (net_ends$y0 + net_ends$y1)
  W <- X$window
  mid_pt_ppp <- ppp(mid_pts_x, mid_pts_y, window=W)
  return(mid_pt_ppp)
}

####################################################################################
arrange_bias_se <- function(bias_vals, se_vals){
  stopifnot(is.matrix(bias_vals) && is.matrix(se_vals))
  stopifnot(dim(bias_vals) == dim(se_vals))
  
  nRow <- nrow(bias_vals)
  nCol <- ncol(bias_vals)
  
  bias_vals <- round(bias_vals, digits=3)
  se_vals <- round(se_vals, digits=3)
  new_table <- matrix(0, nrow=nRow, ncol=nCol)
  for(j in 1:nCol){
    new_table[,j] <- paste0(bias_vals[,j], " ", paste0("(",se_vals[,j],")"))
  }
  return(new_table)
}
#########################################################################################
# Get unified color breaks for a list of linim objects
#########################################################################################
get_unified_ribbon_breaks <- function(linimList = covLinImages_trans[1:5]){
  stopifnot(!missing(linimList))
  n <- length(linimList)
  linplt <- plot(linimList[[1]], main="")
  nbrk <- length(attr(linplt, "stuff")$breaks)
  
  brks <- vector(mode="numeric", length = n*nbrk)
  
  for(i in 1:n){
    plt <- plot(linimList[[i]], main="")
    brk <- attr(plt, "stuff")$breaks
    sid <- (i-1)*nbrk + 1
    eid <- i*nbrk
    brks[sid:eid] <- brk
  }
  
  brksCut <- levels(cut(brks, breaks = (nbrk-1)))
  brk_res <- vector(mode="numeric", length = nbrk)
  
  for(i in seq_along(brksCut)){
    brk <- brksCut[i]
    brkParts <- str_split(brk, ",")
    brk2 <- brkParts[[1]][2]
    if(i == 1){
      brk1 <- brkParts[[1]][1]
      brk_res[i] <- as.numeric(str_replace(brk1, "\\(", ""))
    }
    brk_res[i+1] <- as.numeric(str_replace(brk2, "\\]", ""))
  }
  return(brk_res)
}

#################################################################################
compute_lambda <- function(env=parent.frame()){
  beta0 <<- seq(from=10, to=5, length=5)
  expr <-  expression((2.5)*0.00005 *exp(beta0[1] * v1 + beta0[2] * v2 + beta0[3] * v3 + beta0[4] * v4 + beta0[5] * v5))
  
  Lambda <- eval(expr,env$covLinImages_trans,env)
  return(Lambda)
}

#######################################################################################
simulate_lasso2 <- function(inten, ND){
  net.lpp <- rpoislpp(lambda=inten, drop=T)
  Q <- linequad(net.lpp, nd=ND)
  X <- Q$data
  P <- union.quad(Q)
  form1 <- as.formula(paste("~", paste0("v", 1:nvar, collapse="+")))
  
  prep <- mpl.engine(Q, X, P, trend = form1, interaction=Poisson(), covariates = newenv1$covLinImages_trans, preponly = TRUE)
  
  glmData <- prep$glmdata
  y_response <- glmData$.mpl.Y
  wts <- glmData$.mpl.W
  x_cov <- data.matrix(glmData[3L:(nvar + 2L)])
  
  ## Deviance loss for Cross-Validation
  glm_dev <- cv.glmnet(x_cov, y_response, family = "poisson", weights=wts, type.measure = type_measure, alpha=1)
  
  retval <- list()
  retval$npts <- npoints(net.lpp)
  retval$glm_dev <- glm_dev
  
  
  return(retval)
}
########################################################################################################
compute_bias_se_lpp_models2 <- function(Lamb, 
                                        nsimulate=100,
                                        lambtimes=1L, 
                                        nvar=5,
                                        ND=c(1000,2500,5000), 
                                        EPS,
                                        n_d_model=1L){
  ## All the checks
  stopifnot((length(nsimulate) == 1) && (nsimulate > 4))
  stopifnot((length(nvar) == 1) && (nvar > 0))
  # stopifnot(is.vector(regcoeff) && (length(regcoeff) == nvar))
  # stopifnot((length(lambtimes) == 1) && (lambtimes > 0))
  stopifnot(is.vector(ND))
  stopifnot((!missing(EPS)) && (is.vector(EPS)))
  
  # n_d_model <- 3L # Number of discretised models (poisson, cloglog, and logit) considered. 
  # There are 3 * length(eps) discretised models and length(ND) lppm models
  neps <- length(EPS) # For each EPS value, all 3 discretised models will be fitted
  
  n_models <- ((neps + 1) * n_d_model) + length(ND) # Total number of models
  
  coef_estimates <- vector("list", length=nsimulate)
  npts <- rep(0, length=nsimulate)
  
  for(i in 1:nsimulate){
    newenv1 <- new.env()
    #callObj <- expression(covariate.as.linim(X=net, Psill=0.25, Model='Exp', Range=100, Dummy=TRUE, Beta=1, Nmax=20))
    #covLinImages <- covariate_lpp(callObj, nsim=nvar)
    
    #newenv1$covLinImages_trans <- lapply(covLinImages, FUN = function(z){z/10})
    #newenv1$beta0 <- regcoeff
    #newenv1$lambtimes <- lambtimes
    
    #expr <-  expression(lambtimes*0.00005 *exp(beta0[1] * v1 + beta0[2] * v2 + beta0[3] * v3 + beta0[4] * v4 + beta0[5] * v5))
    #Lamb <- eval(expr, newenv1$covLinImages_trans, newenv1)
    net.lpp <- rpoislpp(lambda=Lamb, drop=T)
    npts[i] <- npoints(net.lpp)
    
    all_model_coefs <- matrix(0, nrow=nvar, ncol=n_models)
    
    ## Coefficient estimates of all 3 discretised models using the original network, 
    ## rendering number of lixels equal to number of segments.
    all_model_coefs[, 1:n_d_model] <- as.matrix(get_coef_of_discretised_models(net.lpp, newenv1$covLinImages_trans))
    
    # Fit all the discretised models using lixellated network.
    for(j in seq_along(EPS)){
      net.lpp.eps <- lixellate(net.lpp, eps=EPS[j])
      model_coefs_eps <- get_coef_of_discretised_models(net.lpp.eps, newenv1$covLinImages_trans)
      all_model_coefs[ ,(n_d_model + (j-1)*n_d_model + 1) : (n_d_model + j*n_d_model)] <- as.matrix(model_coefs_eps)
    }
    # 3+neps
    ## Coefficient estimates of the lppm model
    model_coefs2 <- get_coef_of_lppm(net.lpp, newenv1$covLinImages_trans, ND)
    all_model_coefs[ ,(n_d_model*(neps+1)+1):(n_models)] <- as.matrix(model_coefs2)
    
    coef_estimates[[i]] <- all_model_coefs
  }
  
  reg_coeff <- matrix(rep(regcoeff,n_models), nrow=nvar, ncol=n_models)
  # Estimate the bias
  coef_mean <- Reduce('+', coef_estimates)/nsimulate
  
  coef_bias <- lapply(coef_estimates, function(z, z0){z-z0}, z0=reg_coeff)
  
  coef_bias <- Reduce('+', coef_bias)/nsimulate
  # Estimate se
  coef_centerd_sq <- lapply(coef_estimates, function(z, z0){(z-z0)^2}, z0=coef_mean)
  coef_se <- sqrt(Reduce('+', coef_centerd_sq)/nsimulate)
  npts <- round(mean(npts))
  return(list(bias=coef_bias, se=coef_se, npts=npts))
} 
##############################################################################################
get_bias_se_utils <- function(data){
  data_mut1 <-   data %>% 
    mutate(bias1=b1-10,bias2=b2-9,bias3=b3-8,bias4=b4-7,bias5=b5-6,
           bias_sq1=bias1^2,bias_sq2=bias2^2,bias_sq3=bias3^2,
           bias_sq4=bias4^2,bias_sq5=bias5^2) %>%
    group_by(Method, lixel, Model) %>%
    summarise(avg_bias1=mean(bias1),avg_bias2=mean(bias2),avg_bias3=mean(bias3),
              avg_bias4 = mean(bias4),avg_bias5 = mean(bias5),se1=sqrt(mean(bias_sq1)),
              se2=sqrt(mean(bias_sq2)),se3=sqrt(mean(bias_sq3)),se4=sqrt(mean(bias_sq4)),
              se5=sqrt(mean(bias_sq5))) %>% rowwise() %>% 
    mutate(b1 = paste0(formatC(avg_bias1, digits=3, format="f"), " (", formatC(se1, digits=3, format="f"), ")"),
           b2 = paste0(formatC(avg_bias2, digits=3, format="f"), " (", formatC(se2, digits=3, format="f"), ")"),
           b3 = paste0(formatC(avg_bias3, digits=3, format="f"), " (", formatC(se3, digits=3, format="f"), ")"),
           b4 = paste0(formatC(avg_bias4, digits=3, format="f"), " (", formatC(se4, digits=3, format="f"), ")"),
           b5 = paste0(formatC(avg_bias5, digits=3, format="f"), " (", formatC(se5, digits=3, format="f"), ")")) %>%
    select(Method, lixel, Model, b1:b5)
  
  
  data_mut2 <-   data %>% 
    mutate(bias6=b6,bias7=b7,bias8=b8,bias9=b9,bias10=b10,
           bias_sq6=bias6^2,bias_sq7=bias7^2,bias_sq8=bias8^2,
           bias_sq9=bias9^2,bias_sq10=bias10^2) %>%
    group_by(Method, lixel, Model) %>%
    summarise(avg_bias6=mean(bias6),avg_bias7=mean(bias7),avg_bias8=mean(bias8),
              avg_bias9 = mean(bias9),avg_bias10 = mean(bias10),se6=sqrt(mean(bias_sq6)),
              se7=sqrt(mean(bias_sq7)),se8=sqrt(mean(bias_sq8)),se9=sqrt(mean(bias_sq9)),
              se10=sqrt(mean(bias_sq10))) %>% rowwise() %>% 
    mutate(b6 = paste0(formatC(avg_bias6, digits=3, format="f"), " (", formatC(se6, digits=3, format="f"), ")"),
           b7 = paste0(formatC(avg_bias7, digits=3, format="f"), " (", formatC(se7, digits=3, format="f"), ")"),
           b8 = paste0(formatC(avg_bias8, digits=3, format="f"), " (", formatC(se8, digits=3, format="f"), ")"),
           b9 = paste0(formatC(avg_bias9, digits=3, format="f"), " (", formatC(se9, digits=3, format="f"), ")"),
           b10 = paste0(formatC(avg_bias10, digits=3, format="f"), " (", formatC(se10, digits=3, format="f"), ")")) %>%
    select(Method, lixel, Model, b6:b10)
  
  retlist <- list(bias_se_b1_to_b5 = data_mut1, bias_se_b6_to_b10 = data_mut2)
  return(retlist)
}
###########################################################################################################
get_coef_inclusion_prop <- function(reg_coef_all){
  reg_coef_all1 <- reg_coef_all %>% rowwise() %>% 
    mutate(Model = paste0(Model, ifelse(is.na(lixel),"",lixel))) %>%
    select(-one_of(c("lixel")))
  reg_coef_all2 <- reg_coef_all1 %>% mutate(b1 = ifelse(abs(b1) <= 0.0001, 0L, 1L),
                                            b2 = ifelse(abs(b2) <= 0.0001, 0L, 1L),
                                            b3 = ifelse(abs(b3) <= 0.0001, 0L, 1L),
                                            b4 = ifelse(abs(b4) <= 0.0001, 0L, 1L),
                                            b5 = ifelse(abs(b5) <= 0.0001, 0L, 1L),
                                            b6 = ifelse(abs(b6) <= 0.0001, 0L, 1L),
                                            b7 = ifelse(abs(b7) <= 0.0001, 0L, 1L),
                                            b8 = ifelse(abs(b8) <= 0.0001, 0L, 1L),
                                            b9 = ifelse(abs(b9) <= 0.0001, 0L, 1L),
                                            b10 = ifelse(abs(b10) <= 0.0001, 0L, 1L))
  reg_coef_all3 <- ungroup(reg_coef_all2) %>% group_by(Model, TuningCoef) %>%
    summarize(b1 = mean(b1), b2 = mean(b2), b3 = mean(b3), b4 = mean(b4),
              b5 = mean(b5), b6 = mean(b6), b7 = mean(b7), b8 = mean(b8),
              b9 = mean(b9), b10 = mean(b10))
  
  reg_coef_all3$Model <- factor(reg_coef_all3$Model, levels = c("logit1810",
                                                                "poisson1810",
                                                                "logit3370",
                                                                "poisson3370",
                                                                "logit6481",
                                                                "poisson6481",
                                                                "lppm1000",
                                                                "lppm2500",
                                                                "lppm5000"))
  reg_coef_all3$TuningCoef <- factor(reg_coef_all3$TuningCoef,
                                     levels = c("Min", "Avg", "1se"))
  reg_coef_prop0 <- arrange(reg_coef_all3, Model, TuningCoef)
  return(reg_coef_prop0)
}
############################################################################################################
match_spdlim <- function(road.no, strt.slk, end.slk, matching.data, default.spd=110){
  ret_spd <- default.spd
  road.data <- subset(matching.data, Road_No==road.no)
  n <- nrow(road.data)
  if(n == 0){
    return(ret_spd)
  }
  road.data <- road.data[!duplicated(road.data),c(2,3,4)]
  road.data <- road.data[order(road.data$START_SLK, road.data$END_SLK), ]
  START_SLK <- road.data$START_SLK
  END_SLK <-  road.data$END_SLK
  SPLI_SPEED <- road.data$SPLI_SPEED
  
    min.strt.slk <- min(START_SLK)
    if(strt.slk < min.strt.slk){
      ret_spd <- SPLI_SPEED[1]
    }else{
      max.end.slk <- max(END_SLK)
      if(end.slk > max.end.slk){
        ret_spd <- tail(SPLI_SPEED, n=1)
      }else{
        cond <- (strt.slk <= START_SLK)
        if(any(cond)){
          which.id <- which(cond)[1]
          if(end.slk < START_SLK[which.id]){
            ret_spd <- SPLI_SPEED[which.id -1]
          }else{
            ret_spd <- SPLI_SPEED[which.id]
          }
        }else{
          ret_spd <- tail(SPLI_SPEED, n=1)
        }
      }
    }
    return(ret_spd)
}
###########################################################################################################
scale_road_covariates <- function(covar.data){
  # Dummy coding of the factors
  all_road_vars <- c("SPLI_SPEED","HOAL_CURVE","TOTAL_PAVE","TOTAL_SEAL","TRAFFICABL", 
                     "SHOULDER_S","KERB_L","KERB_R","NO_OF_LANE","FLOODWAY","BRIDGE")
  covars_scaled <- covar.data[all_road_vars]
  covars_scaled$SHOULDER_S <- factor(as.integer(covar.data$SHOULDER_S) -1)
  covars_scaled$KERB_L <- factor(as.integer(covar.data$KERB_L) -1)
  covars_scaled$KERB_R <- factor(as.integer(covar.data$KERB_R) -1)
  covars_scaled$FLOODWAY <- factor(as.integer(covar.data$FLOODWAY) -1)
  covars_scaled$BRIDGE <- factor(as.integer(covar.data$BRIDGE) -1)
  # Scaling (normalization) of the quantitative variables
  covars_scaled$SPLI_SPEED <- (covar.data$SPLI_SPEED - min(covar.data$SPLI_SPEED))/(max(covar.data$SPLI_SPEED) - min(covar.data$SPLI_SPEED))
  covars_scaled$HOAL_CURVE <- (covar.data$HOAL_CURVE - min(covar.data$HOAL_CURVE))/(max(covar.data$HOAL_CURVE) - min(covar.data$HOAL_CURVE))
  covars_scaled$TOTAL_PAVE <- (covar.data$TOTAL_PAVE- min(covar.data$TOTAL_PAVE))/(max(covar.data$TOTAL_PAVE) - min(covar.data$TOTAL_PAVE))
  covars_scaled$TOTAL_SEAL <- (covar.data$TOTAL_SEAL- min(covar.data$TOTAL_SEAL))/(max(covar.data$TOTAL_SEAL) - min(covar.data$TOTAL_SEAL))
  covars_scaled$TRAFFICABL <- (covar.data$TRAFFICABL- min(covar.data$TRAFFICABL))/(max(covar.data$TRAFFICABL) - min(covar.data$TRAFFICABL))
  covars_scaled$NO_OF_LANE <- (covar.data$NO_OF_LANE- min(covar.data$NO_OF_LANE))/(max(covar.data$NO_OF_LANE) - min(covar.data$NO_OF_LANE))
  return(covars_scaled)
}
###############################################################################################################
get_derived_covars <- function(covars_wa_massive_df){
  covars_wa_massive_df$SPLI_SPEED2 <- (covars_wa_massive_df$SPLI_SPEED)^2
  covars_wa_massive_df$HOAL_CURVE2 <- (covars_wa_massive_df$HOAL_CURVE)^2
  covars_wa_massive_df$TOTAL_PAVE2 <- (covars_wa_massive_df$TOTAL_PAVE)^2
  covars_wa_massive_df$TOTAL_SEAL2 <- (covars_wa_massive_df$TOTAL_SEAL)^2
  covars_wa_massive_df$TRAFFICABL2 <- (covars_wa_massive_df$TRAFFICABL)^2
  covars_wa_massive_df$NO_OF_LANE2 <- (covars_wa_massive_df$NO_OF_LANE)^2
  
  covars_wa_massive_df$SHOULDER_S <- as.integer(covars_wa_massive_df$SHOULDER_S)-1
  covars_wa_massive_df$KERB_L <- as.integer(covars_wa_massive_df$KERB_L)-1
  covars_wa_massive_df$KERB_R <- as.integer(covars_wa_massive_df$KERB_R)-1
  covars_wa_massive_df$FLOODWAY <- as.integer(covars_wa_massive_df$FLOODWAY)-1
  covars_wa_massive_df$BRIDGE <- as.integer(covars_wa_massive_df$BRIDGE)-1
  
  # Interaction between KERB_L and others
  covars_wa_massive_df$SPLI_SPEEDxKERB_L <- (covars_wa_massive_df$SPLI_SPEED)*(covars_wa_massive_df$KERB_L)
  covars_wa_massive_df$HOAL_CURVExKERB_L <- (covars_wa_massive_df$HOAL_CURVE)*(covars_wa_massive_df$KERB_L)
  covars_wa_massive_df$TOTAL_PAVExKERB_L <- (covars_wa_massive_df$TOTAL_PAVE)*(covars_wa_massive_df$KERB_L)
  covars_wa_massive_df$TOTAL_SEALxKERB_L <- (covars_wa_massive_df$TOTAL_SEAL)*(covars_wa_massive_df$KERB_L)
  covars_wa_massive_df$TRAFFICABLxKERB_L <- (covars_wa_massive_df$TRAFFICABL)*(covars_wa_massive_df$KERB_L)
  covars_wa_massive_df$NO_OF_LANExKERB_L <- (covars_wa_massive_df$NO_OF_LANE)*(covars_wa_massive_df$KERB_L)
  
  covars_wa_massive_df$SHOULDER_SxKERB_L <- (covars_wa_massive_df$SHOULDER_S)*(covars_wa_massive_df$KERB_L)
  covars_wa_massive_df$KERB_RxKERB_L <- (covars_wa_massive_df$KERB_R)*(covars_wa_massive_df$KERB_L)
  covars_wa_massive_df$FLOODWAYxKERB_L <- (covars_wa_massive_df$FLOODWAY)*(covars_wa_massive_df$KERB_L)
  covars_wa_massive_df$BRIDGExKERB_L <- (covars_wa_massive_df$BRIDGE)*(covars_wa_massive_df$KERB_L)
  
  # Interaction between SHOULDER_S and others
  covars_wa_massive_df$SPLI_SPEEDxSHOULDER_S <- (covars_wa_massive_df$SPLI_SPEED)*(covars_wa_massive_df$SHOULDER_S)
  covars_wa_massive_df$HOAL_CURVExSHOULDER_S <- (covars_wa_massive_df$HOAL_CURVE)*(covars_wa_massive_df$SHOULDER_S)
  covars_wa_massive_df$TOTAL_PAVExSHOULDER_S <- (covars_wa_massive_df$TOTAL_PAVE)*(covars_wa_massive_df$SHOULDER_S)
  covars_wa_massive_df$TOTAL_SEALxSHOULDER_S <- (covars_wa_massive_df$TOTAL_SEAL)*(covars_wa_massive_df$SHOULDER_S)
  covars_wa_massive_df$TRAFFICABLxSHOULDER_S <- (covars_wa_massive_df$TRAFFICABL)*(covars_wa_massive_df$SHOULDER_S)
  covars_wa_massive_df$NO_OF_LANExSHOULDER_S <- (covars_wa_massive_df$NO_OF_LANE)*(covars_wa_massive_df$SHOULDER_S)
  
  covars_wa_massive_df$KERB_RxSHOULDER_S <- (covars_wa_massive_df$KERB_R)*(covars_wa_massive_df$SHOULDER_S)
  covars_wa_massive_df$FLOODWAYxSHOULDER_S <- (covars_wa_massive_df$FLOODWAY)*(covars_wa_massive_df$SHOULDER_S)
  covars_wa_massive_df$BRIDGExSHOULDER_S <- (covars_wa_massive_df$BRIDGE)*(covars_wa_massive_df$SHOULDER_S)
  
  # Interaction between KERB_R and others
  covars_wa_massive_df$SPLI_SPEEDxKERB_R <- (covars_wa_massive_df$SPLI_SPEED)*(covars_wa_massive_df$KERB_R)
  covars_wa_massive_df$HOAL_CURVExKERB_R <- (covars_wa_massive_df$HOAL_CURVE)*(covars_wa_massive_df$KERB_R)
  covars_wa_massive_df$TOTAL_PAVExKERB_R <- (covars_wa_massive_df$TOTAL_PAVE)*(covars_wa_massive_df$KERB_R)
  covars_wa_massive_df$TOTAL_SEALxKERB_R <- (covars_wa_massive_df$TOTAL_SEAL)*(covars_wa_massive_df$KERB_R)
  covars_wa_massive_df$TRAFFICABLxKERB_R <- (covars_wa_massive_df$TRAFFICABL)*(covars_wa_massive_df$KERB_R)
  covars_wa_massive_df$NO_OF_LANExKERB_R <- (covars_wa_massive_df$NO_OF_LANE)*(covars_wa_massive_df$KERB_R)
  
  
  covars_wa_massive_df$FLOODWAYxKERB_R <- (covars_wa_massive_df$FLOODWAY)*(covars_wa_massive_df$KERB_R)
  covars_wa_massive_df$BRIDGExKERB_R <- (covars_wa_massive_df$BRIDGE)*(covars_wa_massive_df$KERB_R)
  
  return(covars_wa_massive_df)
}
###############################################################################################################
get_derived_covars_for_marked_pattern <- function(X){
  X0 <- X[-1] # without marks column
  X0 <- get_derived_covars(X0)
  names0 <- names(X0)
  names1 <- paste0("MARKSx", names0)
  marks <- as.integer(factor(X$marks, levels = c("Low", "High"), labels = c(0, 1))) -1
  #########################################################################################
  X2 <- as.data.frame(apply(X0, 2, function(z){z*marks}))
  names(X2) <- names1
  X0$MARKS <- marks
  X_return <- cbind(X0, X2)
  return(X_return)
}
###################################################








