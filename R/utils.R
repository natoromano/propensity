###############################################################################
# (Script)
# Utility functions

# Author: Yen Low
###############################################################################

require(xlsx)


# INSTALL NEW PACKAGE
installnewpackage <- function(reqpackages) {
  newpackages=reqpackages[!(reqpackages %in% installed.packages()[, "Package"])]
  if(length(newpackages)>0) {
    install.packages(newpackages)
  }
}


# LASSO HELPERS
coeffAtlambda <- function(cvfit,lambda="lambda.1se") {
  coeff_interpretation <- coef(cvfit, s=lambda)
  coeff_interpretmat <- as.matrix(coeff_interpretation)
  coeff_interpretmat[-(coeff_interpretation@i+1),] <- NA # replace . with NA
  coeff <- as.vector(coeff_interpretmat)
  names(coeff) <- rownames(coeff_interpretmat)
  return(coeff)
}


genCI <- function(xmat, yvec, ntrials=100, family="binomial", lambda,
                  alpha=1, penalty.factor=1){
  # loop through 1 to ntrials
  modnet = foreach(j=1:ntrials,.packages="glmnet") %dopar% {
    if(is.vector(yvec)){  # for most models
      mod_ID <- sample(1:length(yvec),replace=T)
      yvec_subset <- yvec[mod_ID]
    } else if (ncol(yvec)==2) {  # for cox model
      mod_ID <- sample(1:nrow(yvec),replace=T)
      yvec_subset <- yvec[mod_ID,]
    } else stop("yvec must be a vector OR matrix/Surv obj of 2 columns")
    
    if (!is.null(lambda)) {
      glmnet(xmat[mod_ID,],yvec_subset,thresh=1e-5,
             alpha=alpha,penalty.factor=penalty.factor,lambda=lambda,
             family=family,standardize=F)
    } else {
      cv.glmnet(xmat[mod_ID,], yvec_subset, thresh=1e-5, nfold=3,
                alpha=alpha, penalty.factor=penalty.factor, nlambda=20,
                family=family, standardize=F, parallel=T)
    }
  } # end of outer j foreach loop
  
  if(!is.null(lambda)) {
    temp <- sapply(modnet, coef)
  } else {
    temp <- sapply(modnet, function(x) coef(x, s=x$lambda.1se))      
  }
  coeffarray <- temp[[1]]
  for(p in 2:length(temp)) coeffarray <- cBind(coeffarray, temp[[p]])
  
  list(coeffarray=coeffarray, type="glmnet")
}


setCL <- function(CIobj, CL=c(0.025,0.5,0.975), maxNA=0.5, verbose=FALSE) {
  betaCIlim <- t(apply(CIobj$coeffarray, 1, quantile, probs=CL, na.rm=T,
                       names=F))
  colnames(betaCIlim) <- c("lowlim", "median", "upplim")
  beta_nonZero <- !apply(betaCIlim, 1, function(x) x[1]==0 & x[3]==0)
  fracNA <- (rowSums(is.na(CIobj$coeffarray))/ncol(CIobj$coeffarray)) < maxNA
  beta_nonZero <- ((beta_nonZero+fracNA)==2)
  return(list(betaCIlim=betaCIlim,beta_nonZero=beta_nonZero))
}
