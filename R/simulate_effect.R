# Generates treatments and outcome on OSIM2 data.
# Parameters:
#    covariate: covariate matrix
#    seed: seed number
#    conf.case: 1, 2, or 3 for confounding strength
#    TE.type: 0, 1, 2 for zero effect, homogeneous, non-linear heterogeneous
simulateTreatmentEffect = function (covariate, seed, conf.case, TE.type) {

  library(logitnorm)
  colnames(covariate)[2] <- 'YOB'
  X <- as.matrix(covariate[ ,-c(1, 432)]) # remove PERSON_ID
  X[, 1] <- scale(X[, 1])   # scale (center and standardise) YOB
  rm(covariate)
  
  #----------------------------------------------------------------------
  # Input Parameters
  #----------------------------------------------------------------------
  N = as.numeric(nrow(X))   # total sample size
  Nvar = as.numeric(ncol(X)) # Number of covariates in total
  Nc = 5 # Number of confounders
  No = 5 # Number of outcome-only predictors
  Ne = 5 # Number of exposure-only predictors
  Ntau = 5 # Number of TE modifiers only
  Ncov = Nc + No + Ne + Ntau
  Nnoise = Nvar - Ncov

  ## Confounding magnitude
  # conf.case = 1   # low confounding
  # conf.case = 2   # moderate confounding
  # conf.case = 3   # high confounding
  
  ## Type of Treatment Effect
  # TE.type=0 # Zero effect
  # TE.type=1 # Marginal effect of 1
  # TE.type=2 # non-linear HTE

  # seed = 12345
  set.seed(seed)
  
  #----------------------------------------------------------------------
  #  Identify 19 most dense binary features
  #----------------------------------------------------------------------
  # tb <- apply(X[,2:Nvar], 2, function(j) {min(table(j))/N} )
  # prev <- cbind(idx = 2:Nvar, tb)
  # prev <- prev[order(prev[,'tb'], decreasing = T), ]
  # idx.dense <- prev[1:19,'idx']  
  # idx.cov <- c(1, idx.dense)   
  # Ncov = 20 (YOB + 19 most dense binary variables as covariate)
  idx.cov <- c(1, 693, 2, 106, 503, 559, 129, 474, 690, 205, 692, 286, 596, 
               665, 648, 631, 210, 658, 421, 607)
  Xc   = X[,idx.cov[1:5]]
  Xe   = X[,idx.cov[6:10]]
  Xo   = X[,idx.cov[11:15]]
  Xtau = X[,idx.cov[16:20]]

  #----------------------------------------------------------------------
  #  Simulate exposure - based on Austin (2015)
  #----------------------------------------------------------------------
  ## Define coefficients
  k = c(1, 2, 5)  # control confounding magnitude
  coeffs = c(1.05, 1.10, 1.20, 1.25, 1.50, 1.75, 2, 1.5, 1.25, 1.1)  
  alpha = c( log(k[conf.case]*coeffs[1:5]), log(k[2]*coeffs[6:10]) )
  # exp(alpha)
  
  ## Compute propensity scores
  e = te = rep(0, N)
  Xtemp = cbind(Xc, Xe)
  slopes = sapply(1:N, function(i) alpha %*% Xtemp[i,])
  # par(mfrow = c(1,2)); hist(slopes); hist(invlogit(slopes))
  slopes.mean = mean(slopes)
  a = c(1.79, 2.064, 2.8)  # Determined iteratively to have mean(e)=0.2
  alpha0 = -a[conf.case] - slopes.mean
  # mean(invlogit(slopes + alpha0))
  t.e = sapply(1:N, function(i) { alpha0 + alpha %*% Xtemp[i,] })
  e = invlogit(t.e)
  # par(mfrow = c(1,2)); hist(t.e); hist(e); mean(t.e); mean(e)
  
  ## Generate a Bernoulli distribution with probability ps
  Z = rbinom(N , 1, e) 
  # round(table(Z)/N*100, digits=1)
  covNames.Z <- colnames(Xtemp)

  #----------------------------------------------------------------------
  #  Simmulate treatment effect 
  #----------------------------------------------------------------------
  tau = t.tau = rep(0, N)
  covNames.tau <- NULL
  if (TE.type==1) { 
    beta = 3 * log(c(2, 1.75, 1.50, 1.25, 1.10, 1.05, 1.50, 1.75, 2, 1.25))
    Xtemp = cbind(Xc, Xtau)
    slopes = sapply(1:N, function(i) beta %*% Xtemp[i,])
    # par(mfrow = c(1,2)); hist(slopes); hist(invlogit(slopes))
    slopes.mean = mean(slopes)
    beta0 = 1.642 - slopes.mean 
    t.tau = sapply(1:N, function(i) { beta0 + beta %*% Xtemp[i,] })   
    min = -2; max = 2
    tau = min + (max - min) * invlogit(t.tau)
    # par(mfrow = c(1,2)); hist(t.tau); hist(tau); mean(tau)
    covNames.tau <- colnames(Xtemp)
  }
  
  #----------------------------------------------------------------------
  #  Simmulate outcome
  #----------------------------------------------------------------------
  ## Define coefficients
  gamma = c(2, 1.75, 1.50, 1.25, 1.10, 1.05, 1.50, 1.75, 2, 1.25)
  
  ## Compute probabilities
  pZ0 = pZ1 = t.pZ0 = t.po = po = TE =  rep(0, N)
  Xtemp = cbind(Xc, Xo)
  slopes = sapply(1:N, function(i) gamma %*% Xtemp[i,])
  # par(mfrow = c(1,2)); hist(slopes); hist(invlogit(slopes))
  slopes.mean = mean(slopes)
  gamma0 = -2.9 - slopes.mean
  # mean(invlogit(slopes + gamma0))
  t.pZ0 = sapply(1:N, function(i) { gamma0 + gamma %*% Xtemp[i,] })
  pZ0 = invlogit(t.pZ0)
  mean(pZ0)
  pZ1 = invlogit(t.pZ0 + tau)
  TE = pZ1 - pZ0   # TE on binary outcome is defined as risk difference
  # par(mfrow=c(2,2)); hist(pZ0); hist(pZ1); hist(TE); mean(TE)   
  # Compression effect (see Green and Kern, 2012)
  t.po = t.pZ0 + tau * Z
  po = invlogit(t.po)
  # par(mfrow = c(1,2)); hist(t.po); hist(po); mean(t.po); mean(po)
  
  ## Generate a Bernoulli distribution with probability ps
  Y = rbinom(N , 1, po) 
  # round(table(Y)/N*100, digits=1)
  
  covNames.Y <- colnames(Xtemp)
  
  #----------------------------------------------------------------------
  # SAVE DATA
  #----------------------------------------------------------------------
  simData=data.frame(X, Z, Y, e, po, tau, TE)
  # head(simData, 20)
  covNames <- list(Z = covNames.Z, Y = covNames.Y, tau = covNames.tau)
  return(list(simData = simData, covNames = covNames))
  
}
