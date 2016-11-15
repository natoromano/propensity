# This script uses the function in simulate_effect generate
# the dataset and split it into training and test sets.

rm(list=ls())
# library(Matching)
# options(java.parameters = "-Xmx8g")
# library(bartMachine)
# set_bart_machine_num_cores(1)
source("simulate_effect.R")
covariate <- read.csv('/labs/shahlab/OSIM2/covariate.csv', check.names =  T)


#-----------------------------------------------------------------------------
# Simulate and format data
#-----------------------------------------------------------------------------
## Confounding magnitude
# conf.case = 1   # low
 conf.case = 2   # moderate
# conf.case = 3   # high

## Type of Treatment Effect
# TE.type = 0 # Zero effect
TE.type = 1 # marginal effect of 1
# TE.type = 2 # non-linear HTE


## Simulate data 
seed = as.numeric(sample(1:1000, 1))
lSim <- simulateTreatmentEffect(covariate, seed, conf.case, TE.type)
simData <- lSim$simData
covNames <- lSim$covNames
# head(simData)
Nvar = 698
N = as.numeric(nrow(simData))
round(table(simData$Z)/N*100, digits=1)
round(table(simData$Y)/N*100, digits=1)
par(mfrow=c(1,2)); hist(simData$tau); hist(simData$TE)
# tau is the coefficient on the linear scale; TE is on probability scale


## Sample size
# s = 1   # original sample size
s = 20000/0.8/N


## Define sample size and randomly partition original records into training and test records
Nsamp = ceiling(N * s)
Ntrain <- ceiling(0.8 * Nsamp)     # number of training records
Ntest <- Nsamp - Ntrain
idx.train <- sample(1:N, Ntrain, replace=F)
idx.test <- sample(c(1:N)[-idx.train], Ntest, replace=F)
dtrain <- simData[idx.train, ]
dtest <- simData[idx.test, ]

write.table(dtrain, file = 'dtrain.txt',
            sep="\t", row.name=F, col.name=T)
write.table(dtest, file = 'dtest.txt',
            sep="\t", row.name=F, col.name=T)

rm(lSim, simData)


#-----------------------------------------------------------------------------------------------------
# Naive ATE estimate to quantify confounding and heterogeneity bias
#-----------------------------------------------------------------------------------------------------
ATEnaive = mean(dtrain$po[dtrain$Z==1]) - mean(dtrain$po[dtrain$Z==0])
ATEtrue = mean(dtrain$TE)
ATE_relBias_naive <- ifelse( TE.type==0,
                             (ATEnaive - ATEtrue) * 100,
                             (ATEnaive - ATEtrue) / ATEtrue * 100)


#-----------------------------------------------------------------------------------------------------
# ATE and ITE with Causal Forest (DSF or PF)
#-----------------------------------------------------------------------------------------------------
## Run python Causal Forest algorithm
# cat(sprintf("Fitting forest\n"))
# cmd <- "python dsf.py"
cmd <- "python pf.py"
system.time(  system(cmd) )
# B = 1000


## Extract estimates and variances
TE_hat_train <- read.table('TE_hat_train.txt')[,1]
TE_hat_test <- read.table('TE_hat_test.txt')[,1]
TE_var_train <- read.table('TE_var_train.txt')[,1]
TE_var_test <- read.table('TE_var_test.txt')[,1]

par(mfrow=c(1,2))
hist(dtrain$TE)
hist(TE_hat_train)


## In-sample ATE
ATE_true = mean(dtrain$TE)
ATE_hat = mean(TE_hat_train, na.rm = T)
ATE_se_1 = sqrt( sum(TE_var_train)/Ntrain^2 )
ATE_95CI = c( ATE_hat - 1.96 * ATE_se_1, ATE_hat + 1.96 * ATE_se_1 )
ATE_relBias_1 = c( (ATE_true - ATE_95CI[2]), (ATE_true - ATE_hat), (ATE_true - ATE_95CI[1]) ) * 100


## Out-of-sample ITE
TE_true = dtest$TE
TE_bias = TE_true - TE_hat_test
TE_rmse_1 = sqrt( mean(TE_bias^2) )
TE_se_1 = sqrt(TE_var_test)
CI_lo = TE_hat_test - 1.96 * TE_se_1
CI_up = TE_hat_test + 1.96 * TE_se_1
# cbind(CI_lo, CI_up)
TE_covRate_1 = sum(CI_lo < TE_true & TE_true < CI_up) / Ntest

print(list(ATE_bias, TE_covRate, TE_rmse))


#-----------------------------------------------------------------------------------------------------
# Save results
#-----------------------------------------------------------------------------------------------------
# out <- list(seed=seed, conf.case=conf.case, TE.type=TE.type, Ntrain=Ntrain,
#             ATE_relBias=ATE_relBias, TE_rmse=TE_rmse, TE_covRate=TE_covRate)
# 
# return(out)
# }
