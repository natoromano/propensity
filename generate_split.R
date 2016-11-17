# This script uses the function in simulate_effect to generate
# the dataset and split it into training, validation and test sets.

rm(list=ls())
setwd("/home/naromano/propensity/")
source("simulate_effect.R")
covariate <- read.csv('/labs/shahlab/OSIM2/covariate.csv', check.names=T)


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
Ntrain <- ceiling(0.55 * Nsamp)
Nval1 <- ceiling(0.25 * Nsamp)
Nval2 <- ceiling(0.12 * Nsamp)
Ntest <- Nsamp - Ntrain - Nval1 - Nval2
idx.train <- sample(1:N, Ntrain, replace=F)
idx.val1 <- sample(c(1:N)[-idx.train], Nval1, replace=F)
idx.val2 <- sample(c(1:N)[-idx.train][-idx.val1], Nval2, replace=F)
idx.test <- sample(c(1:N)[-idx.train][-idx.val1][-idx.val2], Ntest, replace=F)
dtrain <- simData[idx.train, ]
dval1 <- simData[idx.val1, ]
dval2 <- simData[idx.val2, ]
dtest <- simData[idx.test, ]

write.table(dtrain, file='/scratch/users/naromano/OSIM2/dtrain.txt',
            sep="\t", row.name=F, col.name=T)
write.table(dval1, file='/scratch/users/naromano/OSIM2/dval1.txt',
            sep="\t", row.name=F, col.name=T)
write.table(dval2, file='/scratch/users/naromano/OSIM2/dval2.txt',
            sep="\t", row.name=F, col.name=T)
write.table(dtest, file='/scratch/users/naromano/OSIM2/dtest.txt',
            sep="\t", row.name=F, col.name=T)

rm(lSim, simData)

