###############################################################################
# (Script)
# Apply baseline methods to OSIM2 simulated data.
#
# Nathanael Romano
###############################################################################
rm(list=ls())
setwd("/home/naromano/propensity")
source("R/utils.R")

############# IMPORTS #############
installnewpackage(
  c("Matching", "Epi", "glmnet", "randomForest", 
    "vegan", "FNN", "Matrix", "doParallel", "foreach")
)

require(Matching)
require(Epi)  # clogistic
require(glmnet)
require(randomForest)
require(vegan)
require(FNN)
require(Matrix)

# parallel core computing for lassoMV bootstrapping
require(doParallel)
require(foreach)

# specify number of cores req
cl <- makeCluster(10) 
registerDoParallel(cl)

############# SET-UP #############
# required functions
source("R/matching_functions.R")

# load data
dtrain <- read.table("/scratch/users/naromano/OSIM2/dval2.txt", header=T)

# define covariate names
id <- "id"
exposed <- "Z"
outcome <- "Y"
effect <- "TE"
# empirical variables may be matched depending on HD-PS
empvariables <- setdiff(
  colnames(dtrain), c("Z", "Y", "e", "po", "tau", "TE")
)
empvariables_num <- c("YOB")
empvariables_cat <- setdiff(empvariables, empvariables_num)
form = as.formula(
  paste(exposed, paste(empvariables, collapse=" + "), sep=" ~ ")
)

# simpler formula for direct matching
empvariables_subset = c("YOB", "GENDER", 
                        empvariables[sample(3:length(empvariables), size=100)])
form_subset = as.formula(
  paste(exposed, paste(empvariables_subset, collapse=" + "), sep=" ~ ")
)

# prepare result arrays
startAll <- proc.time()[1:3]
Nmethods <- 3
smdmatArray <- array(dim=c(Nmethods + 1, length(empvariables)))
time <- array(dim=c(Nmethods, 3))
dimnames(time)[[1]] <- c("naive", "PSM", "similarity")
dimnames(time)[[2]] <- c("user", "system", "elapsed")
matchedID <- list()

dtrain$id <- 1:nrow(dtrain)
nminor <- min(table(dtrain[, exposed]))
Nexposed <- sum(dtrain[, exposed])
Noutcomes <- sum(dtrain[, outcome])

############# PSM #############
print("############ PSM ############")
start = proc.time()[1:3]

# have to compute dircectly as calling from matching_function bugs
m.out <- matchit(formula=form, data=dtrain, method="nearest", 
                 distance="linear.logit", caliper=0.5, m.order="random")
matchedData <- match.data(m.out)
lassoResults <- extractResults(matchedData, exposed, outcome, verbose=FALSE,
                               method="DA", fmod=NULL)
lassoResults[[1]]["bias"] <- lassoResults[[1]]["coeff"] - 
  mean(matchedData[, effect])
lassoResults[[1]]["MSE"] <- mean(
  (lassoResults[[1]]["coeff"] - matchedData[, effect]) ** 2
)

end <- proc.time()[1:3]
time["PSM",] <- end-start

############# SIMILARITY #############
print("############ similarity ############")
start <- proc.time()[1:3]

matchedData <- match(data=dtrain, formula=form_subset, representation=NULL, 
                     fmod=NULL, id=id, exposed=exposed, outcome=outcome, 
                     method="similarity")
similarityResults <- extractResults(matchedData, exposed, outcome, 
                                    verbose=FALSE, method="DA", fmod=NULL)
similarityResults[[1]]["bias"] <- similarityResults[[1]]["coeff"] - 
                                    mean(matchedData[, effect])
similarityResults[[1]]["MSE"] <- mean(
  (similarityResults[[1]]["coeff"] - matchedData[, effect]) ** 2
)

end <- proc.time()[1:3]
time["similarity",] <- end-start


####### NAIVE ########
print("############ naive model ############")
start <- proc.time()[1:3]

naiveResults <- extractResults(
  dtrain, exposed, outcome, verbose=FALSE,
  method="naive", fmod=NULL
)
naiveResults[[1]]["bias"] <- naiveResults[[1]]["coeff"] - 
                                mean(matchedData[, effect])
naiveResults[[1]]["MSE"] <- mean(
  (naiveResults[[1]]["coeff"] - matchedData[, effect]) ** 2
)

end <- proc.time()[1:3]
time["naive",] <- end-start

# consolidate matchedIDs
matchedID <- list(naiveResults[[2]], lassoResults[[2]], similarityResults[[2]])
names(matchedID) <- dimnames(time)[[1]]

# check if baseline variables are balanced
beforeMatching <- matchStats(numvar=empvariables_num, catvar=empvariables_cat,
                             treatment=exposed, data=dtrain,
                             outXlsx=NULL, verbose=FALSE)
afterMatching <- list()
for(j in c(3, 2, 1)){
  # naive req special handling (shouldn't be matched)
  if (j==1) {
    afterMatching[[j]] <- afterMatching[[3]]
  } else {
    afterMatching[[j]] <- matchStats(numvar=empvariables_num,
                                     catvar=empvariables_cat,
                                     treatment=exposed,
                                     data=dtrain[matchedID[[j]],],
                                     outXlsx=NULL, verbose=FALSE)
  }
}
names(afterMatching) <- names(matchedID)

# extract SMD matrix
# initialize smdmat vector with beforeMatching
smdmat <- extractSmd(beforeMatching)
# append smdmat vector with afterMatching
smdmat <- cbind(smdmat, sapply(afterMatching, extractSmd))
colnames(smdmat)[1] <- "before"
smdmat <- cbind(smdmat)
smdmat[, "naive"] <- NA

# consolidate results
resultsmat <- as.data.frame(rbind(naiveResults[[1]],
                                  lassoResults[[1]],
                                  similarityResults[[1]]))
rownames(resultsmat) <- dimnames(time)[[1]]
resultsArray <- as.matrix(resultsmat)
smdmatArray <- t(smdmat[empvariables,])

endAll <- proc.time()[1:3]
endAll-startAll

# remember to terminate the cores when done
stopCluster(cl)

dimnames(smdmatArray) <- list(colnames(smdmat), empvariables)
save(resultsArray, time, matchedID, Noutcomes,
     Nexposed, smdmatArray, file="baseline_results.RData")

closeAllConnections()
