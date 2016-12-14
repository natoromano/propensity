###############################################################################
# (Script)
# Apply baseline methods to OSIM2 simulated data.
#
# Nathanael Romano
###############################################################################
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

rm(list=ls())
setwd("/home/naromano/propensity")
source("R/utils.R")

############# IMPORTS #############
installnewpackage(c("rJava", "Matching", "Epi", "glmnet", "randomForest", 
                    "vegan", "FNN", "Matrix", "doParallel", "foreach"))

require(rJava)
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

# initialize the Java subsystem. include a jdbc object
.jinit(classpath="/home/naromano/software/pharmacoepi.jar", 
       force.init=TRUE, parameters="-Xmx12gm");
.jclassPath()

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
                        empvariables[sample(3:length(empvariables), size=30)])
form_subset = as.formula(
  paste(exposed, paste(empvariables_subset, collapse=" + "), sep=" ~ ")
)

# prepare result arrays
startAll <- proc.time()[1:3]
Nmethods <- 3 # 4 baselines: unadjusted, hdps, psm, euclidean matching
resultsArray <- array(dim=c(Nmethods, 20))  
pmatArray <- array(dim=c(Nmethods + 1, length(empvariables)))
smdmatArray <- pmatArray
time <- array(dim=c(Nmethods, 3))
dimnames(time)[[1]] <- c("naive", "PSM", "euclidean")
dimnames(time)[[2]] <- c("user", "system", "elapsed")
dimnames(resultsArray)[[1]] <- dimnames(time)[[1]]
matchedID <- list()

dtrain$id <- 1:nrow(dtrain)
nminor <- min(table(dtrain[, exposed]))
trueCoeff <- mean(dtrain$tau)
trueOR <- exp(trueCoeff)
Nexposed <- sum(dtrain[, exposed])
Noutcomes <- sum(dtrain[, outcome])

############# PSM #############
print("############ PSM ############")
start = proc.time()[1:3]

# matchedData <- match(data=dtrain, formula=form, representation=NULL, fmod=NULL, 
#                     id=id, exposed=exposed, outcome=outcome, method="PSM")
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

if (FALSE) {
  xmat = as.matrix(dtrain[, empvariables])
  xmat_ctrl <- xmat[dtrain[, exposed]==0,]
  xmat_trted <- xmat[dtrain[, exposed]==1,]
  rownames(xmat_ctrl) <- dtrain[dtrain[, exposed]==0, id]
  rownames(xmat_trted) <- dtrain[dtrain[, exposed]==1, id]
  rm(xmat)
  
  runSimPS <- function(method="euclidean", caliper=0.7, nsd=3,
                       algorithm="kd_tree") {
    print("STARTING TO MATCH")
    matchedSim <- matchByDist(xmat_ctrl, xmat_trted, method=method,
                              k_neighbors=5, caliper=caliper, nsd=nsd,
                              algorithm=algorithm)
    print("ENDED MATCHING")
    temp <- cbind(c(names(matchedSim$matchedset), matchedSim$matchedset), 
                  set_num=rep(1:length(matchedSim$matchedset), 2))
    mode(temp) <- "numeric"
    rownames(temp) <- NULL
    colnames(temp) <- c("pat_id", "set_num")
    print("ABOUT TO MERGE")
    matcheddata <- merge(data, temp, by.x=id, by.y="pat_id", all.y=T, sort=F)
    print("ENDED MERGING")
    simResults <- extractResults(matchedData,
                                 exposed, outcome, verbose=FALSE, 
                                 method="clogit", fmod=NULL)
    return(simResults)
  }
  
  start <- proc.time()[1:3]
  euclideanResults <- runSimPS(method="euclidean", caliper=0.8,
                               nsd=3, algorithm="kd_tree")
  euclideanResults[[1]] <- c(n0=euclideanResults[[1]]["n0"], 
                             n1=euclideanResults[[1]]["n1"],
                             coeff=euclideanResults[[1]]["coeff_matched"],
                             se=euclideanResults[[1]]["se_matched"],
                             smd=euclideanResults[[1]]["SMD"])
}

matchedData <- match(data=dtrain, formula=form_subset, representation=NULL, 
                     fmod=NULL, id=id, exposed=exposed, outcome=outcome, 
                     method="similarity")
euclideanResults <- extractResults(matchedData, exposed, outcome, verbose=FALSE,
                                   method="DA", fmod=NULL)
euclideanResults[[1]]["bias"] <- euclideanResults[[1]]["coeff"] 
- mean(matchedData[, effect])
euclideanResults[[1]]["MSE"] <- mean(
  (euclideanResults[[1]]["coeff"] - matchedData[, effect]) ** 2
)

end <- proc.time()[1:3]
time["euclidean",] <- end-start


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
matchedID <- list(naiveResults[[2]], lassoResults[[2]],euclideanResults[[2]])
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

# extract smd matrix
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
                                  euclideanResults[[1]]))
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
