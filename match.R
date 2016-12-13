###############################################################################
# (Script)
# Matches on all dimensions's representations.
# Uses validation set 2.
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

baseDir <- getwd();
tempDir <- paste(baseDir, "tmp", sep="/") # output directory
patientfile <- paste(baseDir, "patients.txt", sep="/")
dir.create(tempDir)

# initialize the Java subsystem. include a jdbc object
.jinit(classpath="/home/naromano/software/pharmacoepi.jar", 
       force.init=TRUE, parameters="-Xmx12gm");
.jclassPath()

# required functions
source("R/matching_functions.R")

# load data
dval <- read.table("/scratch/users/naromano/OSIM2/dval2.txt", header=T)

# define covariate names
id <- "id"
exposed <- "Z"
outcome <- "Y"
empvariables <- setdiff(colnames(dval), c("Z", "Y", "e", "po", "tau", "TE"))
empvariables_num <- c("YOB")
empvariables_cat <- setdiff(empvariables, empvariables_num)

desiredOrder <- c("n1", "ORmatched", "ORlow_matched", "ORupp_matched",
                  "coeff_matched", "se_matched", "bias_matched", "inCI_matched",
                  "n0", "ORadj", "ORlow_adj", "ORupp_adj", "coeff_adj",
                  "se_adj", "bias_adj", "inCI_adj", "SMD", "KS", "KL", "Cstat")

# prepare result arrays
startAll <- proc.time()[1:3]
Nmethods <- 9
resultsArray <- array(dim=c(Nmethods, 20))  
pmatArray <- array(dim=c(Nmethods + 1, length(empvariables)))
smdmatArray <- pmatArray
time <- array(dim=c(Nmethods, 3))
dimnames(time)[[1]] <- 2:10
dimnames(time)[[2]] <- c("user", "system", "elapsed")
dimnames(resultsArray)[[1]] <- dimnames(time)[[1]]
matchedID <- list()

dval$id <- 1:nrow(dval)
nminor <- min(table(dval[, exposed]))
trueCoeff <- mean(dval$tau)
trueOR <- exp(trueCoeff)
Nexposed <- sum(dval[, exposed])
Noutcomes <- sum(dval[, outcome])

results <- c()
matchedID <- list()

for (dim in 2:10) {
  name = paste("layerval2_", dim, ".txt", sep="")
  xmat = read.table(paste("/scratch/users/naromano/OSIM2/", name, sep=""), 
                    header=F)
  xmat_ctrl <- xmat[dval[, exposed]==0,]
  xmat <- as.data.frame(scale(xmat))
  xmat_ctrl <- xmat[dval[, exposed]==0,]
  xmat_trted <- xmat[dval[, exposed]==1,]
  rownames(xmat_ctrl) <- dval[dval[, exposed]==0, id]
  rownames(xmat_trted) <- dval[dval[, exposed]==1, id]
  
  runSimPS <- function(method="euclidean", caliper=0.7, nsd=3,
                       algorithm="kd_tree"){
    matchedSim <- matchByDist(xmat_ctrl, xmat_trted, method=method,
                              k_neighbors=5, caliper=caliper, nsd=nsd,
                              algorithm=algorithm)
    simResults <- extractResults(ps=matchedSim, exposurevec=NULL,
                                 data=dval, fmod=NULL, id=id, exposed=exposed,
                                 outcome=outcome, logitFlag=F,
                                 outfile=NULL, verbose=FALSE)
    return(simResults)
  }
  
  start <- proc.time()[1:3]
  euclideanResults <- runSimPS(method="euclidean", caliper=0.8,
                               nsd=0.25, algorithm="brute")
  results[dim-1] <- euclideanResults
  matchedID[[dim-1]] <- euclideanResults[[2]]
  
  end <- proc.time()[1:3]
  time[dim-1,] <- end-start

}

names(matchedID) <- dimnames(time)[[1]]

# check if baseline variables are balanced
beforeMatching <- matchStats(numvar=empvariables_num, catvar=empvariables_cat,
                             treatment=exposed, data=dval,
                             outXlsx=NULL, verbose=FALSE)
afterMatching <- list()
for(j in 1:9){
    afterMatching[[j]] <- matchStats(numvar=empvariables_num,
                                     catvar=empvariables_cat,
                                     treatment=exposed,
                                     data=dval[matchedID[[j]],],
                                     outXlsx=NULL, verbose=FALSE)
}
names(afterMatching) <- names(matchedID)

# extract pval matrix
# initialize pmat vector with beforeMatching
pmat <- extractPval(beforeMatching)

# append pmat vector with afterMatching
pmat <- cbind(pmat, sapply(afterMatching, extractPval))
colnames(pmat)[1] <- "before"

# extract smd matrix
# initialize smdmat vector with beforeMatching
smdmat <- extractSmd(beforeMatching)
# append smdmat vector with afterMatching
smdmat <- cbind(smdmat, sapply(afterMatching, extractSmd))
colnames(smdmat)[1] <- "before"
smdmat <- cbind(smdmat)

# consolidate results
resultsmat <- results[1][[1]]
for (i in 2:9) {
  resultsmat <- rbind(resultsmat, results[i][[1]])
}
resultsmat <- as.data.frame(resultsmat)
rownames(resultsmat) <- dimnames(time)[[1]]
resultsmat$bias_adj <- resultsmat$coeff_adj-trueCoeff
resultsmat$bias_matched <- resultsmat$coeff_matched-trueCoeff
resultsmat$inCI_adj <- (resultsmat$ORlow_adj<=trueOR &
                          resultsmat$ORupp_adj>=trueOR)
resultsmat$inCI_matched <- (resultsmat$ORlow_matched<=trueOR & 
                              resultsmat$ORupp_matched>=trueOR)
resultsArray <- as.matrix(resultsmat[, desiredOrder])
pmatArray <- t(pmat[empvariables,])
smdmatArray <- t(smdmat[empvariables,])

endAll <- proc.time()[1:3]
endAll-startAll

# remember to terminate the cores when done
stopCluster(cl)

dimnames(resultsArray) <- list(rownames(resultsmat), desiredOrder)
dimnames(pmatArray) <- list(colnames(pmat), empvariables)
dimnames(smdmatArray) <- list(colnames(smdmat), empvariables)
save(resultsArray, time, matchedID, Noutcomes,
     Nexposed, trueCoeff, pmatArray, smdmatArray, file="results.RData")
