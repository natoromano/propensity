###############################################################################
# (Script)
# Apply baseline methods to OSIM2 simulated data.
#
# Nathanael Romano
###############################################################################
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

rm(list=ls())
setwd("/home/naromano/propensity")
source("utils.R")

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
source("matching_functions.R")

# load data
dtrain <- read.table("/scratch/users/naromano/OSIM2/dtrain.txt", header=T)

# define covariate names
id <- "id"
exposed <- "Z"
outcome <- "Y"
# empirical variables may be matched depending on HD-PS
empvariables <- setdiff(colnames(dtrain),
                        c("Z", "Y", "e", "po", "tau", "TE"))
empvariables_num <- c("YOB")
empvariables_cat <- setdiff(empvariables, empvariables_num)

desiredOrder <- c("n1", "ORmatched", "ORlow_matched", "ORupp_matched",
                  "coeff_matched", "se_matched", "bias_matched", "inCI_matched",
                  "n0", "ORadj", "ORlow_adj", "ORupp_adj", "coeff_adj",
                  "se_adj", "bias_adj", "inCI_adj", "SMD", "KS", "KL", "Cstat")

# prepare result arrays
startAll <- proc.time()[1:3]
Nmethods <- 3 # 4 baselines: unadjusted, hdps, psm, euclidean matching
resultsArray <- array(dim=c(Nmethods, 20))  
pmatArray <- array(dim=c(Nmethods + 1, length(empvariables)))
smdmatArray <- pmatArray
time <- array(dim=c(Nmethods, 3))
dimnames(time)[[1]] <- c("unadjusted", "PSM", "euclidean")
dimnames(time)[[2]] <- c("user", "system", "elapsed")
dimnames(resultsArray)[[1]] <- dimnames(time)[[1]]
matchedID <- list()

dtrain$id <- 1:nrow(dtrain)
nminor <- min(table(dtrain[, exposed]))
trueCoeff <- mean(dtrain$TE)
trueOR <- exp(trueCoeff)
Nexposed <- sum(dtrain[, exposed])
Noutcomes <- sum(dtrain[, outcome])

############# hdPS #############
if (FALSE) {
  # skipping for now as too many dimensions
  # TODO: only keep top 100 correlated dimensions w/ Z and Y?
  print("############ hdPS ############")
  start <- proc.time()[1:3]
  
  # prepare data into format required by hdps
  # read data file as single string (to use addPatientsFromBuffer)
  patientheader <- paste(id, exposed, outcome, sep="\t")
  patientstring <- paste(paste(dtrain[, id], dtrain[, exposed], 
                               dtrain[, outcome],
                               sep="\t"), collapse="\n")
  # for reading in whole file as single string
  # patientstring=readChar(patientfile, file.info(patientfile)$size)
  # patientstring=gsub("\"","",patientstring)
  # patients=dataset[,c("pid_org",exposed,outcome[1])] #get patient id, exp, out
  # write.table(patients,file=patientfile,sep="\t",
  # col.names=T,row.names=F,na="") # output data to patient.txt;
  datainstring <- paste(patientheader, patientstring, sep="\n")
  variables <- empvariables_cat
  dimdata <- dtrain[, c(id, variables)]
  
  ### INVOKE pharmacoepi.jar ### (handles cat variables only)
  hdpsobj <- hdps(datainstring, dimdata, outDir=tempDir, Nmostfreq=10, k=10,
                  stratifyDim=FALSE, outfile="output_cohort.txt",
                  FullOutput=TRUE, verbose=F, ZeroCellCorrection=F)
  hdpsobj$selectedvariables # (see Fig2, Schneeweiss 2008)
  wantedvar <- hdpsobj$selectedvariables[grep("1Once$",
                                              hdpsobj$selectedvariables)]
  
  var_corExposed <- empvariables_num[abs(cor(dtrain[, exposed],
                                             dtrain[, empvariables_num]))>0.05]
  var_corOutcome <- empvariables_num[abs(cor(dtrain[, outcome], # include in PS
                                             dtrain[, empvariables_num]))>0.05]
  IV <- setdiff(var_corExposed, var_corOutcome) # exclude from PS
  
  # Estimate the PS (force numerical variables into PS model)
  dataPS <- cbind(dtrain[,c(id, exposed, var_corOutcome)],
                  hdpsobj$hdpsdata[,wantedvar])
  hdPsMod <- glm(paste(exposed, "~ . -", id), data=dataPS, family="binomial")
  names(hdPsMod$fitted.values) <- as.character(dataPS[, id])
  summary(hdPsMod$fitted.values)
  
  hdResults <- extractResults(ps=hdPsMod$fitted, exposurevec=hdPsMod$y,
                              fmod=NULL, data=dtrain, id=id, exposed=exposed,
                              outcome=outcome, logitFlag=TRUE,
                              outfile=NULL, verbose=FALSE)
  
  end <- proc.time()[1:3]
  time["hdPS",] <- end-start
}

############# lassoPS #############
print("############ lassoPS ############")
start <- proc.time()[1:3]

xmat <- as.matrix(dtrain[, empvariables])
# penalty.factor=c(0, rep(1, ncol(xmat)-1)) # set 0 to force variable into model
# tune lambda for glmnet using 5-fold CV
lassoPsmod <- cv.glmnet(xmat, dtrain[, exposed], alpha=1, family="binomial",
                        standardize=F, nfold=5)
bestlambda <- lassoPsmod$lambda.1se
lassoPSBeta <- coeffAtlambda(lassoPsmod)[-1]  # exclude intercept

# get estimated ps in logit form
psl <- unlogit(as.numeric(predict(lassoPsmod, xmat, s=lassoPsmod$lambda.1se)))
names(psl) <- dtrain[, id]

lassoResults <- extractResults(ps=psl, exposurevec=dtrain[, exposed], fmod=NULL,
                               data=dtrain, id=id, exposed=exposed,
                               outcome=outcome, logitFlag=TRUE,
                               outfile=NULL, verbose=FALSE)

end <- proc.time()[1:3]
time["PSM",] <- end-start


############# SIMILARITY #############
print("############ similarity ############")

xmat_ctrl <- xmat[dtrain[, exposed]==0,]
xmat_trted <- xmat[dtrain[, exposed]==1,]
rownames(xmat_ctrl) <- dtrain[dtrain[, exposed]==0, id]
rownames(xmat_trted) <- dtrain[dtrain[, exposed]==1, id]

runSimPS <- function(method="euclidean", caliper=0.7, nsd=3,
                     algorithm="kd_tree"){
  matchedSim <- matchByDist(xmat_ctrl, xmat_trted, method=method,
                            k_neighbors=5, caliper=caliper, nsd=nsd,
                            algorithm=algorithm)
  simResults <- extractResults(ps=matchedSim, exposurevec=NULL,
                               data=dtrain, fmod=NULL, id=id, exposed=exposed,
                               outcome=outcome, logitFlag=F,
                               outfile=NULL,verbose=FALSE)
  return(simResults)
}

start <- proc.time()[1:3]
euclideanResults <- runSimPS(method="euclidean", caliper=0.8,
                               nsd=3, algorithm="kd_tree")
end <- proc.time()[1:3]
time["euclidean",] <- end-start


####### UNADJUSTED ########
print("############ lassoMV model ############")
start <- proc.time()[1:3]

xmat <- Matrix(as.matrix(dtrain[, c("YOB", exposed)]), 
               sparse=T)
penalty.factor <- c(0, 1)
lassoUnAdjmod <- cv.glmnet(xmat, dtrain[, outcome], alpha=1, family="binomial",
                           standardize=F, nfold=5, 
                           penalty.factor=penalty.factor)

coeff <- na.omit(sort(coeffAtlambda(lassoUnAdjmod)))

lassomodbootstrap <- genCI(xmat, dtrain[,outcome], ntrials=100,
                           lambda=lassoUnAdjmod$lambda.1se, alpha=1,
                           penalty.factor=penalty.factor)
lassomodCI <- setCL(lassomodbootstrap)

ORCInonZero <- as.data.frame(exp(lassomodCI$betaCIlim[lassomodCI$beta_nonZero,]))
Smd <- smd(dtrain[,c(exposed,outcome)], exposed=exposed, variable=outcome,
           verbose=FALSE, categorical=TRUE)
se_adj <- (lassomodCI$betaCIlim[exposed,3] - 
             lassomodCI$betaCIlim[exposed,1])/2/1.96
results <- c(n0=nrow(dtrain), n1=NA,
             KS=NA, KL=NA, Cstat=NA, SMD=Smd,
             ORadj=ORCInonZero[exposed,2],
             ORlow_adj=ORCInonZero[exposed,1],
             ORupp_adj=ORCInonZero[exposed,3],
             coeff_adj=lassomodCI$betaCIlim[exposed,2],
             se_adj=se_adj,
             rep(NA,5))
names(results) <- colnames(lassoResults[[1]])
lassoUnAdjResults <- list(results=results, matchedID=NULL)

end <- proc.time()[1:3]
time["unadjusted",] <- end-start

# consolidate matchedIDs
matchedID <- list(lassoUnAdjResults[[2]], lassoResults[[2]],
                  euclideanResults[[2]])
names(matchedID) <- dimnames(time)[[1]]

# check if baseline variables are balanced
beforeMatching <- matchStats(numvar=empvariables_num, catvar=empvariables_cat,
                             treatment=exposed, data=dtrain,
                             outXlsx=NULL, verbose=FALSE)
afterMatching <- list()
for(j in c(3, 2, 1)){
  # unadjusted req special handling (shouldn't be matched)
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

# extract pval matrix
# initialize pmat vector with beforeMatching
pmat <- extractPval(beforeMatching)

# append pmat vector with afterMatching
pmat <- cbind(pmat, sapply(afterMatching, extractPval))
colnames(pmat)[1] <- "before"
pmat[, "unadjusted"] <- NA

# extract smd matrix
# initialize smdmat vector with beforeMatching
smdmat <- extractSmd(beforeMatching)
# append smdmat vector with afterMatching
smdmat <- cbind(smdmat, sapply(afterMatching, extractSmd))
colnames(smdmat)[1] <- "before"
smdmat <- cbind(smdmat)
smdmat[, "unadjusted"] <- NA

# consolidate results
resultsmat <- as.data.frame(rbind(lassoUnAdjResults[[1]],
                                  lassoResults[[1]],
                                  euclideanResults[[1]]))
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
save(resultsArray, time, lassoPSBeta, bestlambda, matchedID, Noutcomes,
     Nexposed, trueCoeff, pmatArray, smdmatArray, file="baseline_results.RData")
