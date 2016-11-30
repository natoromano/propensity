############# TUNE CALIPER FOR SIMILARITY METHODS #############
# required functions
source("matching_functions.R")

calipers <- c(-1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

# change for various datasets
dval <- read.table("/scratch/users/naromano/OSIM2/dval2.txt", header=T)
resultsArrayCal <- array(dim=c(1, 10, length(calipers)))
dval$id <- 1:nrow(dval)

empvariables <- setdiff(colnames(dval),
                        c("Z", "Y", "e", "po", "tau", "TE"))
xmat <- as.matrix(dval[, empvariables])

desiredOrder <- c("n1", "ORmatched", "ORlow_matched", "ORupp_matched",
                  "coeff_matched", "se_matched", "bias_matched", "inCI_matched",
                  "n0", "ORadj", "ORlow_adj", "ORupp_adj", "coeff_adj",
                  "se_adj", "bias_adj", "inCI_adj", "SMD", "KS", "KL", "Cstat")

trueCoeff <- mean(dval$TE)
trueOR <- exp(trueCoeff)
Noutcomes <- sum(dval$Y)
id <- "id"
exposed <- "Z"
outcome <- "Y"

xmat_ctrl <- xmat[dval[, exposed]==0,]
xmat_trted <- xmat[dval[, exposed]==1,]
rownames(xmat_ctrl) <- dval[dval[, exposed]==0, id]
rownames(xmat_trted) <- dval[dval[, exposed]==1, id]

for (p in 1:length(calipers)) {
  # compute matched sets
  matchedEuclidean <- matchByDist(xmat_ctrl, xmat_trted, method="euclidean", 
                                  k_neighbors=5, caliper=calipers[p])
  # extract results
  euclideanResults <- extractResults(ps=matchedEuclidean, exposurevec=NULL,
                                     data=dval, fmod=NULL, id=id, exposed=exposed,
                                     outcome=outcome, logitFlag=F, outfile=NULL,
                                     verbose=FALSE)
  
  # consolidate results
  resultsmat <- as.data.frame(rbind(euclideanResults[[1]]))
  rownames(resultsmat) <- c("euclidean")
  resultsmat$bias_matched <- sum(resultsmat$coeff_matched, -trueCoeff)
  resultsmat$inCI_matched <- (resultsmat$ORlow_matched<=trueOR &
                                resultsmat$ORupp_matched>=trueOR)
  resultsArrayCal[,,p] <- as.matrix(resultsmat[, desiredOrder[c(1:9,17)]])
}

dimnames(resultsArrayCal)[1:3] <- list(rownames(resultsmat),
                                       desiredOrder[c(1:9,17)],
                                       calipers)
save(resultsArrayCal, trueOR, trueCoeff, file="tuneCalSim.RData")
