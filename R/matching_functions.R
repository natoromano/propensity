###############################################################################
# (Functions)
# Functions used for matching.
#
# Authors: Nathanael Romano / Yen Low
###############################################################################

installnewpackage(
  c("rJava", "flexmix", "ROCR", "Epi", "reshape2", "scales", "MatchIt")
)

require(flexmix)
require(rJava)
require(ROCR)
require(Epi)
require(reshape)
require(MatchIt)

logit <- function(x) log(x/(1-x))
unlogit <- function(x) 1/(1+exp(-x))


############# MAIN MATCHING AND INFERENCE #############
match <- function(data, formula, representation=NULL, fmod=NULL, 
                  id=id, exposed=exposed, outcome=outcome, method="PSM") {

  if (method == "PSM") {
    m.out <- matchit(formula=formula, data=data, method="nearest", 
                     distance="linear.logit", caliper=0.5, m.order="random")
    m.data <- match.data(m.out)
  } else if (method == "representation") {
    m.out <- matchit(formula=formula, data=representation, method="nearest", 
                     distance="mahalanobis", m.order="random")
    m.data <- match.data(m.out)
    
    m.data <- cbind(m.data, data[rownames(m.data), ])
    m.data <- m.data[, c(colnames(data), "distance")]
  } else if (method == "similarity") {
    m.out <- matchit(formula=formula, data=data, method="nearest", 
                     distance="mahalanobis", m.order="random")
    m.data <- match.data(m.out)
  }
  
  return(m.data)
}


extractResults <- function(matchedData, exposed, outcome, verbose=FALSE,
                           method="clogit", fmod=NULL) {
  if (is.null(fmod)) {
    fmod <- paste(outcome, exposed, sep=" ~ ")
  } 
  
  if (method == "clogit") {
    tryobj3 <- try(clogistic(fmod, strata="set_num", data=matchedData))
    if (class(tryobj3)!="try-error") model <- tryobj3 else model <- NULL
  } else if (method == "DA") {
    model <- glm(fmod, data=matchedData)
  } else if (method == "naive") {
    coeff = mean(matchedData[matchedData[, "Z"]==1, "Y"]) -
            mean(matchedData[matchedData[, "Z"]==0, "Y"])
    model <- NULL
  }

  # get SMD
  Smd <- smd(matchedData, exposed=exposed, variable=outcome, categorical=TRUE)
  
  if (!is.null(model)) {  # if clogit is possible
    coeff <- coef(model)[exposed]
    se_matched <- sqrt(diag(model$var)[exposed])
    print(se_matched)
    # extract results of interest
    results <- c(coeff=coeff[[1]], se=se_matched, smd=Smd[[1]])
  } else {
    results <- c(coeff=coeff, se=NA, smd=Smd[[1]])
  }
  
  list(results=results, matchedID=matchedData[, id])
}


############# HDPS #############
hdps <- function(datainstring, dimdata, outDir, Nmostfreq=100, k=500,
                 stratifyDim=FALSE, outfile="output_cohort.txt",
                 FullOutput=T, verbose=T, ZeroCellCorrection=F) {
  flush.console()
  
  # Instantiate an hd-PS object
  if (exists("x")) rm(x)
  x <- .jnew("org.drugepi.hdps.Hdps");
  
  # Set hd-PS parameters
  .jfield(x, "tempDirectory") = outDir
  .jfield(x, "fullOutputFilename") = outfile # in tempDirectory
  .jfield(x, "topN") = as.integer(Nmostfreq) # N most freq variables
  .jfield(x, "k") = as.integer(k);  # k var most assoc with outcome (by RR)
  .jfield(x, "inferServiceIntensityVars") = as.integer(stratifyDim)
  outputfile = paste(.jcall(x, "S", "getTempDirectory"), 
                     "/", .jcall(x,"S","getFullOutputFilename"), sep="");
  
  if (FullOutput==TRUE) {
    .jfield(x, "doFullOutput") <- as.integer(1)
  } else if(FullOutput==FALSE) {
    .jfield(x, "doSparseOutput") <- as.integer(1)
  } else {
    print("Set FullOutput=T for full output (default) or F for sparse output")
  }
  # ZeroCellCorrection recommended if < 150 outcome events (Rassen 2011)
  if(ZeroCellCorrection) {
    .jfield(x, "useOutcomeZeroCellCorrection") <- as.integer(1)
  }
  
  # Get and display the object version
  if (verbose==TRUE) {
    print(.jcall(x, "S", "getDescription"))
    print(.jcall(x, "S", "getVersion"))
    cat("Outcome type:", .jcall(x, "S", "getOutcomeType"), "\n")
    cat("FullOutputMode =", .jcall(x, "I", "getDoFullOutput"), "\n")
    cat("Outputfile =", outputfile, "\n")
    cat("Stratify Dimensions? AKA InferServiceIntensityVars",
        .jcall(x,"I","getInferServiceIntensityVars"), "\n")
  }
  
  # addPatients from file does not work; 
  # use addPatientsFromBuffer instead (buffer in string)
  # .jcall(x, "V", "addPatients", paste(baseDir, "patients.txt", sep=""));
  .jcall(x, "V", "addPatientsFromBuffer", datainstring);
  
  # addDimension is required; 
  # At least 1 dimension (2 columns: id, dimension; 3rd column date is optional)
  # create dimension files (one per empvariables) 
  # as needed for pharmacoepi.jar (addDimensions)
  dimfile <- paste(outDir, "dim.txt", sep="/")
  
  for (j in 2:ncol(dimdata)) {
    write.table(dimdata[, c(1, j)], file=dimfile, sep="\t", col.names=T,
                row.names=F, na="", quote=F)
    tryobj <- try(.jcall(x, "V", "addDimension", colnames(dimdata)[j], dimfile), 
                  silent=F)
  }
  
  # check number of dimensions
  #	try(.jcall(x, "I", "getNumDimensions"), silent=T)
  cat("Number of variables considered: ", try(.jcall(x, "I", 
                                                     "getNumDimensions"),
                                              silent=F),"\n")
  
  # run the hd-PS algorithm
  print("Running hd-PS. Please wait...")
  .jcall(x, "V", "run")
  
  # read in the output cohort file generated by PE toolbox
  tryobj2 <- try(read.table(outputfile, header=TRUE, fill=TRUE));
  if (class(tryobj2)!="try-error") {
    hdpsdata <- tryobj2	
    # get the names of the selected variables; prepend hdpsVars$ to each
    selectedvariables <- colnames(hdpsdata)[-1]
    list(jobj=x, hdpsdata=hdpsdata, selectedvariables=selectedvariables)
  } else {
    list(jobj=x, hdpsdata=NULL, selectedvariables=NULL)
  }
}


############# SMD ############# 
smd <- function(data, exposed=exposed, variable=outcome, categorical=FALSE) {

  # Computes standardized means difference
  if (categorical==FALSE) {  # continuous variables
    ngrps <- table(data[, exposed], useNA="ifany")
    means <- tapply(data[, variable], data[,exposed], mean, na.rm=T)
    sds <- tapply(data[, variable], data[,exposed], sd, na.rm=T) # sample sd
    spooled <- (ngrps[1]-1)*sds[1]^2 + (ngrps[2]-2)*sds[2]^2
    spooled <- sqrt (spooled / (ngrps[1]+ngrps[2]-2))
    smdval <- (means[2] - means[1]) / spooled
  } else {  # categorical variables
    ns <- table(data[data[, exposed] %in% c(0, 1), c(variable, exposed)], 
                useNA="ifany")
    pctTab <- prop.table(ns, 2)
    sds <- t(apply(pctTab, 1, function(x) c(x[1]*(1-x[1]), x[2]*(1-x[2]))))
    spooled <- apply(sds, 1, function(x) sqrt(sum(x)/2))
    smdval <- ((pctTab[, 2]-pctTab[, 1]) / spooled)[2]
  }

  return(smdval)
}


############# GET OR ############# 
getOR <- function(model, verbose=TRUE){
  OR <- exp(model$coefficients[exposed])
  ORint <- suppressMessages(exp(confint(model)[exposed,]))
  if (verbose==TRUE) {
    cat(OR,"[",ORint,"]\n")
  }
  list(OR=OR, ORint=ORint)
}


############# SIMILARITY ############# 
similarity <- function(x, y=NULL, method="jaccard"){
  # Calculate similarities
  xymat <- tcrossprod(x, y)
  x2sum <- diag(tcrossprod(x, x))
  if(!is.null(y)) y2sum <- diag(tcrossprod(y, y)) else y2sum <- x2sum
  s <- matrix(nrow=length(y2sum), ncol=nrow(x))
  
  if (method=="jaccard") {
    for(i in 1:nrow(x)) for(j in 1:nrow(y)) s[j,i] <- xymat[i,j]/(x2sum[i]+y2sum[j]-xymat[i,j])
  } else if(method=="dice") {
    for(i in 1:nrow(x)) for(j in 1:nrow(y)) s[j,i] <- 2*xymat[i,j]/(x2sum[i]+y2sum[j])
  } else if(method=="cosine") {
    for(i in 1:nrow(x)) for(j in 1:nrow(y)) s[j,i] <- xymat[i,j]/sqrt(x2sum[i]*y2sum[j])
  } else if(method %in% c("pearson","spearman","kendall")) {
    s <- cor(t(y), t(x), method=method)
  } else {
    stop(paste("Error: enter a valid type of similarity (jaccard, dice",
               ",cosine, pearson, spearman, kendall)!"))
  }
  rownames(s) <- rownames(y)
  colnames(s) <- rownames(x)
  return(s)
}


############# NEAREST NEIGHBOR ############# 
NNmat <- function(simmat, k_neighbors=NULL, type="sim") {
  # Create nearest neighbor object
  # Inputs: 
  # - simmat: S matrix returned by "similarity"
  # - k_neighbors: k
  if (is.null(k_neighbors)) k_neighbors <- ncol(simmat)
  if (type=="sim") decFlag <- TRUE else decFlag <- FALSE  
  NNdist <- t(apply(simmat, 1, sort, decreasing=decFlag))[, 1:k_neighbors]
  NNID <- t(apply(simmat, 1, order, decreasing=decFlag))[, 1:k_neighbors]
  NN <- matrix(colnames(simmat)[NNID], nrow=nrow(simmat), ncol=k_neighbors,
               byrow=FALSE)[, 1:k_neighbors]
  list(NN=NN, NNdist=NNdist, NNID=NNID)
}


############# DUPLICATE NN HANDLING ############# 
uniqueNN <- function(relevantNN, NNnames) {
  # Loop to handle duplicate neighbors. Get next unrepeated neighbor.
  # check for duplicated first NN
  dupNN <- names(relevantNN[duplicated(relevantNN[, 1]), 1])
  
  for (i in 1:length(dupNN)) {
    # print(i)
    j <- 2
    temp <- relevantNN[dupNN[i], j]  # get 2nd nearest neighbor
    # print(temp)
    while (!is.na(temp)) {
      # print("temp is NOT na")
      # check jth NN has not been repeated among first NN
      if (!(temp %in% relevantNN[, 1])) {
        relevantNN[dupNN[i], 1]=temp  # replace jth NN as first NN
        # cat(j," as NN\n")
        break  # end case
      } else {  # jth NN among first NN. Thus, get (j+1)th NN
        # cat(j," among 1NN\n")
        j <- j+1
        if (j<=ncol(relevantNN)) {
          temp <- relevantNN[dupNN[i], j]
        } else {
          temp <- NA
          break
        }
      }
    }
    if(is.na(temp)){  #if next NN is NA
      relevantNN <- relevantNN[!(rownames(relevantNN) %in% dupNN[i]),] 
      # remove entire row
      # print("temp is na. Drop Row")
    }
  }
  
  # check that relevantNN has no duplicates
  if (sum(duplicated(relevantNN[,1]))>0) print("Error: NN is duplicated!")
  
  # check that relevantNN has no NA 
  relevantNN <- relevantNN[!is.na(relevantNN[, 1]), ]
  cat("Patients with unique NN: ", nrow(relevantNN), "\n")
  
  # map rowID numbers to xmat_ctrl names
  vec <- NNnames[relevantNN[, 1]]
  names(vec) <- rownames(relevantNN)
  return(vec)
}


############# DISTANCE MATCHING ############# 
matchByDist <- function(xmat_ctrl, xmat_trted, method="jaccard",
                        k_neighbors=10, caliper=NULL, nsd=3, 
                        algorithm="kd_tree") {
  # Putting the previous three functions together...
  flush.console()
  print(paste("xmat_ctrl and xmat_trted should be auto-scaled and ",
              "have patient IDs as row names."))
  
  if (method %in% c("jaccard", "dice", "cosine", "pearson", 
                    "spearman", "kendall")) {
    # checked that the seeds are the minor class
    # create similarity matrix
    if (nrow(xmat_ctrl)>=nrow(xmat_trted)) {
      simmat <- similarity(xmat_ctrl[sample(nrow(xmat_ctrl)),], 
                           xmat_trted[sample(nrow(xmat_trted)),], method)
    }else{
      simmat <- similarity(xmat_trted[sample(nrow(xmat_trted)),], 
                           xmat_ctrl[sample(nrow(xmat_ctrl)),], method)
    }
    NNobj <- NNmat(simmat, k_neighbors <- k_neighbors, type="sim")
    
    if (is.null(caliper)) {  # if caliper not specified, manually set it
      if (method %in% c("jaccard")) {
        caliper <- 0.6
      } else if(method %in% c("dice", "cosine", "pearson",
                              "spearman", "kendall")) {
        caliper <- 0.7
      }
    }
    
  } else if(method=="euclidean") {
    if (nrow(xmat_ctrl)>=nrow(xmat_trted)) {
      fNNobj <- get.knnx(xmat_ctrl, query=xmat_trted, k=k_neighbors,
                         algorithm=algorithm)
    } else {
      fNNobj <- get.knnx(xmat_trted, query=xmat_ctrl, k=k_neighbors,
                         algorithm=algorithm)
    }
    NNobj <- list(NN=NULL, NNdist=fNNobj$nn.dist, NNID=fNNobj$nn.index)
    if (is.null(nsd)) nsd <- 3
    if (is.numeric(nsd)) {
      caliper <- mean(NNobj$NNdist[, 1]) + nsd*sd(NNobj$NNdist[, 1]) 
    } else {
      stop("Error: nsd must be a number. Set 3 for 3 std deviations.")
    }
  }
  
  cat("Method: ", method, "\n")
  cat("Caliper: ", caliper, "\n")
  
  relevantNN <- NNobj$NNID
  # (note that this is in the original row ID of the ctrl group
  # Either NN or NNID could do.)
  
  cat("Patients with NN before applying caliper: ", nrow(relevantNN), "\n")
  if (!is.null(rownames(xmat_trted)) & !is.null(rownames(xmat_ctrl))) {
    if (nrow(xmat_ctrl)>=nrow(xmat_trted)) {
      rownames(relevantNN) <- rownames(xmat_trted)
      majorNames <- rownames(xmat_ctrl)
    } else {
      rownames(relevantNN) <- rownames(xmat_ctrl)
      majorNames <- rownames(xmat_trted)
    }
  } else stop("Error: xmat_trted or xmat_ctrl must have rownames!")
  if (method=="euclidean") {
    relevantNN[NNobj$NNdist>caliper] <- NA
  } else {
    relevantNN[NNobj$NNdist<caliper] <- NA
  }
  relevantNN <- relevantNN[!is.na(relevantNN[,1]),] # drop row if first NN is NA
  cat("Patients with NN after applying caliper: ", nrow(relevantNN), "\n")
  
  tryobj2=try(uniqueNN(relevantNN, NNnames=majorNames))
  if (class(tryobj2)!="try-error") vec <- tryobj2 else vec <- NULL
  list(matchedset=vec, NNobj=NNobj)
}

###### SMD EXTRACTION ######
extractSmd <- function(MatchingObj){
  # initialize smd vector with pvalues of numerical variables
  if (!is.null(MatchingObj$numvar[, "smd"])) {
    smd <- as.numeric(MatchingObj$numvar[, "smd"])
  } else {
    smd <- c()
  }
  names(smd) <- c(rownames(MatchingObj$numvar))
  # append with pvalues of categorical variables
  smd <- c(smd, sapply(MatchingObj$catvar, function(x) x[["smd"]][[1]]))
  return(smd)
}


# MATCH STATS
matchStats <- function(numvar, catvar, treatment, data, 
                       labelsNumvar=NULL, ordervar=NULL, groups=NULL, 
                       outXlsx=NULL, verbose=TRUE) {
  numvardataframe <- matrix(NA, ncol=8, nrow=length(numvar))
  colnames(numvardataframe) <- c("label", "", "meanX", "sdX", "meanY", "sdY",
                                 "p.value", "smd")
  for(j in 1:length(numvar)) {
    if(verbose==TRUE) cat("\n", "\n", numvar[j], "\n")
    x <- data[data[, treatment]==0, numvar[j]]
    y <- data[data[, treatment]==1, numvar[j]]
    ttest <- try(t.test(x, y, alternative="two.sided"))
    smdval <- smd(data, exposed=treatment, variable=numvar[j], categorical=F)
    
    if (class(ttest) != "try-error"){
      if(verbose==TRUE) print(ttest)
      # put in data.frame
      numvardataframe[j, 3:8] <-c(ttest$estimate[1], sd(x,na.rm=T),
                                  ttest$estimate[2], sd(x,na.rm=T),
                                  ttest$p.value, smdval)
    } else {  # if ttest is invalid
      numvardataframe[j,3:8] <- c(rep(NA,5), smdval)
    }
  }
  # manipulate data.frame
  rownames(numvardataframe) <- numvar
  if(is.null(labelsNumvar)) labelsNumvar <- numvar
  numvardataframe[,1] <- labelsNumvar
  if(verbose==TRUE) print(numvardataframe)
  
  # categorical variables
  catvarlist <- list()
  for(i in 1:length(catvar)) {
    if(verbose==TRUE) cat("\n", "\n", catvar[i], "\n")
    tab <- table(data[, catvar[i]], data[,treatment], useNA="ifany")
    pctTab <- prop.table(tab,2)
    mat <- cbind(tab[,1], pctTab[,1], tab[,2], pctTab[,2])
    colnames(mat) <- c("nX", "%X", "nY", "%Y")
    if(sum(is.na(data[,catvar[i]]))>0) rownames(mat)[nrow(mat)] <- "Unknown"
    
    # calc overall independence pvalue from chisq or fisher's test
    if(min(tab)>=5) {
      pvalue <- chisq.test(tab)$p.value
    } else pvalue <- fisher.test(tab, workspace=2e8)$p.value
    
    # calc pvalue for each level/row proportions 
    ngroups <- colSums(tab)
    smdval <- smd(data, exposed=treatment, variable=catvar[i], categorical=T)
    pval.indiv <- c()
    for(j in 1:nrow(tab)){
      pval.indiv[j] <- prop.test(x=tab[j,], n=ngroups)$p.value
    }
    catvarlist[[catvar[i]]] <- list(p.value=c(pvalue, pval.indiv),
                                    mat=as.data.frame(mat), smd=smdval)
    if(verbose==TRUE) print(catvarlist[[i]])
  }
  
  if(!is.null(outXlsx)) {
    if(is.null(groups)) groups=names(table(data[,treatment]))
    outwb <- createWorkbook()
    sheet <- createSheet(outwb, sheetName = "patientChar")
    cb <- CellBlock(sheet, 1, 1, 1000, 10)
    titleStyle <- CellStyle(outwb, font=Font(outwb, isBold=TRUE),
                            alignment=Alignment(h="ALIGN_CENTER", wrapText=T))
    CB.setRowData(cb,c("","",groups[1],"",groups[2]),1,rowStyle=titleStyle)	
    CB.setRowData(cb,c("","","N or mean","% or SD","N or mean","% or SD",
                       "P value","SMD"),2,rowStyle=titleStyle)	
    row <- 3 # min 1
    if(is.null(ordervar)) ordervar=c(numvar,catvar)
    for(k in 1:length(ordervar)) {
      if(ordervar[k] %in% numvar) {  # for numvar
        CB.setRowData(cb,numvardataframe[ordervar[k],], row, showNA=F)
        row <- row+1
      } else if(ordervar[k] %in% catvar) {   # for catvar
        
        if(nrow(catvarlist[[ordervar[k]]]$mat)>2) {  # 3 levels or more
          # fill in overall category info
          newrow <- createRow(sheet, rowIndex=row)
          cellsHeader <- createCell(newrow, colIndex=1:8)
          setCellValue(cellsHeader[[1,1]], ordervar[k]) # fill in category name
          setCellValue(cellsHeader[[1,7]], catvarlist[[ordervar[k]]]$p.value[1])
          # fill in indiv level info
          # CB.setMatrixData(cb,catvarlist[[ordervar[k]]]$mat,row,1)
          addDataFrame(cbind(catvarlist[[ordervar[k]]]$mat,
                             catvarlist[[ordervar[k]]]$p.value[-1],
                             catvarlist[[ordervar[k]]]$smd),
                       sheet,col.names=F, row.names=T, startRow=row+1,
                       startColumn=2)
          row <- row+nrow(catvarlist[[ordervar[k]]]$mat) + 1
        } else if(nrow(catvarlist[[ordervar[k]]]$mat)==2) {
          CB.setRowData(cb,c(ordervar[k] ,"", catvarlist[[ordervar[k]]]$mat[2,],
                             catvarlist[[ordervar[k]]]$p.value[1],
                             catvarlist[[ordervar[k]]]$smd[2]),row)
          row <- row+1
        }
      }
    }  # end of for loop 
    saveWorkbook(outwb, file=outXlsx)
  }
  list(numvar=numvardataframe,catvar=catvarlist)
}
