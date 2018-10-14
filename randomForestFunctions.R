predResultsAllTrees <- function(tree.predictions, testSet){
  aggregateRes <- predictResults(tree.predictions$aggregate, testSet)
  treeRes <- apply(tree.predictions$individual, 2, predictResults, 
                   testSet = testSet)
  toReturn <- list(aggregateRes, treeRes)
  names(toReturn) <- c("aggregate", "individual")
  return(toReturn)
}

predictResults <- function(predictions, testSet){
  posIndex <- which(testSet$Target == 1)
  truePosFalseNeg <- table(predictions[posIndex])
  trueNegFalsePos <- table(predictions[-posIndex])
  noTruePos <- truePosFalseNeg["1"]
  noTrueNeg <- trueNegFalsePos["0"]
  noFalseNeg <- truePosFalseNeg["0"]
  noFalsePos <- trueNegFalsePos["1"]
  if(is.na(noTruePos)) noTruePos <- 0
  if(is.na(noTrueNeg)) noTrueNeg <- 0
  if(is.na(noFalsePos)) noFalsePos <- 0
  if(is.na(noFalseNeg)) noFalseNeg <- 0
  testN <- nrow(testSet)
  testAcuracy <- (noTruePos + noTrueNeg) / testN
  precision <- noTruePos / (noTruePos + noFalsePos)
  recall <- noTruePos / (noTruePos + noFalseNeg)
  fscore <- (2*precision*recall) / (precision + recall)
  negPrec <- noTrueNeg / (noTrueNeg + noFalseNeg)
  specificity <- noTrueNeg / (noTrueNeg + noFalsePos)
  negFScore <- (2*negPrec*specificity) / (negPrec + specificity)
  results <- c(testAcuracy, precision, recall, fscore, negPrec, specificity, 
               negFScore, noTrueNeg, noFalseNeg, noTruePos, noFalsePos)
  results[is.nan(results)] <- 1.0
  names(results) <- c("Accuracy", "Precision", "Recall", "F-score",
                      "Negative Precision", "specificity", "Negative F-score",
                      "# True Negatives", "# False Negatives", 
                      "# True Positives", "# False Positives")
  return(results)
}

predStatistics <- function(rfTest) {
  noTruePos <- rfTest$confusion["1","1"]
  noFalsePos <- rfTest$confusion["0","1"]
  noTrueNeg <- rfTest$confusion["0","0"]
  noFalseNeg <- rfTest$confusion["1","0"]
  testN <- sum(rfTest$confusion[,1:2])
  testAcuracy <- (noTruePos + noTrueNeg) / testN
  precision <- noTruePos / (noTruePos + noFalsePos)
  recall <- noTruePos / (noTruePos + noFalseNeg)
  fscore <- (2*precision*recall) / (precision + recall)
  negPrec <- noTrueNeg / (noTrueNeg + noFalseNeg)
  specificity <- noTrueNeg / (noTrueNeg + noFalsePos)
  negFScore <- (2*negPrec*specificity) / (negPrec + specificity)
  results <- c(testAcuracy, precision, recall, fscore, negPrec, specificity, 
               negFScore, noTrueNeg, noFalseNeg, noTruePos, noFalsePos)
  names(results) <- c("Accuracy", "Precision", "Recall", "F-score",
                      "Negative Precision", "specificity", "Negative F-score",
                      "# True Negatives", "# False Negatives", 
                      "# True Positives", "# False Positives")
  return(results)
}

# balance targets and non targets by only sampling from the bigger set
# an amout equal to the smaller set 
# i.e for 300 targs and 100 non targs, include all non targs but sample
# only 100 targs
balanceFeatMat <- function(featsList, seed){
  balancedFeatsList <- featsList
  set.seed(seed)
  for(i in 1:length(featsList)){
    targCount <- table(featsList[[i]]$Target)
    nonTRows <- which(featsList[[i]]$Target == "0")
    targRows <- which(featsList[[i]]$Target == "1")
    if(targCount["0"] != 0 & targCount["1"] != 0){
      if(targCount["0"] < targCount["1"]){
        balancedFeatsList[[i]] <- 
          featsList[[i]][c(targRows[sample.int(targCount["0"])], nonTRows),]
      } else {
        balancedFeatsList[[i]] <- 
          featsList[[i]][c(nonTRows[sample.int(targCount["1"])], targRows),]
      }
    }
  }
  return(balancedFeatsList)
}

cvTrainingSet <- function(featsList, trainNo){
  trainListList <- featsList[-trainNo]
  tList <- trainListList[[1]]
  predFeats <- colnames(tList)[!colnames(tList) %in% c("GenBank.Accession")]
  for(k in 2:length(trainListList)){
    tList <- rbind(tList, trainListList[[k]])
  }
  return(tList[,predFeats])
}

# Now the 10 x 11 fold CV
# Rewritten for memory efficiency
executeCV <- function(featsList, reps = 10){
  predFeats <- colnames(featsList[[1]])
  predFeats <- predFeats[!predFeats %in% c("Target", "GenBank.Accession")]
  resultsMatrix <- matrix(nrow = reps*length(featsList), ncol = 11)
  mdaList <- list(length=length(featsList)*reps)
  runNo <- 1
  runNames <- c()
  for(i in 1:reps){
    balancedFeatsList <- balanceFeatMat(featsList, runNo)
    for(j in 1:length(balancedFeatsList)){
      runName <- paste0("Run ", i, ": ", names(balancedFeatsList[j]))
      runNames <- c(runNames, runName)
      cat(runName, "\n")
      trainSetCV <- cvTrainingSet(balancedFeatsList, j)
      testSetCV <- balancedFeatsList[[j]]
      y <- trainSetCV$Target
      x <- trainSetCV[,predFeats]
      xtest <- testSetCV[,predFeats]
      ytest <- testSetCV$Target
      set.seed(runNo)
      miRFCV <- randomForest(x = x, y = y, xtest = xtest, ytest = ytest, 
                             keep.forest = TRUE, ntree = 251)
      predResultCV <- predStatistics(miRFCV$test)
      resultsMatrix[runNo,] <- predResultCV
      mdaList[[runNo]] <- measureDA(miRFCV, testSetCV, predFeats, runNo)
      runNo <- runNo + 1
    }
  }
  rownames(resultsMatrix) <- runNames
  colnames(resultsMatrix) <- names(predResultCV)
  names(mdaList) <- runNames
  toReturn <- list(resultsMatrix, mdaList)
  names(toReturn) <- c("resultsMatrix", "mdaList")
  return(toReturn)
}

measureDA <- function(rf, testSet, predFeats, seed){
  daList <- vector("list", length = length(predFeats))
  names(daList) <- predFeats
  colNo <- 1
  for(nam in predFeats){
    permTestSet <- testSet
    set.seed(seed + colNo)
    permTestSet[,nam] <- testSet[shuffle(nrow(testSet)), nam]
    permPred <- predict(rf, permTestSet, predict.all = TRUE)
    resultsTrees <- predResultsAllTrees(permPred, permTestSet)
    if(!"0" %in% testSet$Target) rf$test$err.rate[,2] <- 1
    if(!"1" %in% testSet$Target) rf$test$err.rate[,3] <- 1
    accu <- 1 - rf$test$err.rate
    accuPerm <- t(resultsTrees$individual[c("Accuracy", "specificity", "Recall"), ])
    da <- accu - accuPerm
    daList[[nam]] <- da
    colNo <- colNo + 1
  }
  return(daList)
}

cvVariableImportance <- function(mdaList, nameTestSet, reps){
  numTestSet <- length(nameTestSet)
  predFeats <- names(mdaList[[1]])
  ntree <- nrow(mdaList[[1]][[1]])
  testSetsDAs <- sapply(1:numTestSet, function(i){
    mdaListTi <- mdaList[seq.int(i,numTestSet*reps,numTestSet)]
    featDA <- sapply(predFeats, function(nam){
      da <- sapply(mdaListTi, function(x) x[[nam]][,1])
      mDA <- mean(da)
      stdevDA <- sd(da)
      if(stdevDA == 0) mDANorm <- mDA
      else mDANorm <- mDA / stdevDA
      return(mDANorm)
    })
    return(featDA)
  })
  testSetsDNTAs <- sapply(1:numTestSet, function(i){
    mdaListTi <- mdaList[seq.int(i,numTestSet*reps,numTestSet)]
    featDNTA <- sapply(predFeats, function(nam){
      dnta <- sapply(mdaListTi, function(x) x[[nam]][,2])
      mDNTA <- mean(dnta)
      stdevDNTA <- sd(dnta)
      if(stdevDNTA == 0) mDNTANorm <- mDNTA
      else mDNTANorm <- mDNTA / stdevDNTA
      return(mDNTANorm)
    })
    return(featDNTA)
  })
  testSetsDTAs <- sapply(1:numTestSet, function(i){
    mdaListTi <- mdaList[seq.int(i,numTestSet*reps,numTestSet)]
    featDTA <- sapply(predFeats, function(nam){
      dta <- sapply(mdaListTi, function(x) x[[nam]][,3])
      mDTA <- mean(dta)
      stdevDTA <- sd(dta)
      if(stdevDTA == 0) mDTANorm <- mDTA
      else mDTANorm <- mDTA / stdevDTA
      return(mDTANorm)
    })
    return(featDTA)
  })
  colnames(testSetsDTAs) <- nameTestSet
  colnames(testSetsDAs) <- nameTestSet
  colnames(testSetsDNTAs) <- nameTestSet
  rownames(testSetsDNTAs) <- predFeats
  rownames(testSetsDTAs) <- predFeats
  rownames(testSetsDAs) <- predFeats
  return(list(testSetsDAs = testSetsDAs, testSetsDNTAs=testSetsDNTAs,
              testSetsDTAs = testSetsDTAs))
}

varImpPlot2 <- function (imp, sort = TRUE, n.var = min(30, nrow(imp)), 
          main = deparse(substitute(imp)), xlab=NULL, ...) 
{
  if (ncol(imp) > 2) 
    imp <- imp[, -(1:(ncol(imp) - 2))]
  nmeas <- ncol(imp)
  if (nmeas > 1) {
    op <- par(mfrow = c(1, 2), mar = c(4, 5, 4, 1), 
              mgp = c(2, 0.8, 0), oma = c(0, 0, 2, 0), no.readonly = TRUE)
    on.exit(par(op))
  }
  for (i in 1:nmeas) {
    ord <- if (sort) 
      rev(order(imp[, i], decreasing = TRUE)[1:n.var])
    else 1:n.var
    xmin <- if (colnames(imp)[i] %in% c("IncNodePurity", 
                                        "MeanDecreaseGini")) 
      0
    else min(imp[ord, i])
    xlab <- if (is.null(xlab)) colnames(imp)
      else xlab
    dotchart(imp[ord, i], xlab = xlab[i], ylab = "", 
             main = if (nmeas == 1) 
               main
             else NULL, xlim = c(xmin, max(imp[, i])), ...)
  }
  if (nmeas > 1) 
    mtext(outer = TRUE, side = 3, text = main, cex = 1.2)
  invisible(imp)
}


