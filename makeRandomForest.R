library(randomForest)
library(permute)

trainSetRows <- sample.int(nrow(fullFeatMatNoAcc), 
                           round((4/5)*nrow(fullFeatMatNoAcc)))
testSetRows <- which(!1:nrow(fullFeatMatNoAcc) %in% trainSetRows)
trainingSet <- fullFeatMatNoAcc[trainSetRows, ]
testSet <- fullFeatMatNoAcc[testSetRows, ]

predictResults <- function(predictions, testSet){
  truePos <- testSet[testSet$Target == 1 & predictions == 1, ]
  trueNeg <- testSet[testSet$Target == 0 & predictions == 0, ]
  falseNeg <- testSet[testSet$Target == 1 & predictions == 0, ]
  falsePos <- testSet[testSet$Target == 0 & predictions == 1, ]
  noTruePos <- nrow(truePos)
  noFalsePos <- nrow(falsePos)
  noTrueNeg <- nrow(trueNeg)
  noFalseNeg <- nrow(falseNeg)
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
  names(results) <- c("Accuracy", "Precision", "Recall", "F-score",
                      "Negative Precision", "specificity", "Negative F-score",
                      "# True Negatives", "# False Negatives", 
                      "# True Positives", "# False Positives")
  return(results)
}

trainStart <- Sys.time()
miTargRF <- randomForest(Target ~ ., data = trainingSet)
Sys.time() - trainStart
predictions <- predict(miTargRF, testSet)
predStats <- predictResults(predictions, testSet)

trainStart <- Sys.time()
miTargRF.420 <- randomForest(Target ~ ., data = trainingSet, mtry = 420)
Sys.time() - trainStart
predictions.420 <- predict(miTargRF.420, testSet)
predStats.420 <- predictResults(predictions.420, testSet)

# Test on miRNA not used in the training set
predict.199a <- predict(miTargRF, featMat.199a)
results.199a <- predictResults(predict.199a, featMat.199a)
predict.199a.420 <- predict(miTargRF.420, featMat.199a)
results.199a.420 <- predictResults(predict.199a.420, featMat.199a)

# Retrain with feature set more general to all miRNA:mRNA pairs
trainStart <- Sys.time()
miTargRF.trimFeat <- randomForest(Target ~ ., data = trimTrainingSet)
Sys.time() - trainStart
predTrim <- predict(miTargRF.trimFeat, trimTestSet)
trimResults <- predictResults(predTrim, trimTestSet)
predTrim.199a <- predict(miTargRF.trimFeat, trim199aSet)
trimResults.199a <- predictResults(predTrim.199a, trim199aSet)

# 10 X 11 fold CV divided by miRNA

# Full feature sets with all miRNA:RNA pairs
featsList <- list(featMat.1289, featMat.155, featMat.199a, featMat.424, 
                  featMat.4536, featMat.548d, featMat.584, featMatmir23b, 
                  featMatmir27a, featMatmir3118, featMatmir4307)
names(featsList) <- c("miR-1289", "miR-155", "miR-199a", "miR-424",
                      "miR-4536", "miR-548d", "miR-584", "miR-23b",
                      "miR-27a", "miR-3118", "miR-4307")

# Trail test and training sets with just 5 pairs for each miRNA
# and only 5 features
targCol <- which(colnames(featsList[[1]]) == "Target")
sampleCols <- sample.int(ncol(featsList[[1]]), 5)
sampleCols <- sampleCols[sampleCols != targCol]
sampleCols <- c(sampleCols, targCol)
trialFeatsList <- lapply(
  featsList, function(x) x[sample.int(nrow(x), 5),sampleCols])

# Trim feature sets with balanced targets and non targets

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

# Trim
trimFeatsList <- featsList
for(i in 1:length(trimFeatsList)){
  newFeats <- createNewFeatures(trimFeatsList[[i]])
  combFeats <- cbind(trimFeatsList[[i]], newFeats)
  trimFeatsList[[i]] <- removeBadFeats(combFeats)
}

# Now the training sets
cvTrainingSet <- function(featsList, trainNo){
  trainListList <- featsList[-trainNo]
  tList <- trainListList[[1]]
  for(k in 2:length(trainListList)){
    tList <- rbind(tList, trainListList[[k]])
  }
  return(tList[,-1])
}

# Now the 10 x 11 fold CV
# Rewritten for memory efficiency
executeCV <- function(featsList){
  predFeats <- colnames(featsList[[1]])
  predFeats <- predFeats[!predFeats %in% c("Target", "GenBank.Accession")]
  resultsMatrix <- matrix(nrow = 10*length(featsList), ncol = 11)
  mdaMatrix <- matrix(ncol = 10*length(featsList), nrow = length(predFeats))
  mdfMatrix <- matrix(ncol = 10*length(featsList), nrow = length(predFeats))
  mdnfMatrix <- matrix(ncol = 10*length(featsList), nrow = length(predFeats))
  rownames(mdaMatrix) <- predFeats
  rownames(mdfMatrix) <- predFeats
  rownames(mdnfMatrix) <- predFeats
  runNo <- 1
  runNames <- c()
  for(i in 1:10){
    balancedFeatsList <- balanceFeatMat(featsList, runNo)
    for(j in 1:length(balancedFeatsList)){
      runName <- paste0("Run ", i, ": ", names(balancedFeatsList[j]))
      runNames <- c(runNames, runName)
      cat(runName, "\n")
      trainSetCV <- cvTrainingSet(balancedFeatsList, j)
      testSetCV <- balancedFeatsList[[j]]
      set.seed(runNo)
      miRFCV <- randomForest(Target ~ ., data = trainSetCV)
      predCV <- predict(miRFCV, testSetCV)
      predResultCV <- predictResults(predCV, testSetCV)
      resultsMatrix[runNo,] <- predResultCV
      das <- measureDA(miRFCV, testSetCV, predResultCV, runNo)
      mdaMatrix[,runNo] <- das[1,]
      mdfMatrix[,runNo] <- das[2,]
      mdnfMatrix[,runNo] <- das[3,]
      runNo <- runNo + 1
    }
  }
  rownames(resultsMatrix) <- runNames
  colnames(resultsMatrix) <- names(predResultCV)
  colnames(mdaMatrix) <- runNames
  colnames(mdfMatrix) <- runNames
  colnames(mdnfMatrix) <- runNames
  toReturn <- list(resultsMatrix, mdaMatrix, mdfMatrix, mdnfMatrix)
  names(toReturn) <- c("resultsMatrix", "mdaMatrix", "mdfMatrix", "mdnfMatrix")
  return(list(resultsMatrix, mdaMatrix, mdfMatrix, mdnfMatrix))
}

measureDA <- function(rf, testSet, results, seed){
  targCol <- which(colnames(testSet) == "Target")
  accCol <- which(colnames(testSet) == "GenBank.Accession")
  exclude <- c(accCol, targCol)
  daMat <- matrix(nrow = 3, ncol = length(colnames(testSet)[-exclude]))
  colNo <- 1
  for(nam in colnames(testSet)[-exclude]){
    permTestSet <- testSet
    set.seed(seed)
    permTestSet[,nam] <- testSet[shuffle(nrow(testSet)), nam]
    permPred <- predict(rf, permTestSet)
    permResults <- predictResults(permPred, permTestSet)
    da <- results["Accuracy"] - permResults["Accuracy"]
    df <- results["F-score"] - permResults["F-score"]
    dnf <- results["Negative F-score"] - permResults["Negative F-score"]
    daMat[,colNo] <- c(da, df, dnf)
    colNo <- colNo + 1
  }
  return(daMat)
}

cv.10.11 <- executeCV(trimFeatsList) 


# Random Forest trained on just miR-548d data
new548d <- createNewFeatures(featMat.548d)
comb548d <- cbind(featMat.548d, new548d)
trim548d <- removeBadFeats(comb548d)

set.seed(456)
trainRows <- sample.int(nrow(trim548d), floor((4 * nrow(trim548d))/5))
train548d <- trim548d[trainRows,-1]
table(train548d$Target)
test548d <- trim548d[-trainRows,-1]
table(test548d$Target)
set.seed(567)
rf.548d <- randomForest(Target~., data = train548d, importance = TRUE)
pred.548d <- predict(rf.548d, test548d)
results.548d <- predictResults(pred.548d, test548d)
imp548d <- importance(rf.548d, type = 1)
