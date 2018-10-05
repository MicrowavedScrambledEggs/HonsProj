library(randomForest)

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
  specificity <- noTruePos / (noTruePos + noFalsePos)
  sensitivity <- noTruePos / (noTruePos + noFalseNeg)
  results <- c(testAcuracy, specificity, sensitivity, noTrueNeg, noFalseNeg,
               noTruePos, noFalsePos)
  names(results) <- c("Accuracy", "Specificity", "Sensitivity", 
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
trainList <- list(length = length(featsList))
for(j in 1:length(featsList)){
  trainListList <- featsList[-j]
  tList <- trainListList[[1]]
  for(k in 2:length(trainListList)){
    tList <- rbind(tList, trainListList[[k]])
  }
  trainList[[j]] <- tList[,-1]
}

# Trail test and training sets with just one pair for each miRNA
trialFeatsList <- lapply(featsList, function(x) x[1,])
trialTrainList <- list(length = length(trialFeatsList))
for(j in 1:length(trialFeatsList)){
  trainListList <- trialFeatsList[-j]
  tList <- trainListList[[1]]
  for(k in 2:length(trainListList)){
    tList <- rbind(tList, trainListList[[k]])
  }
  trialTrainList[[j]] <- tList[,-1]
}

# Trial results
resultsMatrix <- matrix(nrow = 10*length(featsList), ncol = 7)
cvForests <- list(length = 10*length(featsList))
runNo <- 1
runNames <- c()
for(i in 1:10){
  for(j in 1:length(trialFeatsList)){
    runName <- paste0("Run ", i, ": ", names(trialFeatsList[j]))
    runNames <- c(runNames, runName)
    cat(runName, "\n")
    trainSetCV <- trialTrainList[[j]]
    testSetCV <- trialFeatsList[[j]]
    miRFCV <- randomForest(Target ~ ., data = trainSetCV)
    predCV <- predict(miRFCV, testSetCV)
    predResultCV <- predictResults(predCV, testSetCV)
    resultsMatrix[runNo,] <- predResultCV
    cvForests[[runNo]] <- miRFCV
    runNo <- runNo + 1
  }
}
rownames(resultsMatrix) <- runNames
colnames(resultsMatrix) <- names(predResultCV)
names(cvForests) <- runNames

# Trim feature sets with balanced targets and non targets
# balance targets and non targets
balancedFeatsList <- featsList
set.seed(234)
for(i in 1:length(trimFeatsList)){
  targCount <- table(featsList[[i]]$Target)
  nonTRows <- which(featsList[[i]]$Target == "0")
  targRows <- which(featsList[[i]]$Target == "1")
  if(targCount["0"] < targCount["1"]){
    balancedFeatsList[[i]] <- 
      featsList[[i]][c(targRows[sample.int(targCount["0"])], nonTRows),]
  } else {
    balancedFeatsList[[i]] <- 
      featsList[[i]][c(nonTRows[sample.int(targCount["1"])], targRows),]
  }
}

# Trim
trimFeatsList <- balancedFeatsList
for(i in 1:length(trimFeatsList)){
  newFeats <- createNewFeatures(trimFeatsList[[i]])
  combFeats <- cbind(trimFeatsList[[i]], newFeats)
  trimFeatsList[[i]] <- removeBadFeats(combFeats)
}

# Now the training sets
trimTrainList <- list(length = length(trimFeatsList))
for(j in 1:length(trimFeatsList)){
  trainListList <- trimFeatsList[-j]
  tList <- trainListList[[1]]
  for(k in 2:length(trainListList)){
    tList <- rbind(tList, trainListList[[k]])
  }
  trimTrainList[[j]] <- tList[,-1]
}

# Now the 10 x 11 fold CV
resultsMatrix <- matrix(nrow = 10*length(trimFeatsList), ncol = 7)
cvForests <- list(length = 10*length(trimFeatsList))
runNo <- 1
runNames <- c()
set.seed(345)
for(i in 1:10){
  for(j in 1:length(trimFeatsList)){
    runName <- paste0("Run ", i, ": ", names(trimFeatsList[j]))
    runNames <- c(runNames, runName)
    cat(runName, "\n")
    trainSetCV <- trimTrainList[[j]]
    testSetCV <- trimFeatsList[[j]]
    miRFCV <- randomForest(Target ~ ., data = trainSetCV)
    predCV <- predict(miRFCV, testSetCV)
    predResultCV <- predictResults(predCV, testSetCV)
    resultsMatrix[runNo,] <- predResultCV
    cvForests[[runNo]] <- miRFCV
    runNo <- runNo + 1
  }
}
rownames(resultsMatrix) <- runNames
colnames(resultsMatrix) <- names(predResultCV)
names(cvForests) <- runNames