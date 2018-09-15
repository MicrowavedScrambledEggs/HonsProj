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

