library(randomForest)

trainSetRows <- sample.int(nrow(fullFeatMatNoAcc), 
                           round((4/5)*nrow(fullFeatMatNoAcc)))
testSetRows <- which(!1:nrow(fullFeatMatNoAcc) %in% trainSetRows)
trainingSet <- fullFeatMatNoAcc[trainSetRows, ]
testSet <- fullFeatMatNoAcc[testSetRows, ]

trainStart <- Sys.time()
miTargRF <- randomForest(Target ~ ., data = trainingSet)
Sys.time() - trainStart
predictions <- predict(miTargRF, testSet)
truePos <- testSet[testSet$Target == 1 & predictions == 1, ]
trueNeg <- testSet[testSet$Target == 0 & predictions == 0, ]
falseNeg <- testSet[testSet$Target == 1 & predictions == 0, ]
falsePos <- testSet[testSet$Target == 0 & predictions == 1, ]
noTruePos <- nrow(truePos)
noFalsePos <- nrow(falsePos)
noTrueNeg <- nrow(trueNeg)
noFalseNeg <- nrow(falseNeg)
# Check that I've got numbers that make sense
testN <- nrow(testSet)
noTrueNeg + noTruePos + noFalsePos + noFalseNeg == testN
testAcuracy <- (noTruePos + noTrueNeg) / testN
specificity <- noTruePos / (noTruePos + noFalsePos)
sensitivity <- noTruePos / (noTruePos + noFalseNeg)

trainStart <- Sys.time()
miTargRF.420 <- randomForest(Target ~ ., data = trainingSet, mtry = 420)
Sys.time() - trainStart
predictions.420 <- predict(miTargRF.420, testSet)
truePos.420 <- testSet[testSet$Target == 1 & predictions.420 == 1, ]
trueNeg.420 <- testSet[testSet$Target == 0 & predictions.420 == 0, ]
falseNeg.420 <- testSet[testSet$Target == 1 & predictions.420 == 0, ]
falsePos.420 <- testSet[testSet$Target == 0 & predictions.420 == 1, ]
noTruePos.420 <- nrow(truePos.420)
noFalsePos.420 <- nrow(falsePos.420)
noTrueNeg.420 <- nrow(trueNeg.420)
noFalseNeg.420 <- nrow(falseNeg.420)

testN.420 <- nrow(testSet)
testAcuracy.420 <- (noTruePos.420 + noTrueNeg.420) / testN.420
specf.420 <- noTruePos.420 / (noTruePos.420 + noFalsePos.420)
senst.420 <- noTruePos.420 / (noTruePos.420 + noFalseNeg.420)