library(randomForest)

trainSetRows <- sample.int(nrow(fullFeatMatNoAcc), 
                           round((7/8)*nrow(fullFeatMatNoAcc)))
testSetRows <- which(!1:nrow(fullFeatMatNoAcc) %in% trainSetRows)
trainingSet <- fullFeatMatNoAcc[trainSetRows, ]
testSet <- fullFeatMatNoAcc[testSetRows, ]

miTargRF <- randomForest(Target ~ ., data = trainingSet)
