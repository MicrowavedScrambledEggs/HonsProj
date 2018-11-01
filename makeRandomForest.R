library(randomForest)
library(permute)
source("randomForestFunctions.R")

trainSetRows <- sample.int(nrow(fullFeatMatNoAcc), 
                           round((4/5)*nrow(fullFeatMatNoAcc)))
testSetRows <- which(!1:nrow(fullFeatMatNoAcc) %in% trainSetRows)
trainingSet <- fullFeatMatNoAcc[trainSetRows, ]
testSet <- fullFeatMatNoAcc[testSetRows, ]



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

# 10 X 10 fold CV divided by miRNA

# Full feature sets with all miRNA:RNA pairs
featsList <- list(featMat.1289, featMat.155, featMat.199a, featMat.424, 
                  featMat.4536, featMat.548d, featMatmir23b, featMat.10a,
                  featMat.182, featMat.30e)
names(featsList) <- c("miR-1289", "miR-155", "miR-199a", "miR-424",
                      "miR-4536", "miR-548d", "miR-23b",
                      "miR-10a", "miR-182", "miR-30e")

# Trail test and training sets with just 5 pairs for each miRNA
# and only 5 features
targCol <- which(colnames(featsList[[1]]) == "Target")
sampleCols <- sample.int(ncol(featsList[[1]]), 5)
sampleCols <- sampleCols[sampleCols != targCol]
sampleCols <- c(sampleCols, targCol)
trialFeatsList <- lapply(
  featsList, function(x) x[sample.int(nrow(x), 5),sampleCols])

# Trim feature sets with balanced targets and non targets
# Trim
trimFeatsList <- featsList
for(i in 1:length(trimFeatsList)){
  newFeats <- createNewFeatures(trimFeatsList[[i]])
  combFeats <- cbind(trimFeatsList[[i]], newFeats)
  trimFeatsList[[i]] <- removeBadFeats(combFeats)
}

cv.10.11 <- executeCV(trimFeatsList, reps=10) 
# Above CV was done when 584 was still part of the feats list,
# stricter definition of sufficiently expressed had not been applied,
# and extra 548d and 1289 data hadn't been added yet
cv.10.10 <- executeCV(trimFeatsList, reps=10)
cv.10.10.balance <- executeCV(trimFeatsList, reps=10)
cv.10.10.fullTest <- executeCV(trimFeatsList, testList = trimFeatsList, reps = 10)


# Looks like miR-548d pairs just get all classed as non-targets
# Could be:
#    - miR-548d data is crap
#    - Every other dataset is crap
#    - miR-548d is that much of a special snowflake 
#        - no GU wobble in seed

# Random Forest trained on just miR-548d data
trim548d <- trimFeatsList$`miR-548d`

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

# CV without miR-548d data
cv.wout.548d <- executeCV(trimFeatsList[-6], 1)

# Training and testing on 548d gets 70% accuracy, 
# CV perfomance appears unaffected by the removal of 548d, c4 features are 
# important in all the CV rfs, but 548d does not have G or U in it's seed.
# From this I conclude that the other rfs have learned to classify all miRNA
# with no c4 seed sites as non-targets.

# Test if we got the correct 30e sequence
featMat.30e <- convertMissingValues(featMat.30e)
new.30e <- createNewFeatures(featMat.30e)
comb.30e <- cbind(featMat.30e, new.30e)
trim.30e <- removeBadFeats(comb.30e)

set.seed(678)
train.30e.row <- sample.int(nrow(trim.30e), round((4/5)*nrow(trim.30e)))
train.30e <- trim.30e[train.30e.row, -1]
test.30e <- trim.30e[-train.30e.row, ]

rf.30e <- randomForest(Target~., data = train.30e, importance = TRUE)
pred.30e <- predict(rf.30e, test.30e)
results.30e <- predictResults(pred.30e, test.30e)

# now try 3p
miR30e.3p <- "CUUUCAGUCGGAUGUUUACAGC"
s.rS.fM.30e.3p <- createFeatureMatrix(miR30e.3p, accToTarg.30e[,1], 
                                   accToTarg.30e[,2])
freeEngFeats.30e.3p <- siteDuplexFreeEnergy(miR30e.3p, s.rS.fM.30e.3p[[2]],
                                         s.rS.fM.30e.3p[[1]])
freeEngFeats.30e.3p$GenBank.Accession <- names(s.rS.fM.30e.3p[[2]])
featMat.30e.3p <- merge(s.rS.fM.30e.3p[[3]], freeEngFeats.30e.3p, by = "GenBank.Accession")
featMat.30e.3p <- convertMissingValues(featMat.30e.3p)

new.30e.3p <- createNewFeatures(featMat.30e.3p)
comb.30e.3p <- cbind(featMat.30e.3p, new.30e.3p)
trim.30e.3p <- removeBadFeats(comb.30e.3p)

set.seed(678)
train.30e.3p.row <- sample.int(nrow(trim.30e.3p), round((4/5)*nrow(trim.30e.3p)))
train.30e.3p <- trim.30e.3p[train.30e.3p.row, -1]
test.30e.3p <- trim.30e.3p[-train.30e.3p.row, ]

rf.30e.3p <- randomForest(Target~., data = train.30e.3p, importance = TRUE)
pred.30e.3p <- predict(rf.30e.3p, test.30e.3p)
results.30e.3p <- predictResults(pred.30e.3p, test.30e.3p)

# Find good ntree
trim155 <- trimFeatsList$`miR-155`

set.seed(456)
trainRows <- sample.int(nrow(trim155), floor((4 * nrow(trim155))/5))
train155 <- trim155[trainRows,-1]
table(train155$Target)
test155 <- trim155[-trainRows,-1]
table(test155$Target)
set.seed(567)
targCol <- which(colnames(test155) == "Target")
rf.155.501 <- randomForest(x=train155[,-targCol], y=train155$Target, 
                       xtest = test155[,-targCol], ytest = test155$Target,
                       ntree = 501)
rf.155.501
rf.155.1001 <- randomForest(x=train155[,-targCol], y=train155$Target, 
                           xtest = test155[,-targCol], ytest = test155$Target,
                           ntree = 1001)
rf.155.1001
rf.155.251 <- randomForest(x=train155[,-targCol], y=train155$Target, 
                           xtest = test155[,-targCol], ytest = test155$Target,
                           ntree = 251)
rf.155.251

# Train and test an RF sampling from all the data
set.seed(654)
balancedFeatsList <- balanceFeatMat(trimFeatsList, 654)
train.all <- data.frame()
test.all <- data.frame()
for(miName in names(balancedFeatsList)){
  trimSet <- balancedFeatsList[[miName]]
  trainRows <- sample.int(nrow(trimSet), floor((4 * nrow(trimSet))/5))
  trainSet <- trimSet[trainRows,-1]
  testSet <- trimSet[-trainRows,-1]
  train.all <- rbind(train.all, trainSet)
  test.all <- rbind(test.all, testSet)
}

targCol <- which(colnames(test.all) == "Target")
all.rf <- randomForest(
  x = train.all[,-targCol], y = train.all$Target, xtest = test.all[,-targCol],
  ytest = test.all$Target, ntree = 251, importance = TRUE)
predStatistics(all.rf$test)
varImpPlot(all.rf)
varImpPlot2(importance(all.rf)[,1:2], xlab = c("Mean Decrease in Specificity", "Mean Decrease in Sensitivity"), 
            main = "Top 30 Important Features as Measured by Specificity and Sensitivity")

new10a <- createNewFeatures(featMat.10a)
comb10a <- cbind(featMat.10a, new10a)
trim10a <- removeBadFeats(comb10a)

set.seed(456)
trainRows <- sample.int(nrow(trim10a), floor((4 * nrow(trim10a))/5))
train10a <- trim10a[trainRows,-1]
table(train10a$Target)
test10a <- trim10a[-trainRows,-1]
table(test10a$Target)
set.seed(567)
rf.10a.251 <- randomForest(x=train10a[,-targCol], y=train10a$Target, 
                           xtest = test10a[,-targCol], ytest = test10a$Target,
                           ntree = 251)

# Results tables for balanced CV
getStats <- function(i, cvResults){
  meanAcc <- mean(cvResults$resultsMatrix[seq(i,100,10),1])
  sdAcc <- sd(cvResults$resultsMatrix[seq(i,100,10),1])
  meanSens <- mean(cvResults$resultsMatrix[seq(i,100,10),3])
  meanSpec <- mean(cvResults$resultsMatrix[seq(i,100,10),6])
  return(c("Mean Accuracy"=meanAcc, "Std Dev Accuracy"=sdAcc, "Mean Specificity"=meanSpec,
           "Mean Sensitivity"=meanSens))
}

# Feature importance

# rbind the tables for each rep for each miRNA for each feature together
featImpTabs <- sapply(names(featsList), function(miNam) {
  i <- which(names(featsList) == miNam)
  featTabLists <- cv.10.10.balance$mdaList[seq(i,100,10)]
  featImpTab <- sapply(names(featTabLists[[1]]), function(featNam){
    featTab <- do.call("rbind", lapply(featTabLists, function(x) x[[featNam]]))
    col.sd <- apply(featTab, 2, sd)
    col.sd <- col.sd / sqrt(210)
    return(colMeans(featTab) / col.sd)
  })
  return(list(t(featImpTab)))
})

featImpTabsUnorm <- sapply(names(featsList), function(miNam) {
  i <- which(names(featsList) == miNam)
  featTabLists <- cv.10.10.balance$mdaList[seq(i,100,10)]
  featImpTab <- sapply(names(featTabLists[[1]]), function(featNam){
    featTab <- do.call("rbind", lapply(featTabLists, function(x) x[[featNam]]))
    return(colMeans(featTab))
  })
  return(list(t(featImpTab)))
})

featImpTabsSTDER <- sapply(names(featsList), function(miNam) {
  i <- which(names(featsList) == miNam)
  featTabLists <- cv.10.10.balance$mdaList[seq(i,100,10)]
  featImpTab <- sapply(names(featTabLists[[1]]), function(featNam){
    featTab <- do.call("rbind", lapply(featTabLists, function(x) x[[featNam]]))
    col.sd <- apply(featTab, 2, sd)
    col.sd <- col.sd / sqrt(210)
    return(col.sd)
  })
  return(list(t(featImpTab)))
})

# Table for each miRNA that have a row for each feature and a column each for mda, mdspec, mdsens

# plotting feat importance
all.imp <- importance(all.rf)
all.imp <- cbind(all.imp[,3], all.imp[,1:2])
varImpPlot3(all.imp, featImpTabs, xlab = c("Mean Decrease in Accuracy", "Mean Decrease in Specificity", 
                                           "Mean Decrease in Sensitivity"))

# Checking that my way of calulating mda is the same as the packages way
all.rf.wfb <- randomForest(
  x = train.all[,-targCol], y = train.all$Target, xtest = test.all[,-targCol],
  ytest = test.all$Target, ntree = 251, importance = TRUE, keep.forest = TRUE, keep.inbag = T)

names(testDAOOB) <- names(all.imp)

featImpTabOOB <- sapply(names(testDAOOB), function(featNam){
  featTab <- testDAOOB[[featNam]]
  col.sd <- apply(featTab, 2, sd)
  col.sd <- col.sd / sqrt((1/3)*1680)
  return(colMeans(featTab) / col.sd)
})

# Data best suited for this feature set:

bestPerf <- list(featMatmir23b, featMatmir27a, featMatmir3118)
fullBestMat <- do.call("rbind", bestPerf)
fullBestMatTrim <- getTrimFeatSet(fullBestMat)
best.size <- nrow(fullBestMatTrim)
table(fullBestMatTrim$Target)
best.targs <- which(fullBestMatTrim$Target == 1)
best.nonTargs <- (1:best.size)[-best.targs]
set.seed(746)
train.targs <- sample(best.targs, (5/6)*length(best.targs))
train.nonTargs <- sample(best.nonTargs, (5/6)*length(best.nonTargs))
best.train <- fullBestMatTrim[c(train.targs, train.nonTargs),-1]
best.test <- fullBestMatTrim[-c(train.targs, train.nonTargs),-1]

best.rf <- randomForest(
  x = best.train[,-targCol], y = best.train$Target, xtest = best.test[,-targCol],
  ytest = best.test$Target, ntree = 251, importance = TRUE)

# Finding the optimal target vote majority percent
all.rf.60 <- randomForest(
  x = train.all[,-targCol], y = train.all$Target, xtest = test.all[,-targCol],
  ytest = test.all$Target, ntree = 251, cutoff = c(0.4,0.6))

all.rf.55 <- randomForest(
  x = train.all[,-targCol], y = train.all$Target, xtest = test.all[,-targCol],
  ytest = test.all$Target, ntree = 251, cutoff = c(0.45,0.55))

train.1289.balanced <- do.call("rbind",balancedFeatsList[-1])

rf.1289.60 <- randomForest(
  x = train.1289.balanced[,-c(1,targCol+1)], y = train.1289.balanced$Target, xtest = trimFeatsList$`miR-1289`[,-c(1,targCol+1)],
  ytest = trimFeatsList$`miR-1289`$Target, ntree = 251, cutoff = c(0.4,0.6))

rf.1289 <- randomForest(
  x = train.1289.balanced[,-c(1,targCol+1)], y = train.1289.balanced$Target, xtest = trimFeatsList$`miR-1289`[,-c(1,targCol+1)],
  ytest = trimFeatsList$`miR-1289`$Target, ntree = 251)

rf.1289.70 <- randomForest(
  x = train.1289.balanced[,-c(1,targCol+1)], y = train.1289.balanced$Target, xtest = trimFeatsList$`miR-1289`[,-c(1,targCol+1)],
  ytest = trimFeatsList$`miR-1289`$Target, ntree = 251, cutoff = c(0.3,0.7))

rf.1289.65 <- randomForest(
  x = train.1289.balanced[,-c(1,targCol+1)], y = train.1289.balanced$Target, xtest = trimFeatsList$`miR-1289`[,-c(1,targCol+1)],
  ytest = trimFeatsList$`miR-1289`$Target, ntree = 251, cutoff = c(0.35,0.65))

rf.1289.amb.8 <- randomForest(
  x = train.1289.balanced[,-c(1,targCol+1)], y = train.1289.balanced$Target, xtest = trimFeatsList$`miR-1289`[,-c(1,targCol+1)],
  ytest = trimFeatsList$`miR-1289`$Target, ntree = 251, cutoff = c(0.8,0.8))

fullTestStats <- t(sapply(1:10, getStats, cvResults = cv.10.10.fullTest))
rownames(fullTestStats) <- names(trimFeatsList)
fullTestStats <- round(fullTestStats, 3)
ntargsNonTargs <- t(sapply(trimFeatsList, function(x) table(x$Target)))
fullTestStats <- cbind(fullTestStats, ntargsNonTargs)

fullTestFullAcc <- mean(cv.10.10.fullTest$resultsMatrix[,1])
fullTestFullSD <- sd(cv.10.10.fullTest$resultsMatrix[,1])
fullTestFullSens <- mean(cv.10.10.fullTest$resultsMatrix[,3])
fullTestFullSpec <- mean(cv.10.10.fullTest$resultsMatrix[,6])


