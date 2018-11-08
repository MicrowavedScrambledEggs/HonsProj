
createNewFeatures <- function(oldFeatures){
  newFeatures <- t(apply(oldFeatures, 1, convertRow))
  newFeatNames <- c(sapply(
    1:9, function(x) paste0("c", x, "_matchfreq", 1:26)))
  newFeatNames <- c(newFeatNames, sapply(
    1:9, function(x) paste0("c", x, "_GUfreq", 1:26)
  ))
  colnames(newFeatures) <- newFeatNames
  newFeatures <- as.data.frame(newFeatures)
  for(i in 1:(9*26*2)){
    newFeatures[,i] <- as.numeric(levels(newFeatures[,i]))[newFeatures[,i]]
  }
  return(newFeatures)
}

convertRow <- function(featRow) {
  matchFeats <- c()
  guFeats <- c()
  for(i in 1:9){
    for(j in 1:26){
      nucIdenFeat <- paste0("mi", j)
      matchFreq <- 0
      guFreq <- 0
      if(featRow[nucIdenFeat] == 'A'){
        attriName <- paste0("c",i,"_freqT",j)
        matchFreq <- featRow[attriName]
      }
      if(featRow[nucIdenFeat] == 'C'){
        attriName <- paste0("c",i,"_freqG",j)
        matchFreq <- featRow[attriName]
      }
      if(featRow[nucIdenFeat] == 'G'){
        attriName <- paste0("c",i,"_freqC",j)
        matchFreq <- featRow[attriName]
        guName <- paste0("c",i,"_freqT",j)
        guFreq <- featRow[guName]
      }
      if(featRow[nucIdenFeat] == 'T'){
        attriName <- paste0("c",i,"_freqA",j)
        matchFreq <- featRow[attriName]
        guName <- paste0("c",i,"_freqG",j)
        guFreq <- featRow[guName]
      }
      matchFeats <- c(matchFeats, matchFreq)
      guFeats <- c(guFeats, guFreq)
    }
  }
  newFeats <- c(matchFeats, guFeats)
  return(newFeats)
}

# newTrainFeats <- createNewFeatures(trainingSet)
# newTrainingSet <- cbind(trainingSet, newTrainFeats)

# Keep only site identity frequency features for sites 20 to 30 and -5 to 5
badCols <- c(sapply(
  1:9, function(x) paste0("c", x, "_freq", c("A", "T", "G", "C"))))
badColNums <- c(sapply(badCols, function(x) paste0(x, "neg", 6:10)))
badColNums <- c(badColNums, sapply(badCols, function(x) paste0(x, c(6:19, 31:36))))

removeBadFeats <- function(dataSet){
  # Remove miRNA nucleotide identity features
  newDataSet <-  dataSet[,!grepl("mi\\d", colnames(dataSet))]
  # Keep only site identity frequency features for sites 20 to 30 and -5 to 5
  newDataSet <- newDataSet[,!colnames(newDataSet) %in% badColNums]
  return(newDataSet)
}

# trimTrainingSet <- removeBadFeats(newTrainingSet)
# newTestFeats <- createNewFeatures(testSet)
# newTestSet <- cbind(testSet, newTestFeats)
# trimTestSet <- removeBadFeats(newTestSet)
# new199aFeats <- createNewFeatures(featMat.199a)
# new199aSet <- cbind(featMat.199a, new199aFeats)
# trim199aSet <- removeBadFeats(new199aSet)

