library(rentrez)
library(Biostrings)
library(dequer)

# Size of window around position of end of sites in which to measure 
# nucleotide frequencies at each position
dwnStrm <- 35
upStrm <- 10
nucWin <- dwnStrm + upStrm + 1

# Directory where the ViennaRNA executables are. Edit when running on
# windows systems where you can not add the directory to PATH
vRNADir <- ""
if(grepl("H:/", getwd())) vRNADir <- '"C:/Users/bjam575/AppData/Local/ViennaRNA Package/'

# Read in already downloaded sequences and annotations so we don't have
# to query dbs every time
# accSeqRegions <- read.csv("accSeqRegions.csv")


querymRNASequence <- function(accNos){
  # Querry nucleotide
  # Parse the sequences out of the FASTA
  # Parse the codon start and stop out of the feature tables
  
  output <- as.data.frame(matrix(nrow = 0, ncol = 4))
  colnames(output) <- c("GenBank Accession", "Full Sequence", "CDS Start",
                        "CDS Stop")
  
  # Because we're usually dealling with 1e4 to 1.5e4 RNAs we need split
  # the request into chunks
  chunkSize = 100
  for(i in seq(1, length(accNos), chunkSize))
  {
    accChunk <- accNos[i:min(i+chunkSize-1, length(accNos))]
    fasta <- entrez_fetch("nucleotide", id=accChunk, rettype = "fasta")
    seqFastas <- strsplit(fasta, "\n\n")[[1]]
    seqFastaList <- strsplit(seqFastas, "\n")
    seqsStrings <- sapply(seqFastaList, function(x) 
      paste0(x[2:length(x)], collapse = ""))
    featTabs <- entrez_fetch("nucleotide", id=accChunk, rettype = "ft")
    featTabsList <- strsplit(featTabs, "\n\n")[[1]]
    hasCDS <- regexpr("\\d+\\t\\d+\\tCDS", featTabsList)
    cdsStartStops <- regmatches(featTabsList, hasCDS)
    cdsList <- strsplit(cdsStartStops, "\t")
    cdsStarts <- sapply(cdsList, function(x) as.integer(x[1]))
    cdsStops <- sapply(cdsList, function(x) as.integer(x[2]))
    chunk <- as.data.frame(cbind(accChunk, seqsStrings))
    colnames(chunk) <- c("GenBank Accession", "Full Sequence")
    chunk$CDS.Start <- rep(NaN, length(accChunk))
    chunk$CDS.Start[hasCDS != -1] <- cdsStarts
    chunk$CDS.Stop <- rep(NaN, length(accChunk))
    chunk$CDS.Stop[hasCDS != -1] <- cdsStops
    output <- rbind(output, chunk)
  }
  return(output)
}


guWob_seqs <- function(seq){
  # Outputs a DNAStringSet that represent every possible site that
  # involves at least one GU complement (with everywhere else perfectly 
  # WC complementing).
  # seq: sequence of the perfect WC complement. The GU complements are bulit
  # from this
  output <- c()
  seqQ <- queue()
  pushback(seqQ, list(seq, 1))
  while(length(seqQ) > 0){
    sTup <- pop(seqQ)
    sCurrent <- sTup[[1]]
    startPos <- sTup[[2]]
    for(i in c(startPos:length(sCurrent))){
      if(letter(sCurrent, i) == "C"){
        snew <- replaceAt(sCurrent, IRanges(i,i), "T")
        if(i < length(sCurrent)) pushback(seqQ, list(snew, i+1))
        output <- c(output, snew)
      }
      if(letter(sCurrent, i) == "A"){
        snew <- replaceAt(sCurrent, IRanges(i,i), "G")
        if(i < length(sCurrent)) pushback(seqQ, list(snew, i+1))
        output <- c(output, snew)
      }
    }
  }
  return(DNAStringSet(output))
}


misMatchNotGU <- function(seq, start = 1, end = NULL, reference = NULL){
  # Outputs a DNAStringSet representing every site with all bp wc
  # complementing one of the DNAStrings in seq except for one mismatch
  # between positions start and end.
  # The introduced mismatch does not allow for nucleotides that would 
  # make the position a GU complement (as that is considered a different
  # site class).
  # If seq represents sites with GU complements, reference should be a 
  # DNAString of the perfect wc complement. This is used to ensure that
  # 1: for DNAString in seq with only one GU complement the mismatch is
  #    not put in the same place
  # 2: for DNAString in seq with only multiple GU complements the mismatch 
  #    at a GU complement doesn't allow for a wc complement there 
  if(class(seq) != "DNAStringSet"){
    seq <- DNAStringSet(seq) # allows inputting a single DNAString
  }
  output <- c()
  
  # IUPAC letter for "any nucleotide that wont make this a wc or gu pair"
  mismatchMap <- c("Y", "R", "H", "V")
  names(mismatchMap) <- c("A", "C", "G", "T")
  
  # for when seq is of sites with GU pairing
  noGUWC <- c("Y", "R")
  names(noGUWC) <- c("G", "T")
  
  theEndWasNull <- is.null(end) 
  for(j in 1:length(seq)){ 
    scurrent <- seq[[j]]
    if(theEndWasNull) end <- length(scurrent)
    badi <- -1 # site of only GU pair
    alreadyGU <- c()
    if(!is.null(reference)){
      # find the gu pairs
      scomp <- compareStrings(scurrent, reference)
      diffs <- gregexpr("\\?", scomp)[[1]]
      if(length(diffs) == 1) badi <- diffs[1]
      else if(length(diffs) > 1) alreadyGU <- diffs
    }
    for(i in start:end){
      if(i != badi){
        if(i %in% alreadyGU){
          snew <- replaceAt(scurrent, IRanges(i,i), 
                            noGUWC[letter(scurrent, i)])
        } else {
          snew <- replaceAt(scurrent, IRanges(i,i), 
                            mismatchMap[letter(scurrent, i)])
        }
        output <- c(output, snew)
      }
    }
  }
  return(unique(DNAStringSet(output)))
}


nucFreqAroundSites <- function(dnastring, siteEnds)
{
  # Measure the frequency of nucleotides at each position in windows
  # centered on the sites 
  
  # Get the sequence for each window 
  siteAreaSeqs <- DNAStringSet(sapply(
    siteEnds, function(x) 
      subseq(dnastring, start = max(x - dwnStrm, 1), 
             end = min(x + upStrm, length(dnastring)))))
  
  # Because sites may be at the ends of the mRNA, so the windows may
  # extend beyond the bounds of the mRNA sequence. This means sites might
  # vary in length. Therefore we have to pad out window sequences such 
  # that the index of the position opposite miRNA position 1 is the same
  # accross all window sequences
  siteMat <- as.data.frame(
    cbind(sapply(siteAreaSeqs, as.character), siteEnds))
  siteMat$siteEnds <- siteEnds
  siteMat$fiveOff <- 0 - (siteEnds - dwnStrm - 1)
  siteMat$fiveOff[siteMat$fiveOff < 0] <- 0
  siteMat$threeOff <- (siteEnds + upStrm) - nchar(dnastring)
  siteMat$threeOff[siteMat$threeOff < 0] <- 0
  shiftSeq <- apply(siteMat, 1, function(x) 
    paste0(strrep("-", x[3]), x[1], strrep("-", x[4])))
  shiftSeq <- DNAStringSet(shiftSeq)
  
  # Measure frequencies
  nucFreqs <- c(sapply(1:nucWin, function(x) 
    nucleotideFrequencyAt(shiftSeq, c(x))))
  nucFreqs <- nucFreqs / length(siteEnds)
  return(nucFreqs)
}


sitesOfClass <- function(classSearch, nameOfClass, rnaSet, rnaCDS, 
                         centeredClass = NULL) 
{
  # TODO: Document this disgustingly long method
  if(class(classSearch) != "DNAStringSet"){
    classSearch <- DNAStringSet(classSearch) # allows inputting a single DNAString
  }
  classMatches <- sapply(classSearch, vmatchPattern, subject = rnaSet,
                         fixed = "subject")
  classEnds <- sapply(classMatches, endIndex)
  classCounts <- sapply(classMatches, elementNROWS)
  allClassEnds <- vector('list', length(rnaSet))
  allClassCounts <- rep(0, length(rnaSet))
  meanDistFrom3 <- rep(NaN, length(rnaSet))
  varDistFrom3 <- rep(NaN, length(rnaSet))
  medianDistFrom3 <- rep(NaN, length(rnaSet))
  meanDistFromStop <- rep(NaN, length(rnaSet))
  varDistFromStop <- rep(NaN, length(rnaSet))
  medianDistFromStop <- rep(NaN, length(rnaSet))
  freqIn3 <- rep(NaN, length(rnaSet))
  freqInCDS <- rep(NaN, length(rnaSet))
  freqIn5 <- rep(NaN, length(rnaSet))
  
  if(!is.null(centeredClass)){
    freqStart3 <- rep(NaN, length(rnaSet))
    freqStart4 <- rep(NaN, length(rnaSet))
    freqStart5 <- rep(NaN, length(rnaSet))
  }
  
  # for counting the nucleotide frequencies in and around the sites
  siteIndex <- seq(dwnStrm+1, -upStrm, -1)
  siteIndex <- siteIndex[siteIndex != 0]
  nucFreqMat <- matrix(NaN, nrow = length(rnaSet), ncol = nucWin*4)
  
  for(i in 1:length(rnaSet)){
    classEnd <- c()
    classCount <- 0
    for(j in 1:length(classMatches)){
      if(!is.null(classEnds[[i, j]])) classEnd <- c(classEnd, classEnds[[i, j]])
      classCount <- classCount + classCounts[i,j]
    }
    if(classCount > 0) {
      cdsStart <- rnaCDS$CDS.Start[[i]]
      cdsStop <- rnaCDS$CDS.Stop[[i]]
      allClassCounts[i] <- classCount
      allClassEnds[[i]] <- classEnd
      distFrom3 <- length(rnaSet[[i]]) - classEnd
      meanDistFrom3[i] <- mean(distFrom3)
      varDistFrom3[i] <- var(distFrom3)
      medianDistFrom3[i] <- median(distFrom3)
      
      nucFreqMat[i,] <- nucFreqAroundSites(rnaSet[[i]], classEnd)
      
      if(!is.na(cdsStop)){
        distFromStop <- cdsStop - classEnd
        meanDistFromStop[i] <- mean(distFromStop)
        varDistFromStop[i] <- var(distFromStop)
        medianDistFromStop[i] <- median(distFromStop)
        freqIn5[i] <- length(which(classEnd < cdsStart)) / length(classEnd)
        freqInCDS[i] <- length(which(classEnd >= cdsStart 
                                     & classEnd <= cdsStop)) / length(classEnd)
        freqIn3[i] <- length(which(classEnd > cdsStop)) / length(classEnd)
      } else { # Non coding RNA
        meanDistFromStop[i] <- NA
        varDistFromStop[i] <- NA
        medianDistFromStop[i] <- NA
        freqIn5[i] <- NA
        freqInCDS[i] <- NA
        freqIn3[i] <- NA
      }
      
      if(!is.null(centeredClass)){
        freqStart3[i] <- sum(classCounts[i, which(centeredClass == 3)]) / classCount
        freqStart4[i] <- sum(classCounts[i, which(centeredClass == 4)]) / classCount
        freqStart5[i] <- sum(classCounts[i, which(centeredClass == 5)]) / classCount
      }
    } 
  }
  
  colnames(nucFreqMat) <- c(sapply(
    siteIndex, function(x) 
      paste0(nameOfClass, "_freq", c("A", "C", "G", "T"), x)))
  
  outputDF <- as.data.frame(
    cbind(allClassCounts, meanDistFrom3, medianDistFrom3, varDistFrom3,
          meanDistFromStop, medianDistFromStop, varDistFromStop, freqIn5, 
          freqInCDS, freqIn3))
  colnames(outputDF) <- c(
    paste0("n", nameOfClass), 
    paste0(nameOfClass, "_",c("end3'distMean", "end3'distMedian", "end3'distVar", 
                              "stopdistMean", "stopdistMedian", "stopdistVar", 
                              "5'UTRFreq", "CDSFreq", "3'UTRFreq")))
  outputDF <- cbind(outputDF, nucFreqMat)
  
  if(!is.null(centeredClass)){
    cenStart <- cbind(freqStart3, freqStart4, freqStart5)
    colnames(cenStart) <- paste0(nameOfClass, "_", "freqStart", c(3,4,5))
    outputDF <- cbind(outputDF, cenStart)
  }
  
  return(list(allClassEnds, outputDF))
}


createFeatureMatrix <- function(miRNAString, acc, targCol = NULL)
{
  # Fills a feature matrix for either training and testing the classifier
  # or for prediction by a trained classifier
  # Also returns the 3' end of the positions of sites found on each mRNA
  # miRNAString is the sequence of the mature microRNA
  # acc is the GenBank accession numbers for the RNA
  # targCol indicates which RNA are miRNA targets, which is needed for
  # building a feature matrix for training and testing
  
  if(class(miRNAString) %in% c("DNAString", "RNAString")){
    miRNAString <- as.character(miRNAString)
  }
  # Converted to DNAString to match the mRNA sequence type
  # Also because i don't know if wildcard nucleotides work when 
  # using matchPattern with RNAString
  miRNAString <- gsub("U", "T", miRNAString)
  
  miLen <- nchar(miRNAString)
  
  # Create feature values for nucleotide identies of miRNA
  # Assumes max length of miRNA is 26 and deals with varying length of
  # miRNA by using character "X" to represent missing nucleotides
  miNucIdenVec <- strsplit(miRNAString, split = "")[[1]]
  miNucIdenFeats <- rep("X", 26)
  miNucIdenFeats[1:length(miNucIdenVec)] <- miNucIdenVec
  
  miRNAString <- DNAString(miRNAString)
  
  # Collect sequences and CDS start-stop for the accession numbers we have 
  # predownloaded sequences for
  seqAndCDS <- accSeqRegions[
    accSeqRegions$GenBank.Accession %in% acc, ]
  
  # For the accession numbers we do not have predownloaded sequences for:
  # Query Nucleotide using the accession numbers to get a data frame of
  # sequence strings and, for the mRNAs that are protein coding, the CDS
  # start and stop sites
  newAcc <- acc[!(acc %in% seqAndCDS$GenBank.Accession)]
  if(length(newAcc) > 0) {
    newSeqAndCDS <- querymRNASequence(newAcc)
    # Save the new ones for next time
    write.csv(newSeqAndCDS, file = "accSeqRegions.csv", append = TRUE,
              row.names = FALSE)
  
    seqAndCDS <- rbind(seqAndCDS, newSeqAndCDS)
  }
  mRNAStrings <- DNAStringSet(seqAndCDS$Full.Sequence)
  names(mRNAStrings) <- seqAndCDS$GenBank.Accession
  
  # Create search strings for the different classes of sites
  # =========================================
  
  miRNA_revComp <- reverseComplement(miRNAString)
  # Basic seed
  seq_6mer <- subseq(miRNA_revComp, start = (length(miRNA_revComp)-6), 
                     end = (length(miRNA_revComp)-1))
  
  # Preventing overlap in searches
  # get identities of nucleotides at position 1 and 8 of site
  iden8 <- letter(miRNA_revComp, (length(miRNA_revComp)-7))
  iden1 <- letter(miRNA_revComp, (length(miRNA_revComp)))
  # The IUPAC letter for "any nucleotide type except one"
  notMatches <- c("V", "H", "D", "B")
  names(notMatches) <- c("T", "G", "C", "A")
  # The IUPAC letter for "any nucleotide type except A and one other"
  notMatchNotA <- c("S", "Y", "K", "B")
  names(notMatchNotA) <- c("T", "G", "C", "A")
  notMatch8 <- notMatches[iden8]
  notMatch1 <- notMatches[iden1]
  notMorA1 <- notMatchNotA[iden1]
  
  # Class1 search strings
  seq_7merm8 <- xscat(subseq(miRNA_revComp, start = (length(miRNA_revComp)-7), 
                             end = (length(miRNA_revComp)-1)), notMatch1)
  seq_8mer <- subseq(miRNA_revComp, start = (length(miRNA_revComp)-7), 
                     end = length(miRNA_revComp))
  seq_7merA1 <- xscat(notMatch8, subseq(miRNA_revComp, start = (length(miRNA_revComp)-6), 
                                        end = (length(miRNA_revComp)-1)), "A")
  class1Search <- DNAStringSet(list(seq_8mer, seq_7merm8, seq_7merA1))
  
  # Class 2 search strings
  # GU wobble not allowed at position 1 or 8 so concatenate the seach letters
  # for 1 and 8 to GU wobble versions of the core 2-7 seed
  c2GU <- guWob_seqs(seq_6mer)
  c2GU7merm8 <- xscat(letter(miRNA_revComp, length(miRNA_revComp)-7), c2GU, 
                      notMatch1)
  c2GU8mer <- xscat(letter(miRNA_revComp, length(miRNA_revComp)-7), c2GU,
                    letter(miRNA_revComp, length(miRNA_revComp)))
  c2GU7merA1 <- xscat(notMatch8, c2GU, "A")
  c2Full <- DNAStringSet(c(c2GU7merA1, c2GU7merm8, c2GU8mer))
  
  # class 3 search strings
  class3search <- misMatchNotGU(class1Search, start = 2, end = 7)
  
  # class 4 search strings
  c4A1 <- misMatchNotGU(c2GU7merA1, start = 2, end = 7, 
                        reference = seq_7merA1)
  c4m8 <- misMatchNotGU(c2GU7merm8, start = 2, end = 7, 
                        reference = seq_7merm8)
  c4_8mer <- misMatchNotGU(c2GU8mer, start = 2, end = 7, 
                           reference = seq_8mer)
  class4search <- DNAStringSet(c(c4A1, c4m8, c4_8mer))
  
  # class 5 search string
  class5Search <- xscat(notMatch8, seq_6mer, notMorA1)
  
  # class 6 search strings
  # Need to add N wild cards at the end so as site ends returned by search
  # correspond to the position opposite site 1 on the miRNA
  cent3 <- xscat(subseq(miRNA_revComp, start = (length(miRNA_revComp)-12), 
                  end = (length(miRNA_revComp)-2)), "NN")
  cent4 <- xscat(subseq(miRNA_revComp, start = (length(miRNA_revComp)-13), 
                  end = (length(miRNA_revComp)-3)), "NNN")
  cent5 <- xscat(subseq(miRNA_revComp, start = (length(miRNA_revComp)-14), 
                  end = (length(miRNA_revComp)-4)), "NNNN")
  class6Search <- DNAStringSet(list(cent3, cent4, cent5))
  c6SearchCenStart <- c(3,4,5)
  
  # class 7 search strings
  GUcent3 <- guWob_seqs(cent3)
  GUcent4 <- guWob_seqs(cent4)
  GUcent5 <- guWob_seqs(cent5)
  class7Search <- DNAStringSet(c(GUcent3, GUcent4, GUcent5))
  c7SearchCenStart <- c(rep(3, length(GUcent3)),rep(4, length(GUcent4)),
                        rep(5, length(GUcent5)))
  
  # class 8 search strings
  mm_cent3 <- misMatchNotGU(cent3, end = 11)
  mm_cent4 <- misMatchNotGU(cent4, end = 11)
  mm_cent5 <- misMatchNotGU(cent5, end = 11)
  class8Search <- DNAStringSet(c(mm_cent3, mm_cent4, mm_cent5))
  c8SearchCenStart <- c(rep(3, length(mm_cent3)),rep(4, length(mm_cent4)),
                        rep(5, length(mm_cent5)))
  
  # class 9 search strings
  c9cent3 <- misMatchNotGU(GUcent3, end = 11, reference = cent3)
  c9cent4 <- misMatchNotGU(GUcent4, end = 11, reference = cent4)
  c9cent5 <- misMatchNotGU(GUcent5, end = 11, reference = cent5)
  class9Search <- DNAStringSet(c(c9cent3, c9cent4, c9cent5))
  c9SearchCenStart <- c(rep(3, length(c9cent3)),rep(4, length(c9cent4)),
                        rep(5, length(c9cent5)))
  
  # Find sites and measure values for site features 
  sitesClass1 <- sitesOfClass(class1Search, "c1", mRNAStrings, seqAndCDS)
  sitesClass2 <- sitesOfClass(c2Full, "c2", mRNAStrings, seqAndCDS) 
  sitesClass3 <- sitesOfClass(class3search, "c3", mRNAStrings, seqAndCDS) 
  sitesClass4 <- sitesOfClass(class4search, "c4", mRNAStrings, seqAndCDS) 
  sitesClass5 <- sitesOfClass(class5Search, "c5", mRNAStrings, seqAndCDS) 
  sitesClass6 <- sitesOfClass(class6Search, "c6", mRNAStrings, 
                              seqAndCDS, c6SearchCenStart) 
  sitesClass7 <- sitesOfClass(class7Search, "c7", mRNAStrings, 
                              seqAndCDS, c7SearchCenStart) 
  sitesClass8 <- sitesOfClass(class8Search, "c8", mRNAStrings, 
                              seqAndCDS, c8SearchCenStart) 
  sitesClass9 <- sitesOfClass(class9Search, "c9", mRNAStrings, 
                              seqAndCDS, c9SearchCenStart) 
  
  # Combine the tables
  fullFeatTable <- cbind(sitesClass1[[2]], sitesClass2[[2]],
                         sitesClass3[[2]], sitesClass4[[2]],
                         sitesClass5[[2]], sitesClass6[[2]],
                         sitesClass7[[2]], sitesClass8[[2]],
                         sitesClass9[[2]])
  
  # Lists of position ends. Can be used for analysis or extracting extra
  # features, etc
  siteEnds <- t(sapply(
    1:length(acc), function(x) list(sitesClass1[[1]][[x]], sitesClass2[[1]][[x]],
                                    sitesClass3[[1]][[x]], sitesClass4[[1]][[x]],
                                    sitesClass5[[1]][[x]], sitesClass6[[1]][[x]],
                                    sitesClass7[[1]][[x]], sitesClass8[[1]][[x]],
                                    sitesClass9[[1]][[x]])))
  
  # Add collumn for total sites found accross all classes
  fullFeatTable$N <- rowSums(fullFeatTable[,paste0("nc", 1:9)])
  
  # Add collumns for identities of miRNA sequence
  miNucIdenFeatsMat <- t(matrix(miNucIdenFeats, ncol = length(acc),
                                nrow = length(miNucIdenFeats)))
  colnames(miNucIdenFeatsMat) <- paste0("mi", 1:26)
  fullFeatTable <- cbind(fullFeatTable, miNucIdenFeatsMat)
  
  # Add collumns for length of miRNA and mRNA
  fullFeatTable$miLen <- miLen
  fullFeatTable$mRNALen <- sapply(
    as.character(seqAndCDS$Full.Sequence), nchar)
  
  fullFeatTable$GenBank.Accession <- seqAndCDS$GenBank.Accession
  
  if(!is.null(targCol)){
    accToTarget <- cbind(acc, targCol)
    colnames(accToTarget) <- c("GenBank.Accession", "Target")
    
    # So we can supervise learning
    fullFeatTable <- merge(fullFeatTable, accToTarget, 
                           by = "GenBank.Accession")
    # Remove targets where no sites were found (should be rare)
    beforeRemoval <- nrow(fullFeatTable[targCol == 1,])
    fullFeatTable <- fullFeatTable[!(targCol == 1 & fullFeatTable$N == 0), ]
    afterRemoval <- nrow(fullFeatTable[fullFeatTable$Target == 1,])
    numTargRemoved <- beforeRemoval - afterRemoval
    print(paste("Removed", numTargRemoved, 
                "targets with no binding sites found out of", beforeRemoval))
  }
  
  return(list(siteEnds, mRNAStrings, fullFeatTable))
}


siteDuplexFreeEnergy <- function(miRNA, mRNAStrings, siteEnds)
{
  output <- as.data.frame(matrix(nrow = 0, ncol = 10*4))
  featNames <- sapply(
    paste0("c", 1:9), function(x) paste0(
      x, paste0("_dupFreeEng", c("Min", "Mean", "Median", "Var"))))
  featNames <- c(featNames, paste0(
    "dupFreeEng", c("Min", "Mean", "Median", "Var")))

  seqFileString <- paste0(">miRNA\n", miRNA, "\n")
  # Indicating which energy value of a site comes from which RNA
  siteFromRNA <- c()
  # Indicating which class of site the energy value is for
  siteFromClass <- c()

  # Writing the site seq file
  for(i in 1:length(mRNAStrings)) {
    mRNAString <- mRNAStrings[[i]]
    mRNAName <- names(mRNAStrings[[i]])
    for(j in 1:9) {
      classSiteEnds <- siteEnds[i,][[j]]
      if(!is.null(classSiteEnds)){
        for(k in 1:length(classSiteEnds)){
          siteEnd <- classSiteEnds[k]
          siteSeq <- subseq(
            mRNAString, start = max(1, siteEnd - 25),
            end = classSiteEnds[k])
          seqFileString <- paste0(
            seqFileString, ">mRNA ", mRNAName, ", site end ", siteEnd,
            ", class ", j, "\n", as.character(siteSeq), "\n")
          siteFromRNA <- c(siteFromRNA, i)
          siteFromClass <- c(siteFromClass, j)
        }
      }
    }
  }
  write(seqFileString, "siteSeqs.seq")

  # Running RNAup on the file
  rnaUpOutput <- shell(
    paste0(vRNADir, 'RNAup.exe" -b --no_output_file < siteSeqs.seq'),
    intern = TRUE)
  rnaUpOutEng <- rnaUpOutput[seq(3, length(rnaUpOutput), 3)]
  freeEngStr <- regmatches(
    rnaUpOutEng, regexpr("=\\s-?\\d+\\.\\d\\d", rnaUpOutEng))
  freeEngStr <- gsub("= ", "", freeEngStr)
  freeEng <- as.numeric(freeEngStr)

  # Creating the features
  engMat <- cbind(freeEng, siteFromRNA, siteFromClass)
  for(i in 1:length(mRNAStrings)){
    totalSiteEngs <- engMat[engMat[,2] == i, ]
    classFeats <- rep(NaN, 4*9)
    everyFeats <- rep(NaN, 4)
    if(nrow(totalSiteEngs) > 0){
      everyFeats <- c(min(totalSiteEngs[,1]), mean(totalSiteEngs[,1]),
                      median(totalSiteEngs[,1]), var(totalSiteEngs[,1]))
      for(j in 1:9) {
        classEngs <- totalSiteEngs[,3] == j
        if(TRUE %in% classEngs){
          classFeats[4*(j-1)+1] <- min(totalSiteEngs[classEngs,1])
          classFeats[4*(j-1)+2] <- mean(totalSiteEngs[classEngs,1])
          classFeats[4*(j-1)+3] <- median(totalSiteEngs[classEngs,1])
          classFeats[4*(j-1)+4] <- var(totalSiteEngs[classEngs,1])
        }
      }
    }
    allFeats <- c(classFeats, everyFeats)
    allFeats <- matrix(allFeats, nrow = 1, ncol = length(allFeats))
    colnames(allFeats) <- featNames
    output <- rbind(output, allFeats)
  }
  return(output)
}

