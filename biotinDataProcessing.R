library(GEOquery)
library(limma)
library(rentrez)
library(Biostrings)
library(dequer)

cloonan23b27a <- getGEO("GSE40410", GSEMatrix =TRUE, AnnotGPL=TRUE)
#comes up with a "13 parsing failures" warning
featDat <- featureData(cloonan23b27a$GSE40410_series_matrix.txt.gz)
featDatTable <- fData(cloonan23b27a$GSE40410_series_matrix.txt.gz)
assayDat <- assayData(cloonan23b27a$GSE40410_series_matrix.txt.gz)

# Volcano plot for miR-23b. Feeling like a script kidde here
layout <- c("pulldown", "control", "pulldown", "control")
design<-model.matrix(~0+layout)
colnames(design)<-c("control","pulldown")
fit<-lmFit(assayDat$exprs[,c(1:4)], design)
cont.matrix <- makeContrasts(PullDownvsControl=pulldown-control, 
                             levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
full_gene_list<-topTable(fit2, coef=1, number=1000000, sort.by="logFC")
plot(full_gene_list$logFC, -log10(full_gene_list$P.Value), 
     xlab="log2 fold change", ylab="-log10 p-value")

# miR-23b targets defined as below log2FoldChange threshold -0.6 and 
# above -log10pValue threshold of -log10(0.05) 
# Maybe need a better reason for these thresholds other than 
# "because that's what other bpd studies used"
targets <- full_gene_list[full_gene_list$logFC < -0.6 
                          & -log10(full_gene_list$P.Value) > -log10(0.05), ]
# Is everything else really non-targets?
nontargets <- full_gene_list[!row.names(full_gene_list) 
                             %in% row.names(targets), ]
# Get genbank accession numbers
genbankTargetAcc <- featDatTable$`GenBank Accession`[featDatTable$ID 
                                               %in% row.names(targets)]
genbankTargetAcc <- genbankTargetAcc[genbankTargetAcc != ""]
# Querry nucleotide
# Parse the sequences out of the FASTA
getmRNASequence <- function(accNo){
  fasta <- entrez_fetch("nucleotide", id=accNo, rettype = "fasta")
  seqLines <- unlist(strsplit(fasta, "\n"))
  return(paste0(seqLines[2:(length(seqLines)-1)], collapse = ""))
}
targetSequences <- cbind(genbankTargetAcc, sapply(genbankTargetAcc, 
                                                  getmRNASequence))

# Converted to DNAString to match the mRNA sequence type
# Also because i don't know if wildcard nucleotides work when 
# using matchPattern with RNAString
targetRNAStrings <- DNAStringSet(targetRNASeq)
miR23b <- DNAString(gsub("U", "T","AUCACAUUGCCAGGGAUUACC"))
miR23b_revComp <- reverseComplement(miR23b)

# Basic seed
seq_6mer <- subseq(miR23b_revComp, start = (length(miR23b_revComp)-6), 
                   end = (length(miR23b_revComp)-1))

# Preventing overlap in searches
# get identities of nucleotides at position 1 and 8 of site
iden8 <- letter(miR23b_revComp, (length(miR23b_revComp)-7))
iden1 <- letter(miR23b_revComp, (length(miR23b_revComp)))
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
seq_7merm8 <- xscat(subseq(miR23b_revComp, start = (length(miR23b_revComp)-7), 
                     end = (length(miR23b_revComp)-1)), notMatch1)
seq_8mer <- subseq(miR23b_revComp, start = (length(miR23b_revComp)-7), 
                   end = length(miR23b_revComp))
seq_7merA1 <- xscat(notMatch8, subseq(miR23b_revComp, start = (length(miR23b_revComp)-6), 
                     end = (length(miR23b_revComp)-1)), "A")
class1Search <- DNAStringSet(list(seq_8mer, seq_7merm8, seq_7merA1))

# Class 2 search strings

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

# GU wobble not allowed at position 1 or 8 so concatenate the seach letters
# for 1 and 8 to GU wobble versions of the core 2-7 seed
c2GU <- guWob_seqs(seq_6mer)
c2GU7merm8 <- xscat(letter(miR23b_revComp, length(miR23b_revComp)-7), c2GU, 
                    notMatch1)
c2GU8mer <- xscat(letter(miR23b_revComp, length(miR23b_revComp)-7), c2GU,
                  letter(miR23b_revComp, length(miR23b_revComp)))
c2GU7merA1 <- xscat(notMatch8, c2GU, "A")
c2Full <- DNAStringSet(list(c2GU7merA1, c2GU7merm8, c2GU8mer))

# class 3 search strings

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
cent3 <- subseq(miR23b_revComp, start = (length(miR23b_revComp)-12), 
                end = (length(miR23b_revComp)-2))
cent4 <- subseq(miR23b_revComp, start = (length(miR23b_revComp)-13), 
                end = (length(miR23b_revComp)-3))
cent5 <- subseq(miR23b_revComp, start = (length(miR23b_revComp)-14), 
                end = (length(miR23b_revComp)-4))
class6Search <- DNAStringSet(list(cent3, cent4, cent5))

# class 7 search strings
GUcent3 <- guWob_seqs(cent3)
GUcent4 <- guWob_seqs(cent4)
GUcent5 <- guWob_seqs(cent5)
class7Search <- DNAStringSet(c(GUcent3, GUcent4, GUcent5))

# class 8 search strings
class8Search <- misMatchNotGU(class6Search)

# class 9 search strings
c9cent3 <- misMatchNotGU(GUcent3, reference = cent3)
c9cent4 <- misMatchNotGU(GUcent4, reference = cent4)
c9cent5 <- misMatchNotGU(GUcent5, reference = cent5)
class9Search <- DNAStringSet(c(c9cent3, c9cent4, c9cent5))

# counting sites
class1_counts <- sapply(class1Search, vcountPattern, subject = targetRNAStrings,
                        fixed = "subject")
class1_count <- rowSums(class1_counts)
class2_counts <- sapply(c2Full, vcountPattern, subject = targetRNAStrings,
                        fixed = "subject")
class2_count <- rowSums(class2_counts)
class3_counts <- sapply(class3search, vcountPattern, subject = targetRNAStrings,
                        fixed = "subject")
class3_count <- rowSums(class3_counts)
class4_counts <- sapply(class4search, vcountPattern, subject = targetRNAStrings,
                        fixed = "subject")
class4_count <- rowSums(class4_counts)
class5_counts <- vcountPattern(class5Search, targetRNAStrings, fixed = "subject")
class6_counts <- sapply(class6Search, vcountPattern, subject = targetRNAStrings,
                        fixed = "subject")
class6_count <- rowSums(class6_counts)
class7_counts <- sapply(class7Search, vcountPattern, subject = targetRNAStrings,
                        fixed = "subject")
class7_count <- rowSums(class7_counts)
class8_counts <- sapply(class8Search, vcountPattern, subject = targetRNAStrings,
                        fixed = "subject")
class8_count <- rowSums(class8_counts)
class9_counts <- sapply(class9Search, vcountPattern, subject = targetRNAStrings,
                        fixed = "subject")
class9_count <- rowSums(class9_counts)
siteCountMat <- cbind(class1_count, class2_count, class3_count, class4_count,
                      class5_counts, class6_count, class7_count, class8_count, 
                      class9_count)

totalSites <- rowSums(siteCountMat)



