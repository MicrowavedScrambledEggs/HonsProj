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
targetRNASeq <- gsub("T", "U", targetSequences[,2])
targetRNAStrings <- RNAStringSet(targetRNASeq)

miR23b <- RNAString("AUCACAUUGCCAGGGAUUACC")
miR23b_revComp <- reverseComplement(miR23b)

# Class5 search strings
seq_6mer <- subseq(miR23b_revComp, start = (length(miR23b_revComp)-6), 
                   end = (length(miR23b_revComp)-1))

# Class1 and 3 search strings
seq_7merm8 <- subseq(miR23b_revComp, start = (length(miR23b_revComp)-7), 
                     end = (length(miR23b_revComp)-1))
seq_8mer <- subseq(miR23b_revComp, start = (length(miR23b_revComp)-7), 
                   end = length(miR23b_revComp))
seq_7merA1 <- xscat(subseq(miR23b_revComp, start = (length(miR23b_revComp)-6), 
                     end = (length(miR23b_revComp)-1)), "A")
class1Search <- RNAStringSet(list(seq_8mer, seq_7merm8, seq_7merA1))

# Class 2 and 4 search strings
guWob_seqs <- function(seq){
  output <- c()
  seqQ <- queue()
  pushback(seqQ, list(seq, 1))
  while(length(seqQ) > 0){
    sTup <- pop(seqQ)
    sCurrent <- sTup[[1]]
    startPos <- sTup[[2]]
    for(i in c(startPos:length(sCurrent))){
      if(letter(sCurrent, i) == "C"){
        snew <- replaceAt(sCurrent, IRanges(i,i), "U")
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
  return(RNAStringSet(output))
}
c2GU <- guWob_seqs(seq_6mer)
c2GU7merm8 <- xscat(letter(miR23b_revComp, length(miR23b_revComp)-7), c2GU)
c2GU8mer <- xscat(letter(miR23b_revComp, length(miR23b_revComp)-7), c2GU,
                  letter(miR23b_revComp, length(miR23b_revComp)))
c2GU7merA1 <- xscat(c2GU, "A")
c2Full <- RNAStringSet(c(c2GU7merA1, c2GU7merm8, c2GU8mer))

# class 6 and 8 search strings
cent3 <- subseq(miR23b_revComp, start = (length(miR23b_revComp)-12), 
                end = (length(miR23b_revComp)-2))
cent4 <- subseq(miR23b_revComp, start = (length(miR23b_revComp)-13), 
                end = (length(miR23b_revComp)-3))
cent5 <- subseq(miR23b_revComp, start = (length(miR23b_revComp)-14), 
                end = (length(miR23b_revComp)-4))
class6Search <- RNAStringSet(c(cent3, cent4, cent5))

# class 7 and 9 search strings
GUcent3 <- guWob_seqs(cent3)
GUcent4 <- guWob_seqs(cent4)
GUcent5 <- guWob_seqs(cent5)
class7Search <- RNAStringSet(c(GUcent3, GUcent4, GUcent5))

# TO:DO sort out overlaping matches in the same class

class1_counts <- sapply(class1Search, vcountPattern, subject = targetRNAStrings)
class2_counts <- sapply(c2Full, vcountPattern, subject = targetRNAStrings)
class3_counts <- sapply(class1Search, vcountPattern, subject = targetRNAStrings,
                        min.mismatch = 1, max.mismatch = 1)
class4_counts <- sapply(c2Full, vcountPattern, subject = targetRNAStrings,
                        min.mismatch = 1, max.mismatch = 1)
class5_counts <- vcountPattern(seq_6mer, targetRNAStrings)
class6_counts <- sapply(class6Search, vcountPattern, subject = targetRNAStrings)
class7_counts <- sapply(class7Search, vcountPattern, subject = targetRNAStrings)
class8_counts <- sapply(class6Search, vcountPattern, subject = targetRNAStrings,
                        min.mismatch = 1, max.mismatch = 1)
class9_counts <- sapply(class7Search, vcountPattern, subject = targetRNAStrings,
                        min.mismatch = 1, max.mismatch = 1)
siteCountMat <- cbind(class1_counts, class2_counts, class3_counts, class4_counts,
                      class5_counts, class6_counts, class7_counts, class8_counts, 
                      class9_counts)




