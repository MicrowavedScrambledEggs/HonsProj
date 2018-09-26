library(GEOquery)
library(limma)
source("CreateFeatureMatrix.R")

cloonan23b27a <- getGEO("GSE40410", GSEMatrix =TRUE, AnnotGPL=TRUE)
#comes up with a "13 parsing failures" warning
featDat.23b.27a <- featureData(cloonan23b27a$GSE40410_series_matrix.txt.gz)
featDatTable.23b.27a <- fData(cloonan23b27a$GSE40410_series_matrix.txt.gz)
assayDat.23b.27a <- assayData(cloonan23b27a$GSE40410_series_matrix.txt.gz)

measurelogFC <- function(assayDat, layout, cols){
  design<-model.matrix(~0+layout)
  colnames(design)<-c("control","pulldown")
  fit<-lmFit(assayDat$exprs[,cols], design)
  cont.matrix <- makeContrasts(PullDownvsControl=pulldown-control, 
                               levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  full_gene_list<-topTable(fit2, coef=1, number=1000000, sort.by="logFC")
  plot(full_gene_list$logFC, -log10(full_gene_list$P.Value), 
       xlab="log2 fold change", ylab="-log10 p-value")
  return(full_gene_list)
}

miR23b <- "AUCACAUUGCCAGGGAUUACC"
# Volcano plot for miR-23b.
layout <- c("pulldown", "control", "pulldown", "control")
full_gene_list_23b <- measurelogFC(assayDat.23b.27a, layout, 1:4)

accToTargMat <- function(full_gene_list, featDatTable){
  # targets defined as below log2FoldChange threshold 1 and 
  # below p-value threshold of 0.05
  # This is to ensure statistically significant targets that showed up at least
  # twice as much in the pull down compared to the control
  targets <- full_gene_list[full_gene_list$logFC > 1 
                            & full_gene_list$P.Value < 0.05, ]
  # choose NOn targets that barely showed up in the pull downs
  nontargets <- full_gene_list[full_gene_list$logFC < -1 
                               & full_gene_list$P.Value < 0.05, ]
  # Get genbank accession numbers
  genbankAcc <- featDatTable$`GenBank Accession`[
    row.names(featDatTable) %in% row.names(targets)
    ]
  genbankAcc <- c(genbankAcc, featDatTable$`GenBank Accession`[
    row.names(featDatTable) %in% row.names(nontargets)])
  
  # Create a column for 1 being target and 0 being non target
  targetCol <- c(rep(1, nrow(targets)), rep(0, nrow(nontargets)))
  
  accToTarget <- cbind(genbankAcc, targetCol)
  # Removing repeats and RNA with no acc numbers
  accToTarget <- accToTarget[accToTarget[,1] != "", ]
  accToTarget <- dealWithRepAcc(accToTarget)
  return(accToTarget)
}

dealWithRepAcc <- function(accToTarget){
  accCounts <- table(accToTarget[,1])
  repeatedAcc <- names(accCounts[accCounts > 1])
  toRemove <- c()
  for(acc in repeatedAcc){
    targCount <- table(accToTarget[accToTarget[,1] == acc, 2])
    # Mark the ambigous Accession numbers for removal
    if(length(targCount) > 1){
      toRemove <- c(toRemove, acc)
    }
  }
  accToTarget <- accToTarget[!(accToTarget[,1] %in% toRemove),]
  accToTarget <- unique(accToTarget)
  return(accToTarget)
}

accToTarget_23b <- accToTargMat(full_gene_list_23b, featDatTable.23b.27a)

mir23b_featMat <- createFeatureMatrix(miR23b, accToTarget_23b[,1], 
                                      accToTarget_23b[,2])
freeEngFeats_23b <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]], mir23b_featMat[[1]])

# Processing miR27a data
miR27a <- "UUCACAGUGGCUAAGUUCCGC"
# Columns 5:8 are miR27a
full_gene_list_27a <- measurelogFC(assayDat.23b.27a, layout, 5:8)
accToTarget_27a <- accToTargMat(full_gene_list_27a, featDatTable.23b.27a)

mir27a_featMat <- createFeatureMatrix(miR27a, accToTarget_27a[,1], 
                                      accToTarget_27a[,2])
freeEngFeats_27a <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]], mir27a_featMat[[1]])

freeEngFeats_23b$GenBank.Accession <- names(mir23b_featMat[[2]])
freeEngFeats_27a$GenBank.Accession <- names(mir27a_featMat[[2]])
# Have to merge instead of cbind as rows are not necesarily in the same order
featMatmir23b <- merge(mir23b_featMat[[3]], freeEngFeats_23b, 
                       by = "GenBank.Accession")
featMatmir27a <- merge(mir27a_featMat[[3]], freeEngFeats_27a,
                       by = "GenBank.Accession")

convertMissingValues <- function(featMat){
  # Dealing with missing values
  # NaN for mean, median and variance of distances of sites (of each class) 
  # from 3' end and from stop codon, and frequencies of nucleotides at positions
  # around sites, represent that there were no sites of that class
  #   - For nucleotide frequencies, set NaN to 0
  #   - For distance from 3' end set NaN to -1
  #   - For distance from stop codon set NaN to -max double
  #   - For siteFreeEng set NaN to max double
  # NA for variance of distances of sites (of each class) represent that
  # there was only one site of that class
  #   - For variance, set NA to 0
  # NA for mean meadian and variance of distances of sites from stop codon
  # and for frequencies of sites in regions of mRNA represent the RNA is 
  # not mRNA
  #   - For mean and median set NA to max double
  #   - For variance set NA to 0
  #   - For region frequencies set NA to 0
  
  # NaN nucleotide frequencies
  nucFreqCols <- colnames(featMat)[
    grepl("c\\d_freq", colnames(featMat))]
  featMat[,nucFreqCols] <- t(
    apply(featMat[,nucFreqCols], 1, 
          function(x) {x[is.nan(x)] <- 0; return(x)}))
  # NaN 3' end distances
  endDistCols <- colnames(featMat)[
    grepl("c\\d_end", colnames(featMat))]
  featMat[,endDistCols] <- t(
    apply(featMat[,endDistCols], 1, 
          function(x) {x[is.nan(x)] <- -1; return(x)}))
  # NaN stop codon distances
  stopDistCols <- colnames(featMat)[
    grepl("c\\d_stop", colnames(featMat))]
  featMat[,stopDistCols] <- t(
    apply(featMat[,stopDistCols], 1, 
          function(x) {x[is.nan(x)] <- -.Machine$double.xmax; return(x)}))
  # NaN duplex free energy
  dupFreeEngCols <- colnames(featMat)[
    grepl("dup", colnames(featMat))]
  featMat[,dupFreeEngCols] <- t(
    apply(featMat[,dupFreeEngCols], 1, 
          function(x) {x[is.nan(x)] <- .Machine$double.xmax; return(x)}))
  # NA variance
  varCols <- colnames(featMat)[
    grepl("Var", colnames(featMat))]
  featMat[,varCols] <- t(
    apply(featMat[,varCols], 1, 
          function(x) {x[is.na(x)] <- 0; return(x)}))
  # NA stop codon distances
  featMat[,stopDistCols] <- t(
    apply(featMat[,stopDistCols], 1, 
          function(x) {x[is.na(x)] <- .Machine$double.xmax; return(x)}))
  # NA region frequencies
  regionCols <- colnames(featMat)[
    grepl("c\\d_(5'UTR|CDS|3'UTR)", colnames(featMat))]
  featMat[,regionCols] <- t(
    apply(featMat[,regionCols], 1, 
          function(x) {x[is.na(x)] <- 0; return(x)}))
  
  # Random Forest does not like column names with symbols in them
  colnames(featMat) <- gsub("'", "", colnames(featMat))
  colnames(featMat) <- gsub("-", "neg", colnames(featMat))
  featMat$Target <- as.factor(featMat$Target)
  
  # Make sure categorical features have all factor levels
  miIdenCols <- colnames(featMat)[grepl("mi\\d", colnames(featMat))]
  for(colu in miIdenCols)
    featMat[,colu] <- factor(featMat[,colu], levels = c("A", "C", "G", "T", "X"))
  
  return(featMat)
}
featMatmir23b <- convertMissingValues(featMatmir23b)
featMatmir27a <- convertMissingValues(featMatmir27a)

# MORE DATA
cloonan3118_4307 <- getGEO("GSE40409", GSEMatrix =TRUE, AnnotGPL=TRUE)
featDat3118_4307 <- featureData(cloonan3118_4307$GSE40409_series_matrix.txt.gz)
featDatTable3118_4307 <- fData(cloonan3118_4307$GSE40409_series_matrix.txt.gz)
assayDat3118_4307 <- assayData(cloonan3118_4307$GSE40409_series_matrix.txt.gz)
phenoDat3118_4307 <- phenoData(cloonan3118_4307$GSE40409_series_matrix.txt.gz)
pDatTable3118_4307 <- pData(cloonan3118_4307$GSE40409_series_matrix.txt.gz)

miR3118 <- "UGUGACUGCAUUAUGAAAAUUCU"
layout <- c("pulldown", "control", "pulldown", "control", "pulldown", "control")
cols.3118 <- c(1,2,5,6,9,10)
full_gene_list_3118 <- measurelogFC(assayDat3118_4307, layout, cols.3118)

miR4307 <- "AAUGUUUUUUCCUGUUUCC"
# Columns 3,4,7,8,11,12 are miR-4307
cols.4307 <- c(3,4,7,8,11,12)
full_gene_list_4307 <- measurelogFC(assayDat3118_4307, layout, cols.4307)

accToTarget_3118 <- accToTargMat(full_gene_list_3118, featDatTable3118_4307)
mir3118_featMat <- createFeatureMatrix(miR3118, accToTarget_3118[,1], 
                                       accToTarget_3118[,2])

# miR4307 targets and features
accToTarget_4307 <- accToTargMat(full_gene_list_4307, featDatTable3118_4307)
mir4307_featMat <- createFeatureMatrix(miR4307, accToTarget_4307[,1], 
                                       accToTarget_4307[,2])

freeEngFeats_3118 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]], mir3118_featMat[[1]])

freeEngFeats_4307 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]], mir4307_featMat[[1]])

freeEngFeats_3118$GenBank.Accession <- names(mir3118_featMat[[2]])
freeEngFeats_4307$GenBank.Accession <- names(mir4307_featMat[[2]])

featMatmir3118 <- merge(mir3118_featMat[[3]], freeEngFeats_3118, 
                        by = "GenBank.Accession")
featMatmir3118 <- convertMissingValues(featMatmir3118)
featMatmir4307 <- merge(mir4307_featMat[[3]], freeEngFeats_4307, 
                        by = "GenBank.Accession")
featMatmir4307 <- convertMissingValues(featMatmir4307)

fullFeatMat <- rbind(featMatmir23b, featMatmir27a, featMatmir3118, 
                     featMatmir4307)
fullFeatMatNoAcc <- fullFeatMat[, !colnames(fullFeatMat) %in% c("GenBank.Accession")]

# More data to test performance on miRNA not included in training set

cloonan.199a.424 <- getGEO("GSE40406", GSEMatrix =TRUE, AnnotGPL=TRUE)
featDat.199a.424 <- featureData(cloonan.199a.424$GSE40406_series_matrix.txt.gz)
featDatTable.199a.424 <- fData(cloonan.199a.424$GSE40406_series_matrix.txt.gz)
assayDat.199a.424 <- assayData(cloonan.199a.424$GSE40406_series_matrix.txt.gz)
phenoDat.199a.424 <- phenoData(cloonan.199a.424$GSE40406_series_matrix.txt.gz)
pDatTable.199a.424 <- pData(cloonan.199a.424$GSE40406_series_matrix.txt.gz)
# check experiment design
pDatTable.199a.424[, c("title", "pull-down:ch1", "transfection molecule:ch1", "data_processing")]

layout.199a.424 <- c("pulldown", "pulldown", "pulldown","control","control","control")
cols.199a <- 1:6
fgl_199a <- measurelogFC(assayDat.199a.424, layout.199a.424, cols.199a)
accToTarget.199a <- accToTargMat(fgl_199a, featDatTable.199a.424)
miR199a3p <- "ACAGUAGUCUGCACAUUGGUUA"
sites.rnaSeq.featMat.199a <- createFeatureMatrix(miR199a3p, accToTarget.199a[,1],
                                                 accToTarget.199a[,2])
freeEngFeats.199a <- siteDuplexFreeEnergy(
  miR199a3p, sites.rnaSeq.featMat.199a[[2]], sites.rnaSeq.featMat.199a[[1]])
freeEngFeats.199a$GenBank.Accession <- names(sites.rnaSeq.featMat.199a[[2]])
featMat.199a <- merge(sites.rnaSeq.featMat.199a[[3]], freeEngFeats.199a,
                      by = "GenBank.Accession")
featMat.199a <- convertMissingValues(featMat.199a)
featMat.199a.noAcc <- featMat.199a[
  , !colnames(featMat.199a) %in% c("GenBank.Accession")]

cols.424 <- 7:12
fgl_424 <- measurelogFC(assayDat.199a.424, layout.199a.424, cols.424)
accToTarget.424 <- accToTargMat(fgl_424, featDatTable.199a.424)
miR424.3p <- "CAAAACGUGAGGCGCUGCUAU"
sites.rnaSeq.featMat.424 <- createFeatureMatrix(miR424.3p, accToTarget.424[,1],
                                                 accToTarget.424[,2])
freeEngFeats.424 <- siteDuplexFreeEnergy(
  miR424.3p, sites.rnaSeq.featMat.424[[2]], sites.rnaSeq.featMat.424[[1]])
freeEngFeats.424$GenBank.Accession <- names(sites.rnaSeq.featMat.424[[2]])
featMat.424 <- merge(sites.rnaSeq.featMat.424[[3]], freeEngFeats.424,
                      by = "GenBank.Accession")
featMat.424 <- convertMissingValues(featMat.424)

# RAW DATA



