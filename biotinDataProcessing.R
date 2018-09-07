library(GEOquery)
library(limma)
source("CreateFeatureMatrix.R")

cloonan23b27a <- getGEO("GSE40410", GSEMatrix =TRUE, AnnotGPL=TRUE)
#comes up with a "13 parsing failures" warning
featDat <- featureData(cloonan23b27a$GSE40410_series_matrix.txt.gz)
featDatTable <- fData(cloonan23b27a$GSE40410_series_matrix.txt.gz)
assayDat <- assayData(cloonan23b27a$GSE40410_series_matrix.txt.gz)

miR23b <- "AUCACAUUGCCAGGGAUUACC"
# Volcano plot for miR-23b. Feeling like a script kidde here
layout <- c("pulldown", "control", "pulldown", "control")
design<-model.matrix(~0+layout)
colnames(design)<-c("control","pulldown")
# Columns 1:4 are miR23b
fit_23b<-lmFit(assayDat$exprs[,c(1:4)], design)
cont.matrix.23b <- makeContrasts(PullDownvsControl=pulldown-control, 
                             levels=design)
fit2_23b <- contrasts.fit(fit_23b, cont.matrix.23b)
fit2_23b <- eBayes(fit2_23b)
full_gene_list_23b<-topTable(fit2_23b, coef=1, number=1000000, sort.by="logFC")
plot(full_gene_list_23b$logFC, -log10(full_gene_list_23b$P.Value), 
     xlab="log2 fold change", ylab="-log10 p-value")

# miR-23b targets defined as below log2FoldChange threshold -1 and 
# above -log10pValue threshold of -log10(0.05) 
# Maybe need a better reason for these thresholds other than 
# "because that's what other bpd studies used"
targets_23b <- full_gene_list_23b[full_gene_list_23b$logFC < -1 
                          & -log10(full_gene_list_23b$P.Value) > -log10(0.05), ]
# Is everything else really non-targets?
nontargets_23b <- full_gene_list_23b[!row.names(full_gene_list_23b) 
                             %in% row.names(targets_23b), ]
# Get genbank accession numbers
genbankAcc_23b <- featDatTable$`GenBank Accession`
# Create a column for 1 being target and 0 being non target
targetCol_23b <- rep(0, nrow(full_gene_list_23b))
targetCol_23b[row.names(full_gene_list_23b) %in% row.names(targets_23b)] <- 1

accToTarget_23b <- cbind(genbankAcc_23b, targetCol_23b)
# Removing repeats and RNA with no acc numbers
accToTarget_23b <- accToTarget_23b[accToTarget_23b[,1] != "", ]

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
accToTarget_23b <- dealWithRepAcc(accToTarget_23b)

mir23b_featMat <- createFeatureMatrix(miR23b, accToTarget_23b[,1], 
                                      accToTarget_23b[,2])
freeEngFeats_23b1 <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]][1:1000], mir23b_featMat[[1]][1:1000,])
freeEngFeats_23b2 <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]][1001:2000], mir23b_featMat[[1]][1001:2000,])
freeEngFeats_23b3 <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]][2001:3000], mir23b_featMat[[1]][2001:3000,])
freeEngFeats_23b4 <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]][3001:4000], mir23b_featMat[[1]][3001:4000,])
freeEngFeats_23b5 <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]][4001:5000], mir23b_featMat[[1]][4001:5000,])
freeEngFeats_23b6 <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]][5001:6000], mir23b_featMat[[1]][5001:6000,])
freeEngFeats_23b7 <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]][6001:7000], mir23b_featMat[[1]][6001:7000,])
freeEngFeats_23b8 <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]][7001:8000], mir23b_featMat[[1]][7001:8000,])
freeEngFeats_23b9 <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]][8001:9000], mir23b_featMat[[1]][8001:9000,])
freeEngFeats_23b10 <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]][9001:10000], mir23b_featMat[[1]][9001:10000,])
freeEngFeats_23b11 <- siteDuplexFreeEnergy(
  miR23b, mir23b_featMat[[2]][10001:length(mir23b_featMat[[2]])], 
  mir23b_featMat[[1]][10001:length(mir23b_featMat[[2]]),])

# Processing miR27a data
miR27a <- "UUCACAGUGGCUAAGUUCCGC"
# Columns 5:8 are miR27a
fit_27a<-lmFit(assayDat$exprs[,c(5:8)], design)
cont.matrix.27a <- makeContrasts(PullDownvsControl=pulldown-control, 
                                 levels=design)
fit2_27a <- contrasts.fit(fit_27a, cont.matrix.27a)
fit2_27a <- eBayes(fit2_27a)
full_gene_list_27a<-topTable(fit2_27a, coef=1, number=1000000, sort.by="logFC")
plot(full_gene_list_27a$logFC, -log10(full_gene_list_27a$P.Value), 
     xlab="log2 fold change", ylab="-log10 p-value")

# miR-27a targets defined as below log2FoldChange threshold -1 and 
# above -log10pValue threshold of -log10(0.05) 
targets_27a <- full_gene_list_27a[full_gene_list_27a$logFC < -1 
                                  & -log10(full_gene_list_27a$P.Value) > -log10(0.05), ]
# Is everything else really non-targets?
nontargets_27a <- full_gene_list_27a[!row.names(full_gene_list_27a) 
                                     %in% row.names(targets_27a), ]
# Get genbank accession numbers
genbankAcc_27a <- featDatTable$`GenBank Accession`
# Create a column for 1 being target and 0 being non target
targetCol_27a <- rep(0, nrow(full_gene_list_27a))
targetCol_27a[row.names(full_gene_list_27a) %in% row.names(targets_27a)] <- 1

accToTarget_27a <- cbind(genbankAcc_27a, targetCol_27a)
# Removing repeats and RNA with no acc numbers
accToTarget_27a <- accToTarget_27a[accToTarget_27a[,1] != "", ]
accToTarget_27a <- dealWithRepAcc(accToTarget_27a)

mir27a_featMat <- createFeatureMatrix(miR27a, accToTarget_27a[,1], 
                                      accToTarget_27a[,2])
freeEngFeats_27a1 <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]][1:1000], mir27a_featMat[[1]][1:1000,])
freeEngFeats_27a2 <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]][1001:2000], mir27a_featMat[[1]][1001:2000,])
freeEngFeats_27a3 <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]][2001:3000], mir27a_featMat[[1]][2001:3000,])
freeEngFeats_27a4 <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]][3001:4000], mir27a_featMat[[1]][3001:4000,])
freeEngFeats_27a5 <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]][4001:5000], mir27a_featMat[[1]][4001:5000,])
freeEngFeats_27a6 <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]][5001:6000], mir27a_featMat[[1]][5001:6000,])
freeEngFeats_27a7 <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]][6001:7000], mir27a_featMat[[1]][6001:7000,])
freeEngFeats_27a8 <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]][7001:8000], mir27a_featMat[[1]][7001:8000,])
freeEngFeats_27a9 <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]][8001:9000], mir27a_featMat[[1]][8001:9000,])
freeEngFeats_27a10 <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]][9001:10000], mir27a_featMat[[1]][9001:10000,])
freeEngFeats_27a11 <- siteDuplexFreeEnergy(
  miR27a, mir27a_featMat[[2]][10001:length(mir27a_featMat[[2]])], 
  mir27a_featMat[[1]][10001:length(mir27a_featMat[[2]]),])

# Combining the tables
freeEngFeats_23b <- rbind(freeEngFeats_23b1, freeEngFeats_23b2, 
                          freeEngFeats_23b3, freeEngFeats_23b4,
                          freeEngFeats_23b5, freeEngFeats_23b6,
                          freeEngFeats_23b7, freeEngFeats_23b8,
                          freeEngFeats_23b9, freeEngFeats_23b10,
                          freeEngFeats_23b11)
freeEngFeats_27a <- rbind(freeEngFeats_27a1, freeEngFeats_27a2, 
                          freeEngFeats_27a3, freeEngFeats_27a4,
                          freeEngFeats_27a5, freeEngFeats_27a6,
                          freeEngFeats_27a7, freeEngFeats_27a8,
                          freeEngFeats_27a9, freeEngFeats_27a10,
                          freeEngFeats_27a11)
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
  
  return(featMat)
}
featMatmir23b <- convertMissingValues(featMatmir23b)
featMatmir27a <- convertMissingValues(featMatmir27a)

# MORE DATA
cloonan23b27a <- getGEO("GSE40410", GSEMatrix =TRUE, AnnotGPL=TRUE)


