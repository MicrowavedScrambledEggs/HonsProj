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

# miR-23b targets defined as below log2FoldChange threshold 1 and 
# above -log10pValue threshold of -log10(0.05) 
# Maybe need a better reason for these thresholds other than 
# "because that's what other bpd studies used"
targets_23b <- full_gene_list_23b[full_gene_list_23b$logFC < 1 
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
targets_27a <- full_gene_list_27a[full_gene_list_27a$logFC < 1 
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
cloonan3118_4307 <- getGEO("GSE40409", GSEMatrix =TRUE, AnnotGPL=TRUE)
featDat3118_4307 <- featureData(cloonan3118_4307$GSE40409_series_matrix.txt.gz)
featDatTable3118_4307 <- fData(cloonan3118_4307$GSE40409_series_matrix.txt.gz)
assayDat3118_4307 <- assayData(cloonan3118_4307$GSE40409_series_matrix.txt.gz)
phenoDat3118_4307 <- phenoData(cloonan3118_4307$GSE40409_series_matrix.txt.gz)
pDatTable3118_4307 <- pData(cloonan3118_4307$GSE40409_series_matrix.txt.gz)

miR3118 <- "UGUGACUGCAUUAUGAAAAUUCU"
# Volcano plot for miR-23b. Feeling like a script kidde here
layout <- c("pulldown", "control", "pulldown", "control", "pulldown", "control")
design<-model.matrix(~0+layout)
colnames(design)<-c("control","pulldown")
# Columns 1,2,5,6,9,10 are miR-3118
fit_3118<-lmFit(assayDat3118_4307$exprs[,c(1,2,5,6,9,10)], design)
cont.matrix.3118 <- makeContrasts(PullDownvsControl=pulldown-control, 
                                 levels=design)
fit2_3118 <- contrasts.fit(fit_3118, cont.matrix.3118)
fit2_3118 <- eBayes(fit2_3118)
full_gene_list_3118<-topTable(fit2_3118, coef=1, number=1000000, sort.by="logFC")
plot(full_gene_list_3118$logFC, -log10(full_gene_list_3118$P.Value), 
     xlab="log2 fold change", ylab="-log10 p-value")

miR4307 <- "AAUGUUUUUUCCUGUUUCC"
# Columns 3,4,7,8,11,12 are miR-4307
fit_4307<-lmFit(assayDat3118_4307$exprs[,c(3,4,7,8,11,12)], design)
cont.matrix.4307 <- makeContrasts(PullDownvsControl=pulldown-control, 
                                  levels=design)
fit2_4307 <- contrasts.fit(fit_4307, cont.matrix.4307)
fit2_4307 <- eBayes(fit2_4307)
full_gene_list_4307<-topTable(fit2_4307, coef=1, number=1000000, sort.by="logFC")
plot(full_gene_list_4307$logFC, -log10(full_gene_list_4307$P.Value), 
     xlab="log2 fold change", ylab="-log10 p-value")

targets_3118 <- full_gene_list_3118[full_gene_list_3118$logFC < 1 
                                  & -log10(full_gene_list_3118$P.Value) > -log10(0.05), ]
# Is everything else really non-targets?
nontargets_3118 <- full_gene_list_3118[!row.names(full_gene_list_3118) 
                                     %in% row.names(targets_3118), ]
# Get genbank accession numbers
genbankAcc_3118 <- featDatTable3118_4307$`GenBank Accession`
# Create a column for 1 being target and 0 being non target
targetCol_3118 <- rep(0, nrow(full_gene_list_3118))
targetCol_3118[row.names(full_gene_list_3118) %in% row.names(targets_3118)] <- 1

accToTarget_3118 <- cbind(genbankAcc_3118, targetCol_3118)
# Removing repeats and RNA with no acc numbers
accToTarget_3118 <- accToTarget_3118[accToTarget_3118[,1] != "", ]
accToTarget_3118 <- dealWithRepAcc(accToTarget_3118)

mir3118_featMat <- createFeatureMatrix(miR3118, accToTarget_3118[,1], 
                                      accToTarget_3118[,2])

# miR4307 targets and features
targets_4307 <- full_gene_list_4307[full_gene_list_4307$logFC < 1 
                                    & -log10(full_gene_list_4307$P.Value) > -log10(0.05), ]
# Is everything else really non-targets?
nontargets_4307 <- full_gene_list_4307[!row.names(full_gene_list_4307) 
                                       %in% row.names(targets_4307), ]
# Get genbank accession numbers
genbankAcc_4307 <- featDatTable3118_4307$`GenBank Accession`
# Create a column for 1 being target and 0 being non target
targetCol_4307 <- rep(0, nrow(full_gene_list_4307))
targetCol_4307[row.names(full_gene_list_4307) %in% row.names(targets_4307)] <- 1

accToTarget_4307 <- cbind(genbankAcc_4307, targetCol_4307)
# Removing repeats and RNA with no acc numbers
accToTarget_4307 <- accToTarget_4307[accToTarget_4307[,1] != "", ]
accToTarget_4307 <- dealWithRepAcc(accToTarget_4307)

mir4307_featMat <- createFeatureMatrix(miR4307, accToTarget_4307[,1], 
                                       accToTarget_4307[,2])

freeEngFeats_3118.1 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][1:100], mir3118_featMat[[1]][1:100,])
freeEngFeats_3118.2 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][101:300], mir3118_featMat[[1]][101:300,])
freeEngFeats_3118.3 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][301:600], mir3118_featMat[[1]][301:600,])
freeEngFeats_3118.4 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][601:1000], mir3118_featMat[[1]][601:1000,])
freeEngFeats_3118.5 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][1001:1500], mir3118_featMat[[1]][1001:1500,])
freeEngFeats_3118.6 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][1501:2100], mir3118_featMat[[1]][1501:2100,])
freeEngFeats_3118.7 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][2101:2800], mir3118_featMat[[1]][2101:2800,])
freeEngFeats_3118.8 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][2801:3600], mir3118_featMat[[1]][2801:3600,])
freeEngFeats_3118.9 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][3601:4500], mir3118_featMat[[1]][3601:4500,])
freeEngFeats_3118.10 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][4501:5500], mir3118_featMat[[1]][4501:5500,])
freeEngFeats_3118.11 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][5501:6500], mir3118_featMat[[1]][5501:6500,])
freeEngFeats_3118.12 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][6501:7500], mir3118_featMat[[1]][6501:7500,])
freeEngFeats_3118.13 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][7501:8500], mir3118_featMat[[1]][7501:8500,])
freeEngFeats_3118.14 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][8501:9500], mir3118_featMat[[1]][8501:9500,])
freeEngFeats_3118.15 <- siteDuplexFreeEnergy(
  miR3118, mir3118_featMat[[2]][9501:length(mir3118_featMat[[2]])], 
  mir3118_featMat[[1]][9501:length(mir3118_featMat[[2]]),])

freeEngFeats_3118 <- rbind(freeEngFeats_3118.1, freeEngFeats_3118.2, 
                           freeEngFeats_3118.3, freeEngFeats_3118.4,
                           freeEngFeats_3118.5, freeEngFeats_3118.6,
                           freeEngFeats_3118.7, freeEngFeats_3118.8,
                           freeEngFeats_3118.9, freeEngFeats_3118.10,
                           freeEngFeats_3118.11, freeEngFeats_3118.12,
                           freeEngFeats_3118.13, freeEngFeats_3118.14,
                           freeEngFeats_3118.15)

freeEngFeats_4307.1 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]][1:1000], mir4307_featMat[[1]][1:1000,])
freeEngFeats_4307.2 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]][1001:2000], mir4307_featMat[[1]][1001:2000,])
freeEngFeats_4307.3 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]][2001:3000], mir4307_featMat[[1]][2001:3000,])
freeEngFeats_4307.4 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]][3001:4000], mir4307_featMat[[1]][3001:4000,])
freeEngFeats_4307.5 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]][4001:5000], mir4307_featMat[[1]][4001:5000,])
freeEngFeats_4307.6 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]][5001:6000], mir4307_featMat[[1]][5001:6000,])
freeEngFeats_4307.7 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]][6001:7000], mir4307_featMat[[1]][6001:7000,])
freeEngFeats_4307.8 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]][7001:8000], mir4307_featMat[[1]][7001:8000,])
freeEngFeats_4307.9 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]][8001:9000], mir4307_featMat[[1]][8001:9000,])
freeEngFeats_4307.10 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]][9001:10000], mir4307_featMat[[1]][9001:10000,])
freeEngFeats_4307.11 <- siteDuplexFreeEnergy(
  miR4307, mir4307_featMat[[2]][10001:length(mir4307_featMat[[2]])], 
  mir4307_featMat[[1]][10001:length(mir4307_featMat[[2]]),])

freeEngFeats_4307 <- rbind(freeEngFeats_4307.1, freeEngFeats_4307.2, 
                           freeEngFeats_4307.3, freeEngFeats_4307.4,
                           freeEngFeats_4307.5, freeEngFeats_4307.6,
                           freeEngFeats_4307.7, freeEngFeats_4307.8,
                           freeEngFeats_4307.9, freeEngFeats_4307.10,
                           freeEngFeats_4307.11)

freeEngFeats_3118$GenBank.Accession <- names(mir3118_featMat[[2]])
freeEngFeats_4307$GenBank.Accession <- names(mir4307_featMat[[2]])

featMatmir3118 <- merge(mir3118_matrix, freeEngFeats_3118, 
                        by = "GenBank.Accession")
featMatmir3118 <- convertMissingValues(featMatmir3118)
featMatmir4307 <- merge(mir4307_matrix, freeEngFeats_4307, 
                        by = "GenBank.Accession")
featMatmir4307 <- convertMissingValues(featMatmir4307)

fullFeatMat <- rbind(featMatmir23b, featMatmir27a, featMatmir3118, 
                     featMatmir4307)
fullFeatMatNoAcc <- fullFeatMat[, !colnames(fullFeatMat) %in% c("GenBank.Accession")]
