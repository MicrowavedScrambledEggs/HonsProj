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
fit<-lmFit(assayDat$exprs[,c(1:4)], design)
cont.matrix <- makeContrasts(PullDownvsControl=pulldown-control, 
                             levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
full_gene_list<-topTable(fit2, coef=1, number=1000000, sort.by="logFC")
plot(full_gene_list$logFC, -log10(full_gene_list$P.Value), 
     xlab="log2 fold change", ylab="-log10 p-value")

# miR-23b targets defined as below log2FoldChange threshold -1 and 
# above -log10pValue threshold of -log10(0.05) 
# Maybe need a better reason for these thresholds other than 
# "because that's what other bpd studies used"
targets <- full_gene_list[full_gene_list$logFC < -1 
                          & -log10(full_gene_list$P.Value) > -log10(0.05), ]
# Is everything else really non-targets?
nontargets <- full_gene_list[!row.names(full_gene_list) 
                             %in% row.names(targets), ]
# Get genbank accession numbers
genbankAcc <- featDatTable$`GenBank Accession`
# Create a column for 1 being target and 0 being non target
targetCol <- rep(0, nrow(full_gene_list))
targetCol[row.names(full_gene_list) %in% row.names(targets)] <- 1

accToTarget <- cbind(genbankAcc, targetCol)
# Removing repeats and RNA with no acc numbers
accToTarget <- accToTarget[accToTarget[,1] != "", ]
accToTarget <- unique(accToTarget)

mir23b_featMat <- createFeatureMatrix(miR23b, accToTarget[,1], accToTarget[,2])
freeEngFeats <- siteDuplexFreeEnergy(miR23b, mir23b_featMat[[2]], 
                                     mir23b_featMat[[1]])

# startTime <- Sys.time()
# tst <- siteDuplexFreeEnergy(miR23b, test_featMat[[2]], 
#                      test_featMat[[1]])
# Sys.time() - startTime

# Combining the tables
featMatmir23b <- cbind(mir23b_featMat[[3]], freeEngFeats)

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
nucFreqCols <- colnames(featMatmir23b)[
  grepl("c\\d_freq", colnames(featMatmir23b))]
featMatmir23b[,nucFreqCols] <- t(
  apply(featMatmir23b[,nucFreqCols], 1, 
        function(x) {x[is.nan(x)] <- 0; return(x)}))
# NaN 3' end distances
endDistCols <- colnames(featMatmir23b)[
  grepl("c\\d_end", colnames(featMatmir23b))]
featMatmir23b[,endDistCols] <- t(
  apply(featMatmir23b[,endDistCols], 1, 
        function(x) {x[is.nan(x)] <- -1; return(x)}))
# NaN stop codon distances
stopDistCols <- colnames(featMatmir23b)[
  grepl("c\\d_stop", colnames(featMatmir23b))]
featMatmir23b[,stopDistCols] <- t(
  apply(featMatmir23b[,stopDistCols], 1, 
        function(x) {x[is.nan(x)] <- -.Machine$double.xmax; return(x)}))
# NaN duplex free energy
dupFreeEngCols <- colnames(featMatmir23b)[
  grepl("dup", colnames(featMatmir23b))]
featMatmir23b[,dupFreeEngCols] <- t(
  apply(featMatmir23b[,dupFreeEngCols], 1, 
        function(x) {x[is.nan(x)] <- .Machine$double.xmax; return(x)}))
# NA variance
varCols <- colnames(featMatmir23b)[
  grepl("Var", colnames(featMatmir23b))]
featMatmir23b[,varCols] <- t(
  apply(featMatmir23b[,varCols], 1, 
        function(x) {x[is.na(x)] <- 0; return(x)}))
# NA stop codon distances
featMatmir23b[,stopDistCols] <- t(
  apply(featMatmir23b[,stopDistCols], 1, 
        function(x) {x[is.na(x)] <- .Machine$double.xmax; return(x)}))
# NA region frequencies
regionCols <- colnames(featMatmir23b)[
  grepl("c\\d_(5'UTR|CDS|3'UTR)", colnames(featMatmir23b))]
featMatmir23b[,regionCols] <- t(
  apply(featMatmir23b[,regionCols], 1, 
        function(x) {x[is.na(x)] <- 0; return(x)}))

# Random Forest does not like column names with symbols in them
colnames(featMatmir23b) <- gsub("'", "", colnames(featMatmir23b))
colnames(featMatmir23b) <- gsub("-", "neg", colnames(featMatmir23b))
featMatmir23b$Target <- as.factor(featMatmir23b$Target)

library(randomForest)
miTarget_randomForest <- randomForest(Target ~ ., data = featMatmir23b)
