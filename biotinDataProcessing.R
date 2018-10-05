library(GEOquery)
library(limma)
source("CreateFeatureMatrix.R")

cloonan23b27a <- getGEO("GSE40410", GSEMatrix =TRUE, AnnotGPL=TRUE)
#comes up with a "13 parsing failures" warning
featDat.23b.27a <- featureData(cloonan23b27a$GSE40410_series_matrix.txt.gz)
featDatTable.23b.27a <- fData(cloonan23b27a$GSE40410_series_matrix.txt.gz)
assayDat.23b.27a <- assayData(cloonan23b27a$GSE40410_series_matrix.txt.gz)

measurelogFC <- function(assayDat, layout, cols){
  if(class(assayDat) == "environment"){
    assayDat <- assayDat$exprs
  }
  design<-model.matrix(~0+layout)
  colnames(design)<-c("control","pulldown")
  fit<-lmFit(assayDat[,cols], design)
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

accToTargMatRAW <- function(full_gene_list){
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
  genbankAcc <- full_gene_list$SEARCH_KEY[
    row.names(full_gene_list) %in% row.names(targets)
  ]
  genbankAcc <- c(genbankAcc, full_gene_list$SEARCH_KEY[
    row.names(full_gene_list) %in% row.names(nontargets)
  ])
  
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

# miR-96
geo.96 <- getGEO("GSE117103", GSEMatrix =TRUE, AnnotGPL=TRUE)
featDat.96 <- featureData(geo.96$GSE117103_series_matrix.txt.gz)
featDatTable.96 <- fData(geo.96$GSE117103_series_matrix.txt.gz)
assayDat.96 <- assayData(geo.96$GSE117103_series_matrix.txt.gz)
phenoDat.96 <- phenoData(geo.96$GSE117103_series_matrix.txt.gz)
pDatTable.96 <- pData(geo.96$GSE117103_series_matrix.txt.gz)
# Note these were trasfected with 30nM

layout.96 <- c("pulldown", "control", "pulldown", "control", "pulldown", "control")
cols.96 <- c(2,4,6,8,10,12)
fgl.96 <- measurelogFC(assayDat.96, layout.96, cols.96)
# That is a weird volcano plot


# RAW DATA
# miR 30e and miR 4536
# miR-4536
dat.dir <- "Data/30e-4536/"
CAKI.4536.array.labels <- readTargets(file = paste0(dat.dir,"Targets_CAKI_4536.txt"))
CAKI.4536.all.data<-read.ilmn(
  files=paste0(dat.dir, "sample_CAKI_4536.txt"), 
  annotation = c("TargetID", "SEARCH_KEY", "SYMBOL"),
  ctrlfiles=paste0(dat.dir, "control_CAKI_4536.txt"),probeid="PROBE_ID", 
  other.columns="Detection", sep="\t", verbose=TRUE
)
#examine the array controls
boxplot(CAKI.4536.all.data$E~CAKI.4536.all.data$genes$Status, log="y", las=2)

boxplot(log2(CAKI.4536.all.data$E), range=0, ylab="log 2 intensity", las=2)
#background correct and normalize, using the ERCC genes as the negative control
y.4536 <- neqc(CAKI.4536.all.data, negctrl="ERCC")
# Remove probes not expressed in at least 3 arrays
# (as 3 is the size of the pulldown/control groups of arrays)
expressed <- rowSums(y.4536$other$Detection < 0.05) >= 2
y.4536 <- y.4536[expressed,]
plotMDS(y.4536, labels = CAKI.4536.array.labels$Type)
plotMDS(y.4536, labels = CAKI.4536.array.labels$Sample.Group)
# control sample A doesn't look right
# recorrect without the bad sample A controls
y.4536 <- neqc(CAKI.4536.all.data[,-1], negctrl = "ERCC")
plotMDS(y.4536, labels = CAKI.4536.array.labels$Type[-1])
plotMDS(y.4536, labels = CAKI.4536.array.labels$Sample.Group[-1])
# Reremove probes not expressed in at least 2 arrays
# (as 2 is the smallest size of the pulldown/control groups of arrays)
expressed <- rowSums(y.4536$other$Detection < 0.05) >= 2
y.4536 <- y.4536[expressed,]

layout.CAKI.4536 <- CAKI.4536.array.labels$Type[-1]
cols.CAKI.4536 <- (1:5)
fgl.CAKI.4536 <- measurelogFC(y.4536, layout.CAKI.4536, cols.CAKI.4536) 
accToTarg.4536 <- accToTargMatRAW(fgl.CAKI.4536)
table(accToTarg.4536[,2])
# High target to non-target ratio. Should I use all targets when training?

#5p
miR4536.5p <- "UGUGGUAGAUAUAUGCACGAU"
sites.rnaSeq.featMat.4536 <- createFeatureMatrix(miR4536.5p, accToTarg.4536[,1],
                                                accToTarg.4536[,2])
freeEngFeats.4536 <- siteDuplexFreeEnergy(
  miR4536.5p, sites.rnaSeq.featMat.4536[[2]], sites.rnaSeq.featMat.4536[[1]])
freeEngFeats.4536$GenBank.Accession <- names(sites.rnaSeq.featMat.4536[[2]])
featMat.4536 <- merge(sites.rnaSeq.featMat.4536[[3]], freeEngFeats.4536,
                     by = "GenBank.Accession")
featMat.4536 <- convertMissingValues(featMat.4536)

# miR-30e
SKOV1.30e.array.labels <- readTargets(file = paste0(dat.dir,"Targets_RNA_DNA_30e.txt"))
SKOV1.30e.all.data<-read.ilmn(
  files=paste0(dat.dir, "sample_RNA_DNA_30e.txt"), 
  annotation = c("TargetID", "SEARCH_KEY", "SYMBOL"),
  ctrlfiles=paste0(dat.dir, "control_RNA_DNA_30e.txt"),probeid="PROBE_ID", 
  other.columns="Detection", sep="\t", verbose=TRUE
)
boxplot(log2(y.30e$E),range=0,ylab="log2 intensity", las=2)
#examine the array controls
boxplot(SKOV1.30e.all.data$E~SKOV1.30e.all.data$genes$Status, log="y", las=2)
#background correct and normalize, using the ERCC genes as the negative control
y.30e <- neqc(SKOV1.30e.all.data, negctrl="ERCC")
# Reremove probes not expressed in at least 4 arrays
# (as 4 is the smallest size of the pulldown/control groups of arrays)
expressed <- rowSums(y.30e$other$Detection < 0.05) >= 4
y.30e <- y.30e[expressed,]
plotMDS(y.30e, labels = SKOV1.30e.array.labels$Type)
plotMDS(y.30e, labels = SKOV1.30e.array.labels$Sample.Group)
plotMDS(y.30e, labels = SKOV1.30e.array.labels$Molecule)

layout.SKOV1.30e <- SKOV1.30e.array.labels$Type
cols.SKOV1.30e <- 1:10
fgl.SKOV1.30e <- measurelogFC(y.30e, layout.SKOV1.30e, cols.SKOV1.30e) 
accToTarg.30e <- accToTargMatRAW(fgl.SKOV1.30e)
table(accToTarg.30e[,2])
# Also large amount of targets to non targets,

# miR-155-5p and miR-584-5p
dat.dir <- "Data/miRs-155-584_SKOV3_PANC1/"
array.labels.155.584 <- readTargets(file = paste0(dat.dir,"Targets.txt"))
all.data.155.584<-read.ilmn(
  files=paste0(dat.dir, "june_sample_probe_profile.txt"), 
  annotation = c("TargetID", "SEARCH_KEY", "SYMBOL"),
  ctrlfiles=paste0(dat.dir, "june_control_probe_profile.txt"),probeid="PROBE_ID", 
  other.columns="Detection", sep="\t", verbose=TRUE
)
boxplot(all.data.155.584$E~all.data.155.584$genes$Status, las=2)
boxplot(log2(all.data.155.584$E),range=0,ylab="log2 intensity", las=2)
y.155.584 <- neqc(all.data.155.584, negctrl = "NEGATIVE")
boxplot(log2(y.155.584$E),range=0,ylab="log2 intensity", las=2)
plotMDS(y.155.584, labels = array.labels.155.584$Type)
pancCols.155 <- c(4,5,6,10,11,12)
skovCols.155 <- c(1,2,3,7,8,9)
pancCols.584 <- c(16,17,18,22,23,24)
skovCols.584 <- c(13,14,15,19,20,21)
# look at panc cell type
plotMDS(y.155.584[,c(pancCols.584, pancCols.155)], 
        labels = array.labels.155.584$Type[c(pancCols.584, pancCols.155)])
# Ooh a dodgy pull down
plotMDS(y.155.584[,c(pancCols.584, pancCols.155)], 
        labels = array.labels.155.584$Sample.Group[c(pancCols.584, pancCols.155)])
# 200162410075_L is dodgy (col 24)
pancCols.584 <- c(16,17,18,22,23)
plotMDS(y.155.584[,c(pancCols.584, pancCols.155)], 
        labels = array.labels.155.584$Type[c(pancCols.584, pancCols.155)])
# look at panc 155
plotMDS(y.155.584[,c(pancCols.155)], 
        labels = array.labels.155.584$Sample.Group[c(pancCols.155)])
plotMDS(y.155.584[,c(pancCols.155)], 
        labels = array.labels.155.584$Type[c(pancCols.155)])
# Controls clustering to gether so probably legit

# look at panc 584
plotMDS(y.155.584[,c(pancCols.584)], 
        labels = array.labels.155.584$Sample.Group[c(pancCols.584)])
plotMDS(y.155.584[,c(pancCols.584)], 
        labels = array.labels.155.584$Type[c(pancCols.584)])
# hmm me no likey

# look at skov
plotMDS(y.155.584[,c(skovCols.584, skovCols.155)], 
        labels = array.labels.155.584$Type[c(skovCols.584, skovCols.155)])
# 200162410075_G is dodgy (col 19). Though is it...? Maybe 75_H is cos it's on Controls
skovCols.584 <- c(13,14,15,19,21)
plotMDS(y.155.584[,c(skovCols.584, skovCols.155)], 
        labels = array.labels.155.584$miRNA[c(skovCols.584, skovCols.155)])
# look at skov 155
plotMDS(y.155.584[,c(skovCols.155)], 
        labels = array.labels.155.584$Type[c(skovCols.155)])
# hmmm ... Maybe try without I and C
plotMDS(y.155.584[,c(1,2,7,8)], 
        labels = array.labels.155.584$Type[c(1,2,7,8)])
# look at skov 584
plotMDS(y.155.584[,c(skovCols.584)], 
        labels = array.labels.155.584$Type[c(skovCols.584)])
skovCols.155 <- c(1,2,7,8)

# Use SKOV arrays
# exclude probes not expressed (pvalue >0.05) in at least 4 (num pulldowns) arrays
# 4 is the smallest size of the pulldown/control groups of arrays
expressed <- rowSums(
  y.155.584$other$Detection[,c(skovCols.155, skovCols.584)] < 0.05) >= 4
y.155.584 <- y.155.584[expressed,]
# miR155
layout.155 <- array.labels.155.584$Type[skovCols.155]
fgl.155 <- measurelogFC(y.155.584, layout.155, skovCols.155)
accToTarg.155 <- accToTargMatRAW(fgl.155)
table(accToTarg.155[,2])
# small sample!

# miR-584
layout.584 <- array.labels.155.584$Type[skovCols.584]
fgl.584 <- measurelogFC(y.155.584, layout.584, skovCols.584)
accToTarg.584 <- accToTargMatRAW(fgl.584)
table(accToTarg.584[,2])
# relatively ballanced for once

# Get features!

# Could not find which arm the 30e biotinylated sequence is from
# 5p is the dominant one according to: 
#   http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0000749
# so I'll use that for now
miR30e.5p <- "UGUAAACAUCCUUGACUGGAAG"
s.rS.fM.30e <- createFeatureMatrix(miR30e.5p, accToTarg.30e[,1], 
                                   accToTarg.30e[,2])
freeEngFeats.30e <- siteDuplexFreeEnergy(miR30e.5p, s.rS.fM.30e[[2]],
                                         s.rS.fM.30e[[1]])
freeEngFeats.30e$GenBank.Accession <- names(s.rS.fM.30e[[2]])
featMat.30e <- merge(s.rS.fM.30e[[3]], freeEngFeats.30e, by = "GenBank.Accession")
# There were two targets where no binding sites were found
# So maybe the 30e used was from the 3p branch?

# 155 feats
miR155.5p <- "UUAAUGCUAAUCGUGAUAGGGGUU"
s.rS.fM.155 <- createFeatureMatrix(miR155.5p, accToTarg.155[,1], 
                                   accToTarg.155[,2])
freeEngFeats.155 <- siteDuplexFreeEnergy(miR155.5p, s.rS.fM.155[[2]],
                                         s.rS.fM.155[[1]])
freeEngFeats.155$GenBank.Accession <- names(s.rS.fM.155[[2]])
featMat.155 <- merge(s.rS.fM.155[[3]], freeEngFeats.155, by = "GenBank.Accession")
featMat.155 <- convertMissingValues(featMat.155)

# 584 feats
miR584.5p <- "UUAUGGUUUGCCUGGGACUGAG"
s.rS.fM.584 <- createFeatureMatrix(miR584.5p, accToTarg.584[,1], 
                                   accToTarg.584[,2])
freeEngFeats.584 <- siteDuplexFreeEnergy(miR584.5p, s.rS.fM.584[[2]],
                                          s.rS.fM.584[[1]])
freeEngFeats.584$GenBank.Accession <- names(s.rS.fM.584[[2]])
featMat.584 <- merge(s.rS.fM.584[[3]], freeEngFeats.584, by = "GenBank.Accession")
featMat.584 <- convertMissingValues(featMat.584)

# miR-548d & miR-1289
dat.dir <- "Data/548d-1289/"
array.labels.548d.1289 <- readTargets(file = paste0(dat.dir,"Targets.txt"))
all.data.548d.1289<-read.ilmn(
  files=paste0(dat.dir, "sample probe profile.txt"), 
  annotation = c("TargetID", "SEARCH_KEY", "SYMBOL"),
  ctrlfiles=paste0(dat.dir, "control probe profile.txt"),probeid="PROBE_ID", 
  other.columns="Detection", sep="\t", verbose=TRUE
)
boxplot(all.data.548d.1289$E~all.data.548d.1289$genes$Status, las=2)
boxplot(log2(all.data.548d.1289$E)~all.data.548d.1289$genes$Status, las=2)
# A fair bit of exrpession in the negative controls at the levels at which 
# the bulk of the regular genes are expressed
# normalise with ERCC instead
y.548d.1289 <- neqc(all.data.548d.1289, negctrl = "ERCC")
boxplot(log2(y.548d.1289$E), las = 2)

# Check MDS plots
skovCols.548d <- 1:6
cakiCols.548d <- 7:12
skovCols.1289 <- 13:18
cakiCols.1289 <- 19:24
plotMDS(y.548d.1289, labels = array.labels.548d.1289$Type)
plotMDS(y.548d.1289[,skovCols.548d], 
        labels = array.labels.548d.1289$Type[skovCols.548d])
plotMDS(y.548d.1289[,cakiCols.548d], 
        labels = array.labels.548d.1289$Type[cakiCols.548d])
plotMDS(y.548d.1289[,skovCols.1289], 
        labels = array.labels.548d.1289$Type[skovCols.1289])
plotMDS(y.548d.1289[,cakiCols.1289], 
        labels = array.labels.548d.1289$Type[cakiCols.1289])
# So beautiful :')

layout.548d.skov <- array.labels.548d.1289$Type[skovCols.548d]
fgl.548d.skov <- measurelogFC(y.548d.1289, layout.548d.skov, skovCols.548d)
accToTarg.548d.skov <- accToTargMatRAW(fgl.548d.skov)
table(accToTarg.548d.skov[,2])
# That's alotta data!
# Wait whoops forgot to filter out lowly expressed genes
expressed <- rowSums(y.548d.1289$other$Detection[,skovCols.548d] < 0.05) >= 3
y.548d.skov <- y.548d.1289[expressed,skovCols.548d]
fgl.548d.skov <- measurelogFC(y.548d.skov, layout.548d.skov, 1:6)
accToTarg.548d.skov <- accToTargMatRAW(fgl.548d.skov)
table(accToTarg.548d.skov[,2])
# ah didn't change number of targets and non-targets anyway

# Filter out lowly expressed genes for 1289
expressed <- rowSums(y.548d.1289$other$Detection[,skovCols.1289] < 0.05) >= 3
y.1289.skov <- y.548d.1289[expressed, skovCols.1289]

layout.1289.skov <- array.labels.548d.1289$Type[skovCols.1289]
fgl.1289.skov <- measurelogFC(y.1289.skov, layout.1289.skov, 1:6)
accToTarg.1289.skov <- accToTargMatRAW(fgl.1289.skov)
table(accToTarg.1289.skov[,2])
#nice

# miR-548d feats
miR548d.3p <- "CAAAAACCACAGUUUCUUUUGC" 
s.rS.fM.548d <- createFeatureMatrix(miR548d.3p, accToTarg.548d.skov[,1],
                                    accToTarg.548d.skov[,2])
fEF.548d <- siteDuplexFreeEnergy(miR548d.3p, s.rS.fM.548d[[2]],
                                 s.rS.fM.548d[[1]])
fEF.548d$GenBank.Accession <- names(s.rS.fM.548d[[2]])
featMat.548d <- merge(s.rS.fM.548d[[3]], fEF.548d, by = "GenBank.Accession")
featMat.548d <- convertMissingValues(featMat.548d)

# miR-1289 feats
miR1289 <- "UGGAGUCCAGGAAUCUGCAUUUU"
s.rS.fM.1289 <- createFeatureMatrix(miR1289, accToTarg.1289.skov[,1],
                                    accToTarg.1289.skov[,2])
fEF.1289 <- siteDuplexFreeEnergy(miR1289, s.rS.fM.1289[[2]],
                                 s.rS.fM.1289[[1]])
fEF.1289$GenBank.Accession <- names(s.rS.fM.1289[[2]])
featMat.1289 <- merge(s.rS.fM.1289[[3]], fEF.1289, by = "GenBank.Accession")
featMat.1289 <- convertMissingValues(featMat.1289)
