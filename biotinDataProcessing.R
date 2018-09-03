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
