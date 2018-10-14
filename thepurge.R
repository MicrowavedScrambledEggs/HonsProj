# Applying the stricter definition of sufficiently expressed to the raw data, 
# then filtering out from the datasets the RNA that do not meet this defiiniton

# miR-4536
y.4536 <- neqc(CAKI.4536.all.data, negctrl="ERCC")
# Remove probes not expressed in all arrays of at least one group of control 
# or pulldown
expressed.cont <- rowSums(y.4536[,1:3]$other$Detection < 0.05) >= 3
expressed.pull <- rowSums(y.4536[,4:6]$other$Detection < 0.05) >= 3
expressed <- expressed.cont | expressed.pull
y.4536 <- y.4536[expressed,]
plotMDS(y.4536, labels = CAKI.4536.array.labels$Type)
plotMDS(y.4536, labels = CAKI.4536.array.labels$Sample.Group)
# control sample A looks like it's miles away
# recorrect without the bad sample A controls
y.4536 <- neqc(CAKI.4536.all.data[,-1], negctrl = "ERCC")
expressed.cont <- rowSums(y.4536[,1:2]$other$Detection < 0.05) >= 2
expressed <- expressed.cont | expressed.pull
y.4536 <- y.4536[expressed,]
plotMDS(y.4536, labels = CAKI.4536.array.labels$Type[-1])
plotMDS(y.4536, labels = CAKI.4536.array.labels$Sample.Group[-1])

# Save the old target set so know which RNA to remove
accToTarg.4536.old <- accToTarg.4536
# Remeasure logFC
layout.CAKI.4536 <- CAKI.4536.array.labels$Type[-1]
cols.CAKI.4536 <- (1:5)
fgl.CAKI.4536 <- measurelogFC(y.4536, layout.CAKI.4536, cols.CAKI.4536) 
accToTarg.4536 <- accToTargMatRAW(fgl.CAKI.4536)
# look for new RNA (should not be any?)
new.4536 <- accToTarg.4536[
  !(accToTarg.4536[,1] %in% accToTarg.4536.old[,1]), ]
new.4536
# There's two. Probably not worth extracting features for
# Find the RNA to purge from featMat
purge.4536 <- accToTarg.4536.old[
  !(accToTarg.4536.old[,1] %in% accToTarg.4536[,1]),1]
featMat.4536 <- featMat.4536[!(featMat.4536$GenBank.Accession %in% purge.4536),]

# miR 155 i got more data so will have to train a bunch of new feats anyway
# miR 584 I don't trust the data so I'll probably won't use anymore

# miR 548d I accidentally wrote over the old accToTarg.548d.SKOV
# But I can use the genbank from the featmat
# See if there are any new RNA
missing.548d.SKOV <- accToTarg.548d.skov[
  !(accToTarg.548d.skov[,1] %in% featMat.548d$GenBank.Accession), ]
# Ok there's a bunch. Probably targets that had no detected sites. Will find features
# for just in case, and I'll be having to do for new RNA from CAKI anyway
# Find the RNA to purge
purge.548d <- featMat.548d$GenBank.Accession[
  !(featMat.548d$GenBank.Accession %in% accToTarg.548d.skov[,1])] 
featMat.548d <- featMat.548d[!(featMat.548d$GenBank.Accession %in% purge.548d),]

# miR-1289
accToTarg.1289.skov.old <- accToTarg.1289.skov
# look for new RNA
new.1289 <- accToTarg.1289.skov[
  !(accToTarg.1289.skov[,1] %in% accToTarg.1289.skov.old[,1]), ]
new.1289
# only 3, so not worth training feats for
purge.1289 <- accToTarg.1289.skov.old[
  !(accToTarg.1289.skov.old[,1] %in% accToTarg.1289.skov[,1]),1
]
featMat.1289 <- featMat.1289[!(featMat.1289$GenBank.Accession %in% purge.1289), ]  

# Learn new feats for new pairs
s.rS.fM.155 <- createFeatureMatrix(miR155.5p, accToTarg.155[,1], 
                                   accToTarg.155[,2])
freeEngFeats.155 <- siteDuplexFreeEnergy(miR155.5p, s.rS.fM.155[[2]],
                                         s.rS.fM.155[[1]])
freeEngFeats.155$GenBank.Accession <- names(s.rS.fM.155[[2]])
featMat.155 <- merge(s.rS.fM.155[[3]], freeEngFeats.155, by = "GenBank.Accession")
featMat.155 <- convertMissingValues(featMat.155)

# new 548d pair feats
new.548d <- accToTarg.548d[!(accToTarg.548d[,1] %in% featMat.548d$GenBank.Accession), ]
s.rS.fM.548d.new <- createFeatureMatrix(miR548d.3p, new.548d[,1], new.548d[,2])
fEF.548d.new <- siteDuplexFreeEnergy(miR548d.3p, s.rS.fM.548d.new[[2]],
                                 s.rS.fM.548d.new[[1]])
fEF.548d.new$GenBank.Accession <- names(s.rS.fM.548d.new[[2]])
featMat.548d.new <- merge(s.rS.fM.548d.new[[3]], fEF.548d.new, by = "GenBank.Accession")
featMat.548d.new <- convertMissingValues(featMat.548d.new)
featMat.548d <- rbind(featMat.548d, featMat.548d.new)

# new 1289 pair feats
new.1289 <- accToTarg.1289[!(accToTarg.1289[,1] %in% featMat.1289$GenBank.Accession), ]
s.rS.fM.1289.new <- createFeatureMatrix(miR1289, new.1289[,1], new.1289[,2])
fEF.1289.new <- siteDuplexFreeEnergy(miR1289, s.rS.fM.1289.new[[2]],
                                     s.rS.fM.1289.new[[1]])
fEF.1289.new$GenBank.Accession <- names(s.rS.fM.1289.new[[2]])
featMat.1289.new <- merge(s.rS.fM.1289.new[[3]], fEF.1289.new, by = "GenBank.Accession")
featMat.1289.new <- convertMissingValues(featMat.1289.new)
featMat.1289 <- rbind(featMat.1289, featMat.1289.new)

