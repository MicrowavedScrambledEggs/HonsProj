library(Biostrings)
library(seqTools)
library(seqinr)

miR23b <- "AUCACAUUGCCAGGGAUUACC"
miR27a <- "UUCACAGUGGCUAAGUUCCGC"
miR3118 <- "UGUGACUGCAUUAUGAAAAUUCU"
miR4307 <- "AAUGUUUUUUCCUGUUUCC"
miR199a3p <- "ACAGUAGUCUGCACAUUGGUUA"
miR424.3p <- "CAAAACGUGAGGCGCUGCUAU"
miR4536.5p <- "UGUGGUAGAUAUAUGCACGAU"
miR30e.5p <- "UGUAAACAUCCUUGACUGGAAG"
miR155.5p <- "UUAAUGCUAAUCGUGAUAGGGGUU"
miR584.5p <- "UUAUGGUUUGCCUGGGACUGAG"
miR548d.3p <- "CAAAAACCACAGUUUCUUUUGC"
miR1289 <- "UGGAGUCCAGGAAUCUGCAUUUU"
miR10a <- "UACCCUGUAGAUCCGAAUUUGUG"
miR182 <- "UUUGGCAAUGGUAGAACUCACACU"

all.Mature <- readRNAStringSet("Data/mature.fa/mature.fa", "fasta")
all.Mature.hsa <- all.Mature[which(grepl("hsa-",names(all.Mature)))]

par(mfrow = c(3,1), mar = c(3,4,2,1))

hsa.Con <- t(consensusMatrix(all.Mature.hsa, baseOnly = TRUE, as.prob = TRUE))
matplot(hsa.Con[1:25,-5], type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
        ylim = c(0.0, 0.6))
legend(legend = colnames(hsa.Con)[-5],"top",col=1:4, lty=1:4, lwd=2, horiz = T)

all.Con <- t(consensusMatrix(all.Mature, baseOnly = TRUE, as.prob = TRUE))
matplot(all.Con[1:25,-5], type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
        ylim = c(0.0, 0.6))
legend(legend = colnames(all.Con)[-5],"top",col=1:4, lty=1:4, lwd=2)

miRNA.used <- c(miR10a, miR1289, miR155.5p, miR182, miR23b, miR27a, miR199a3p,
                miR30e.5p, miR3118, miR424.3p, miR4307, miR4536.5p, miR548d.3p,
                miR584.5p)
names(miRNA.used) <- c("miR10a", "miR1289", "miR155.5p", 'miR182', 'miR23b', "miR27a", 'miR199a3p',
                       'miR30e.5p', 'miR3118', 'miR424.3p', 'miR4307', 'miR4536.5p', 'miR548d.3p',
                       'miR584.5p')
miRNA.used <- RNAStringSet(miRNA.used)

miRNA.CV <- c(miR10a,  miR1289, miR155.5p, miR182, miR23b, miR199a3p,
              miR30e.5p, miR424.3p, miR4536.5p, miR548d.3p)
miRNA.CV <- RNAStringSet(miRNA.CV)

cv.Con <- t(consensusMatrix(miRNA.CV, baseOnly = TRUE, as.prob = TRUE))
matplot(cv.Con[1:22,-5], type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
        xlim = c(1, 25))
legend(legend = colnames(cv.Con)[-5],"topright",col=1:4, lty=1:4, lwd=2)

used.Con <- t(consensusMatrix(miRNA.used, baseOnly = TRUE, as.prob = TRUE))
matplot(used.Con[1:22,-5], type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
        main = "Nucleotide frequency at each position for all miRNA used")
legend(legend = colnames(used.Con)[-5],"topright",col=1:4, lty=1:4, lwd=2)

avTable <- rbind(colMeans(cv.Con[1:22,]), colMeans(hsa.Con[1:22,]))
avTable <- rbind(avTable, colMeans(all.Con[1:22,]))

avTableSeed <- rbind(colMeans(cv.Con[2:8,]), colMeans(hsa.Con[2:8,]))
avTableSeed <- rbind(avTableSeed, colMeans(all.Con[2:8,]))

avTable10P <- rbind(colMeans(cv.Con[10:22,]), colMeans(hsa.Con[10:22,]))
avTable10P <- rbind(avTable10P, colMeans(all.Con[10:22,]))

avTable1 <- rbind(cv.Con[1,], hsa.Con[1,])
avTable1 <- rbind(avTable1, all.Con[1,])

write.csv(round(avTable, 3), "nucAv.csv")
write.csv(round(avTableSeed, 3), "nucAvSeed.csv")
write.csv(round(avTable10P, 3), "nucAvNotSeed.csv")
write.csv(round(avTable1, 3), "nucAv1st.csv")

# Looks like to get close to the distribution, can not go lower than 40, so as to get the 2.5%
# differences in frequencies seen accross the positions.

# sample 10 sets of 40 hsa-miRs and plot their freq distribution agains the human freq dist
for(i in 1:10){
  samplmiRNA <- all.Mature.hsa[sample.int(length(all.Mature.hsa),400)]
  sampl.con <- t(consensusMatrix(samplmiRNA, baseOnly = TRUE, as.prob = TRUE))
  comp.con <- cbind(hsa.Con[1:22,1:4],sampl.con[1:22,1:4])
  colnames(comp.con) <- c("A HSA", "C HSA", "G HSA", "U HSA", "A Sample", "C Sample", "G Sample", "U Sample")
  matplot(comp.con, type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
          main = "Nucleotide frequency: Full human vs Sample", col = c(1:4,1:4),
          lty = c(1,1,1,1,3,3,3,3))
  legend(legend = colnames(comp.con),"top",col = c(1:4,1:4), lty = c(1,1,1,1,3,3,3,3), lwd=2, horiz = T, cex = 0.5)
}

# Multiple alignment stuff

seed.seqs <- sapply(miRNA.used, subseq, start = 1, end = 8)
centered.seqs <- sapply(miRNA.used, subseq, start = 3, end = 15)
names(seed.seqs) <- paste0(names(miRNA.used), "Seed")
names(centered.seqs) <- paste0(names(miRNA.used), "Centered")
to.File <- c(sapply(seed.seqs, as.character), sapply(centered.seqs, as.character))
to.File <- as.list(to.File)
names(to.File) <- c(names(seed.seqs), names(centered.seqs))
write.fasta(to.File, names = names(to.File), file.out = "seedAndCentered.fa")

miRNA.names.split <- strsplit(names(all.Mature), " ")
miRNA.Species <- sapply(miRNA.names.split, function(x) paste(x[3], x[4]))
table(miRNA.Species)
