library(Biostrings)
library(seqTools)
library(seqinr)

all.Mature <- readRNAStringSet("Data/mature.fa/mature.fa", "fasta")
all.Mature.hsa <- all.Mature[which(grepl("hsa-",names(all.Mature)))]

hsa.Con <- t(consensusMatrix(all.Mature.hsa, baseOnly = TRUE, as.prob = TRUE))
matplot(hsa.Con[1:25,-5], type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
        main = "Nucleotide frequency at each position for all known human miRNA")
legend(legend = colnames(hsa.Con)[-5],"top",col=1:4, lty=1:4, lwd=2, horiz = T)

all.Con <- t(consensusMatrix(all.Mature, baseOnly = TRUE, as.prob = TRUE))
matplot(all.Con[1:25,-5], type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
        main = "Nucleotide frequency at each position for all known miRNA")
legend(legend = colnames(all.Con)[-5],"topright",col=1:4, lty=1:4, lwd=2)

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

used.Con <- t(consensusMatrix(miRNA.used, baseOnly = TRUE, as.prob = TRUE))
matplot(used.Con[1:22,-5], type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
        main = "Nucleotide frequency at each position for all miRNA used")
legend(legend = colnames(used.Con)[-5],"topright",col=1:4, lty=1:4, lwd=2)

cv.Con <- t(consensusMatrix(miRNA.CV, baseOnly = TRUE, as.prob = TRUE))
matplot(cv.Con[1:22,-5], type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
        main = "Nucleotide frequency at each position for all miRNA used in 10 fold CV")
legend(legend = colnames(cv.Con)[-5],"topright",col=1:4, lty=1:4, lwd=2)

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
