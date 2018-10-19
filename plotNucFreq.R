library(Biostrings)
library(seqTools)

all.Mature <- readRNAStringSet("Data/mature.fa/mature.fa", "fasta")
all.Mature.hsa <- all.Mature[which(grepl("hsa-",names(all.Mature)))]

hsa.Con <- t(consensusMatrix(all.Mature.hsa, baseOnly = TRUE, as.prob = TRUE))
matplot(hsa.Con[1:25,-5], type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
        main = "Nucleotide frequency at each position for all known human miRNA")
legend(legend = colnames(hsa.Con)[-5],"topright",col=1:4, lty=1:4, lwd=2)

all.Con <- t(consensusMatrix(all.Mature, baseOnly = TRUE, as.prob = TRUE))
matplot(all.Con[1:25,-5], type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
        main = "Nucleotide frequency at each position for all known miRNA")
legend(legend = colnames(all.Con)[-5],"topright",col=1:4, lty=1:4, lwd=2)

miRNA.used <- c(miR10a, miR1289, miR155.5p, miR182, miR23b, miR27a, miR199a3p,
                miR30e.5p, miR3118, miR424.3p, miR4307, miR4536.5p, miR548d.3p,
                miR584.5p)
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
  samplmiRNA <- all.Mature.hsa[sample.int(length(all.Mature.hsa),40)]
  sampl.con <- t(consensusMatrix(samplmiRNA, baseOnly = TRUE, as.prob = TRUE))
  matplot(used.Con[1:22,-5], type="l", lwd=2, xlab="Sequence Position", ylab= "Base frequency",
          main = "Nucleotide frequency at each position for all miRNA used", col = c(1:4,1:4))
  legend(legend = colnames(used.Con)[-5],"topright",col=1:4, lty=1:4, lwd=2)
}



