
plot.targ.nontarg.diff <- function(fullGeneList){
  cuts <- seq(0.0,1.0,0.1)
  diffs <- c()
  for(cut in cuts) {
    targs <- nrow(fullGeneList[fullGeneList$logFC > cut & fullGeneList$P.Value < 0.05, ])
    nontargs <- nrow(fullGeneList[fullGeneList$logFC < -cut & fullGeneList$P.Value < 0.05, ])
    diffs <- c(diffs, targs - nontargs)
  }
  plot(cuts, diffs, xlab = "Symetric logFC cutoff", ylab = "#targets - #nontargets")
}

plot.targ.nontarg.diff(full_gene_list_3118)
plot.targ.nontarg.diff(full_gene_list_4307)



