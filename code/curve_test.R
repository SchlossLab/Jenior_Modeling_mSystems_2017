

# Load dependencies
deps <- c('vegan', 'igraph', 'ggplot2', 'shape', 'wesanderson', 'matrixStats');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}


growth <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/formatted_growth.tsv', sep='\t', header=TRUE, row.names=1)
growth <- as.data.frame(t(growth))

y_gluc_y_aa <- cbind(growth$B3, growth$B4, growth$B5) - growth$B2
y_gluc_y_aa[y_gluc_y_aa < 0] <- 0
y_gluc_n_aa <- cbind(growth$C3, growth$C4, growth$C5) - growth$C2
y_gluc_n_aa[y_gluc_n_aa < 0] <- 0
n_gluc_y_aa <- cbind(growth$D3, growth$D4, growth$D5) - growth$D2
n_gluc_y_aa[n_gluc_y_aa < 0] <- 0
bhi <- cbind(growth$F3, growth$F4, growth$F5) - growth$F2
bhi[bhi < 0] <- 0
rm(growth)

y_gluc_y_aa_median <- apply(y_gluc_y_aa, 1, median)
y_gluc_y_aa_sd <- rowSds(y_gluc_y_aa)
y_gluc_n_aa_median <- apply(y_gluc_n_aa, 1, median)
y_gluc_n_aa_sd <- rowSds(y_gluc_n_aa)
n_gluc_y_aa_median <- apply(n_gluc_y_aa, 1, median)
n_gluc_y_aa_sd <- rowSds(n_gluc_y_aa)
bhi_median <- apply(bhi, 1, median)
bhi_sd <- rowSds(bhi)
rm(y_gluc_y_aa, y_gluc_n_aa, n_gluc_y_aa, bhi)



pdf(file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/growth_test.pdf', width=9, height=8)
par(mar=c(7,7,1.5,2), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE, mgp=c(4,2,0))
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,50), ylim=c(-0.03,1.0), lwd=2, pch=15, xlab='Time (hours)', ylab=expression(OD[600]), cex=2.3)
abline(h=seq(0,0.9,0.1), lty=3, col='gray68') # adding gridlines
abline(v=seq(1,50,2), lty=3, col='gray68') # adding gridlines
axis(1, at=seq(1,49,4), labels=seq(0,24,2), tck=-0.018)
axis(2, at=seq(0.0,1.0,0.2), labels=c('0.0','0.2','0.4','0.6','0.8','1.0'), tck=-0.018)
lines(y_gluc_y_aa_median, type='o', col='black', lwd=2.5, pch=16, cex=3)
lines(y_gluc_n_aa_median, type='o', col='blue', lwd=2.5, pch=16, cex=3)
lines(n_gluc_y_aa_median, type='o', col='red', lwd=2.5, pch=16, cex=3)
lines(bhi_median, type='o', col='green', lwd=2.5, pch=16, cex=3)
legend('topleft', legend=c('y_gluc_y_aa','y_gluc_n_aa','n_gluc_y_aa_median','bhi'), 
       col=c('black','blue','red','green'), pch=16, cex=1.7, pt.cex=2, bg='white', lwd=2.5)

dev.off()





