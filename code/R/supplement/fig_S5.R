
# Start with blank slate
rm(list=ls())
gc()

# Load dependency
if ('wesanderson' %in% installed.packages()[,"Package"] == FALSE){
  install.packages(as.character('wesanderson'), quiet=TRUE);
}
library('wesanderson', verbose=FALSE, character.only=TRUE)


# Select files
metabolome <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/metabolomics.scaled_intensities.tsv'
metadata <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'

# Read in data
metabolome <- read.delim(metabolome, sep='\t', header=T, row.names=1)
metabolome <- metabolome[, !colnames(metabolome) %in% c('CefC5M2','StrepC4M1')] # Remove possible contamination
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination

# Merge metabolomics with metadata
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$type <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome <- as.data.frame(t(metabolome))
metabolome <- merge(metadata, metabolome, by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
rm(metadata)

# Subset metabolites - mock vs infected comparison
glycine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('glycine')))]
glycine_untreated <- subset(glycine, abx == 'none')
glycine_untreated$abx <- NULL
colnames(glycine_untreated) <- c('infection', 'substrate')
glycine_untreated$infection <- factor(glycine_untreated$infection, levels=c('mock','630'))
glycine_cef <- subset(glycine, abx == 'cefoperazone')
glycine_cef$abx <- NULL
colnames(glycine_cef) <- c('infection', 'substrate')
glycine_cef$infection <- factor(glycine_cef$infection, levels=c('mock','630'))
glycine_strep <- subset(glycine, abx == 'streptomycin')
glycine_strep$abx <- NULL
colnames(glycine_strep) <- c('infection', 'substrate')
glycine_strep$infection <- factor(glycine_strep$infection, levels=c('mock','630'))
glycine_clinda <- subset(glycine, abx == 'clindamycin')
glycine_clinda$abx <- NULL
colnames(glycine_clinda) <- c('infection', 'substrate')
glycine_clinda$infection <- factor(glycine_clinda$infection, levels=c('mock','630'))
glycine_gf <- subset(glycine, abx == 'germfree')
glycine_gf$abx <- NULL
colnames(glycine_gf) <- c('infection', 'substrate')
glycine_gf$infection <- factor(glycine_gf$infection, levels=c('mock','630'))
rm(glycine)
hydroxyproline <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('trans-4-hydroxyproline')))]
hydroxyproline_untreated <- subset(hydroxyproline, abx == 'none')
hydroxyproline_untreated$abx <- NULL
colnames(hydroxyproline_untreated) <- c('infection', 'substrate')
hydroxyproline_untreated$infection <- factor(hydroxyproline_untreated$infection, levels=c('mock','630'))
hydroxyproline_cef <- subset(hydroxyproline, abx == 'cefoperazone')
hydroxyproline_cef$abx <- NULL
colnames(hydroxyproline_cef) <- c('infection', 'substrate')
hydroxyproline_cef$infection <- factor(hydroxyproline_cef$infection, levels=c('mock','630'))
hydroxyproline_strep <- subset(hydroxyproline, abx == 'streptomycin')
hydroxyproline_strep$abx <- NULL
colnames(hydroxyproline_strep) <- c('infection', 'substrate')
hydroxyproline_strep$infection <- factor(hydroxyproline_strep$infection, levels=c('mock','630'))
hydroxyproline_clinda <- subset(hydroxyproline, abx == 'clindamycin')
hydroxyproline_clinda$abx <- NULL
colnames(hydroxyproline_clinda) <- c('infection', 'substrate')
hydroxyproline_clinda$infection <- factor(hydroxyproline_clinda$infection, levels=c('mock','630'))
hydroxyproline_gf <- subset(hydroxyproline, abx == 'germfree')
hydroxyproline_gf$abx <- NULL
colnames(hydroxyproline_gf) <- c('infection', 'substrate')
hydroxyproline_gf$infection <- factor(hydroxyproline_gf$infection, levels=c('mock','630'))
rm(hydroxyproline)
leucine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('leucine')))]
leucine_untreated <- subset(leucine, abx == 'none')
leucine_untreated$abx <- NULL
colnames(leucine_untreated) <- c('infection', 'substrate')
leucine_untreated$infection <- factor(leucine_untreated$infection, levels=c('mock','630'))
leucine_cef <- subset(leucine, abx == 'cefoperazone')
leucine_cef$abx <- NULL
colnames(leucine_cef) <- c('infection', 'substrate')
leucine_cef$infection <- factor(leucine_cef$infection, levels=c('mock','630'))
leucine_strep <- subset(leucine, abx == 'streptomycin')
leucine_strep$abx <- NULL
colnames(leucine_strep) <- c('infection', 'substrate')
leucine_strep$infection <- factor(leucine_strep$infection, levels=c('mock','630'))
leucine_clinda <- subset(leucine, abx == 'clindamycin')
leucine_clinda$abx <- NULL
colnames(leucine_clinda) <- c('infection', 'substrate')
leucine_clinda$infection <- factor(leucine_clinda$infection, levels=c('mock','630'))
leucine_gf <- subset(leucine, abx == 'germfree')
leucine_gf$abx <- NULL
colnames(leucine_gf) <- c('infection', 'substrate')
leucine_gf$infection <- factor(leucine_gf$infection, levels=c('mock','630'))
rm(leucine)
isoleucine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('isoleucine')))]
isoleucine_untreated <- subset(isoleucine, abx == 'none')
isoleucine_untreated$abx <- NULL
colnames(isoleucine_untreated) <- c('infection', 'substrate')
isoleucine_untreated$infection <- factor(isoleucine_untreated$infection, levels=c('mock','630'))
isoleucine_cef <- subset(isoleucine, abx == 'cefoperazone')
isoleucine_cef$abx <- NULL
colnames(isoleucine_cef) <- c('infection', 'substrate')
isoleucine_cef$infection <- factor(isoleucine_cef$infection, levels=c('mock','630'))
isoleucine_strep <- subset(isoleucine, abx == 'streptomycin')
isoleucine_strep$abx <- NULL
colnames(isoleucine_strep) <- c('infection', 'substrate')
isoleucine_strep$infection <- factor(isoleucine_strep$infection, levels=c('mock','630'))
isoleucine_clinda <- subset(isoleucine, abx == 'clindamycin')
isoleucine_clinda$abx <- NULL
colnames(isoleucine_clinda) <- c('infection', 'substrate')
isoleucine_clinda$infection <- factor(isoleucine_clinda$infection, levels=c('mock','630'))
isoleucine_gf <- subset(isoleucine, abx == 'germfree')
isoleucine_gf$abx <- NULL
colnames(isoleucine_gf) <- c('infection', 'substrate')
isoleucine_gf$infection <- factor(isoleucine_gf$infection, levels=c('mock','630'))
rm(isoleucine)
rm(metabolome)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Untreated vs Mock
p.adjust(c(wilcox.test(subset(glycine_untreated, infection=='mock')[,2], subset(glycine_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_untreated, infection=='mock')[,2], subset(glycine_strep, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_untreated, infection=='mock')[,2], subset(glycine_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_untreated, infection=='mock')[,2], subset(glycine_cef, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_untreated, infection=='mock')[,2], subset(glycine_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_untreated, infection=='mock')[,2], subset(glycine_clinda, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_untreated, infection=='mock')[,2], subset(glycine_gf, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_untreated, infection=='mock')[,2], subset(glycine_gf, infection=='630')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(hydroxyproline_untreated, infection=='mock')[,2], subset(hydroxyproline_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_untreated, infection=='mock')[,2], subset(hydroxyproline_strep, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_untreated, infection=='mock')[,2], subset(hydroxyproline_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_untreated, infection=='mock')[,2], subset(hydroxyproline_cef, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_untreated, infection=='mock')[,2], subset(hydroxyproline_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_untreated, infection=='mock')[,2], subset(hydroxyproline_clinda, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_untreated, infection=='mock')[,2], subset(hydroxyproline_gf, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_untreated, infection=='mock')[,2], subset(hydroxyproline_gf, infection=='630')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(isoleucine_untreated, infection=='mock')[,2], subset(isoleucine_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(isoleucine_untreated, infection=='mock')[,2], subset(isoleucine_strep, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(isoleucine_untreated, infection=='mock')[,2], subset(isoleucine_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(isoleucine_untreated, infection=='mock')[,2], subset(isoleucine_cef, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(isoleucine_untreated, infection=='mock')[,2], subset(isoleucine_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(isoleucine_untreated, infection=='mock')[,2], subset(isoleucine_clinda, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(isoleucine_untreated, infection=='mock')[,2], subset(isoleucine_gf, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(isoleucine_untreated, infection=='mock')[,2], subset(isoleucine_gf, infection=='630')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(leucine_untreated, infection=='mock')[,2], subset(leucine_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(leucine_untreated, infection=='mock')[,2], subset(leucine_strep, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(leucine_untreated, infection=='mock')[,2], subset(leucine_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(leucine_untreated, infection=='mock')[,2], subset(leucine_cef, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(leucine_untreated, infection=='mock')[,2], subset(leucine_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(leucine_untreated, infection=='mock')[,2], subset(leucine_clinda, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(leucine_untreated, infection=='mock')[,2], subset(leucine_gf, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(leucine_untreated, infection=='mock')[,2], subset(leucine_gf, infection=='630')[,2], exact=F)$p.value), method='BH')

# Mock vs Infected
p.adjust(c(wilcox.test(subset(glycine_strep, infection=='630')[,2], subset(glycine_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_cef, infection=='630')[,2], subset(glycine_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_clinda, infection=='630')[,2], subset(glycine_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_gf, infection=='630')[,2], subset(glycine_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(hydroxyproline_strep, infection=='630')[,2], subset(hydroxyproline_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_cef, infection=='630')[,2], subset(hydroxyproline_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_clinda, infection=='630')[,2], subset(hydroxyproline_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_gf, infection=='630')[,2], subset(hydroxyproline_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(isoleucine_strep, infection=='630')[,2], subset(isoleucine_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(isoleucine_cef, infection=='630')[,2], subset(isoleucine_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(isoleucine_clinda, infection=='630')[,2], subset(isoleucine_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(isoleucine_gf, infection=='630')[,2], subset(isoleucine_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(leucine_strep, infection=='630')[,2], subset(leucine_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(leucine_cef, infection=='630')[,2], subset(leucine_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(leucine_clinda, infection=='630')[,2], subset(leucine_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(leucine_gf, infection=='630')[,2], subset(leucine_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S5.pdf'
pdf(file=plot_file, width=7, height=9)
layout(matrix(c(1,
                2,
                3,
                4), nrow=4, ncol=1, byrow=TRUE))
par(mgp=c(2.3,0.7,0), xpd=FALSE, las=1, mar=c(3,4,1,1))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Leucine
stripchart(substrate~infection, data=leucine_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15)
stripchart(substrate~infection, data=leucine_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=leucine_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=leucine_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=leucine_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
axis(side=2, at=seq(0,4,1), labels=c('0','1','2','3','4'))
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.3, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.3, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.8)
mtext('A', side=2, line=2, las=2, adj=1, padj=-5, cex=1.3)
legend('topright', legend='Leucine', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(leucine_untreated[,2]),
              median(subset(leucine_strep, infection=='mock')[,2]), median(subset(leucine_strep, infection=='630')[,2]),
              median(subset(leucine_cef, infection=='mock')[,2]), median(subset(leucine_cef, infection=='630')[,2]),
              median(subset(leucine_clinda, infection=='mock')[,2]), median(subset(leucine_clinda, infection=='630')[,2]),
              median(subset(leucine_gf, infection=='mock')[,2]), median(subset(leucine_gf, infection=='630')[,2])), 
         y1=c(median(leucine_untreated[,2]),
              median(subset(leucine_strep, infection=='mock')[,2]), median(subset(leucine_strep, infection=='630')[,2]),
              median(subset(leucine_cef, infection=='mock')[,2]), median(subset(leucine_cef, infection=='630')[,2]),
              median(subset(leucine_clinda, infection=='mock')[,2]), median(subset(leucine_clinda, infection=='630')[,2]),
              median(subset(leucine_gf, infection=='mock')[,2]), median(subset(leucine_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=c(6,9,12), y0=c(2.9,2.9,2), x1=c(7,10,13), y1=c(2.9,2.9,2), lwd=2)
text(x=c(6.5,9.5,12.5), y=c(3.2,3.2,2.3), '*', font=2, cex=1.5)
legend('topleft', legend='Donor', pt.cex=0, bty='n')
mtext(rep('*',7), side=3, adj=c(0.21,0.28,
                                0.5,
                                0.645,0.715,
                                0.863,0.933), padj=0.4, font=2, cex=1.3, col='gray40') # Untreated vs Mock significance

#-----------------#

# Isoleucine
stripchart(substrate~infection, data=isoleucine_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15)
stripchart(substrate~infection, data=isoleucine_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=isoleucine_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=isoleucine_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=isoleucine_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
axis(side=2, at=seq(0,4,1), labels=c('0','1','2','3','4'))
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.3, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.3, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.8)
mtext('B', side=2, line=2, las=2, adj=1, padj=-5, cex=1.3)
legend('topright', legend='Isoleucine', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(isoleucine_untreated[,2]),
              median(subset(isoleucine_strep, infection=='mock')[,2]), median(subset(isoleucine_strep, infection=='630')[,2]),
              median(subset(isoleucine_cef, infection=='mock')[,2]), median(subset(isoleucine_cef, infection=='630')[,2]),
              median(subset(isoleucine_clinda, infection=='mock')[,2]), median(subset(isoleucine_clinda, infection=='630')[,2]),
              median(subset(isoleucine_gf, infection=='mock')[,2]), median(subset(isoleucine_gf, infection=='630')[,2])), 
         y1=c(median(isoleucine_untreated[,2]),
              median(subset(isoleucine_strep, infection=='mock')[,2]), median(subset(isoleucine_strep, infection=='630')[,2]),
              median(subset(isoleucine_cef, infection=='mock')[,2]), median(subset(isoleucine_cef, infection=='630')[,2]),
              median(subset(isoleucine_clinda, infection=='mock')[,2]), median(subset(isoleucine_clinda, infection=='630')[,2]),
              median(subset(isoleucine_gf, infection=='mock')[,2]), median(subset(isoleucine_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=c(6,9,12), y0=c(2.9,2.9,2), x1=c(7,10,13), y1=c(2.9,2.9,2), lwd=2)
text(x=c(6.5,9.5,12.5), y=c(3.2,3.2,2.3), '*', font=2, cex=1.5)
legend('topleft', legend='Donor', pt.cex=0, bty='n')
mtext(rep('*',6), side=3, adj=c(0.28,
                                0.5,
                                0.645,0.715,
                                0.863,0.933), padj=0.4, font=2, cex=1.3, col='gray40') # Untreated vs Mock significance

#-----------------#

# Hydroxyproline
stripchart(substrate~infection, data=hydroxyproline_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,5), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15)
stripchart(substrate~infection, data=hydroxyproline_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,5), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=hydroxyproline_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,5), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=hydroxyproline_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,5), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=hydroxyproline_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,5), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
axis(side=2, at=seq(0,5,1), labels=c('0','1','2','3','4','5'))
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.3, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.3, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.8)
mtext('C', side=2, line=2, las=2, adj=1, padj=-5, cex=1.3)
legend('topright', legend='Hydroxyproline', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(hydroxyproline_untreated[,2]),
              median(subset(hydroxyproline_strep, infection=='mock')[,2]), median(subset(hydroxyproline_strep, infection=='630')[,2]),
              median(subset(hydroxyproline_cef, infection=='mock')[,2]), median(subset(hydroxyproline_cef, infection=='630')[,2]),
              median(subset(hydroxyproline_clinda, infection=='mock')[,2]), median(subset(hydroxyproline_clinda, infection=='630')[,2]),
              median(subset(hydroxyproline_gf, infection=='mock')[,2]), median(subset(hydroxyproline_gf, infection=='630')[,2])), 
         y1=c(median(hydroxyproline_untreated[,2]),
              median(subset(hydroxyproline_strep, infection=='mock')[,2]), median(subset(hydroxyproline_strep, infection=='630')[,2]),
              median(subset(hydroxyproline_cef, infection=='mock')[,2]), median(subset(hydroxyproline_cef, infection=='630')[,2]),
              median(subset(hydroxyproline_clinda, infection=='mock')[,2]), median(subset(hydroxyproline_clinda, infection=='630')[,2]),
              median(subset(hydroxyproline_gf, infection=='mock')[,2]), median(subset(hydroxyproline_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=c(3,6,9,12), y0=c(4.1,2.9,2.9,2), x1=c(4,7,10,13), y1=c(4.1,2.9,2.9,2), lwd=2)
text(x=c(3.5,6.5,9.5,12.5), y=c(4.4,3.2,3.2,2.3), '*', font=2, cex=1.5)
legend('topleft', legend='Acceptor', pt.cex=0, bty='n')
mtext(rep('*',5), side=3, adj=c(0.21,
                                0.43,
                                0.645,
                                0.863,0.933), padj=0.4, font=2, cex=1.3, col='gray40') # Untreated vs Mock significance

#-----------------#

# Glycine
stripchart(substrate~infection, data=glycine_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,3.5), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15)
stripchart(substrate~infection, data=glycine_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,3.5), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=glycine_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,3.5), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=glycine_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,3.5), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=glycine_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,3.5), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
axis(side=2, at=seq(0,3.5,0.5), labels=c('0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'))
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.3, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.3, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.8)
mtext('D', side=2, line=2, las=2, adj=1, padj=-5, cex=1.3)
legend('topright', legend='Glycine', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(glycine_untreated[,2]),
              median(subset(glycine_strep, infection=='mock')[,2]), median(subset(glycine_strep, infection=='630')[,2]),
              median(subset(glycine_cef, infection=='mock')[,2]), median(subset(glycine_cef, infection=='630')[,2]),
              median(subset(glycine_clinda, infection=='mock')[,2]), median(subset(glycine_clinda, infection=='630')[,2]),
              median(subset(glycine_gf, infection=='mock')[,2]), median(subset(glycine_gf, infection=='630')[,2])), 
         y1=c(median(glycine_untreated[,2]),
              median(subset(glycine_strep, infection=='mock')[,2]), median(subset(glycine_strep, infection=='630')[,2]),
              median(subset(glycine_cef, infection=='mock')[,2]), median(subset(glycine_cef, infection=='630')[,2]),
              median(subset(glycine_clinda, infection=='mock')[,2]), median(subset(glycine_clinda, infection=='630')[,2]),
              median(subset(glycine_gf, infection=='mock')[,2]), median(subset(glycine_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=c(6), y0=c(2.8), x1=c(7), y1=c(2.8), lwd=2)
text(x=c(6.5), y=c(3), '*', font=2, cex=1.5)
legend('topleft', legend='Acceptor', pt.cex=0, bty='n')
mtext(rep('*',4), side=3, adj=c(0.43,0.5,
                                0.863,0.933), padj=0.4, font=2, cex=1.3, col='gray40') # Untreated vs Mock significance

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
detach('package:wesanderson', character.only = TRUE)
rm(list=ls())
gc()

