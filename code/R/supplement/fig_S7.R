
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

# Subset metabolomics - untreated vs mock infected comparison
mock_metabolome <- subset(metabolome, infection == 'mock')
mock_metabolome$infection <- NULL
glycine_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('glycine')))]
colnames(glycine_mock) <- c('abx', 'substrate')
glycine_mock$abx <- factor(glycine_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
hydroxyproline_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('trans-4-hydroxyproline')))]
colnames(hydroxyproline_mock) <- c('abx', 'substrate')
hydroxyproline_mock$abx <- factor(hydroxyproline_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))

# Subset metabolites - mock vs infected comparison
glycine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('glycine')))]
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

#-------------------------------------------------------------------------------------------------------------------------------------#

# glycine - Untreated vs Mock
p.adjust(c(wilcox.test(subset(glycine_mock, abx=='none')[,2], subset(glycine_mock, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_mock, abx=='none')[,2], subset(glycine_mock, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_mock, abx=='none')[,2], subset(glycine_mock, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_mock, abx=='none')[,2], subset(glycine_mock, abx=='germfree')[,2], exact=F)$p.value), method='BH')
# hydroxyproline - Untreated vs Mock
p.adjust(c(wilcox.test(subset(hydroxyproline_mock, abx=='none')[,2], subset(hydroxyproline_mock, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_mock, abx=='none')[,2], subset(hydroxyproline_mock, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_mock, abx=='none')[,2], subset(hydroxyproline_mock, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_mock, abx=='none')[,2], subset(hydroxyproline_mock, abx=='germfree')[,2], exact=F)$p.value), method='BH')

# glycine - Mock vs Infected
p.adjust(c(wilcox.test(subset(glycine_strep, infection=='630')[,2], subset(glycine_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_cef, infection=='630')[,2], subset(glycine_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_clinda, infection=='630')[,2], subset(glycine_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine_gf, infection=='630')[,2], subset(glycine_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
# hydroxyproline - Mock vs Infected
p.adjust(c(wilcox.test(subset(hydroxyproline_strep, infection=='630')[,2], subset(hydroxyproline_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_cef, infection=='630')[,2], subset(hydroxyproline_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_clinda, infection=='630')[,2], subset(hydroxyproline_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(hydroxyproline_gf, infection=='630')[,2], subset(hydroxyproline_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up mock infection data for plotting
glycine_mock <- subset(glycine_mock, abx == 'none')
hydroxyproline_mock <- subset(hydroxyproline_mock, abx == 'none')

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S7.pdf'
pdf(file=plot_file, width=7, height=5)
layout(matrix(c(1,
                2), nrow=2, ncol=1, byrow=TRUE))
par(mgp=c(2.3,0.7,0), xpd=FALSE, las=1, mar=c(2,4,1.5,1))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Glycine
stripchart(substrate~abx, data=glycine_mock, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,3.5), xlim=c(0.5,13.5),
           cex=1.2, ylab='Scaled Intesity', method='jitter', jitter=0.15)
stripchart(substrate~infection, data=glycine_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,3.5), xlim=c(0.5,13.5),
           cex=1.2, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=glycine_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,3.5), xlim=c(0.5,13.5),
           cex=1.2, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=glycine_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,3.5), xlim=c(0.5,13.5),
           cex=1.2, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=glycine_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,3.5), xlim=c(0.5,13.5),
           cex=1.2, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
axis(side=2, at=seq(0,3.5,0.5), labels=c('0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5'))
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.3, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.3, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.8)
mtext('A', side=2, line=2, las=2, adj=2, padj=-5, cex=1.3)
legend('topright', legend='Glycine', pt.cex=0, bty='n')
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(glycine_mock[,2]),
              median(subset(glycine_strep, infection=='mock')[,2]), median(subset(glycine_strep, infection=='630')[,2]),
              median(subset(glycine_cef, infection=='mock')[,2]), median(subset(glycine_cef, infection=='630')[,2]),
              median(subset(glycine_clinda, infection=='mock')[,2]), median(subset(glycine_clinda, infection=='630')[,2]),
              median(subset(glycine_gf, infection=='mock')[,2]), median(subset(glycine_gf, infection=='630')[,2])), 
         y1=c(median(glycine_mock[,2]),
              median(subset(glycine_strep, infection=='mock')[,2]), median(subset(glycine_strep, infection=='630')[,2]),
              median(subset(glycine_cef, infection=='mock')[,2]), median(subset(glycine_cef, infection=='630')[,2]),
              median(subset(glycine_clinda, infection=='mock')[,2]), median(subset(glycine_clinda, infection=='630')[,2]),
              median(subset(glycine_gf, infection=='mock')[,2]), median(subset(glycine_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=c(6), y0=c(2.8), x1=c(7), y1=c(2.8), lwd=2)
text(x=c(6.5), y=c(3), '*', font=2, cex=1.5)
mtext(c('*','*'), side=3, adj=c(0.43,0.863), padj=0.5, font=2, cex=1.7, col='gray40') # Untreated vs Mock significance

#-----------------#

# Hydroxyproline
stripchart(substrate~abx, data=hydroxyproline_mock, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,5), xlim=c(0.5,13.5),
           cex=1.2, ylab='Scaled Intesity', method='jitter', jitter=0.15)
stripchart(substrate~infection, data=hydroxyproline_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,5), xlim=c(0.5,13.5),
           cex=1.2, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=hydroxyproline_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,5), xlim=c(0.5,13.5),
           cex=1.2, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=hydroxyproline_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,5), xlim=c(0.5,13.5),
           cex=1.2, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
stripchart(substrate~infection, data=hydroxyproline_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,5), xlim=c(0.5,13.5),
           cex=1.2, ylab='Scaled Intensity', method='jitter', jitter=0.15, add=TRUE)
axis(side=2, at=seq(0,5,1), labels=c('0','1','2','3','4','5'))
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.3, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.3, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.8)
mtext('B', side=2, line=2, las=2, adj=2, padj=-5, cex=1.3)
legend('topright', legend='Hydroxyproline', pt.cex=0, bty='n')
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(hydroxyproline_mock[,2]),
              median(subset(hydroxyproline_strep, infection=='mock')[,2]), median(subset(hydroxyproline_strep, infection=='630')[,2]),
              median(subset(hydroxyproline_cef, infection=='mock')[,2]), median(subset(hydroxyproline_cef, infection=='630')[,2]),
              median(subset(hydroxyproline_clinda, infection=='mock')[,2]), median(subset(hydroxyproline_clinda, infection=='630')[,2]),
              median(subset(hydroxyproline_gf, infection=='mock')[,2]), median(subset(hydroxyproline_gf, infection=='630')[,2])), 
         y1=c(median(hydroxyproline_mock[,2]),
              median(subset(hydroxyproline_strep, infection=='mock')[,2]), median(subset(hydroxyproline_strep, infection=='630')[,2]),
              median(subset(hydroxyproline_cef, infection=='mock')[,2]), median(subset(hydroxyproline_cef, infection=='630')[,2]),
              median(subset(hydroxyproline_clinda, infection=='mock')[,2]), median(subset(hydroxyproline_clinda, infection=='630')[,2]),
              median(subset(hydroxyproline_gf, infection=='mock')[,2]), median(subset(hydroxyproline_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=c(3,6,9,12), y0=c(4.1,2.9,2.9,2), x1=c(4,7,10,13), y1=c(4.1,2.9,2.9,2), lwd=2)
text(x=c(3.5,6.5,9.5,12.5), y=c(4.4,3.2,3.2,2.3), '*', font=2, cex=1.5)
mtext(c('*','*','*','*'), side=3, adj=c(0.21,0.43,0.645,0.863), padj=0.5, font=2, cex=1.7, col='gray40') # Untreated vs Mock significance

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
detach('package:wesanderson', character.only = TRUE)
rm(list=ls())
gc()
