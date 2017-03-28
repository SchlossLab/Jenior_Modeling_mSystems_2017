
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

# Subset metabolites
acetylglucosamine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-acetylglucosamine/N-acetylgalactosamine')))]
acetylglucosamine_untreated <- subset(acetylglucosamine, abx == 'none')
acetylglucosamine_untreated$abx <- NULL
colnames(acetylglucosamine_untreated) <- c('infection', 'substrate')
acetylglucosamine_untreated$infection <- factor(acetylglucosamine_untreated$infection, levels=c('mock','630'))
acetylglucosamine_cef <- subset(acetylglucosamine, abx == 'cefoperazone')
acetylglucosamine_cef$abx <- NULL
colnames(acetylglucosamine_cef) <- c('infection', 'substrate')
acetylglucosamine_cef$infection <- factor(acetylglucosamine_cef$infection, levels=c('mock','630'))
acetylglucosamine_strep <- subset(acetylglucosamine, abx == 'streptomycin')
acetylglucosamine_strep$abx <- NULL
colnames(acetylglucosamine_strep) <- c('infection', 'substrate')
acetylglucosamine_strep$infection <- factor(acetylglucosamine_strep$infection, levels=c('mock','630'))
acetylglucosamine_clinda <- subset(acetylglucosamine, abx == 'clindamycin')
acetylglucosamine_clinda$abx <- NULL
colnames(acetylglucosamine_clinda) <- c('infection', 'substrate')
acetylglucosamine_clinda$infection <- factor(acetylglucosamine_clinda$infection, levels=c('mock','630'))
acetylglucosamine_gf <- subset(acetylglucosamine, abx == 'germfree')
acetylglucosamine_gf$abx <- NULL
colnames(acetylglucosamine_gf) <- c('infection', 'substrate')
acetylglucosamine_gf$infection <- factor(acetylglucosamine_gf$infection, levels=c('mock','630'))
rm(acetylglucosamine)
proline <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('proline')))]
proline_untreated <- subset(proline, abx == 'none')
proline_untreated$abx <- NULL
colnames(proline_untreated) <- c('infection', 'substrate')
proline_untreated$infection <- factor(proline_untreated$infection, levels=c('mock','630'))
proline_cef <- subset(proline, abx == 'cefoperazone')
proline_cef$abx <- NULL
colnames(proline_cef) <- c('infection', 'substrate')
proline_cef$infection <- factor(proline_cef$infection, levels=c('mock','630'))
proline_strep <- subset(proline, abx == 'streptomycin')
proline_strep$abx <- NULL
colnames(proline_strep) <- c('infection', 'substrate')
proline_strep$infection <- factor(proline_strep$infection, levels=c('mock','630'))
proline_clinda <- subset(proline, abx == 'clindamycin')
proline_clinda$abx <- NULL
colnames(proline_clinda) <- c('infection', 'substrate')
proline_clinda$infection <- factor(proline_clinda$infection, levels=c('mock','630'))
proline_gf <- subset(proline, abx == 'germfree')
proline_gf$abx <- NULL
colnames(proline_gf) <- c('infection', 'substrate')
proline_gf$infection <- factor(proline_gf$infection, levels=c('mock','630'))
rm(proline)
mannitolsorbitol <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('mannitol/sorbitol')))]
mannitolsorbitol_untreated <- subset(mannitolsorbitol, abx == 'none')
mannitolsorbitol_untreated$abx <- NULL
colnames(mannitolsorbitol_untreated) <- c('infection', 'substrate')
mannitolsorbitol_untreated$infection <- factor(mannitolsorbitol_untreated$infection, levels=c('mock','630'))
mannitolsorbitol_cef <- subset(mannitolsorbitol, abx == 'cefoperazone')
mannitolsorbitol_cef$abx <- NULL
colnames(mannitolsorbitol_cef) <- c('infection', 'substrate')
mannitolsorbitol_cef$infection <- factor(mannitolsorbitol_cef$infection, levels=c('mock','630'))
mannitolsorbitol_strep <- subset(mannitolsorbitol, abx == 'streptomycin')
mannitolsorbitol_strep$abx <- NULL
colnames(mannitolsorbitol_strep) <- c('infection', 'substrate')
mannitolsorbitol_strep$infection <- factor(mannitolsorbitol_strep$infection, levels=c('mock','630'))
mannitolsorbitol_clinda <- subset(mannitolsorbitol, abx == 'clindamycin')
mannitolsorbitol_clinda$abx <- NULL
colnames(mannitolsorbitol_clinda) <- c('infection', 'substrate')
mannitolsorbitol_clinda$infection <- factor(mannitolsorbitol_clinda$infection, levels=c('mock','630'))
mannitolsorbitol_gf <- subset(mannitolsorbitol, abx == 'germfree')
mannitolsorbitol_gf$abx <- NULL
colnames(mannitolsorbitol_gf) <- c('infection', 'substrate')
mannitolsorbitol_gf$infection <- factor(mannitolsorbitol_gf$infection, levels=c('mock','630'))
rm(mannitolsorbitol)
acetylneuraminate <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-acetylneuraminate')))]
acetylneuraminate_untreated <- subset(acetylneuraminate, abx == 'none')
acetylneuraminate_untreated$abx <- NULL
colnames(acetylneuraminate_untreated) <- c('infection', 'substrate')
acetylneuraminate_untreated$infection <- factor(acetylneuraminate_untreated$infection, levels=c('mock','630'))
acetylneuraminate_cef <- subset(acetylneuraminate, abx == 'cefoperazone')
acetylneuraminate_cef$abx <- NULL
colnames(acetylneuraminate_cef) <- c('infection', 'substrate')
acetylneuraminate_cef$infection <- factor(acetylneuraminate_cef$infection, levels=c('mock','630'))
acetylneuraminate_strep <- subset(acetylneuraminate, abx == 'streptomycin')
acetylneuraminate_strep$abx <- NULL
colnames(acetylneuraminate_strep) <- c('infection', 'substrate')
acetylneuraminate_strep$infection <- factor(acetylneuraminate_strep$infection, levels=c('mock','630'))
acetylneuraminate_clinda <- subset(acetylneuraminate, abx == 'clindamycin')
acetylneuraminate_clinda$abx <- NULL
colnames(acetylneuraminate_clinda) <- c('infection', 'substrate')
acetylneuraminate_clinda$infection <- factor(acetylneuraminate_clinda$infection, levels=c('mock','630'))
acetylneuraminate_gf <- subset(acetylneuraminate, abx == 'germfree')
acetylneuraminate_gf$abx <- NULL
colnames(acetylneuraminate_gf) <- c('infection', 'substrate')
acetylneuraminate_gf$infection <- factor(acetylneuraminate_gf$infection, levels=c('mock','630'))
rm(acetylneuraminate)
succinate <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('succinate')))]
succinate_untreated <- subset(succinate, abx == 'none')
succinate_untreated$abx <- NULL
colnames(succinate_untreated) <- c('infection', 'substrate')
succinate_untreated$infection <- factor(succinate_untreated$infection, levels=c('mock','630'))
succinate_strep <- subset(succinate, abx == 'streptomycin')
succinate_strep$abx <- NULL
colnames(succinate_strep) <- c('infection', 'substrate')
succinate_strep$infection <- factor(succinate_strep$infection, levels=c('mock','630'))
succinate_cef <- subset(succinate, abx == 'cefoperazone')
succinate_cef$abx <- NULL
colnames(succinate_cef) <- c('infection', 'substrate')
succinate_cef$infection <- factor(succinate_cef$infection, levels=c('mock','630'))
succinate_clinda <- subset(succinate, abx == 'clindamycin')
succinate_clinda$abx <- NULL
colnames(succinate_clinda) <- c('infection', 'substrate')
succinate_clinda$infection <- factor(succinate_clinda$infection, levels=c('mock','630'))
succinate_gf <- subset(succinate, abx == 'germfree')
succinate_gf$abx <- NULL
colnames(succinate_gf) <- c('infection', 'substrate')
succinate_gf$infection <- factor(succinate_gf$infection, levels=c('mock','630'))
rm(succinate)
rm(metabolome)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Calculate significant differences

# Untreated vs Everything
p.adjust(c(wilcox.test(subset(acetylglucosamine_untreated, infection=='mock')[,2], subset(acetylglucosamine_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_untreated, infection=='mock')[,2], subset(acetylglucosamine_strep, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_untreated, infection=='mock')[,2], subset(acetylglucosamine_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_untreated, infection=='mock')[,2], subset(acetylglucosamine_cef, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_untreated, infection=='mock')[,2], subset(acetylglucosamine_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_untreated, infection=='mock')[,2], subset(acetylglucosamine_clinda, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_untreated, infection=='mock')[,2], subset(acetylglucosamine_gf, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_untreated, infection=='mock')[,2], subset(acetylglucosamine_gf, infection=='630')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(proline_untreated, infection=='mock')[,2], subset(proline_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_untreated, infection=='mock')[,2], subset(proline_strep, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_untreated, infection=='mock')[,2], subset(proline_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_untreated, infection=='mock')[,2], subset(proline_cef, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_untreated, infection=='mock')[,2], subset(proline_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_untreated, infection=='mock')[,2], subset(proline_clinda, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_untreated, infection=='mock')[,2], subset(proline_gf, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_untreated, infection=='mock')[,2], subset(proline_gf, infection=='630')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(mannitolsorbitol_untreated, infection=='mock')[,2], subset(mannitolsorbitol_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_untreated, infection=='mock')[,2], subset(mannitolsorbitol_strep, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_untreated, infection=='mock')[,2], subset(mannitolsorbitol_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_untreated, infection=='mock')[,2], subset(mannitolsorbitol_cef, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_untreated, infection=='mock')[,2], subset(mannitolsorbitol_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_untreated, infection=='mock')[,2], subset(mannitolsorbitol_clinda, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_untreated, infection=='mock')[,2], subset(mannitolsorbitol_gf, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_untreated, infection=='mock')[,2], subset(mannitolsorbitol_gf, infection=='630')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(succinate_untreated, infection=='mock')[,2], subset(succinate_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(succinate_untreated, infection=='mock')[,2], subset(succinate_strep, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(succinate_untreated, infection=='mock')[,2], subset(succinate_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(succinate_untreated, infection=='mock')[,2], subset(succinate_cef, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(succinate_untreated, infection=='mock')[,2], subset(succinate_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(succinate_untreated, infection=='mock')[,2], subset(succinate_clinda, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(succinate_untreated, infection=='mock')[,2], subset(succinate_gf, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(succinate_untreated, infection=='mock')[,2], subset(succinate_gf, infection=='630')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(acetylneuraminate_untreated, infection=='mock')[,2], subset(acetylneuraminate_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_untreated, infection=='mock')[,2], subset(acetylneuraminate_strep, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_untreated, infection=='mock')[,2], subset(acetylneuraminate_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_untreated, infection=='mock')[,2], subset(acetylneuraminate_cef, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_untreated, infection=='mock')[,2], subset(acetylneuraminate_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_untreated, infection=='mock')[,2], subset(acetylneuraminate_clinda, infection=='630')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_untreated, infection=='mock')[,2], subset(acetylneuraminate_gf, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_untreated, infection=='mock')[,2], subset(acetylneuraminate_gf, infection=='630')[,2], exact=F)$p.value), method='BH')

#------------------#

# Mock vs Infected
p.adjust(c(wilcox.test(subset(acetylglucosamine_strep, infection=='630')[,2], subset(acetylglucosamine_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_cef, infection=='630')[,2], subset(acetylglucosamine_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_clinda, infection=='630')[,2], subset(acetylglucosamine_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_gf, infection=='630')[,2], subset(acetylglucosamine_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(proline_strep, infection=='630')[,2], subset(proline_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_cef, infection=='630')[,2], subset(proline_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_clinda, infection=='630')[,2], subset(proline_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_gf, infection=='630')[,2], subset(proline_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(mannitolsorbitol_strep, infection=='630')[,2], subset(mannitolsorbitol_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_cef, infection=='630')[,2], subset(mannitolsorbitol_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_clinda, infection=='630')[,2], subset(mannitolsorbitol_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_gf, infection=='630')[,2], subset(mannitolsorbitol_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(acetylneuraminate_strep, infection=='630')[,2], subset(acetylneuraminate_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_cef, infection=='630')[,2], subset(acetylneuraminate_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_clinda, infection=='630')[,2], subset(acetylneuraminate_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_gf, infection=='630')[,2], subset(acetylneuraminate_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
p.adjust(c(wilcox.test(subset(succinate_strep, infection=='630')[,2], subset(succinate_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(succinate_cef, infection=='630')[,2], subset(succinate_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(succinate_clinda, infection=='630')[,2], subset(succinate_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(succinate_gf, infection=='630')[,2], subset(succinate_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_5.pdf'
pdf(file=plot_file, width=6, height=9)
layout(matrix(c(1,
                2,
                3,
                4,
                5), nrow=5, ncol=1, byrow=TRUE))

par(mar=c(3,5,1.5,1), xpd=FALSE, las=1, mgp=c(3,0.7,0))

#------------------#

# N-acetylglucosamine
stripchart(substrate~infection, data=acetylglucosamine_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,14), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(substrate~infection, data=acetylglucosamine_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,14), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylglucosamine_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,14), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylglucosamine_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,14), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylglucosamine_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,14), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
axis(side=2, at=seq(0,14,3.5), labels=c('0.0','3.5','7.0','10.5','14.0'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.5, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-Germfree'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.8)
mtext('A', side=2, line=2, las=2, adj=2, padj=-3.5, cex=1.3)
legend('topright', legend='GlcNAc/GalNAc', pt.cex=0, cex=1.1, bty='n')
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(acetylglucosamine_untreated[,2]),
              median(subset(acetylglucosamine_strep, infection=='mock')[,2]), median(subset(acetylglucosamine_strep, infection=='630')[,2]),
              median(subset(acetylglucosamine_cef, infection=='mock')[,2]), median(subset(acetylglucosamine_cef, infection=='630')[,2]),
              median(subset(acetylglucosamine_clinda, infection=='mock')[,2]), median(subset(acetylglucosamine_clinda, infection=='630')[,2]),
              median(subset(acetylglucosamine_gf, infection=='mock')[,2]), median(subset(acetylglucosamine_gf, infection=='630')[,2])), 
         y1=c(median(acetylglucosamine_untreated[,2]),
              median(subset(acetylglucosamine_strep, infection=='mock')[,2]), median(subset(acetylglucosamine_strep, infection=='630')[,2]),
              median(subset(acetylglucosamine_cef, infection=='mock')[,2]), median(subset(acetylglucosamine_cef, infection=='630')[,2]),
              median(subset(acetylglucosamine_clinda, infection=='mock')[,2]), median(subset(acetylglucosamine_clinda, infection=='630')[,2]),
              median(subset(acetylglucosamine_gf, infection=='mock')[,2]), median(subset(acetylglucosamine_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=12, y0=2.5, x1=13, y1=2.5, lwd=2)
text(x=12.5, y=3.5, '*', font=2, cex=2.5)
mtext(rep('*',8), side=3, adj=c(0.21,0.28,
                                0.43,0.5,
                                0.645,0.715,
                                0.863,0.933), padj=0.4, font=2, cex=1.3, col='gray40') # Untreated vs Mock significance

#------------------#

# Proline
stripchart(substrate~infection, data=proline_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(substrate~infection, data=proline_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=proline_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=proline_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=proline_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
axis(side=2, at=c(0:4), labels=c('0.0','1.0','2.0','3.0', '4.0'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.5, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-Germfree'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.8)
mtext('B', side=2, line=2, las=2, adj=2, padj=-3.5, cex=1.3)
legend('topright', legend='Proline', pt.cex=0, bty='n', cex=1.1)
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(proline_untreated[,2]),
              median(subset(proline_strep, infection=='mock')[,2]), median(subset(proline_strep, infection=='630')[,2]),
              median(subset(proline_cef, infection=='mock')[,2]), median(subset(proline_cef, infection=='630')[,2]),
              median(subset(proline_clinda, infection=='mock')[,2]), median(subset(proline_clinda, infection=='630')[,2]),
              median(subset(proline_gf, infection=='mock')[,2]), median(subset(proline_gf, infection=='630')[,2])), 
         y1=c(median(proline_untreated[,2]),
              median(subset(proline_strep, infection=='mock')[,2]), median(subset(proline_strep, infection=='630')[,2]),
              median(subset(proline_cef, infection=='mock')[,2]), median(subset(proline_cef, infection=='630')[,2]),
              median(subset(proline_clinda, infection=='mock')[,2]), median(subset(proline_clinda, infection=='630')[,2]),
              median(subset(proline_gf, infection=='mock')[,2]), median(subset(proline_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=c(3,6,9,12), y0=c(2.9,2.8,2.5,2.9), x1=c(4,7,10,13), y1=c(2.9,2.8,2.5,2.9), lwd=2)
text(x=c(3.5,6.5,9.5,12.5), y=c(3.2,3.1,2.8,3.2), '*', font=2, cex=2.5)
mtext(rep('*',7), side=3, adj=c(0.21,0.28,
                                0.43,0.5,
                                0.645,
                                0.863,0.933), padj=0.4, font=2, cex=1.3, col='gray40') # Untreated vs Mock significance

#------------------#

# Mannitol / Sorbitol
stripchart(substrate~infection, data=mannitolsorbitol_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,70), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(substrate~infection, data=mannitolsorbitol_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,70), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=mannitolsorbitol_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,70), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=mannitolsorbitol_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,70), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=mannitolsorbitol_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,70), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
axis(side=2, at=seq(0,70,17.5), labels=c('0.0','17.5','35','52.5','70.0'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.5, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-Germfree'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.8)
mtext('C', side=2, line=2, las=2, adj=2, padj=-3.5, cex=1.3)
legend('topright', legend='Mannitol/Sorbitol', pt.cex=0, cex=1.1, bty='n')
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(mannitolsorbitol_untreated[,2]),
              median(subset(mannitolsorbitol_strep, infection=='mock')[,2]), median(subset(mannitolsorbitol_strep, infection=='630')[,2]),
              median(subset(mannitolsorbitol_cef, infection=='mock')[,2]), median(subset(mannitolsorbitol_cef, infection=='630')[,2]),
              median(subset(mannitolsorbitol_clinda, infection=='mock')[,2]), median(subset(mannitolsorbitol_clinda, infection=='630')[,2]),
              median(subset(mannitolsorbitol_gf, infection=='mock')[,2]), median(subset(mannitolsorbitol_gf, infection=='630')[,2])), 
         y1=c(median(mannitolsorbitol_untreated[,2]),
              median(subset(mannitolsorbitol_strep, infection=='mock')[,2]), median(subset(mannitolsorbitol_strep, infection=='630')[,2]),
              median(subset(mannitolsorbitol_cef, infection=='mock')[,2]), median(subset(mannitolsorbitol_cef, infection=='630')[,2]),
              median(subset(mannitolsorbitol_clinda, infection=='mock')[,2]), median(subset(mannitolsorbitol_clinda, infection=='630')[,2]),
              median(subset(mannitolsorbitol_gf, infection=='mock')[,2]), median(subset(mannitolsorbitol_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=12, y0=45, x1=13, y1=45, lwd=2)
text(x=12.5, y=50, '*', font=2, cex=2.5)
mtext(rep('*',4), side=3, adj=c(0.43,0.5,
                                0.863,0.933), padj=0.4, font=2, cex=1.3, col='gray40') # Untreated vs Mock significance

#------------------#

# Succinate
stripchart(substrate~infection, data=succinate_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,40), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(substrate~infection, data=succinate_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,40), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=succinate_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,40), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=succinate_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,40), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=succinate_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,40), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
axis(side=2, at=seq(0,40,10), labels=c('0','10','20','30','40'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.5, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-Germfree'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.8)
mtext('D', side=2, line=2, las=2, adj=2, padj=-3.5, cex=1.3)
legend('topright', legend='Succinate', pt.cex=0, cex=1.1, bty='n')
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(succinate_untreated[,2]),
              median(subset(succinate_strep, infection=='mock')[,2]), median(subset(succinate_strep, infection=='630')[,2]),
              median(subset(succinate_cef, infection=='mock')[,2]), median(subset(succinate_cef, infection=='630')[,2]),
              median(subset(succinate_clinda, infection=='mock')[,2]), median(subset(succinate_clinda, infection=='630')[,2]),
              median(subset(succinate_gf, infection=='mock')[,2]), median(subset(succinate_gf, infection=='630')[,2])), 
         y1=c(median(succinate_untreated[,2]),
              median(subset(succinate_strep, infection=='mock')[,2]), median(subset(succinate_strep, infection=='630')[,2]),
              median(subset(succinate_cef, infection=='mock')[,2]), median(subset(succinate_cef, infection=='630')[,2]),
              median(subset(succinate_clinda, infection=='mock')[,2]), median(subset(succinate_clinda, infection=='630')[,2]),
              median(subset(succinate_gf, infection=='mock')[,2]), median(subset(succinate_gf, infection=='630')[,2])),
         lwd=3)
mtext(rep('*',7), side=3, adj=c(0.21,
                                0.43,0.5,
                                0.645,0.715,
                                0.863,0.933), padj=0.4, font=2, cex=1.3, col='gray40') # Untreated vs Mock significance

#------------------#

# N-acetylneuraminate
stripchart(substrate~infection, data=acetylneuraminate_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(substrate~infection, data=acetylneuraminate_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylneuraminate_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylneuraminate_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylneuraminate_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
axis(side=2, at=c(0:4), labels=c('0.0','1.0','2.0','3.0', '4.0'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.5, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-Germfree'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.8)
mtext('E', side=2, line=2, las=2, adj=2, padj=-3.5, cex=1.3)
legend('topright', legend='Neu5Ac', pt.cex=0, bty='n', cex=1.1)
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(acetylneuraminate_untreated[,2]),
              median(subset(acetylneuraminate_strep, infection=='mock')[,2]), median(subset(acetylneuraminate_strep, infection=='630')[,2]),
              median(subset(acetylneuraminate_cef, infection=='mock')[,2]), median(subset(acetylneuraminate_cef, infection=='630')[,2]),
              median(subset(acetylneuraminate_clinda, infection=='mock')[,2]), median(subset(acetylneuraminate_clinda, infection=='630')[,2]),
              median(subset(acetylneuraminate_gf, infection=='mock')[,2]), median(subset(acetylneuraminate_gf, infection=='630')[,2])), 
         y1=c(median(acetylneuraminate_untreated[,2]),
              median(subset(acetylneuraminate_strep, infection=='mock')[,2]), median(subset(acetylneuraminate_strep, infection=='630')[,2]),
              median(subset(acetylneuraminate_cef, infection=='mock')[,2]), median(subset(acetylneuraminate_cef, infection=='630')[,2]),
              median(subset(acetylneuraminate_clinda, infection=='mock')[,2]), median(subset(acetylneuraminate_clinda, infection=='630')[,2]),
              median(subset(acetylneuraminate_gf, infection=='mock')[,2]), median(subset(acetylneuraminate_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=c(6,12), y0=c(2.3,2.3), x1=c(7,13), y1=c(2.3,2.3), lwd=2)
text(x=c(6.5,12.5), y=c(2.6,2.6), '*', font=2, cex=2.5)
mtext(rep('*',4), side=3, adj=c(0.43,
                                0.715,
                                0.863,0.933), padj=0.4, font=2, cex=1.3, col='gray40') # Untreated vs Mock significance

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
detach('package:wesanderson', character.only = TRUE)
rm(list=ls())
gc()

