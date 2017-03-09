
# Start with blank slate
rm(list=ls())
gc()

# Load dependency
if ('wesanderson' %in% installed.packages()[,"Package"] == FALSE){
  install.packages(as.character('wesanderson'), quiet=TRUE);
}
library('wesanderson', verbose=FALSE, character.only=TRUE)


# Select files
scfa <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/cef_acetate_630.txt'
metabolome <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/metabolomics.tsv'
metadata <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'

# Read in data
scfa <- read.delim(scfa, sep='\t', header=T)
metabolome <- read.delim(metabolome, sep='\t', header=T, row.names=1)
metabolome <- metabolome[, !colnames(metabolome) %in% c('CefC5M2','StrepC4M1')] # Remove possible contamination
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination

# Format the data
scfa$group <- factor(scfa$group, levels=c('infected','mock'))
scfa$acetate <- as.numeric(as.character(scfa$acetate))

# Subset for stats
mock <- as.numeric(scfa[scfa$group == 'mock', 2])
infected <- as.numeric(scfa[scfa$group == 'infected', 2])

# Subset untargeted metabolomics
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
acetylglucosamine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-acetylglucosamine/N-acetylgalactosamine')))]
acetylglucosamine_cef <- subset(acetylglucosamine, abx == 'cefoperazone')
acetylglucosamine_cef$abx <- NULL
colnames(acetylglucosamine_cef) <- c('infection', 'substrate')
acetylglucosamine_strep <- subset(acetylglucosamine, abx == 'streptomycin')
acetylglucosamine_strep$abx <- NULL
colnames(acetylglucosamine_strep) <- c('infection', 'substrate')
acetylglucosamine_clinda <- subset(acetylglucosamine, abx == 'clindamycin')
acetylglucosamine_clinda$abx <- NULL
colnames(acetylglucosamine_clinda) <- c('infection', 'substrate')
acetylglucosamine_gf <- subset(acetylglucosamine, abx == 'germfree')
acetylglucosamine_gf$abx <- NULL
colnames(acetylglucosamine_gf) <- c('infection', 'substrate')
rm(acetylglucosamine)
glycine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('glycine')))]
glycine_cef <- subset(glycine, abx == 'cefoperazone')
glycine_cef$abx <- NULL
colnames(glycine_cef) <- c('infection', 'substrate')
glycine_strep <- subset(glycine, abx == 'streptomycin')
glycine_strep$abx <- NULL
colnames(glycine_strep) <- c('infection', 'substrate')
glycine_clinda <- subset(glycine, abx == 'clindamycin')
glycine_clinda$abx <- NULL
colnames(glycine_clinda) <- c('infection', 'substrate')
glycine_gf <- subset(glycine, abx == 'germfree')
glycine_gf$abx <- NULL
colnames(glycine_gf) <- c('infection', 'substrate')
rm(glycine)
proline <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('proline')))]
proline_cef <- subset(proline, abx == 'cefoperazone')
proline_cef$abx <- NULL
colnames(proline_cef) <- c('infection', 'substrate')
proline_strep <- subset(proline, abx == 'streptomycin')
proline_strep$abx <- NULL
colnames(proline_strep) <- c('infection', 'substrate')
proline_clinda <- subset(proline, abx == 'clindamycin')
proline_clinda$abx <- NULL
colnames(proline_clinda) <- c('infection', 'substrate')
proline_gf <- subset(proline, abx == 'germfree')
proline_gf$abx <- NULL
colnames(proline_gf) <- c('infection', 'substrate')
rm(proline)
mannitolsorbitol <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('mannitol/sorbitol')))]
mannitolsorbitol_cef <- subset(mannitolsorbitol, abx == 'cefoperazone')
mannitolsorbitol_cef$abx <- NULL
colnames(mannitolsorbitol_cef) <- c('infection', 'substrate')
mannitolsorbitol_strep <- subset(mannitolsorbitol, abx == 'streptomycin')
mannitolsorbitol_strep$abx <- NULL
colnames(mannitolsorbitol_strep) <- c('infection', 'substrate')
rm(mannitolsorbitol)
salicin <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('salicylate','glucose')))]
salicin_clinda <- subset(salicin, abx == 'clindamycin')
salicin_clinda$abx <- NULL
colnames(salicin_clinda) <- c('infection', 'salicylate','glucose')
rm(salicin)
acetylneuraminate <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-acetylneuraminate')))]
acetylneuraminate_gf <- subset(acetylneuraminate, abx == 'germfree')
acetylneuraminate_gf$abx <- NULL
colnames(acetylneuraminate_gf) <- c('infection', 'substrate')
rm(acetylneuraminate)
rm(metabolome)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_6.pdf'
pdf(file=plot_file, width=7, height=9)
layout(matrix(c(1,1,
                2,2,
                3,4,
                5,6), nrow=4, ncol=2, byrow=TRUE))
par(las=1, mgp=c(2.3,0.7,0))

#-------------------------------------------------------------------------------------------------------------------------------------#

# N-acetylglucosamine
par(mar=c(3,5,1,1))
stripchart(substrate~infection, data=acetylglucosamine_strep, vertical=T, pch=19, at=c(1,2),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,6), xlim=c(0,12),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
stripchart(substrate~infection, data=acetylglucosamine_cef, vertical=T, pch=19, at=c(4,5),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,6), xlim=c(0,12),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylglucosamine_clinda, vertical=T, pch=19, at=c(7,8),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,6), xlim=c(0,12),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylglucosamine_gf, vertical=T, pch=19, at=c(10,11),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,6), xlim=c(0,12),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
axis(side=2, at=c(0:6), labels=c('0.0','1.0','2.0','3.0','4.0','5.0','6.0'), cex.axis=1.2)
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.9)
mtext(c('+','-','+','-','+','-','+','-'), side=1, 
      at=c(1,2,4,5,7,8,10,11), padj=0.5, cex=1.2)
mtext(c('Streptomycin','Cefoperazone','Clindamycin','Gnotobiotic'), side=1, 
      at=c(1.5,4.5,7.5,10.6), padj=2)
mtext('A', side=2, line=2, las=2, adj=1.7, padj=-5, cex=1.3)
legend('topright', legend='N-Acetylglucosamine / N-Acetylgalactosamine', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,1.6,3.6,4.6,6.6,7.6,9.6,10.6), x1=c(1.4,2.4,4.4,5.4,7.4,8.4,10.4,11.4),
         y0=c(median(subset(acetylglucosamine_strep, infection=='630')[,2]), median(subset(acetylglucosamine_strep, infection=='mock')[,2]),
              median(subset(acetylglucosamine_cef, infection=='630')[,2]), median(subset(acetylglucosamine_cef, infection=='mock')[,2]),
              median(subset(acetylglucosamine_clinda, infection=='630')[,2]), median(subset(acetylglucosamine_clinda, infection=='mock')[,2]),
              median(subset(acetylglucosamine_gf, infection=='630')[,2]), median(subset(acetylglucosamine_gf, infection=='mock')[,2])), 
         y1=c(median(subset(acetylglucosamine_strep, infection=='630')[,2]), median(subset(acetylglucosamine_strep, infection=='mock')[,2]),
              median(subset(acetylglucosamine_cef, infection=='630')[,2]), median(subset(acetylglucosamine_cef, infection=='mock')[,2]),
              median(subset(acetylglucosamine_clinda, infection=='630')[,2]), median(subset(acetylglucosamine_clinda, infection=='mock')[,2]),
              median(subset(acetylglucosamine_gf, infection=='630')[,2]), median(subset(acetylglucosamine_gf, infection=='mock')[,2])),
         lwd=3)
p.adjust(c(wilcox.test(subset(acetylglucosamine_strep, infection=='630')[,2], subset(acetylglucosamine_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_cef, infection=='630')[,2], subset(acetylglucosamine_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_clinda, infection=='630')[,2], subset(acetylglucosamine_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_gf, infection=='630')[,2], subset(acetylglucosamine_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
text(x=10, y=2, '*', font=2, cex=2.5)

#------------------#

# Proline
par(mar=c(3,5,1,1))
stripchart(substrate~infection, data=proline_strep, vertical=T, pch=19, at=c(1,2),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,3), xlim=c(0,12),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
stripchart(substrate~infection, data=proline_cef, vertical=T, pch=19, at=c(4,5),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,3), xlim=c(0,12),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=proline_clinda, vertical=T, pch=19, at=c(7,8),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,3), xlim=c(0,12),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=proline_gf, vertical=T, pch=19, at=c(10,11),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,3), xlim=c(0,12),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
axis(side=2, at=c(0:3), labels=c('0.0','1.0','2.0','3.0'), cex.axis=1.2)
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.9)
mtext(c('+','-','+','-','+','-','+','-'), side=1, 
      at=c(1,2,4,5,7,8,10,11), padj=0.5, cex=1.2)
mtext(c('Streptomycin','Cefoperazone','Clindamycin','Gnotobiotic'), side=1, 
      at=c(1.5,4.5,7.5,10.6), padj=2)
mtext('B', side=2, line=2, las=2, adj=1.7, padj=-5, cex=1.3)
legend('topright', legend='Proline', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,1.6,3.6,4.6,6.6,7.6,9.6,10.6), x1=c(1.4,2.4,4.4,5.4,7.4,8.4,10.4,11.4),
         y0=c(median(subset(proline_strep, infection=='630')[,2]), median(subset(proline_strep, infection=='mock')[,2]),
              median(subset(proline_cef, infection=='630')[,2]), median(subset(proline_cef, infection=='mock')[,2]),
              median(subset(proline_clinda, infection=='630')[,2]), median(subset(proline_clinda, infection=='mock')[,2]),
              median(subset(proline_gf, infection=='630')[,2]), median(subset(proline_gf, infection=='mock')[,2])), 
         y1=c(median(subset(proline_strep, infection=='630')[,2]), median(subset(proline_strep, infection=='mock')[,2]),
              median(subset(proline_cef, infection=='630')[,2]), median(subset(proline_cef, infection=='mock')[,2]),
              median(subset(proline_clinda, infection=='630')[,2]), median(subset(proline_clinda, infection=='mock')[,2]),
              median(subset(proline_gf, infection=='630')[,2]), median(subset(proline_gf, infection=='mock')[,2])),
         lwd=3)
p.adjust(c(wilcox.test(subset(proline_strep, infection=='630')[,2], subset(proline_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_cef, infection=='630')[,2], subset(proline_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_clinda, infection=='630')[,2], subset(proline_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_gf, infection=='630')[,2], subset(proline_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
text(x=c(1,4,7,10), y=c(1,1.5,2,1), '*', font=2, cex=2.5)

#------------------#

# Salicin
par(mar=c(3,5,1,1))
stripchart(salicylate~infection, data=salicin_clinda, vertical=T, pch=19, at=c(1,2),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,7), xlim=c(0.5,2.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
axis(side=2, at=c(0:7), labels=c('0','1.0','2.0','3.0','4.0','5.0','6.0','7.0'), cex.axis=1.2)
mtext('CDI:', side=1, at=0.5, padj=0.5, cex=0.9)
mtext(c('+','-'), side=1, 
      at=c(1,2), padj=0.5, cex=1.2)
mtext(c('Clindamycin'), side=1, 
      at=c(1.5), padj=2)
mtext('C', side=2, line=2, las=2, adj=1.7, padj=-5, cex=1.3)
legend('topright', legend='Salicylate', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,1.6), x1=c(1.4,2.4),
         y0=c(median(subset(salicin_clinda, infection=='630')[,2]), median(subset(salicin_clinda, infection=='mock')[,2])), 
         y1=c(median(subset(salicin_clinda, infection=='630')[,2]), median(subset(salicin_clinda, infection=='mock')[,2])),
         lwd=3)

par(mar=c(3,5,1,1))
stripchart(glucose~infection, data=salicin_clinda, vertical=T, pch=19, at=c(1,2),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,6), xlim=c(0.5,2.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
axis(side=2, at=c(0:7), labels=c('0','1.0','2.0','3.0','4.0','5.0','6.0','7.0'), cex.axis=1.2)
mtext('CDI:', side=1, at=0.5, padj=0.5, cex=0.9)
mtext(c('+','-'), side=1, 
      at=c(1,2), padj=0.5, cex=1.2)
mtext(c('Clindamycin'), side=1, 
      at=c(1.5), padj=2)
mtext('D', side=2, line=2, las=2, adj=1.7, padj=-5, cex=1.3)
legend('topright', legend='Glucose', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,1.6), x1=c(1.4,2.4),
         y0=c(median(subset(salicin_clinda, infection=='630')[,3]), median(subset(salicin_clinda, infection=='mock')[,3])), 
         y1=c(median(subset(salicin_clinda, infection=='630')[,3]), median(subset(salicin_clinda, infection=='mock')[,3])),
         lwd=3)
p.adjust(c(wilcox.test(subset(salicin_clinda, infection=='630')[,2], subset(salicin_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(salicin_clinda, infection=='630')[,3], subset(salicin_clinda, infection=='mock')[,3], exact=F)$p.value), method='BH')

#------------------#

# Mannitol / Sorbitol
par(mar=c(3,5,1,1))
stripchart(substrate~infection, data=mannitolsorbitol_strep, vertical=T, pch=19, at=c(1,2),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,65), xlim=c(0,6),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
stripchart(substrate~infection, data=mannitolsorbitol_cef, vertical=T, pch=19, at=c(4,5),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,65), xlim=c(0,6),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
axis(side=2, at=c(0,10,20,30,40,50,60), labels=c('0','10','20','30','40','50','60'), cex.axis=1.2)
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.9)
mtext(c('+','-','+','-'), side=1, 
      at=c(1,2,4,5), padj=0.5, cex=1.2)
mtext(c('Streptomycin','Cefoperazone'), side=1, 
      at=c(1.5,4.5), padj=2)
mtext('E', side=2, line=2, las=2, adj=1.7, padj=-5, cex=1.3)
legend('topright', legend='Mannitol / Sorbitol', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,1.6,3.6,4.6), x1=c(1.4,2.4,4.4,5.4),
         y0=c(median(subset(mannitolsorbitol_strep, infection=='630')[,2]), median(subset(mannitolsorbitol_strep, infection=='mock')[,2]),
              median(subset(mannitolsorbitol_cef, infection=='630')[,2]), median(subset(mannitolsorbitol_cef, infection=='mock')[,2])), 
         y1=c(median(subset(mannitolsorbitol_strep, infection=='630')[,2]), median(subset(mannitolsorbitol_strep, infection=='mock')[,2]),
              median(subset(mannitolsorbitol_cef, infection=='630')[,2]), median(subset(mannitolsorbitol_cef, infection=='mock')[,2])),
         lwd=3)
p.adjust(c(wilcox.test(subset(mannitolsorbitol_strep, infection=='630')[,2], subset(mannitolsorbitol_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_cef, infection=='630')[,2], subset(mannitolsorbitol_cef, infection=='mock')[,2], exact=F)$p.value), method='BH')

#------------------#

# N-acetylneuraminate
par(mar=c(3,5,1,1))
stripchart(substrate~infection, data=acetylneuraminate_gf, vertical=T, pch=19, at=c(1,2),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,3), xlim=c(0.5,2.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
axis(side=2, at=c(0:3), labels=c('0','1.0','2.0','3.0'), cex.axis=1.2)
mtext('CDI:', side=1, at=0.5, padj=0.5, cex=0.9)
mtext(c('+','-'), side=1, 
      at=c(1,2), padj=0.5, cex=1.2)
mtext(c('Gnotobiotic'), side=1, 
      at=c(1.5), padj=2)
mtext('F', side=2, line=2, las=2, adj=1.7, padj=-5, cex=1.3)
legend('topright', legend='N-Acetylneuraminate', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.7,1.7), x1=c(1.3,2.3),
         y0=c(median(subset(acetylneuraminate_gf, infection=='630')[,2]), median(subset(acetylneuraminate_gf, infection=='mock')[,2])), 
         y1=c(median(subset(acetylneuraminate_gf, infection=='630')[,2]), median(subset(acetylneuraminate_gf, infection=='mock')[,2])),
         lwd=3)
wilcox.test(subset(acetylneuraminate_gf, infection=='630')[,2], subset(acetylneuraminate_gf, infection=='mock')[,2], exact=F)$p.value
text(x=1, y=1.2, '*', font=2, cex=2.5)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
detach('package:wesanderson', character.only = TRUE)
rm(list=ls())
gc()
