
# Start with clear environment
rm(list=ls())
gc()

# Load dependencies
deps <- c('wesanderson', 'plyr');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Select files
metabolome <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/metabolomics.tsv'
metadata <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'

# Read in data
metabolome <- read.delim(metabolome, sep='\t', header=T, row.names=1)
metabolome <- metabolome[, !colnames(metabolome) %in% c('CefC5M2','StrepC4M1')] # Remove possible contamination
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination

# Merge
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


# Subset metabolomics
metabolome <- subset(metabolome, infection == 'mock')
metabolome$infection <- NULL
acetylglucosamine <- metabolome[, c(1,which(colnames(metabolome) %in% c('N-acetylglucosamine/N-acetylgalactosamine')))]
colnames(acetylglucosamine) <- c('abx', 'substrate')
acetylglucosamine$abx <- factor(acetylglucosamine$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
glycine <- metabolome[, c(1,which(colnames(metabolome) %in% c('glycine')))]
colnames(glycine) <- c('abx', 'substrate')
glycine$abx <- factor(glycine$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
proline <- metabolome[, c(1,which(colnames(metabolome) %in% c('proline')))]
colnames(proline) <- c('abx', 'substrate')
proline$abx <- factor(proline$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
mannitolsorbitol <- metabolome[, c(1,which(colnames(metabolome) %in% c('mannitol/sorbitol')))]
colnames(mannitolsorbitol) <- c('abx', 'substrate')
mannitolsorbitol$abx <- factor(mannitolsorbitol$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
salicin <- metabolome[, c(1,which(colnames(metabolome) %in% c('salicylate','glucose')))]
salicin$abx <- factor(salicin$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
acetylneuraminate <- metabolome[, c(1,which(colnames(metabolome) %in% c('N-acetylneuraminate')))]
colnames(acetylneuraminate) <- c('abx', 'substrate')
acetylneuraminate$abx <- factor(acetylneuraminate$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
rm(metabolome)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S6.pdf'
select_palette <- c('gray40', wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen')
pdf(file=plot_file, width=14, height=10)
layout(matrix(c(1,2,
                3,4,
                5,6,
                7,8),
              nrow=4, ncol=2, byrow = TRUE))
par(las=1, mgp=c(2.5,0.7,0))

#--------------------------------#

# N-acetylglucosamine
par(mar=c(3,5,1,1))
stripchart(substrate~abx, data=acetylglucosamine, vertical=T, pch=19, lwd=3, 
           xaxt='n', yaxt='n', col=select_palette, ylim=c(0,14),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25)
mtext(c('No Antibiotics\nSPF','Streptomycin\nSPF','Cefoperazone\nSPF','Clindamycin\nSPF','Gnotobiotic\nGF'), side=1, 
      at=c(1:5), padj=1)
axis(side=2, at=seq(0,14,2), labels=c('0.0','2.0','4.0','6.0','8.0','10.0','12.0','14.0'))
mtext('A', side=2, line=2, las=2, adj=2, padj=-4, cex=1.7)
legend('topright', legend='N-acetylglucosamine / N-acetylgalactosamine', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,1.6,2.6,3.6,4.6), x1=c(1.4,2.4,3.4,4.4,5.4),
         y0=c(median(subset(acetylglucosamine, abx=='none')[,2]), 
              median(subset(acetylglucosamine, abx=='streptomycin')[,2]), 
              median(subset(acetylglucosamine, abx=='cefoperazone')[,2]), 
              median(subset(acetylglucosamine, abx=='clindamycin')[,2]), 
              median(subset(acetylglucosamine, abx=='germfree')[,2])), 
         y1=c(median(subset(acetylglucosamine, abx=='none')[,2]), 
                median(subset(acetylglucosamine, abx=='streptomycin')[,2]), 
                median(subset(acetylglucosamine, abx=='cefoperazone')[,2]), 
                median(subset(acetylglucosamine, abx=='clindamycin')[,2]), 
                median(subset(acetylglucosamine, abx=='germfree')[,2])),
         lwd=3)
p.adjust(c(wilcox.test(subset(acetylglucosamine, abx=='none')[,2], subset(acetylglucosamine, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine, abx=='none')[,2], subset(acetylglucosamine, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine, abx=='none')[,2], subset(acetylglucosamine, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine, abx=='none')[,2], subset(acetylglucosamine, abx=='germfree')[,2], exact=F)$p.value), method='BH')
text(x=c(2:5), y=c(7,4,4,4), '*', font=2, cex=3)

#--------------------------------#

# N-acetylneuraminate
par(mar=c(3,5,1,1))
stripchart(substrate~abx, data=acetylneuraminate, vertical=T, pch=19, lwd=3, 
           xaxt='n', yaxt='n', col=select_palette, ylim=c(0,3),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25)
mtext(c('No Antibiotics\nSPF','Streptomycin\nSPF','Cefoperazone\nSPF','Clindamycin\nSPF','Gnotobiotic\nGF'), side=1, 
      at=c(1:5), padj=1)
axis(side=2, at=c(0:3), labels=c('0.0','1.0','2.0','3.0'))
mtext('B', side=2, line=2, las=2, adj=2, padj=-4, cex=1.7)
legend('topright', legend='N-acetylneuraminate', pt.cex=0, bty='n', cex=1.4)
segments(x0=c(0.6,1.6,2.6,3.6,4.6), x1=c(1.4,2.4,3.4,4.4,5.4),
         y0=c(median(subset(acetylneuraminate, abx=='none')[,2]), 
              median(subset(acetylneuraminate, abx=='streptomycin')[,2]), 
              median(subset(acetylneuraminate, abx=='cefoperazone')[,2]), 
              median(subset(acetylneuraminate, abx=='clindamycin')[,2]), 
              median(subset(acetylneuraminate, abx=='germfree')[,2])), 
         y1=c(median(subset(acetylneuraminate, abx=='none')[,2]), 
              median(subset(acetylneuraminate, abx=='streptomycin')[,2]), 
              median(subset(acetylneuraminate, abx=='cefoperazone')[,2]), 
              median(subset(acetylneuraminate, abx=='clindamycin')[,2]), 
              median(subset(acetylneuraminate, abx=='germfree')[,2])),
         lwd=3)
p.adjust(c(wilcox.test(subset(acetylneuraminate, abx=='none')[,2], subset(acetylneuraminate, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate, abx=='none')[,2], subset(acetylneuraminate, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate, abx=='none')[,2], subset(acetylneuraminate, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate, abx=='none')[,2], subset(acetylneuraminate, abx=='germfree')[,2], exact=F)$p.value), method='BH')
text(x=c(3,5), y=2.3, '*', font=2, cex=3)

#--------------------------------#

# Proline
par(mar=c(3,5,1,1))
stripchart(substrate~abx, data=proline, vertical=T, pch=19, lwd=3, 
           xaxt='n', yaxt='n', col=select_palette, ylim=c(0,3), 
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25)
mtext(c('No Antibiotics\nSPF','Streptomycin\nSPF','Cefoperazone\nSPF','Clindamycin\nSPF','Gnotobiotic\nGF'), side=1, 
      at=c(1:5), padj=1)
axis(side=2, at=c(0:3), labels=c('0.0','1.0','2.0','3.0'))
mtext('C', side=2, line=2, las=2, adj=2, padj=-4, cex=1.7)
legend('topright', legend='Proline', pt.cex=0, bty='n', cex=1.4)
segments(x0=c(0.6,1.6,2.6,3.6,4.6), x1=c(1.4,2.4,3.4,4.4,5.4),
         y0=c(median(subset(proline, abx=='none')[,2]), 
              median(subset(proline, abx=='streptomycin')[,2]), 
              median(subset(proline, abx=='cefoperazone')[,2]), 
              median(subset(proline, abx=='clindamycin')[,2]), 
              median(subset(proline, abx=='germfree')[,2])), 
         y1=c(median(subset(proline, abx=='none')[,2]), 
              median(subset(proline, abx=='streptomycin')[,2]), 
              median(subset(proline, abx=='cefoperazone')[,2]), 
              median(subset(proline, abx=='clindamycin')[,2]), 
              median(subset(proline, abx=='germfree')[,2])),
         lwd=3)
p.adjust(c(wilcox.test(subset(proline, abx=='none')[,2], subset(proline, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(proline, abx=='none')[,2], subset(proline, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(proline, abx=='none')[,2], subset(proline, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(proline, abx=='none')[,2], subset(proline, abx=='germfree')[,2], exact=F)$p.value), method='BH')
text(x=c(2:5), y=2.9, '*', font=2, cex=3)

#--------------------------------#

# Glycine
par(mar=c(3,5,1,1))
stripchart(substrate~abx, data=glycine, vertical=T, pch=19, lwd=3, 
           xaxt='n', yaxt='n', col=select_palette, ylim=c(0,3),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25)
mtext(c('No Antibiotics\nSPF','Streptomycin\nSPF','Cefoperazone\nSPF','Clindamycin\nSPF','Gnotobiotic\nGF'), side=1, 
      at=c(1:5), padj=1)
axis(side=2, at=c(0:3), labels=c('0.0','1.0','2.0','3.0'))
mtext('D', side=2, line=2, las=2, adj=2, padj=-4, cex=1.7)
legend('topright', legend='Glycine', pt.cex=0, bty='n', cex=1.4)
segments(x0=c(0.6,1.6,2.6,3.6,4.6), x1=c(1.4,2.4,3.4,4.4,5.4),
         y0=c(median(subset(glycine, abx=='none')[,2]), 
              median(subset(glycine, abx=='streptomycin')[,2]), 
              median(subset(glycine, abx=='cefoperazone')[,2]), 
              median(subset(glycine, abx=='clindamycin')[,2]), 
              median(subset(glycine, abx=='germfree')[,2])), 
         y1=c(median(subset(glycine, abx=='none')[,2]), 
              median(subset(glycine, abx=='streptomycin')[,2]), 
              median(subset(glycine, abx=='cefoperazone')[,2]), 
              median(subset(glycine, abx=='clindamycin')[,2]), 
              median(subset(glycine, abx=='germfree')[,2])),
         lwd=3)
p.adjust(c(wilcox.test(subset(glycine, abx=='none')[,2], subset(glycine, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine, abx=='none')[,2], subset(glycine, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine, abx=='none')[,2], subset(glycine, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(glycine, abx=='none')[,2], subset(glycine, abx=='germfree')[,2], exact=F)$p.value), method='BH')
text(x=c(3,5), y=2.3, '*', font=2, cex=3)

#--------------------------------#

# Salicin 1
par(mar=c(3,5,1,1))
stripchart(salicylate~abx, data=salicin, vertical=T, pch=19, lwd=3, 
           xaxt='n', yaxt='n', col=select_palette, ylim=c(0,5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25)
mtext(c('No Antibiotics\nSPF','Streptomycin\nSPF','Cefoperazone\nSPF','Clindamycin\nSPF','Gnotobiotic\nGF'), side=1, 
      at=c(1:5), padj=1)
axis(side=2, at=c(0:5), labels=c('0.0','1.0','2.0','3.0','4.0','5.0'))
mtext('E', side=2, line=2, las=2, adj=2, padj=-4, cex=1.7)
legend('topright', legend='Salicylate', pt.cex=0, bty='n', cex=1.4)
segments(x0=c(0.6,1.6,2.6,3.6,4.6), x1=c(1.4,2.4,3.4,4.4,5.4),
         y0=c(median(subset(salicin, abx=='none')[,3]), 
              median(subset(salicin, abx=='streptomycin')[,3]), 
              median(subset(salicin, abx=='cefoperazone')[,3]), 
              median(subset(salicin, abx=='clindamycin')[,3]), 
              median(subset(salicin, abx=='germfree')[,3])), 
         y1=c(median(subset(salicin, abx=='none')[,3]), 
              median(subset(salicin, abx=='streptomycin')[,3]), 
              median(subset(salicin, abx=='cefoperazone')[,3]), 
              median(subset(salicin, abx=='clindamycin')[,3]), 
              median(subset(salicin, abx=='germfree')[,3])),
         lwd=3)
p.adjust(c(wilcox.test(subset(salicin, abx=='none')[,3], subset(salicin, abx=='streptomycin')[,3], exact=F)$p.value,
           wilcox.test(subset(salicin, abx=='none')[,3], subset(salicin, abx=='cefoperazone')[,3], exact=F)$p.value,
           wilcox.test(subset(salicin, abx=='none')[,3], subset(salicin, abx=='clindamycin')[,3], exact=F)$p.value,
           wilcox.test(subset(salicin, abx=='none')[,3], subset(salicin, abx=='germfree')[,3], exact=F)$p.value), method='BH')
text(x=c(2:5), y=c(3,3,3,2), '*', font=2, cex=3)

#--------------------------------#

# Salicin 2
par(mar=c(3,5,1,1))
stripchart(glucose~abx, data=salicin, vertical=T, pch=19, lwd=3, 
           xaxt='n', yaxt='n', col=select_palette, ylim=c(0,25),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25)
mtext(c('No Antibiotics\nSPF','Streptomycin\nSPF','Cefoperazone\nSPF','Clindamycin\nSPF','Gnotobiotic\nGF'), side=1, 
      at=c(1:5), padj=1)
axis(side=2, at=seq(0,25,5), labels=c('0.0','5.0','10.0','15.0','20.0','25.0'))
mtext('F', side=2, line=2, las=2, adj=2, padj=-4, cex=1.7)
legend('topright', legend='Glucose', pt.cex=0, bty='n', cex=1.4)
segments(x0=c(0.6,1.6,2.6,3.6,4.6), x1=c(1.4,2.4,3.4,4.4,5.4),
         y0=c(median(subset(salicin, abx=='none')[,2]), 
              median(subset(salicin, abx=='streptomycin')[,2]), 
              median(subset(salicin, abx=='cefoperazone')[,2]), 
              median(subset(salicin, abx=='clindamycin')[,2]), 
              median(subset(salicin, abx=='germfree')[,2])), 
         y1=c(median(subset(salicin, abx=='none')[,2]), 
              median(subset(salicin, abx=='streptomycin')[,2]), 
              median(subset(salicin, abx=='cefoperazone')[,2]), 
              median(subset(salicin, abx=='clindamycin')[,2]), 
              median(subset(salicin, abx=='germfree')[,2])),
         lwd=3)
p.adjust(c(wilcox.test(subset(salicin, abx=='none')[,2], subset(salicin, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(salicin, abx=='none')[,2], subset(salicin, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(salicin, abx=='none')[,2], subset(salicin, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(salicin, abx=='none')[,2], subset(salicin, abx=='germfree')[,2], exact=F)$p.value), method='BH')
text(x=c(2:5), y=5, '*', font=2, cex=3)

#--------------------------------#

# Mannitol + Sorbitol
par(mar=c(4,5,1,1))
stripchart(substrate~abx, data=mannitolsorbitol, vertical=T, pch=19, lwd=3, 
           xaxt='n', yaxt='n', col=select_palette, ylim=c(0,50),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25)
mtext(c('No Antibiotics\nSPF','Streptomycin\nSPF','Cefoperazone\nSPF','Clindamycin\nSPF','Gnotobiotic\nGF'), side=1, 
      at=c(1:5), padj=1)
axis(side=2, at=seq(0,50,10), labels=c('0.0','10.0','20.0','30.0','40.0','50.0'))
mtext('G', side=2, line=2, las=2, adj=1.7, padj=-4, cex=1.7)
legend('topleft', legend='Mannitol / Sorbitol', pt.cex=0, bty='n', cex=1.4)
segments(x0=c(0.6,1.6,2.6,3.6,4.6), x1=c(1.4,2.4,3.4,4.4,5.4),
         y0=c(median(subset(mannitolsorbitol, abx=='none')[,2]), 
              median(subset(mannitolsorbitol, abx=='streptomycin')[,2]), 
              median(subset(mannitolsorbitol, abx=='cefoperazone')[,2]), 
              median(subset(mannitolsorbitol, abx=='clindamycin')[,2]), 
              median(subset(mannitolsorbitol, abx=='germfree')[,2])), 
         y1=c(median(subset(mannitolsorbitol, abx=='none')[,2]), 
              median(subset(mannitolsorbitol, abx=='streptomycin')[,2]), 
              median(subset(mannitolsorbitol, abx=='cefoperazone')[,2]), 
              median(subset(mannitolsorbitol, abx=='clindamycin')[,2]), 
              median(subset(mannitolsorbitol, abx=='germfree')[,2])),
         lwd=3)
p.adjust(c(wilcox.test(subset(mannitolsorbitol, abx=='none')[,2], subset(mannitolsorbitol, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol, abx=='none')[,2], subset(mannitolsorbitol, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol, abx=='none')[,2], subset(mannitolsorbitol, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol, abx=='none')[,2], subset(mannitolsorbitol, abx=='germfree')[,2], exact=F)$p.value), method='BH')
text(x=c(3,5), y=48, '*', font=2, cex=3)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

