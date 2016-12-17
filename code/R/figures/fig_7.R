
# Load dependencies
deps <- c('wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Select files
concentrations <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/ms_substrates.tsv'
metadata <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'

# Read in data
concentrations <- t(read.delim(concentrations, sep='\t', header=T, row.names=1))
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

# Format and merge tables
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$type <- NULL
concentrations <- merge(concentrations, metadata, by='row.names')
rownames(concentrations) <- concentrations$Row.names
concentrations$Row.names <- NULL
concentrations <- subset(concentrations, infection != '630')
concentrations$infection <- NULL

# Sort by treatment group
concentrations$abx <- factor(concentrations$abx, levels=c('none','streptomycin', 'cefoperazone', 'clindamycin', 'germfree'))

# Subset each metabolites - main body
acetylglucosamine <- concentrations[,c(1,14)]
mannitol_sorbitol <- concentrations[,c(3,14)]
galactitol <- concentrations[,c(4,14)]
salicylate <- concentrations[,c(5,14)]
acetylneuriminate <- concentrations[,c(6,14)]

# Subset each metabolites - supplement
proline <- concentrations[,c(2,14)]
ribose <- concentrations[,c(7,14)]
lysine <- concentrations[,c(8,14)]
isoleucine <- concentrations[,c(9,14)]
leucine <- concentrations[,c(10,14)]
methionine <- concentrations[,c(11,14)]
threonine <- concentrations[,c(12,14)]
glycine <- concentrations[,c(13,14)]

# Remove raw data
rm(concentrations, metadata)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Calculate significant differences and correct p-values
acetylglucosamine_p <- p.adjust(c(wilcox.test(subset(acetylglucosamine, acetylglucosamine$abx=='none')[,1], 
                                     subset(acetylglucosamine, acetylglucosamine$abx=='streptomycin')[,1], exact=FALSE)$p.value,
                         wilcox.test(subset(acetylglucosamine, acetylglucosamine$abx=='none')[,1], 
                                     subset(acetylglucosamine, acetylglucosamine$abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                         wilcox.test(subset(acetylglucosamine, acetylglucosamine$abx=='none')[,1], 
                                     subset(acetylglucosamine, acetylglucosamine$abx=='clindamycin')[,1], exact=FALSE)$p.value,
                         wilcox.test(subset(acetylglucosamine, acetylglucosamine$abx=='none')[,1], 
                                     subset(acetylglucosamine, acetylglucosamine$abx=='germfree')[,1], exact=FALSE)$p.value), method='holm')
mannitol_sorbitol_p <- p.adjust(c(wilcox.test(subset(mannitol_sorbitol, abx=='none')[,1], 
                                              subset(mannitol_sorbitol, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                  wilcox.test(subset(mannitol_sorbitol, abx=='none')[,1], 
                                              subset(mannitol_sorbitol, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                  wilcox.test(subset(mannitol_sorbitol, abx=='none')[,1], 
                                              subset(mannitol_sorbitol, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                  wilcox.test(subset(mannitol_sorbitol, abx=='none')[,1], 
                                              subset(mannitol_sorbitol, abx=='germfree')[,1], exact=FALSE)$p.value), method='holm')
galactitol_p <- p.adjust(c(wilcox.test(subset(galactitol, abx=='none')[,1], subset(galactitol, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(galactitol, abx=='none')[,1], subset(galactitol, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(galactitol, abx=='none')[,1], subset(galactitol, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(galactitol, abx=='none')[,1], subset(galactitol, abx=='germfree')[,1], exact=FALSE)$p.value),
                         method='holm')
salicylate_p <- p.adjust(c(wilcox.test(subset(salicylate, abx=='none')[,1], subset(salicylate, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(salicylate, abx=='none')[,1], subset(salicylate, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(salicylate, abx=='none')[,1], subset(salicylate, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(salicylate, abx=='none')[,1], subset(salicylate, abx=='germfree')[,1], exact=FALSE)$p.value),
                         method='holm')
acetylneuriminate_p <- p.adjust(c(wilcox.test(subset(acetylneuriminate, abx=='none')[,1], 
                                              subset(acetylneuriminate, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                  wilcox.test(subset(acetylneuriminate, abx=='none')[,1], 
                                              subset(acetylneuriminate, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                  wilcox.test(subset(acetylneuriminate, abx=='none')[,1], 
                                              subset(acetylneuriminate, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                  wilcox.test(subset(acetylneuriminate, abx=='none')[,1], 
                                              subset(acetylneuriminate, abx=='germfree')[,1], exact=FALSE)$p.value), method='holm')
proline_p <- p.adjust(c(wilcox.test(subset(proline, abx=='none')[,1], subset(proline, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                        wilcox.test(subset(proline, abx=='none')[,1], subset(proline, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                        wilcox.test(subset(proline, abx=='none')[,1], subset(proline, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                        wilcox.test(subset(proline, abx=='none')[,1], subset(proline, abx=='germfree')[,1], exact=FALSE)$p.value), method='holm')
ribose_p <- p.adjust(c(wilcox.test(subset(ribose, abx=='none')[,1], subset(ribose, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                       wilcox.test(subset(ribose, abx=='none')[,1], subset(ribose, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                       wilcox.test(subset(ribose, abx=='none')[,1], subset(ribose, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                       wilcox.test(subset(ribose, abx=='none')[,1], subset(ribose, abx=='germfree')[,1], exact=FALSE)$p.value), method='holm')
lysine_p <- p.adjust(c(wilcox.test(subset(lysine, abx=='none')[,1], subset(lysine, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                       wilcox.test(subset(lysine, abx=='none')[,1], subset(lysine, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                       wilcox.test(subset(lysine, abx=='none')[,1], subset(lysine, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                       wilcox.test(subset(lysine, abx=='none')[,1], subset(lysine, abx=='germfree')[,1], exact=FALSE)$p.value), method='holm')
isoleucine_p <- p.adjust(c(wilcox.test(subset(isoleucine, abx=='none')[,1], subset(isoleucine, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(isoleucine, abx=='none')[,1], subset(isoleucine, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(isoleucine, abx=='none')[,1], subset(isoleucine, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(isoleucine, abx=='none')[,1], subset(isoleucine, abx=='germfree')[,1], exact=FALSE)$p.value),
                         method='holm')
leucine_p <- p.adjust(c(wilcox.test(subset(leucine, abx=='none')[,1], subset(leucine, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                        wilcox.test(subset(leucine, abx=='none')[,1], subset(leucine, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                        wilcox.test(subset(leucine, abx=='none')[,1], subset(leucine, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                        wilcox.test(subset(leucine, abx=='none')[,1], subset(leucine, abx=='germfree')[,1], exact=FALSE)$p.value), method='holm')
methionine_p <- p.adjust(c(wilcox.test(subset(methionine, abx=='none')[,1], subset(methionine, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(methionine, abx=='none')[,1], subset(methionine, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(methionine, abx=='none')[,1], subset(methionine, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                           wilcox.test(subset(methionine, abx=='none')[,1], subset(methionine, abx=='germfree')[,1], exact=FALSE)$p.value),
                         method='holm')
threonine_p <- p.adjust(c(wilcox.test(subset(threonine, abx=='none')[,1], subset(threonine, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                          wilcox.test(subset(threonine, abx=='none')[,1], subset(threonine, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                          wilcox.test(subset(threonine, abx=='none')[,1], subset(threonine, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                          wilcox.test(subset(threonine, abx=='none')[,1], subset(threonine, abx=='germfree')[,1], exact=FALSE)$p.value),
                        method='holm')
glycine_p <- p.adjust(c(wilcox.test(subset(glycine, abx=='none')[,1], subset(glycine, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                          wilcox.test(subset(glycine, abx=='none')[,1], subset(glycine, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                          wilcox.test(subset(glycine, abx=='none')[,1], subset(glycine, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                          wilcox.test(subset(glycine, abx=='none')[,1], subset(glycine, abx=='germfree')[,1], exact=FALSE)$p.value),
                        method='holm')

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_7.pdf'
select_palette <- c('gray', wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen')
pdf(file=plot_file, width=12, height=8)
layout(matrix(c(1,1,2,2,
                3,3,4,4,
                5,5,6,7),
              nrow=3, ncol=4, byrow = TRUE))

#--------------------------------#

# Acetylglucosamine
par(las=1, mar=c(0.3,4,1,1), mgp=c(2.5,0.7,0))
boxplot(Nacetylglucosamine_Nacetylgalactosamine~abx, data=acetylglucosamine, col=select_palette, ylim=c(0,14), whisklty=1, cex.lab=1.1,
           xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('A', side=2, line=2, las=2, adj=1.7, padj=-9, cex=1.1)
legend('topleft', 'N-Acetylglucosamine + N-Acetylgalactosamine', bty='n')
axis(2, at=seq(0,14,3.5), labels=c('0.0','3.5','7.0','10.5','14.0'))
text(c(2,3,4,5), c(7,2.3,1.6,2.3), labels=c('*','*','*','*'), cex=2.5)

#--------------------------------#

# Salicylate
par(las=1, mar=c(4,3.5,0.3,1), mgp=c(2.5,0.7,0))
boxplot(salicylate~abx, data=salicylate, col=select_palette, ylim=c(0,5.2), whisklty=1, cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('B', side=2, line=2, las=2, adj=1.3, padj=-7, cex=1.1)
legend('topleft', 'Salicylate', bty='n')
axis(2, at=seq(0,5.2,1.3), labels=c('0.0','1.3','2.6','3.9','5.2'))
text(c(2,3,4,5), c(2.8,1.8,1.9,0.6), labels=c('*','*','*','*'), cex=2.5)
mtext(c('No Antibiotics\nSPF','Streptomycin\nSPF','Cefoperazone\nSPF','Clindamycin\nSPF','No Antibiotics\nGF'), side=1, at=c(1:5), cex=0.95, padj=1)
mtext('Treatment:', side=1, at=0.14, padj=1.3, cex=0.7)
mtext('Mice:', side=1, at=0.12, padj=3.5, cex=0.7)

#--------------------------------#

# Galactitol
par(las=1, mar=c(0.3,3.5,1,1), mgp=c(2.5,0.7,0))
boxplot(galactitol~abx, data=galactitol, col=select_palette, ylim=c(0,3), whisklty=1, cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('C', side=2, line=2, las=2, adj=1.4, padj=-9, cex=1.1)
legend('topleft', 'Galactitol', bty='n')
axis(2, at=seq(0,3,0.75), labels=c('0.0','0.75','1.5','2.25','3.0'))
text(c(2,3,5), c(1.3,2.7,2.35), labels=c('*','*'), cex=2.5)

#--------------------------------#

# Mannitol - Sorbitol
par(las=1, mar=c(0.3,4,0.3,1), mgp=c(2.5,0.7,0))
boxplot(mannitol_sorbitol~abx, data=mannitol_sorbitol, col=select_palette, ylim=c(0,50), whisklty=1, cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('D', side=2, line=2, las=2, adj=1.7, padj=-9, cex=1.1)
legend('topleft', 'Mannitol + Sorbitol', bty='n')
axis(2, at=seq(0,50,12.5), labels=c('0.0','12.5','25.0','37.5','50.0'))
text(c(3,5), c(47,43), labels=c('*','*'), cex=2.5)


#--------------------------------#

# Acetylneuriminate
par(las=1, mar=c(4,4,0.3,1), mgp=c(2.5,0.7,0))
boxplot(Nacetylneuraminate~abx, data=acetylneuriminate, col=select_palette, ylim=c(0,3), whisklty=1, cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('E', side=2, line=2, las=2, adj=1.7, padj=-7, cex=1.1)
legend('topleft', 'N-Acetylneuriminate', bty='n')
axis(2, at=seq(0,3,0.75), labels=c('0.0','0.75','1.5','2.25','3.0'))
text(c(3,5), c(2.2,2.3), labels=c('*','*'), cex=2.5)
mtext(c('No Antibiotics\nSPF','Streptomycin\nSPF','Cefoperazone\nSPF','Clindamycin\nSPF','No Antibiotics\nGF'), side=1, at=c(1:5), cex=0.95, padj=1)
mtext('Treatment:', side=1, at=0.14, padj=1.3, cex=0.7)
mtext('Mice:', side=1, at=0.12, padj=3.5, cex=0.7)

#--------------------------------#

# Proline
par(las=1, mar=c(4,3.5,0.5,1), mgp=c(2.5,0.7,0))
boxplot(proline~abx, data=proline, col=select_palette, whisklty=1, ylim=c(0,3), cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
legend('topleft', 'Proline', bty='n')
mtext('F', side=2, line=2, las=2, adj=1.3, padj=-7, cex=1.1)
axis(2, at=seq(0,3,0.75), labels=c('0.0','0.75','1.5','2.25','3.0'))
text(c(2,3,4,5), c(3.0,2.8,2.3,2.85), labels=c('*','*','*','*'), cex=2.5)
mtext(c('No Abx\nSPF','Strep\nSPF','Cef\nSPF','Clinda\nSPF','No Abx\nGF'), side=1, at=c(1:5), cex=0.75, padj=1)

# Glycine
par(las=1, mar=c(4,3.5,0.5,1), mgp=c(2.5,0.7,0))
boxplot(glycine~abx, data=glycine, col=select_palette, whisklty=1, ylim=c(0,3), cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
legend('topleft', 'Glycine', bty='n')
mtext('G', side=2, line=2, las=2, adj=1.3, padj=-7, cex=1.1)
axis(2, at=seq(0,3,0.75), labels=c('0.0','0.75','1.5','2.25','3.0'))
text(c(3,5), c(2.5,2.1), labels=c('*','*'), cex=2.5)
mtext(c('No Abx\nSPF','Strep\nSPF','Cef\nSPF','Clinda\nSPF','No Abx\nGF'), side=1, at=c(1:5), cex=0.75, padj=1)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
rm(plot_file, select_palette, acetylglucosamine, mannitol_sorbitol, galactitol, salicylate, acetylneuriminate, 
   proline, ribose, lysine, isoleucine, leucine, methionine, threonine)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

