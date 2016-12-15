
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
acetylglucosamine <- concentrations[,c(1,13)]
mannitol_sorbitol <- concentrations[,c(3,13)]
galactitol <- concentrations[,c(4,13)]
salicylate <- concentrations[,c(5,13)]
acetylneuriminate <- concentrations[,c(6,13)]

# Subset each metabolites - supplement
proline <- concentrations[,c(2,13)]
ribose <- concentrations[,c(7,13)]
lysine <- concentrations[,c(8,13)]
isoleucine <- concentrations[,c(9,13)]
leucine <- concentrations[,c(10,13)]
methionine <- concentrations[,c(11,13)]
threonine <- concentrations[,c(12,13)]

# Remove raw data
rm(concentrations, metadata)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_7.pdf'
select_palette <- c('gray', wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen')
pdf(file=plot_file, width=7, height=10)
layout(matrix(c(1,
                2,
                3,
                4,
                5),
              nrow=5, ncol=1, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# A. Acetylglucosamine
par(las=1, mar=c(0.5,4,1,1), mgp=c(2.5,0.7,0))
boxplot(Nacetylglucosamine_Nacetylgalactosamine~abx, data=acetylglucosamine, col=select_palette, ylim=c(0,14), whisklty=1,
           xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('A', side=2, line=2, las=2, adj=1.7, padj=-6.5, cex=1.1)
legend('topright', legend='N-Acetylglucosamine / N-Acetylgalactosamine', pt.cex=0, bty='n', cex=1.2)
axis(2, at=seq(0,14,3.5), labels=c('0.0','3.5','7.0','10.5','14.0'))
text(c(2,3,4,5), c(7,2.3,1.6,2.3), labels=c('*','*','*','*'), cex=2.5)

#--------------------------------#

# B. Galactitol
par(las=1, mar=c(0.5,4,0.5,1), mgp=c(2.5,0.7,0))
boxplot(galactitol~abx, data=galactitol, col=select_palette, ylim=c(0,3), whisklty=1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('B', side=2, line=2, las=2, adj=1.7, padj=-6.5, cex=1.1)
legend('topright', legend='Galactitol', pt.cex=0, bty='n', cex=1.2)
axis(2, at=seq(0,3,0.75), labels=c('0.0','0.75','1.5','2.25','3.0'))
text(c(2,3,5), c(1.3,2.7,2.35), labels=c('*','*'), cex=2.5)

#--------------------------------#

# C. Mannitol - Sorbitol
par(las=1, mar=c(0.5,4,0.5,1), mgp=c(2.5,0.7,0))
boxplot(mannitol_sorbitol~abx, data=mannitol_sorbitol, col=select_palette, ylim=c(0,50), whisklty=1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('C', side=2, line=2, las=2, adj=1.7, padj=-6.6, cex=1.1)
legend('topright', legend='Mannitol / Sorbitol', pt.cex=0, bty='n', cex=1.2)
axis(2, at=seq(0,50,12.5), labels=c('0.0','12.5','25.0','37.5','50.0'))
text(c(3,5), c(47,43), labels=c('*','*'), cex=2.5)

#--------------------------------#

# D. Salicylate
par(las=1, mar=c(0.5,4,0.5,1), mgp=c(2.5,0.7,0))
boxplot(salicylate~abx, data=salicylate, col=select_palette, ylim=c(0,5), whisklty=1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('D', side=2, line=2, las=2, adj=1.7, padj=-6.6, cex=1.1)
legend('topright', legend='Salicylate', pt.cex=0, bty='n', cex=1.2)
axis(2, at=seq(0,5,1.25), labels=c('0.0','1.25','2.5','3.75','5.0'))
text(c(2,3,4,5), c(2.8,1.8,1.9,0.6), labels=c('*','*','*','*'), cex=2.5)

#--------------------------------#

# E. Acetylneuriminate
par(las=1, mar=c(3,4,0.5,1), mgp=c(2.5,0.7,0))
boxplot(Nacetylneuraminate~abx, data=acetylneuriminate, col=select_palette, ylim=c(0,3), whisklty=1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('E', side=2, line=2, las=2, adj=1.7, padj=-5.5, cex=1.1)
legend('topright', legend='N-Acetylneuriminate', pt.cex=0, bty='n', cex=1.2)
axis(2, at=seq(0,3,0.75), labels=c('0.0','0.75','1.5','2.25','3.0'))
text(c(3,5), c(2.2,2.3), labels=c('*','*'), cex=2.5)

axis(1, at=c(1:5), tick=FALSE, labels=c('No Antibiotics (SPF)','Streptomycin (SPF)','Cefoperazone (SPF)','Clindamycin (SPF)','No Antibiotics (GF)'), cex=1.2)

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
dev.off()
rm(plot_file, select_palette, acetylglucosamine, mannitol_sorbitol, galactitol, salicylate, acetylneuriminate, 
   proline, ribose, lysine, isoleucine, leucine, methionine, threonine)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

