
# Load dependencies
deps <- c('wesanderson','vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,'Package'] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

#--------------------------------------------------------------------------------------------------------------#

# Define variables
cefoperazone_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/select_genes/cefoperazone_630.RNA_reads2select.all.norm.txt'
clindamycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/select_genes/clindamycin_630.RNA_reads2select.all.norm.txt'
streptomycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/select_genes/streptomycin_630.RNA_reads2select.all.norm.txt'
germfree_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/select_genes/germfree.RNA_reads2select.all.norm.txt'

# Open files
cefoperazone <- read.delim(cefoperazone_file, sep='\t', header=FALSE)
colnames(cefoperazone) <- c('gene', 'Cefoperazone')
clindamycin <- read.delim(clindamycin_file, sep='\t', header=FALSE)
colnames(clindamycin) <- c('gene', 'Clindamycin')
streptomycin <- read.delim(streptomycin_file, sep='\t', header=FALSE)
colnames(streptomycin) <- c('gene', 'Streptomycin')
germfree <- read.delim(germfree_file, sep='\t', header=FALSE)
colnames(germfree) <- c('gene', 'Germfree')

# Clean up
rm(cefoperazone_file, clindamycin_file, streptomycin_file, germfree_file)

# Merge tables
combined_mapping <- merge(streptomycin, cefoperazone, by='gene')
combined_mapping <- merge(combined_mapping, clindamycin, by='gene')
combined_mapping <- merge(combined_mapping, germfree, by='gene')
combined_mapping$gene <- gsub('Clostridium_difficile_630\\|','', combined_mapping$gene)
combined_mapping$gene <- gsub('ENA\\|CDT20869\\|CDT20869.1\\|Clostridium_difficile_putative_phage_replication_protein_','', combined_mapping$gene)
combined_mapping$gene <- gsub('_',' ', combined_mapping$gene)
rownames(combined_mapping) <- combined_mapping$gene
combined_mapping$gene <- NULL
rm(cefoperazone, clindamycin, streptomycin, germfree)

#--------------------------------------------------------------------------------------------------------------#

# Define subsets of interest
sigma_keep <- c('CodY', 'CcpA', 'SigH', 'Spo0A', 'PrdR', 'Rex')
paloc_keep <- c('TcdR','TcdC','TcdE','CdtR','TcdA','TcdB')
sporulation_keep <- c('SpoIIAB','SpoIIE','SpoVS','SpoIVA','SpoVFB','SpoVB','SpoVG','CdeC','CotJB2','CotD',
                      'SspA','SspB','SigK','SigF') # variance > 0.01
quorum_keep <- c('LuxS', 'AgrD', 'AgrB')

# Pull of the genes of interest
sigma <- subset(combined_mapping, rownames(combined_mapping) %in% sigma_keep)
paloc <- subset(combined_mapping, rownames(combined_mapping) %in% paloc_keep)
sporulation <- subset(combined_mapping, rownames(combined_mapping) %in% sporulation_keep)
quorum <- subset(combined_mapping, rownames(combined_mapping) %in% quorum_keep)
combined_mapping <- rbind(sigma, paloc, sporulation, quorum)

# Rarefy mappings to be equal within grouping
abx_sub_size <- round(min(colSums(combined_mapping[,1:3])) * 0.9) # 2037
gf_sub_size <- round(sum(combined_mapping[,4]) * 0.9) # 303
cefoperazone_sub <- t(rrarefy(combined_mapping$Cefoperazone, sample=abx_sub_size))
clindamycin_sub <- t(rrarefy(combined_mapping$Clindamycin, sample=abx_sub_size))
streptomycin_sub <- t(rrarefy(combined_mapping$Streptomycin, sample=abx_sub_size))
germfree_sub <- t(rrarefy(combined_mapping$Germfree, sample=gf_sub_size))
for (index in 1:999) {
  cefoperazone_sub <- cbind(cefoperazone_sub, t(rrarefy(combined_mapping$Cefoperazone, sample=abx_sub_size)))
  clindamycin_sub <- cbind(clindamycin_sub, t(rrarefy(combined_mapping$Clindamycin, sample=abx_sub_size)))
  streptomycin_sub <- cbind(streptomycin_sub, t(rrarefy(combined_mapping$Streptomycin, sample=abx_sub_size)))
  germfree_sub <- cbind(germfree_sub, t(rrarefy(combined_mapping$Germfree, sample=gf_sub_size)))
}
combined_mapping$Cefoperazone <- apply(cefoperazone_sub, 1, median)
combined_mapping$Clindamycin <- apply(clindamycin_sub, 1, median)
combined_mapping$Streptomycin <- apply(streptomycin_sub, 1, median)
combined_mapping$Germfree <- apply(germfree_sub, 1, median)
rm(index, abx_sub_size, gf_sub_size, cefoperazone_sub, clindamycin_sub, streptomycin_sub, germfree_sub)

# Convert to relative abundances
totals <- colSums(combined_mapping)
combined_mapping$Streptomycin <- combined_mapping$Streptomycin / totals[1]
combined_mapping$Cefoperazone <- combined_mapping$Cefoperazone / totals[2]
combined_mapping$Clindamycin <- combined_mapping$Clindamycin / totals[3]
combined_mapping$Germfree <- combined_mapping$Germfree / totals[4]
rm(totals)

combined_mapping <- combined_mapping * 10

# Separate the mappings again

# Sigma factors
# Integration of Metabolism and Virulence by Clostridium difficile CodY
# Global transcriptional control by glucose and carbon regulator CcpA in Clostridium difficile.
# Proline-Dependent Regulation of Clostridium difficile Stickland Metabolism
# The Clostridium difficile spo0A Gene Is a Persistence and Transmission Factor
sigma <- subset(combined_mapping, rownames(combined_mapping) %in% sigma_keep)
rownames(sigma) <- c('ccpA', 'codY', 'prdR', 'rex', 'sigH', 'spo0A')
sigma <- t(sigma)

# PaLoc
paloc <- subset(combined_mapping, rownames(combined_mapping) %in% paloc_keep)
rownames(paloc) <- c('cdtR', 'tcdA', 'tcdB', 'tcdC', 'tcdE', 'tcdR')
paloc <- t(paloc)

# Quorum sensing
quorum <- subset(combined_mapping, rownames(combined_mapping) %in% quorum_keep)
rownames(quorum) <- c('agrB', 'agrD', 'luxS')
quorum <- t(quorum)

# Sporulation
# Genes with >0.01 variance included in final analysis
sporulation <- subset(combined_mapping, rownames(combined_mapping) %in% sporulation_keep)
rownames(sporulation) <- c('cdeC','cotD','cotJB2','sigF','sigK','spoIIAB','spoIIE',
                           'spoIVA','spoVB','spoVFB','spoVG','spoVS','sspA','sspB')
sporulation <- t(sporulation)

# Clean up
rm(combined_mapping, sigma_keep, paloc_keep, sporulation_keep, quorum_keep)

#--------------------------------------------------------------------------------------------------------------#

# Set the color palette and plotting environment
select_palette <- c(wes_palette('FantasticFox')[1], wes_palette('FantasticFox')[3], wes_palette('FantasticFox')[5], 'forestgreen')
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_3.pdf'
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
pdf(file=plot_file, width=14, height=10)
layout(matrix(c(1,1,2,2,3,
                1,1,2,2,3,
                4,5,5,5,5,
                4,5,5,5,5),
              nrow=4, ncol=5, byrow = TRUE))

#--------------------------------------------------------------------------------------------------------------#

# Sigma factors
par(las=1, mar=c(4,5.4,1,1), mgp=c(3.9, 1, 0))
barplot(cbind(sigma[,1], sigma[,2], sigma[,6], sigma[,3], sigma[,4], sigma[,5]), 
        col=select_palette, space=c(0,1.5), beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,2.5), cex.lab=1.2)
abline(h=c(0.833, 1.666), lty=2)
barplot(cbind(sigma[,1], sigma[,2], sigma[,6], sigma[,3], sigma[,4], sigma[,5]), 
        col=select_palette, space=c(0,1.5), beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,2.5), add=TRUE, cex.lab=1.2)
box()
axis(side=2, at=c(0,0.833,1.666,2.5), c('0.0%','8.33%','16.67%','25.0%'), tick=TRUE, las=1, cex.axis=1)
#legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, cex=1.2,
#       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,33,5.5), y=par()$usr[3]-0.035*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(c('ccpA', 'codY', 'spo0A', 'prdR', 'rex', 'sigH')), srt=45, adj=1, xpd=TRUE, cex=1.6)
legend('topright', legend='Sigma factors', pt.cex=0, bty='n', cex=1.8)
mtext('a', side=2, line=2, las=2, adj=3.3, padj=-16, cex=1.1, font=2)

# Pathogenicity
par(las=1, mar=c(4,5,1,1), mgp=c(3.9, 1, 0))
barplot(cbind(paloc[,2], paloc[,3], paloc[,4], paloc[,5], paloc[,6], paloc[,1]), 
        col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
                    ylab='Relative Transcript Abundance', ylim=c(0,0.2), cex.lab=1.2)
abline(h=c(.0666,.1333), lty=2)
barplot(cbind(paloc[,2], paloc[,3], paloc[,4], paloc[,5], paloc[,6], paloc[,1]), 
        col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,0.2), add=TRUE, cex.lab=1.2)
box()
axis(side=2, at=c(0,.0666,.1333,0.2), c('0.0%','0.67%','1.33%','2.0%'), tick=TRUE, las=1, cex.axis=1)
#legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, cex=1.2,
#       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,33,5.5), y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(c('tcdA', 'tcdB', 'tcdC', 'tcdE', 'tcdR', 'cdtR')), srt=45, adj=1, xpd=TRUE, cex=1.6)
legend('topright', legend='Pathogenicity', pt.cex=0, bty='n', cex=1.8)
mtext('b', side=2, line=2, las=2, adj=3.3, padj=-16, cex=1.1, font=2)

# Quorum sensing
par(las=1, mar=c(4,5,1,1), mgp=c(3.9, 1, 0))
barplot(quorum, col=select_palette, beside=TRUE, xaxt='n', yaxt='n', 
                    ylab='Relative Transcript Abundance', ylim=c(0,0.25), cex.lab=1.2)
abline(h=c(.0833,.1667), lty=2)
barplot(quorum, col=select_palette, beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,0.2), add=TRUE, cex.lab=1.2)
box()
axis(side=2, at=c(0,.0833,.1667,0.25), c('0.0%','0.83%','1.67%','2.5%'), tick=TRUE, las=1, cex.axis=1)
#legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, cex=1.2,
#       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=c(2.7,8.2,13.7), y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(colnames(quorum)), srt=45, adj=1, xpd=TRUE, cex=1.6)
legend('topright', legend='Quorum sensing', pt.cex=0, bty='n', cex=1.8)
mtext('c', side=2, line=2, las=2, adj=3.3, padj=-16, cex=1.1, font=2)

# Legend plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
legend('center', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germ free'), pt.cex=3.6, cex=2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)

# Sporulation
# Reordered by stage of sporulation
par(las=1, mar=c(6.4,5.4,1,1), mgp=c(3.9, 1, 0))
barplot(cbind(sporulation[,4], sporulation[,6], sporulation[,7], rep(NA,4),
              sporulation[,8], sporulation[,9], sporulation[,13], sporulation[,14], rep(NA,4),
              sporulation[,1], sporulation[,2], sporulation[,3], sporulation[,5], sporulation[,10], sporulation[,11], sporulation[,12]), 
        col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,3), cex.lab=1.2)
abline(h=c(1:2), lty=2)
barplot(cbind(sporulation[,4], sporulation[,6], sporulation[,7], rep(NA,4),
              sporulation[,8], sporulation[,9], sporulation[,13], sporulation[,14], rep(NA,4),
              sporulation[,1], sporulation[,2], sporulation[,3], sporulation[,5], sporulation[,10], sporulation[,11], sporulation[,12]), 
        col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,3), add=TRUE, cex.lab=1.2)
box()
axis(side=2, at=c(0:3), c('0.0%','10.0%','20.0%','30.0%'), tick=TRUE, las=1, cex.axis=1)
#legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, cex=1.2,
#       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,88,5.5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(c('sigF', 'spoIIAB', 'spoIIE', '',
                          'spoIVA', 'spoVB', 'sspA', 'sspB', '',
                          'cdeC', 'cotD', 'cotJB2', 'sigK', 'spoVFB', 'spoVG', 'spoVS')), 
     srt=45, adj=1, xpd=TRUE, cex=1.5)
legend('topright', legend='Sporulation', pt.cex=0, bty='n', cex=1.8)
mtext('d', side=2, line=2, las=2, adj=3.3, padj=-15, cex=1.1, font=2)

# Add groups for sporulation stages
segments(x0=c(2,24,51.5), y0=par()$usr[3]-0.155*(par()$usr[4]-par()$usr[3]), x1=c(16,43.5,87.5), y1=par()$usr[3]-0.155*(par()$usr[4]-par()$usr[3]), lwd=2, xpd=TRUE)
text(x=c(8,17,54.5,77), y=par()$usr[3]-0.185*(par()$usr[4]-par()$usr[3]), labels=c('Stage:','Early','Intermediate','Late'), adj=3, xpd=TRUE, cex=1.5)


dev.off()

#--------------------------------------------------------------------------------------------------------------#

# Clean up
rm(quorum, sigma, sporulation, paloc, plot_file, select_palette, make.italic)
for (dep in deps){
  pkg <- paste('package:', dep,sep='')
   detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

