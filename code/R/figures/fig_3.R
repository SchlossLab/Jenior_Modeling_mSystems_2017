
# Load dependencies
deps <- c('wesanderson','vegan', 'matrixStats', 'plotrix');
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

# Define subsets of interest
sigma_keep <- c('CodY', 'CcpA', 'SigH', 'Spo0A', 'PrdR', 'Rex')
paloc_keep <- c('TcdR','TcdC','TcdE','CdtR','TcdA','TcdB')
sporulation_keep <- c('DpaA', 'SpoIID','SpoIIID','SpoIIAA','SpoIIAB','SpoIIIAA','SpoIIIAB','SpoIIIAC','SpoIIIAD',
                      'SpoIIIAE','SpoIIIAG','SpoIIIAH','SpoIIP','SpoIIGA','SpoIIE','SpoIIR','SpoVAC','SpoVAD',
                      'SpoVAE','SpoIVB2','SpoVS','SpoIV','SpoIVA','SpoVE','SpoVD','SpoVFB','SpoVFA','SpoVB',
                      'SpoVT','SpoVG','CD1492','CD2492','CdeC','CotA','SodA','CotJB2','CotD',
                      'Gpr','SspA','SspB','BclA3', 'SigA2','SigK','SigE','SigF','SigG')
quorum_keep <- c('LuxS', 'AgrD', 'AgrB')

# Pull of the genes of interest
sigma <- subset(combined_mapping, rownames(combined_mapping) %in% sigma_keep)
paloc <- subset(combined_mapping, rownames(combined_mapping) %in% paloc_keep)
sporulation <- subset(combined_mapping, rownames(combined_mapping) %in% sporulation_keep)
quorum <- subset(combined_mapping, rownames(combined_mapping) %in% quorum_keep)
combined_mapping <- rbind(sigma, paloc, sporulation, quorum)

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
# The Key Sigma Factor of Transition Phase, SigH, Controls Sporulation, Metabolism, and Virulence Factor Expression in Clostridium difficile
sigma <- subset(combined_mapping, rownames(combined_mapping) %in% sigma_keep)
rownames(sigma) <- c('ccpA', 'codY', 'prdR', 'rex', 'sigH', 'spo0A')

# PaLoc
paloc <- subset(combined_mapping, rownames(combined_mapping) %in% paloc_keep)
rownames(paloc) <- c('cdtR', 'tcdA', 'tcdB', 'tcdC', 'tcdE', 'tcdR')

# Sporulation
# citation for genes
sporulation <- subset(combined_mapping, rownames(combined_mapping) %in% sporulation_keep)
rownames(sporulation) <- c('bclA3', 'CD1492', 'CD2492', 'cdeC', 'cotA', 
                                   'cotD', 'cotJB2', 'dpaA', 'gpr', 'sigA2','sigK','sigE','sigF','sigG', 'sodA', 'spoIIAA', 'spoIIAB', 'spoIID', 
                                   'spoIIE', 'spoIIGA', 'spoIIIAA', 'spoIIIAB', 'spoIIIAC', 'spoIIIAD', 'spoIIIAE', 
                                   'spoIIIAG', 'spoIIIAH', 'spoIIID', 'spoIIP', 'spoIIR', 'spoIV', 'spoIVA',  
                                   'spoIVB2', 'spoVAC', 'spoVAD', 'spoVAE', 'spoVB', 'spoVD', 'spoVE', 'spoVFA', 
                                   'spoVFB', 'spoVG', 'spoVS', 'spoVT', 'sspA', 'sspB')

# Quorum sensing
quorum <- subset(combined_mapping, rownames(combined_mapping) %in% quorum_keep)
rownames(quorum) <- c('agrB', 'agrD', 'luxS')

# Clean up
rm(combined_mapping, sigma_keep, paloc_keep, sporulation_keep, quorum_keep)

#--------------------------------------------------------------------------------------------------------------#

# Set the color palette and plotting environment
select_palette <- c(wes_palette('FantasticFox')[1], wes_palette('FantasticFox')[3], wes_palette('FantasticFox')[5], 'forestgreen')
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_3.pdf'
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
pdf(file=plot_file, width=14, height=10)
layout(matrix(c(1,1,2,2,
                1,1,2,2,
                3,3,3,3,
                3,3,3,3),
              nrow=4, ncol=4, byrow = TRUE))

#--------------------------------------------------------------------------------------------------------------#

# Sigma factors
par(las=1, mar=c(4,5.4,1,1), mgp=c(3.2, 1, 0))
x_coords <- barplot(t(sigma), col=select_palette, space=c(0,1.5), beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,2.5))
abline(h=1.25, lty=2)
barplot(t(sigma), col=select_palette, space=c(0,1.5), beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,2.5), add=TRUE)
box()
axis(side=2, at=c(0,1.25,2.5), c('','',''), tick=TRUE, las=1)
mtext(text=c('0.0%','12.5%','25.0%'), side=2, at=c(0,1.25,2.5), cex=c(0.7,0.6,0.7), adj=c(1.8,1.4,1.5))
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,36.7,5.5), y=par()$usr[3]-0.035*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(rownames(sigma)), srt=45, adj=1, xpd=TRUE, cex=1.2)
legend('topright', legend='Sigma factors', pt.cex=0, bty='n', cex=1.8)
mtext('A', side=2, line=2, las=2, adj=2.3, padj=-11, cex=1.5)

# Pathogenicity
par(las=1, mar=c(4,5,1,1), mgp=c(3.2, 1, 0))
x_coords <- barplot(t(paloc), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
                    ylab='Relative Transcript Abundance', ylim=c(0,0.2))
abline(h=0.1, lty=2)
barplot(t(paloc), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,0.2), add=TRUE)
box()
axis(side=2, at=c(0,0.1,0.2), c('0.0%','1.0%','2.0%'), tick=TRUE, las=1)
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,36.7,5.5), y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(rownames(paloc)), srt=45, adj=1, xpd=TRUE, cex=1.4)
legend('topright', legend='Pathogenicity', pt.cex=0, bty='n', cex=1.8)
mtext('B', side=2, line=2, las=2, adj=1.9, padj=-11, cex=1.5)

# Sporulation
par(las=1, mar=c(4,5.4,1,1), mgp=c(3.5, 1, 0))
x_coords <- barplot(t(sporulation), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,3))
abline(h=c(1:3), lty=2)
barplot(t(sporulation), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,3), add=TRUE)
box()
axis(side=2, at=c(0:3), c('0.0%','10.0%','20.0%','30.0%'), tick=TRUE, las=1)
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,251.2,5.5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(rownames(sporulation)), srt=45, adj=1, xpd=TRUE, cex=0.9)
legend('topright', legend='Sporulation', pt.cex=0, bty='n', cex=1.8)
mtext('C', side=2, line=2, las=2, adj=2.2, padj=-11, cex=1.5)

dev.off()

#--------------------------------------------------------------------------------------------------------------#

plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_quorum.pdf'
pdf(file=plot_file, width=7, height=7)

# Quorum sensing
par(las=1, mar=c(4,5,1,1), mgp=c(3, 1, 0))
x_coords <- barplot(t(quorum), col=select_palette, beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,1))
abline(h=0.5, lty=2)
barplot(t(quorum), col=select_palette, beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,1), add=TRUE)
box()
axis(side=2, at=c(0:1), c('0.0%','10.0%'), tick=TRUE, las=1, cex=1.7)
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=c(2.7,8.2,13.7), y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(rownames(quorum)), srt=45, adj=1, xpd=TRUE, cex=1.6)
legend('topright', legend='Quorum sensing', pt.cex=0, bty='n', cex=1.8)

dev.off()

#--------------------------------------------------------------------------------------------------------------#

# Clean up
rm(quorum, sigma, sporulation, paloc, plot_file, select_palette, x_coords, make.italic)
for (dep in deps){
  pkg <- paste('package:', dep,sep='')
   detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

