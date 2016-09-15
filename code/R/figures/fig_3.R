
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

# Convert to relative abundances
totals <- colSums(combined_mapping)
combined_mapping$Streptomycin <- combined_mapping$Streptomycin / totals[1]
combined_mapping$Cefoperazone <- combined_mapping$Cefoperazone / totals[2]
combined_mapping$Clindamycin <- combined_mapping$Clindamycin / totals[3]
combined_mapping$Germfree <- combined_mapping$Germfree / totals[4]
rm(totals)

#--------------------------------------------------------------------------------------------------------------#

# Break up the data and calculate stats

# Sigma factors
# Integration of Metabolism and Virulence by Clostridium difficile CodY
# Global transcriptional control by glucose and carbon regulator CcpA in Clostridium difficile.
# Proline-Dependent Regulation of Clostridium difficile Stickland Metabolism
# The Clostridium difficile spo0A Gene Is a Persistence and Transmission Factor
# The Key Sigma Factor of Transition Phase, SigH, Controls Sporulation, Metabolism, and Virulence Factor Expression in Clostridium difficile
sigma_keep <- c('CodY', 'CcpA', 'SigH', 'Spo0A', 'PrdR', 'Rex')
sigma <- subset(combined_mapping, rownames(combined_mapping) %in% sigma_keep)
rownames(sigma) <- c('ccpA', 'codY', 'prdR', 'rex', 'sigH', 'spo0A')

# Iteratively rarefy mappings
#sub_size <- round(min(colSums(sigma[,1:4])) * 0.9)
#cefoperazone <- t(rrarefy(sigma$Cefoperazone, sample=sub_size))
#clindamycin <- t(rrarefy(sigma$Clindamycin, sample=sub_size))
#streptomycin <- t(rrarefy(sigma$Streptomycin, sample=sub_size))
#germfree <- t(rrarefy(sigma$Germfree, sample=sub_size))
#for (index in 1:999) {
#  cefoperazone <- cbind(cefoperazone, t(rrarefy(sigma$Cefoperazone, sample=sub_size)))
#  clindamycin <- cbind(clindamycin, t(rrarefy(sigma$Clindamycin, sample=sub_size)))
#  streptomycin <- cbind(streptomycin, t(rrarefy(sigma$Streptomycin, sample=sub_size)))
#  germfree <- cbind(germfree, t(rrarefy(sigma$Germfree, sample=sub_size)))
#}
#rm(index)
# Log transform data
#cefoperazone <- log10(cefoperazone + 1)
#clindamycin <- log10(clindamycin + 1)
#streptomycin <- log10(streptomycin + 1)
#germfree <- log10(germfree + 1)
# Medians
#sigma$Cefoperazone <- rowMedians(cefoperazone)
#sigma$Clindamycin <- rowMedians(clindamycin)
#sigma$Streptomycin <- rowMedians(streptomycin)
#sigma$Germfree <- rowMedians(germfree)
#sigma_medians <- sigma
# SDs
#sigma$Cefoperazone <- rowSds(cefoperazone)
#sigma$Clindamycin <- rowSds(clindamycin)
#sigma$Streptomycin <- rowSds(streptomycin)
#sigma$Germfree <- rowSds(germfree)
#sigma_sds <- sigma * 1.95
#rownames(sigma_sds) <- c('ccpA', 'codY', 'prdR', 'rex', 'sigH', 'spo0A')

################

# PaLoc
paloc_keep <- c('TcdR','TcdC','TcdE','CdtR','TcdA','TcdB')
paloc <- subset(combined_mapping, rownames(combined_mapping) %in% paloc_keep)
rownames(paloc) <- c('cdtR', 'tcdA', 'tcdB', 'tcdC', 'tcdE', 'tcdR')

# Iteratively rarefy mappings
#sub_size <- round(min(colSums(paloc[,1:4])) * 0.9)
#cefoperazone <- t(rrarefy(paloc$Cefoperazone, sample=sub_size))
#clindamycin <- t(rrarefy(paloc$Clindamycin, sample=sub_size))
#streptomycin <- t(rrarefy(paloc$Streptomycin, sample=sub_size))
#germfree <- t(rrarefy(paloc$Germfree, sample=sub_size))
#for (index in 1:999) {
#  cefoperazone <- cbind(cefoperazone, t(rrarefy(paloc$Cefoperazone, sample=sub_size)))
#  clindamycin <- cbind(clindamycin, t(rrarefy(paloc$Clindamycin, sample=sub_size)))
#  streptomycin <- cbind(streptomycin, t(rrarefy(paloc$Streptomycin, sample=sub_size)))
#  germfree <- cbind(germfree, t(rrarefy(paloc$Germfree, sample=sub_size)))
#}
#rm(index)
# Log transform data
#cefoperazone <- log10(cefoperazone + 1)
#clindamycin <- log10(clindamycin + 1)
#streptomycin <- log10(streptomycin + 1)
#germfree <- log10(germfree + 1)
# Medians
#paloc$Cefoperazone <- rowMedians(cefoperazone)
#paloc$Clindamycin <- rowMedians(clindamycin)
#paloc$Streptomycin <- rowMedians(streptomycin)
#paloc$Germfree <- rowMedians(germfree)
#paloc_medians <- paloc
#rownames(paloc_medians) <- c('cdtR', 'tcdA', 'tcdB', 'tcdC', 'tcdE', 'tcdR')
# SDs
#paloc$Cefoperazone <- rowSds(cefoperazone)
#paloc$Clindamycin <- rowSds(clindamycin)
#paloc$Streptomycin <- rowSds(streptomycin)
#paloc$Germfree <- rowSds(germfree)
#paloc_sds <- paloc * 1.95
#rownames(paloc_sds) <- c('cdtR', 'tcdA', 'tcdB', 'tcdC', 'tcdE', 'tcdR')

################

# Sporulation
# citation for genes
sporulation_keep <- c('DpaA', 'SpoIID','SpoIIID','SpoIIAA','SpoIIAB','SpoIIIAA','SpoIIIAB','SpoIIIAC','SpoIIIAD',
                      'SpoIIIAE','SpoIIIAG','SpoIIIAH','SpoIIP','SpoIIGA','SpoIIE','SpoIIR','SpoVAC','SpoVAD',
                      'SpoVAE','SpoIVB2','SpoVS','SpoIV','SpoIVA','SpoVE','SpoVD','SpoVFB','SpoVFA','SpoVB',
                      'SpoVT','SpoVG','CD1492','CD2492','CdeC','CotA','SodA','CotJB2','CotD',
                      'Gpr','SspA','SspB','BclA3', 'SigA2','SigK','SigE','SigF','SigG')
sporulation <- subset(combined_mapping, rownames(combined_mapping) %in% sporulation_keep)
rownames(sporulation) <- c('bclA3', 'CD1492', 'CD2492', 'cdeC', 'cotA', 
                                   'cotD', 'cotJB2', 'dpaA', 'gpr', 'sigA2','sigK','sigE','sigF','sigG', 'sodA', 'spoIIAA', 'spoIIAB', 'spoIID', 
                                   'spoIIE', 'spoIIGA', 'spoIIIAA', 'spoIIIAB', 'spoIIIAC', 'spoIIIAD', 'spoIIIAE', 
                                   'spoIIIAG', 'spoIIIAH', 'spoIIID', 'spoIIP', 'spoIIR', 'spoIV', 'spoIVA',  
                                   'spoIVB2', 'spoVAC', 'spoVAD', 'spoVAE', 'spoVB', 'spoVD', 'spoVE', 'spoVFA', 
                                   'spoVFB', 'spoVG', 'spoVS', 'spoVT', 'sspA', 'sspB')

# Iteratively rarefy mappings
#sub_size <- round(min(colSums(sporulation[,1:4])) * 0.9)
#cefoperazone <- t(rrarefy(sporulation$Cefoperazone, sample=sub_size))
#clindamycin <- t(rrarefy(sporulation$Clindamycin, sample=sub_size))
#streptomycin <- t(rrarefy(sporulation$Streptomycin, sample=sub_size))
#germfree <- t(rrarefy(sporulation$Germfree, sample=sub_size))
#for (index in 1:999) {
#  cefoperazone <- cbind(cefoperazone, t(rrarefy(sporulation$Cefoperazone, sample=sub_size)))
#  clindamycin <- cbind(clindamycin, t(rrarefy(sporulation$Clindamycin, sample=sub_size)))
#  streptomycin <- cbind(streptomycin, t(rrarefy(sporulation$Streptomycin, sample=sub_size)))
#  germfree <- cbind(germfree, t(rrarefy(sporulation$Germfree, sample=sub_size)))
#}
#rm(index)
# Log transform data
#cefoperazone <- log10(cefoperazone + 1)
#clindamycin <- log10(clindamycin + 1)
#streptomycin <- log10(streptomycin + 1)
#germfree <- log10(germfree + 1)
# Medians
#sporulation$Cefoperazone <- rowMedians(cefoperazone)
#sporulation$Clindamycin <- rowMedians(clindamycin)
#sporulation$Streptomycin <- rowMedians(streptomycin)
#sporulation$Germfree <- rowMedians(germfree)
#sporulation_medians <- sporulation
#rownames(sporulation_medians) <- c('bclA3', 'CD1492', 'CD2492', 'cdeC', 'cotA', 
#                                   'cotD', 'cotJB2', 'dpaA', 'gpr', 'sodA', 'spoIIAA', 'spoIIAB', 'spoIID', 
#                                   'spoIIE', 'spoIIGA', 'spoIIIAA', 'spoIIIAB', 'spoIIIAC', 'spoIIIAD', 'spoIIIAE', 
#                                   'spoIIIAG', 'spoIIIAH', 'spoIIID', 'spoIIP', 'spoIIR', 'spoIV', 'spoIVA',  
#                                   'spoIVB2', 'spoVAC', 'spoVAD', 'spoVAE', 'spoVB', 'spoVD', 'spoVE', 'spoVFA', 
#                                   'spoVFB', 'spoVG', 'spoVS', 'spoVT', 'sspA', 'sspB')
# SDs
#sporulation$Cefoperazone <- rowSds(cefoperazone)
#sporulation$Clindamycin <- rowSds(clindamycin)
#sporulation$Streptomycin <- rowSds(streptomycin)
#sporulation$Germfree <- rowSds(germfree)
#sporulation_sds <- sporulation * 1.95
#rownames(sporulation_sds) <- c('bclA3', 'CD1492', 'CD2492', 'cdeC', 'cotA', 
#                               'cotD', 'cotJB2', 'dpaA', 'gpr', 'sodA', 'spoIIAA', 'spoIIAB', 'spoIID', 
#                               'spoIIE', 'spoIIGA', 'spoIIIAA', 'spoIIIAB', 'spoIIIAC', 'spoIIIAD', 'spoIIIAE', 
#                               'spoIIIAG', 'spoIIIAH', 'spoIIID', 'spoIIP', 'spoIIR', 'spoIV', 'spoIVA',  
#                               'spoIVB2', 'spoVAC', 'spoVAD', 'spoVAE', 'spoVB', 'spoVD', 'spoVE', 'spoVFA', 
#                               'spoVFB', 'spoVG', 'spoVS', 'spoVT', 'sspA', 'sspB')

################

# Quorum sensing
quorum_keep <- c('LuxS', 'AgrD', 'AgrB')
quorum <- subset(combined_mapping, rownames(combined_mapping) %in% quorum_keep)
rownames(quorum) <- c('agrB', 'agrD', 'luxS')

# Iteratively rarefy mappings
#sub_size <- round(min(colSums(quorum[,1:4])) * 0.9)
#cefoperazone <- t(rrarefy(quorum$Cefoperazone, sample=sub_size))
#clindamycin <- t(rrarefy(quorum$Clindamycin, sample=sub_size))
#streptomycin <- t(rrarefy(quorum$Streptomycin, sample=sub_size))
#germfree <- t(rrarefy(quorum$Germfree, sample=sub_size))
#for (index in 1:999) {
#  cefoperazone <- cbind(cefoperazone, t(rrarefy(quorum$Cefoperazone, sample=sub_size)))
#  clindamycin <- cbind(clindamycin, t(rrarefy(quorum$Clindamycin, sample=sub_size)))
#  streptomycin <- cbind(streptomycin, t(rrarefy(quorum$Streptomycin, sample=sub_size)))
#  germfree <- cbind(germfree, t(rrarefy(quorum$Germfree, sample=sub_size)))
#}
#rm(index)
# Log transform data
#cefoperazone <- log10(cefoperazone + 1)
#clindamycin <- log10(clindamycin + 1)
#streptomycin <- log10(streptomycin + 1)
#germfree <- log10(germfree + 1)
# Medians
#quorum$Cefoperazone <- rowMedians(cefoperazone)
#quorum$Clindamycin <- rowMedians(clindamycin)
#quorum$Streptomycin <- rowMedians(streptomycin)
#quorum$Germfree <- rowMedians(germfree)
#quorum_medians <- quorum
#rownames(quorum_medians) <- c('agrB', 'agrD', 'luxS')
# SDs
#quorum$Cefoperazone <- rowSds(cefoperazone)
#quorum$Clindamycin <- rowSds(clindamycin)
#quorum$Streptomycin <- rowSds(streptomycin)
#quorum$Germfree <- rowSds(germfree)
#quorum_sds <- quorum * 1.95
#rownames(quorum_sds) <- c('agrB', 'agrD', 'luxS')

################

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

# A - Sigma factors
par(las=1, mar=c(4,5,1,1), mgp=c(2.5, 1, 0))
x_coords <- barplot(t(sigma), col=select_palette, space=c(0,1.5), beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,4))
abline(h=c(1:3), lty=2)
barplot(t(sigma), col=select_palette, space=c(0,1.5), beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,4), add=TRUE)
box()
labelsY <- c(0, parse(text=paste(rep(10,4), '^', seq(1,4,1), sep='')))
axis(side=2, at=c(0:4), labelsY, tick=TRUE, las=1, cex=1.7)
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,36.7,5.5), y=par()$usr[3]-0.035*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(rownames(sigma)), srt=45, adj=1, xpd=TRUE, cex=1.2)
legend('topright', legend='Sigma factors', pt.cex=0, bty='n', cex=1.8)

#x_coords <- as.data.frame(t(x_coords))
#colnames(x_coords) <- c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree')
#segments(x0=x_coords$Cefoperazone, y0=c(sigma_medians$Cefoperazone+sigma_sds$Cefoperazone), 
#         x1=x_coords$Cefoperazone, y1=c(sigma_medians$Cefoperazone-sigma_sds$Cefoperazone), lwd=1.2)
#segments(x0=x_coords$Clindamycin, y0=c(sigma_medians$Clindamycin+sigma_sds$Clindamycin), 
#         x1=x_coords$Clindamycin, y1=c(sigma_medians$Clindamycin-sigma_sds$Clindamycin), lwd=1.2)
#segments(x0=x_coords$Streptomycin, y0=c(sigma_medians$Streptomycin+sigma_sds$Streptomycin), 
#         x1=x_coords$Streptomycin, y1=c(sigma_medians$Streptomycin-sigma_sds$Streptomycin), lwd=1.2)

mtext('A', side=2, line=2, las=2, adj=1.6, padj=-10, cex=1.5)

################

# B - PaLoc
par(las=1, mar=c(4,5,1,1), mgp=c(2.5, 1, 0))
x_coords <- barplot(t(paloc), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
                    ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,3))
abline(h=c(1:2), lty=2)
barplot(t(paloc), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,3), add=TRUE)
box()
labelsY <- c(0, parse(text=paste(rep(10,3), '^', seq(1,3,1), sep='')))
axis(side=2, at=c(0:3), labelsY, tick=TRUE, las=1, cex=1.7)
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,36.7,5.5), y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(rownames(paloc)), srt=45, adj=1, xpd=TRUE, cex=1.4)
legend('topright', legend='Pathogenicity loci', pt.cex=0, bty='n', cex=1.8)

#x_coords <- as.data.frame(t(x_coords))
#colnames(x_coords) <- c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree')
#segments(x0=x_coords$Cefoperazone, y0=c(paloc_medians$Cefoperazone+paloc_sds$Cefoperazone), 
#         x1=x_coords$Cefoperazone, y1=c(paloc_medians$Cefoperazone-paloc_sds$Cefoperazone), lwd=1.2)
#segments(x0=x_coords$Clindamycin, y0=c(paloc_medians$Clindamycin+paloc_sds$Clindamycin), 
#         x1=x_coords$Clindamycin, y1=c(paloc_medians$Clindamycin-paloc_sds$Clindamycin), lwd=1.2)
#segments(x0=x_coords$Streptomycin, y0=c(paloc_medians$Streptomycin+paloc_sds$Streptomycin), 
#         x1=x_coords$Streptomycin, y1=c(paloc_medians$Streptomycin-paloc_sds$Streptomycin), lwd=1.2)

mtext('B', side=2, line=2, las=2, adj=1.6, padj=-10, cex=1.5)

################

# C - Sporulation
par(las=1, mar=c(4,5,1,1), mgp=c(2.5, 1, 0))
x_coords <- barplot(t(sporulation), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,4))
abline(h=c(1:3), lty=2)
barplot(t(sporulation), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,4), add=TRUE)
box()
labelsY <- c(0, parse(text=paste(rep(10,4), '^', seq(1,4,1), sep='')))
axis(side=2, at=c(0:4), labelsY, tick=TRUE, las=1, cex=1.7)
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,251.2,5.5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(rownames(sporulation)), srt=45, adj=1, xpd=TRUE, cex=0.9)
legend('topright', legend='Sporulation', pt.cex=0, bty='n', cex=1.8)

#x_coords <- as.data.frame(t(x_coords))
#colnames(x_coords) <- c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree')
#segments(x0=x_coords$Cefoperazone, y0=c(sporulation_medians$Cefoperazone+sporulation_sds$Cefoperazone), 
#         x1=x_coords$Cefoperazone, y1=c(sporulation_medians$Cefoperazone-sporulation_sds$Cefoperazone), lwd=1.2)
#segments(x0=x_coords$Clindamycin, y0=c(sporulation_medians$Clindamycin+sporulation_sds$Clindamycin), 
#         x1=x_coords$Clindamycin, y1=c(sporulation_medians$Clindamycin-sporulation_sds$Clindamycin), lwd=1.2)
#segments(x0=x_coords$Streptomycin, y0=c(sporulation_medians$Streptomycin+sporulation_sds$Streptomycin), 
#         x1=x_coords$Streptomycin, y1=c(sporulation_medians$Streptomycin-sporulation_sds$Streptomycin), lwd=1.2)

mtext('C', side=2, line=2, las=2, adj=1.6, padj=-10, cex=1.5)

dev.off()

#--------------------------------------------------------------------------------------------------------------#

plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_quorum.pdf'
pdf(file=plot_file, width=7, height=7)

# Quorum sensing
par(las=1, mar=c(4,5,1,1), mgp=c(2.5, 1, 0))
x_coords <- barplot(t(quorum), col=select_palette, beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,2))
abline(h=1, lty=2)
barplot(t(quorum), col=select_palette, beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,2), add=TRUE)
box()
labelsY <- c(0, parse(text=paste(rep(10,2), '^', seq(1,2,1), sep='')))
axis(side=2, at=c(0:2), labelsY, tick=TRUE, las=1, cex=1.7)
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=c(2.7,8.2,13.7), y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(rownames(quorum)), srt=45, adj=1, xpd=TRUE, cex=1.6)
legend('topright', legend='Quorum sensing', pt.cex=0, bty='n', cex=1.8)

#x_coords <- as.data.frame(t(x_coords))
#colnames(x_coords) <- c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree')
#segments(x0=x_coords$Cefoperazone, y0=c(quorum_medians$Cefoperazone+quorum_sds$Cefoperazone), 
#         x1=x_coords$Cefoperazone, y1=c(quorum_medians$Cefoperazone-quorum_sds$Cefoperazone), lwd=1.2)
#segments(x0=x_coords$Clindamycin, y0=c(quorum_medians$Clindamycin+quorum_sds$Clindamycin), 
#         x1=x_coords$Clindamycin, y1=c(quorum_medians$Clindamycin-quorum_sds$Clindamycin), lwd=1.2)
#segments(x0=x_coords$Streptomycin, y0=c(quorum_medians$Streptomycin+quorum_sds$Streptomycin), 
#         x1=x_coords$Streptomycin, y1=c(quorum_medians$Streptomycin-quorum_sds$Streptomycin), lwd=1.2)

dev.off()

#--------------------------------------------------------------------------------------------------------------#

# Clean up
rm(quorum, sigma, sporulation, paloc, plot_file, select_palette, x_coords, labelsY, make.italic)
for (dep in deps){
  pkg <- paste('package:', dep,sep='')
   detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

