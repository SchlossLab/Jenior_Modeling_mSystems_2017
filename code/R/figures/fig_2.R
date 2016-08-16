
deps <- c('wesanderson','vegan', 'matrixStats', 'plotrix');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
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

#cefoperazone_file <- '/media/mjenior/Jenior\ HD/data/mapping/cdifficile630/select_genes/cefoperazone_630.RNA_reads2select.all.norm.txt'
#clindamycin_file <- '/media/mjenior/Jenior\ HD/data/mapping/cdifficile630/select_genes/clindamycin_630.RNA_reads2select.all.norm.txt'
#streptomycin_file <- '/media/mjenior/Jenior\ HD/data/mapping/cdifficile630/select_genes/streptomycin_630.RNA_reads2select.all.norm.txt'
#germfree_file <- '/media/mjenior/Jenior\ HD/data/mapping/cdifficile630/select_genes/germfree.RNA_reads2select.all.norm.txt'

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
combined_mapping$gene <- gsub("Clostridium_difficile_630\\|","", combined_mapping$gene)
combined_mapping$gene <- gsub('ENA\\|CDT20869\\|CDT20869.1\\|Clostridium_difficile_putative_phage_replication_protein_','', combined_mapping$gene)
combined_mapping$gene <- gsub('_',' ', combined_mapping$gene)
rownames(combined_mapping) <- combined_mapping$gene
combined_mapping$gene <- NULL
rm(cefoperazone, clindamycin, streptomycin, germfree)

#--------------------------------------------------------------------------------------------------------------#

# Break up the data and calculate stats

# Sigma factors
sigma_keep <- c('CodY','CcpA','CdtR','SigH','SigB','SigA1','SigA2','SigE','SigF','SigG','SigK','SigV','FliA',
                'TetR_family|1','TetR_family|2','TetR_family|3','TetR_family|4','TetR_family|5','TetR_family|6',
                'TetR_family|7','TetR_family|8','TetR_family|9','TetR_family|10','TetR_family|11','TetR_family|12',
                'TetR_family|13','TetR_family|14','TetR_family|15','TetR_family|16','Rex','PrdR','Spo0A','DpaA')
sigma <- subset(combined_mapping, rownames(combined_mapping) %in% sigma_keep)
# Iteratively rarefy mappings
sub_size <- round(min(colSums(sigma[,1:3])) * 0.9)
cefoperazone <- t(rrarefy(sigma$Cefoperazone, sample=sub_size))
clindamycin <- t(rrarefy(sigma$Clindamycin, sample=sub_size))
streptomycin <- t(rrarefy(sigma$Streptomycin, sample=sub_size))
for (index in 1:999) {
  cefoperazone <- cbind(cefoperazone, t(rrarefy(sigma$Cefoperazone, sample=sub_size)))
  clindamycin <- cbind(clindamycin, t(rrarefy(sigma$Clindamycin, sample=sub_size)))
  streptomycin <- cbind(streptomycin, t(rrarefy(sigma$Streptomycin, sample=sub_size)))
}
# Log transform data
cefoperazone[cefoperazone == 0] <- 1
cefoperazone <- log10(cefoperazone)
clindamycin[clindamycin == 0] <- 1
clindamycin <- log10(clindamycin)
streptomycin[streptomycin == 0] <- 1
streptomycin <- log10(streptomycin)
# Medians
sigma$Cefoperazone <- rowMedians(cefoperazone)
sigma$Clindamycin <- rowMedians(clindamycin)
sigma$Streptomycin <- rowMedians(streptomycin)
sigma_medians <- sigma
# SDs
sigma$Cefoperazone <- rowSds(cefoperazone)
sigma$Clindamycin <- rowSds(clindamycin)
sigma$Streptomycin <- rowSds(streptomycin)
sigma_sds <- sigma * 1.95
# Clean up
rm(sub_size, cefoperazone, clindamycin, streptomycin, sigma, index)

# PaLoc
paloc_keep <- c('TcdR','TcdC','TcdE','CdtR','TcdA','TcdB')
paloc <- subset(combined_mapping, rownames(combined_mapping) %in% paloc_keep)
# Iteratively rarefy mappings
sub_size <- round(min(colSums(paloc[,1:3])) * 0.9)
cefoperazone <- t(rrarefy(paloc$Cefoperazone, sample=sub_size))
clindamycin <- t(rrarefy(paloc$Clindamycin, sample=sub_size))
streptomycin <- t(rrarefy(paloc$Streptomycin, sample=sub_size))
for (index in 1:999) {
  cefoperazone <- cbind(cefoperazone, t(rrarefy(paloc$Cefoperazone, sample=sub_size)))
  clindamycin <- cbind(clindamycin, t(rrarefy(paloc$Clindamycin, sample=sub_size)))
  streptomycin <- cbind(streptomycin, t(rrarefy(paloc$Streptomycin, sample=sub_size)))
}
# Log transform data
cefoperazone[cefoperazone == 0] <- 1
cefoperazone <- log10(cefoperazone)
clindamycin[clindamycin == 0] <- 1
clindamycin <- log10(clindamycin)
streptomycin[streptomycin == 0] <- 1
streptomycin <- log10(streptomycin)
# Medians
paloc$Cefoperazone <- rowMedians(cefoperazone)
paloc$Clindamycin <- rowMedians(clindamycin)
paloc$Streptomycin <- rowMedians(streptomycin)
paloc_medians <- paloc
# SDs
paloc$Cefoperazone <- rowSds(cefoperazone)
paloc$Clindamycin <- rowSds(clindamycin)
paloc$Streptomycin <- rowSds(streptomycin)
paloc_sds <- paloc * 1.95
# Clean up
rm(sub_size, cefoperazone, clindamycin, streptomycin, paloc, index)

# Sporulation
sporulation_keep <- c('SpoIID','SpoIIID','SpoIIAA','SpoIIAB','SpoIIIAA','SpoIIIAB','SpoIIIAC','SpoIIIAD',
                      'SpoIIIAE','SpoIIIAG','SpoIIIAH','SpoIIP','SpoIIGA','SpoIIE','SpoIIR','SpoVAC','SpoVAD',
                      'SpoVAE','SpoIVB2','SpoIVB','SpoVS','SpoIV','SpoIVA','SpoVE','SpoVD','SpoVFB','SpoVFA','SpoVB',
                      'SpoVT','SpoVG','CD1579','CD1492','CD2492','CotF','CotCB','CdeC','CotA','SodA','CotJB2','CotD',
                      'Gpr','SspA','BclA2','SspB','BclA3')
sporulation <- subset(combined_mapping, rownames(combined_mapping) %in% sporulation_keep)
# Iteratively rarefy mappings
sub_size <- round(min(colSums(sporulation[,1:3])) * 0.9)
cefoperazone <- t(rrarefy(sporulation$Cefoperazone, sample=sub_size))
clindamycin <- t(rrarefy(sporulation$Clindamycin, sample=sub_size))
streptomycin <- t(rrarefy(sporulation$Streptomycin, sample=sub_size))
for (index in 1:999) {
  cefoperazone <- cbind(cefoperazone, t(rrarefy(sporulation$Cefoperazone, sample=sub_size)))
  clindamycin <- cbind(clindamycin, t(rrarefy(sporulation$Clindamycin, sample=sub_size)))
  streptomycin <- cbind(streptomycin, t(rrarefy(sporulation$Streptomycin, sample=sub_size)))
}
# Log transform data
cefoperazone[cefoperazone == 0] <- 1
cefoperazone <- log10(cefoperazone)
clindamycin[clindamycin == 0] <- 1
clindamycin <- log10(clindamycin)
streptomycin[streptomycin == 0] <- 1
streptomycin <- log10(streptomycin)
# Medians
sporulation$Cefoperazone <- rowMedians(cefoperazone)
sporulation$Clindamycin <- rowMedians(clindamycin)
sporulation$Streptomycin <- rowMedians(streptomycin)
sporulation_medians <- sporulation
# SDs
sporulation$Cefoperazone <- rowSds(cefoperazone)
sporulation$Clindamycin <- rowSds(clindamycin)
sporulation$Streptomycin <- rowSds(streptomycin)
sporulation_sds <- sporulation * 1.95
# Clean up
rm(sub_size, cefoperazone, clindamycin, streptomycin, sporulation, index)

# Quorum sensing
quorum_keep <- c('LuxS','AgrD','AgrB')
quorum <- subset(combined_mapping, rownames(combined_mapping) %in% quorum_keep)
# Iteratively rarefy mappings
sub_size <- round(min(colSums(quorum[,1:3])) * 0.9)
cefoperazone <- t(rrarefy(quorum$Cefoperazone, sample=sub_size))
clindamycin <- t(rrarefy(quorum$Clindamycin, sample=sub_size))
streptomycin <- t(rrarefy(quorum$Streptomycin, sample=sub_size))
for (index in 1:999) {
  cefoperazone <- cbind(cefoperazone, t(rrarefy(quorum$Cefoperazone, sample=sub_size)))
  clindamycin <- cbind(clindamycin, t(rrarefy(quorum$Clindamycin, sample=sub_size)))
  streptomycin <- cbind(streptomycin, t(rrarefy(quorum$Streptomycin, sample=sub_size)))
}
# Log transform data
cefoperazone[cefoperazone == 0] <- 1
cefoperazone <- log10(cefoperazone)
clindamycin[clindamycin == 0] <- 1
clindamycin <- log10(clindamycin)
streptomycin[streptomycin == 0] <- 1
streptomycin <- log10(streptomycin)
# Medians
quorum$Cefoperazone <- rowMedians(cefoperazone)
quorum$Clindamycin <- rowMedians(clindamycin)
quorum$Streptomycin <- rowMedians(streptomycin)
quorum_medians <- quorum
# SDs
quorum$Cefoperazone <- rowSds(cefoperazone)
quorum$Clindamycin <- rowSds(clindamycin)
quorum$Streptomycin <- rowSds(streptomycin)
quorum_sds <- quorum * 1.95
# Clean up
rm(sub_size, cefoperazone, clindamycin, streptomycin, quorum, index)

# Clean up
rm(combined_mapping, sigma_keep, paloc_keep, sporulation_keep, quorum_keep)

#--------------------------------------------------------------------------------------------------------------#

# Set the color palette and plotting environment
select_palette <- c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5])
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_2.pdf'
pdf(file=plot_file, width=12, height=14)
layout(matrix(c(1,1,
                2,2,
                3,4),
              nrow=3, ncol=2, byrow = TRUE))

#--------------------------------------------------------------------------------------------------------------#

# A - Sigma factors
par(las=1, mar=c(4,5,1,1), mgp=c(2.5, 1, 0))
x_coords <- barplot(t(sigma_medians), col=select_palette, space=c(0,1.5), beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,3))
abline(h=c(1:2), lty=2)
barplot(t(sigma_medians), col=select_palette, space=c(0,1.5), beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,3), add=TRUE)
box()
labelsY <- c(0, parse(text=paste(rep(10,3), '^', seq(1,3,1), sep='')))
axis(side=2, at=c(0:3), labelsY, tick=TRUE, las=1, cex=1.7)
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,79.2,4.5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(sigma_medians), srt=45, adj=1, xpd=TRUE, cex=1.2)
legend('topright', legend='Sigma factors', pt.cex=0, bty='n', cex=1.8)

x_coords <- as.data.frame(t(x_coords))
colnames(x_coords) <- c('Streptomycin', 'Cefoperazone', 'Clindamycin')
segments(x0=x_coords$Cefoperazone, y0=c(sigma_medians$Cefoperazone+sigma_sds$Cefoperazone), x1=x_coords$Cefoperazone, y1=c(sigma_medians$Cefoperazone-sigma_sds$Cefoperazone), lwd=1.2)
segments(x0=x_coords$Clindamycin, y0=c(sigma_medians$Clindamycin+sigma_sds$Clindamycin), x1=x_coords$Clindamycin, y1=c(sigma_medians$Clindamycin-sigma_sds$Clindamycin), lwd=1.2)
segments(x0=x_coords$Streptomycin, y0=c(sigma_medians$Streptomycin+sigma_sds$Streptomycin), x1=x_coords$Streptomycin, y1=c(sigma_medians$Streptomycin-sigma_sds$Streptomycin), lwd=1.2)

mtext('A', side=2, line=2, las=2, adj=1.6, padj=-10, cex=1.5)

#--------------------------------------------------------------------------------------------------------------#

# B - Sporulation
par(las=1, mar=c(4,5,1,1), mgp=c(2.5, 1, 0))
x_coords <- barplot(t(sporulation_medians), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,3))
abline(h=c(1:2), lty=2)
barplot(t(sporulation_medians), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,3), add=TRUE)
box()
labelsY <- c(0, parse(text=paste(rep(10,3), '^', seq(1,3,1), sep='')))
axis(side=2, at=c(0:3), labelsY, tick=TRUE, las=1, cex=1.7)
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,205.2,4.5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(sporulation_medians), srt=45, adj=1, xpd=TRUE, cex=0.8)
legend('topright', legend='Sporulation', pt.cex=0, bty='n', cex=1.8)

x_coords <- as.data.frame(t(x_coords))
colnames(x_coords) <- c('Streptomycin', 'Cefoperazone', 'Clindamycin')
segments(x0=x_coords$Cefoperazone, y0=c(sporulation_medians$Cefoperazone+sporulation_sds$Cefoperazone), x1=x_coords$Cefoperazone, y1=c(sporulation_medians$Cefoperazone-sporulation_sds$Cefoperazone), lwd=1.2)
segments(x0=x_coords$Clindamycin, y0=c(sporulation_medians$Clindamycin+sporulation_sds$Clindamycin), x1=x_coords$Clindamycin, y1=c(sporulation_medians$Clindamycin-sporulation_sds$Clindamycin), lwd=1.2)
segments(x0=x_coords$Streptomycin, y0=c(sporulation_medians$Streptomycin+sporulation_sds$Streptomycin), x1=x_coords$Streptomycin, y1=c(sporulation_medians$Streptomycin-sporulation_sds$Streptomycin), lwd=1.2)

mtext('B', side=2, line=2, las=2, adj=1.6, padj=-10, cex=1.5)

#--------------------------------------------------------------------------------------------------------------#

# C - PaLoc
par(las=1, mar=c(4,5,1,1), mgp=c(2.5, 1, 0))
x_coords <- barplot(t(paloc_medians), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,2))
abline(h=1, lty=2)
barplot(t(paloc_medians), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,2), add=TRUE)
box()
labelsY <- c(0, parse(text=paste(rep(10,2), '^', seq(1,2,1), sep='')))
axis(side=2, at=c(0:2), labelsY, tick=TRUE, las=1, cex=1.7)
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,29.7,4.5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(paloc_medians), srt=45, adj=1, xpd=TRUE, cex=1.4)
legend('topright', legend='Pathogenicity loci', pt.cex=0, bty='n', cex=1.8)

x_coords <- as.data.frame(t(x_coords))
colnames(x_coords) <- c('Streptomycin', 'Cefoperazone', 'Clindamycin')
segments(x0=x_coords$Cefoperazone, y0=c(paloc_medians$Cefoperazone+paloc_sds$Cefoperazone), x1=x_coords$Cefoperazone, y1=c(paloc_medians$Cefoperazone-paloc_sds$Cefoperazone), lwd=1.2)
segments(x0=x_coords$Clindamycin, y0=c(paloc_medians$Clindamycin+paloc_sds$Clindamycin), x1=x_coords$Clindamycin, y1=c(paloc_medians$Clindamycin-paloc_sds$Clindamycin), lwd=1.2)
segments(x0=x_coords$Streptomycin, y0=c(paloc_medians$Streptomycin+paloc_sds$Streptomycin), x1=x_coords$Streptomycin, y1=c(paloc_medians$Streptomycin-paloc_sds$Streptomycin), lwd=1.2)

mtext('C', side=2, line=2, las=2, adj=1.6, padj=-10, cex=1.5)

#--------------------------------------------------------------------------------------------------------------#

# D - Quorum sensing
par(las=1, mar=c(4,5,1,1), mgp=c(2.5, 1, 0))
x_coords <- barplot(t(quorum_medians), col=select_palette, beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,2))
abline(h=1, lty=2)
barplot(t(quorum_medians), col=select_palette, beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,2), add=TRUE)
box()
labelsY <- c(0, parse(text=paste(rep(10,2), '^', seq(1,2,1), sep='')))
axis(side=2, at=c(0:2), labelsY, tick=TRUE, las=1, cex=1.7)
legend('topleft', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=c(2.7,6.7,10.7), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(quorum_medians), srt=45, adj=1, xpd=TRUE, cex=1.6)
legend('topright', legend='Quorum sensing', pt.cex=0, bty='n', cex=1.8)

x_coords <- as.data.frame(t(x_coords))
colnames(x_coords) <- c('Streptomycin', 'Cefoperazone', 'Clindamycin')
segments(x0=x_coords$Cefoperazone, y0=c(quorum_medians$Cefoperazone+quorum_sds$Cefoperazone), x1=x_coords$Cefoperazone, y1=c(quorum_medians$Cefoperazone-quorum_sds$Cefoperazone), lwd=1.2)
segments(x0=x_coords$Clindamycin, y0=c(quorum_medians$Clindamycin+quorum_sds$Clindamycin), x1=x_coords$Clindamycin, y1=c(quorum_medians$Clindamycin-quorum_sds$Clindamycin), lwd=1.2)
segments(x0=x_coords$Streptomycin, y0=c(quorum_medians$Streptomycin+quorum_sds$Streptomycin), x1=x_coords$Streptomycin, y1=c(quorum_medians$Streptomycin-quorum_sds$Streptomycin), lwd=1.2)

mtext('D', side=2, line=2, las=2, adj=1.6, padj=-10, cex=1.5)

#--------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
rm(quorum_medians, quorum_sds, sigma_medians, sigma_sds, sporulation_medians, sporulation_sds, paloc_medians, paloc_sds, plot_file, select_palette, x_coords, labelsY)
for (dep in deps){
  pkg <- paste('package:', dep,sep='')
   detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

