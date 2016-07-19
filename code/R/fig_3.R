
deps <- c('wesanderson','vegan', 'matrixStats');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

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

#--------------------------------------------------------------------------------------------------------------#

# Merge tables
combined_mapping <- merge(cefoperazone, clindamycin, by='gene')
combined_mapping <- merge(combined_mapping, streptomycin, by='gene')
combined_mapping$gene <- gsub("Clostridium_difficile_630\\|","", combined_mapping$gene)
combined_mapping$gene <- gsub('ENA\\|CDT20869\\|CDT20869.1\\|Clostridium_difficile_putative_phage_replication_protein_','', combined_mapping$gene)
combined_mapping$gene <- gsub('_',' ', combined_mapping$gene)
rownames(combined_mapping) <- combined_mapping$gene
combined_mapping$gene <- NULL
rm(cefoperazone, clindamycin, streptomycin, germfree)

# Iteratively rarefy mappings
sub_size <- round(min(colSums(combined_mapping[,1:3])) * 0.9) # 2668
cefoperazone <- t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size))
clindamycin <- t(rrarefy(combined_mapping$Clindamycin, sample=sub_size))
streptomycin <- t(rrarefy(combined_mapping$Streptomycin, sample=sub_size))
for (index in 1:1000) {
  cefoperazone <- cbind(cefoperazone, t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size)))
  clindamycin <- cbind(clindamycin, t(rrarefy(combined_mapping$Clindamycin, sample=sub_size)))
  streptomycin <- cbind(streptomycin, t(rrarefy(combined_mapping$Streptomycin, sample=sub_size)))
}
# Medians
combined_mapping$Cefoperazone <- rowMedians(cefoperazone)
combined_mapping$Clindamycin <- rowMedians(clindamycin)
combined_mapping$Streptomycin <- rowMedians(streptomycin)
normalized_mapping <- combined_mapping
# SDs
combined_mapping$Cefoperazone <- rowSds(cefoperazone)
combined_mapping$Clindamycin <- rowSds(clindamycin)
combined_mapping$Streptomycin <- rowSds(streptomycin)
sd_mapping <- combined_mapping
# Clean up
rm(sub_size, cefoperazone, clindamycin, streptomycin, combined_mapping, index)

# Log10 transform the data
normalized_mapping[normalized_mapping == 0] <- 1
transformed_mapping <- log10(normalized_mapping)
rm(normalized_mapping)

#--------------------------------------------------------------------------------------------------------------#

# Subset for sigma factors
sigma_keep <- c('CodY','CcpA','CdtR','SigH','SigB','SigA1','SigA2','SigE','SigF','SigG','SigK','SigV','FliA',
           'TetR_family|1','TetR_family|2','TetR_family|3','TetR_family|4','TetR_family|5','TetR_family|6',
           'TetR_family|7','TetR_family|8','TetR_family|9','TetR_family|10','TetR_family|11','TetR_family|12',
           'TetR_family|13','TetR_family|14','TetR_family|15','TetR_family|16','Rex','PrdR','Spo0A','RstA','DpaA')
sigma_medians <- subset(transformed_mapping, rownames(transformed_mapping) %in% sigma_keep)
sigma_sds <- subset(transformed_mapping, rownames(transformed_mapping) %in% sigma_keep)

# Subset for PaLoc
paloc_keep <- c('TcdR','TcdC','TcdE','CdtR','TcdA','TcdB')
paloc_medians <- subset(transformed_mapping, rownames(transformed_mapping) %in% paloc_keep)
paloc_sds <- subset(sd_mapping, rownames(sd_mapping) %in% paloc_keep) * 3

# Subset for sporulation
sporulation_keep <- c('SpoIID','SpoIIID','SpoIIAA','SpoIIAB','SpoIIIAA','SpoIIIAB','SpoIIIAC','SpoIIIAD',
                 'SpoIIIAE','SpoIIIAG','SpoIIIAH','SpoIIP','SpoIIGA','SpoIIE','SpoIIR','SpoVAC','SpoVAD',
                 'SpoVAE','SpoIVB2','SpoIVB','SpoVS','SpoIV','SpoIVA','SpoVE','SpoVD','SpoVFB','SpoVFA','SpoVB',
                 'SpoVT','SpoVG','CD1579','CD1492','CD2492','CotF','CotCB','CdeC','CotA','SodA','CotJB2','CotD',
                 'Gpr','SspA','BclA2','SspB','BclA3')
sporulation_medians <- subset(transformed_mapping, rownames(transformed_mapping) %in% sporulation_keep)
sporulation_sds <- subset(sd_mapping, rownames(sd_mapping) %in% sporulation_keep)

# Quorum sensing
quorum_keep <- c('LuxS','AgrD','AgrB')
quorum_medians <- subset(transformed_mapping, rownames(transformed_mapping) %in% quorum_keep)
quorum_sds <- subset(sd_mapping, rownames(sd_mapping) %in% quorum_keep)

# Clean up
rm(transformed_mapping, sigma_keep, paloc_keep, sporulation_keep, quorum_keep)

#--------------------------------------------------------------------------------------------------------------#

# Set the color palette and plotting environment
select_palette <- c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5])
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_3.pdf'
pdf(file=plot_file, width=12, height=14)
layout(matrix(c(1,1,
                2,2,
                3,4),
              nrow=3, ncol=2, byrow = TRUE))

#--------------------------------------------------------------------------------------------------------------#

# A - Sigma factors
par(las=1, mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
x_coords <- barplot(t(sigma_medians), col=select_palette, space=c(0,1.5), beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,3))
box()
axis(side=2, at=c(1:4), parse(text=paste(rep(10,4), '^', seq(1,4,1), sep='')), tick=TRUE, las=1, cex=1.4)
abline(h=c(1:4), lty=2)
legend('topleft', legend=c('Cefoperazone', 'Clindamycin', 'Streptomycin'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,83.7,4.5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(sigma_medians), srt=45, adj=1, xpd=TRUE, cex=1.2)
legend('topright', legend='Sigma factors', pt.cex=0, bty='n', cex=1.8)

# 95% confidence interval
x_coords <- as.data.frame(t(x_coords))
colnames(x_coords) <- c('Cefoperazone', 'Clindamycin', 'Streptomycin')
segments(x0=x_coords$Cefoperazone, y0=c(sigma_medians$Cefoperazone+sigma_sds$Cefoperazone), x1=x_coords$Cefoperazone, y1=c(sigma_medians$Cefoperazone-sigma_sds$Cefoperazone), lwd=1.2)
segments(x0=x_coords$Clindamycin, y0=c(sigma_medians$Clindamycin+sigma_sds$Clindamycin), x1=x_coords$Clindamycin, y1=c(sigma_medians$Clindamycin-sigma_sds$Clindamycin), lwd=1.2)
segments(x0=x_coords$Streptomycin, y0=c(sigma_medians$Streptomycin+sigma_sds$Streptomycin), x1=x_coords$Streptomycin, y1=c(sigma_medians$Streptomycin-sigma_sds$Streptomycin), lwd=1.2)

#--------------------------------------------------------------------------------------------------------------#

# B - Sporulation
par(las=1, mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
barplot(t(sporulation_medians), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,3))
box()
axis(side=2, at=c(1:4), parse(text=paste(rep(10,4), '^', seq(1,4,1), sep='')), tick=TRUE, las=1, cex=1.4)
abline(h=c(1:4), lty=2)
legend('topleft', legend=c('Cefoperazone', 'Clindamycin', 'Streptomycin'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,205.2,4.5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(sporulation_medians), srt=45, adj=1, xpd=TRUE, cex=0.8)
legend('topright', legend='Sporulation', pt.cex=0, bty='n', cex=1.8)

#--------------------------------------------------------------------------------------------------------------#

# C - PaLoc
par(las=1, mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
barplot(t(paloc_medians), col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,2))
box()
axis(side=2, at=c(1:4), parse(text=paste(rep(10,4), '^', seq(1,4,1), sep='')), tick=TRUE, las=1, cex=1.4)
abline(h=c(1:4), lty=2)
legend('topleft', legend=c('Cefoperazone', 'Clindamycin', 'Streptomycin'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=seq(3.7,29.7,4.5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(paloc_medians), srt=45, adj=1, xpd=TRUE, cex=1.4)
legend('topright', legend='Pathogenicity Locus', pt.cex=0, bty='n', cex=1.8)

#--------------------------------------------------------------------------------------------------------------#

# D - Quorum sensing
par(las=1, mar=c(4,4,1,1), mgp=c(2.5, 1, 0))
barplot(t(quorum_medians), col=select_palette, beside=TRUE, xaxt='n', yaxt='n', 
        ylab=expression(paste('Transcript Abundance (',Log[10],')')), ylim=c(0,2))
box()
axis(side=2, at=c(1:4), parse(text=paste(rep(10,4), '^', seq(1,4,1), sep='')), tick=TRUE, las=1, cex=1.4)
abline(h=c(1:4), lty=2)
legend('topleft', legend=c('Cefoperazone', 'Clindamycin', 'Streptomycin'), pt.cex=2.3, bty='n', cex=1.2,
       pch=22, col='black', pt.bg=select_palette, ncol=1)
text(x=c(2.7,6.7,10.7), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(quorum_medians), srt=45, adj=1, xpd=TRUE, cex=1.6)
legend('topright', legend='Quorum sensing', pt.cex=0, bty='n', cex=1.8)

#--------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
rm(paloc, quorum, sigma, sporulation, plot_file, select_palette)


