
# Load dependencies
deps <- c('wesanderson','vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,'Package'] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
set.seed(42)

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
sigma_keep <- c('SigK', 'SigF', 'CodY', 'CcpA', 'SigH', 'Spo0A', 'PrdR', 'Rex', 'SigE', 'SigG')
paloc_keep <- c('TcdR','TcdC','TcdE','CdtR','TcdA','TcdB')
sporulation_keep <- c('SpoIIAB','SpoIIE','SpoVS','SpoIVA','SpoVFB','SpoVB','SpoVG','CdeC','CotJB2','CotD',
                      'SspA','SspB') # variance > 0.01
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
streptomycin_sub <- t(rrarefy(combined_mapping$Streptomycin, sample=abx_sub_size))
streptomycin_sub <- streptomycin_sub / colSums(streptomycin_sub)
cefoperazone_sub <- t(rrarefy(combined_mapping$Cefoperazone, sample=abx_sub_size))
cefoperazone_sub <- cefoperazone_sub / colSums(cefoperazone_sub)
clindamycin_sub <- t(rrarefy(combined_mapping$Clindamycin, sample=abx_sub_size))
clindamycin_sub <- clindamycin_sub / colSums(clindamycin_sub)
germfree_sub <- t(rrarefy(combined_mapping$Germfree, sample=gf_sub_size))
germfree_sub <- germfree_sub / colSums(germfree_sub)
for (index in 1:999) {
  temp <- t(rrarefy(combined_mapping$Streptomycin, sample=abx_sub_size))
  temp <- temp / colSums(temp)
  streptomycin_sub <- cbind(streptomycin_sub, temp)
  temp <- t(rrarefy(combined_mapping$Cefoperazone, sample=abx_sub_size))
  temp <- temp / colSums(temp)
  cefoperazone_sub <- cbind(cefoperazone_sub, temp)
  temp <- t(rrarefy(combined_mapping$Clindamycin, sample=abx_sub_size))
  temp <- temp / colSums(temp)
  clindamycin_sub <- cbind(clindamycin_sub, temp)
  temp <- t(rrarefy(combined_mapping$Germfree, sample=gf_sub_size))
  temp <- temp / colSums(temp)
  germfree_sub <- cbind(germfree_sub, temp)
}

combined_mapping$Streptomycin <- apply(streptomycin_sub, 1, median)
combined_mapping$Cefoperazone <- apply(cefoperazone_sub, 1, median)
combined_mapping$Clindamycin <- apply(clindamycin_sub, 1, median)
combined_mapping$Germfree <- apply(germfree_sub, 1, median)

combined_mapping_sd <- cbind(apply(cefoperazone_sub, 1, sd), apply(clindamycin_sub, 1, sd), apply(streptomycin_sub, 1, sd), apply(germfree_sub, 1, sd))
colnames(combined_mapping_sd) <- c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Germfree')
rownames(combined_mapping_sd) <- rownames(combined_mapping)

rm(index, abx_sub_size, gf_sub_size, cefoperazone_sub, clindamycin_sub, streptomycin_sub, germfree_sub, temp)

# Convert to percentages
combined_mapping <- combined_mapping * 100
combined_mapping_sd <- combined_mapping_sd * 100
# Separate the mappings again

# Sigma factors
# Integration of Metabolism and Virulence by Clostridium difficile CodY
# Global transcriptional control by glucose and carbon regulator CcpA in Clostridium difficile.
# Proline-Dependent Regulation of Clostridium difficile Stickland Metabolism
# The Clostridium difficile spo0A Gene Is a Persistence and Transmission Factor
sigma <- subset(combined_mapping, rownames(combined_mapping) %in% sigma_keep)
rownames(sigma) <- c('ccpA', 'codY', 'prdR', 'rex', 'sigE', 'sigF', 'sigG', 'sigH', 'sigK', 'spo0A')
sigma <- t(sigma)
sigma <-cbind(sigma[,10], sigma[,1], sigma[,2], sigma[,3], sigma[,4], sigma[,5], sigma[,6], sigma[,7], sigma[,8], sigma[,9])
sigma[sigma == 0] <- NA
sigma_sd <- subset(combined_mapping_sd, rownames(combined_mapping_sd) %in% sigma_keep)
rownames(sigma_sd) <- c('ccpA', 'codY', 'prdR', 'rex', 'sigE', 'sigF', 'sigG', 'sigH', 'sigK', 'spo0A')
sigma_sd <- t(sigma_sd)
sigma_sd <-cbind(sigma_sd[,10], sigma_sd[,1], sigma_sd[,2], sigma_sd[,3], sigma_sd[,4], sigma_sd[,5], sigma_sd[,6], sigma_sd[,7], sigma_sd[,8], sigma_sd[,9])

# Pathogenicity
paloc <- subset(combined_mapping, rownames(combined_mapping) %in% paloc_keep)
rownames(paloc) <- c('cdtR', 'tcdA', 'tcdB', 'tcdC', 'tcdE', 'tcdR')
paloc <- t(paloc)
paloc[paloc == 0] <- NA
paloc_sd <- subset(combined_mapping_sd, rownames(combined_mapping_sd) %in% paloc_keep)
rownames(paloc_sd) <- c('cdtR', 'tcdA', 'tcdB', 'tcdC', 'tcdE', 'tcdR')
paloc_sd <- t(paloc_sd)

# Quorum sensing
quorum <- subset(combined_mapping, rownames(combined_mapping) %in% quorum_keep)
rownames(quorum) <- c('agrB', 'agrD', 'luxS')
quorum <- t(quorum)
quorum[quorum == 0] <- NA
quorum_sd <- subset(combined_mapping_sd, rownames(combined_mapping_sd) %in% quorum_keep)
rownames(quorum_sd) <- c('agrB', 'agrD', 'luxS')
quorum_sd <- t(quorum_sd)

# Sporulation
# Genes with >0.01 variance included in final analysis
sporulation <- subset(combined_mapping, rownames(combined_mapping) %in% sporulation_keep)
rownames(sporulation) <- c('cdeC','cotD','cotJB2','spoIIAB','spoIIE',
                           'spoIVA','spoVB','spoVFB','spoVG','spoVS','sspA','sspB')
sporulation <- t(sporulation)
sporulation <- cbind(sporulation[,4], sporulation[,5], sporulation[,9], sporulation[,10], rep(NA,4),
      sporulation[,6], sporulation[,7], rep(NA,4),
      sporulation[,1], sporulation[,2], sporulation[,3], sporulation[,8], sporulation[,11], sporulation[,12])
colnames(sporulation) <- c('spoIIAB', 'spoIIE', 'spoVG', 'spoVS', '', 'spoIVA', 'spoVB', '', 'cdeC', 'cotD', 'cotJB2', 'spoVFB', 'sspA', 'sspB')
sporulation[sporulation == 0] <- NA
sporulation_sd <- subset(combined_mapping_sd, rownames(combined_mapping_sd) %in% sporulation_keep)
sporulation_sd <- t(sporulation_sd)
sporulation_sd <- cbind(sporulation_sd[,4], sporulation_sd[,5], sporulation_sd[,9], sporulation_sd[,10], rep(NA,4),
                        sporulation_sd[,6], sporulation_sd[,7], rep(NA,4),
                        sporulation_sd[,1], sporulation_sd[,2], sporulation_sd[,3], sporulation_sd[,8], sporulation_sd[,11], sporulation_sd[,12])
colnames(sporulation_sd) <- c('spoIIAB', 'spoIIE', 'spoVG', 'spoVS', '', 'spoIVA', 'spoVB', '', 'cdeC', 'cotD', 'cotJB2', 'spoVFB', 'sspA', 'sspB')

# Clean up
rm(combined_mapping, combined_mapping_sd, sigma_keep, paloc_keep, sporulation_keep, quorum_keep)

#--------------------------------------------------------------------------------------------------------------#

# Set the color palette and plotting environment
select_palette <- c(wes_palette('FantasticFox')[1], wes_palette('FantasticFox')[3], wes_palette('FantasticFox')[5], 'forestgreen')
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_3.pdf'
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
pdf(file=plot_file, width=14, height=10)
layout(matrix(c(1,2,2,2,3,
                1,2,2,2,3,
                4,4,5,5,5,
                4,4,5,5,5),
              nrow=4, ncol=5, byrow = TRUE))

#--------------------------------------------------------------------------------------------------------------#

# Legend plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
legend('center', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germ free'), pt.cex=3.6, cex=2.1,
       pch=22, col='black', pt.bg=select_palette, ncol=1, bty='n')
      
# Sporulation
par(las=1, mar=c(7,5,1,1), mgp=c(3.9, 1, 0))
x_coords <- (barplot(sporulation, col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
                     ylab='Relative Transcript Abundance', ylim=c(0,30), cex.lab=1.5))
abline(h=c(10,20), lty=2)
barplot(sporulation, col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,30), add=TRUE, cex.lab=1.5)
box()
axis(side=2, at=c(0,10,20,30), c('0%','10%','20%','30%'), tick=TRUE, las=1, cex.axis=1.3)
text(x=seq(3.7,77,5.5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(c('spoIIAB', 'spoIIE', 'spoVG', 'spoVS', '',
                          'spoIVA', 'spoVB', '',
                          'cdeC', 'cotD', 'cotJB2', 'spoVFB', 'sspA', 'sspB')), 
     srt=45, adj=1, xpd=TRUE, cex=1.5)
segments(x0=x_coords, y0=sporulation+sporulation_sd, x1=x_coords, y1=sporulation-sporulation_sd, lwd=1.5)
segments(x0=x_coords-0.2, y0=sporulation+sporulation_sd, x1=x_coords+0.2, y1=sporulation+sporulation_sd, lwd=1.5)
segments(x0=x_coords-0.2, y0=sporulation-sporulation_sd, x1=x_coords+0.2, y1=sporulation-sporulation_sd, lwd=1.5)
legend('topright', legend='Sporulation', pt.cex=0, bty='n', cex=1.8)
segments(x0=c(1.5,29,46), y0=par()$usr[3]-0.155*(par()$usr[4]-par()$usr[3]), 
         x1=c(22,38,76.5), y1=par()$usr[3]-0.155*(par()$usr[4]-par()$usr[3]), lwd=2, xpd=TRUE)
text(x=c(22,58.5,70), y=par()$usr[3]-0.185*(par()$usr[4]-par()$usr[3]), 
     labels=c('Early','Intermediate','Late'), adj=3, xpd=TRUE, cex=1.5)
mtext('a', side=2, line=2, las=2, adj=3.3, padj=-12, cex=1.3, font=2)

# Quorum sensing
par(las=1, mar=c(4.5,5,1,1), mgp=c(3.9, 1, 0))
x_coords <- (barplot(quorum, col=select_palette, beside=TRUE, xaxt='n', yaxt='n', 
                     ylab='Relative Transcript Abundance', ylim=c(0,2.7), cex.lab=1.5))
abline(h=c(0.9,1.8), lty=2)
barplot(quorum, col=select_palette, beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,2.7), add=TRUE, cex.lab=1.5)
box()
axis(side=2, at=c(0,0.9,1.8,2.7), c('0%','0.9%','1.8%','2.7%'), tick=TRUE, las=1, cex.axis=1.3)
text(x=c(2.7,8.2,13.7), y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(colnames(quorum)), srt=45, adj=1, xpd=TRUE, cex=1.6)
segments(x0=x_coords, y0=quorum+quorum_sd, x1=x_coords, y1=quorum-quorum_sd, lwd=1.5)
segments(x0=x_coords-0.2, y0=quorum+quorum_sd, x1=x_coords+0.2, y1=quorum+quorum_sd, lwd=1.5)
segments(x0=x_coords-0.2, y0=quorum-quorum_sd, x1=x_coords+0.2, y1=quorum-quorum_sd, lwd=1.5)
legend('topright', legend='Quorum sensing', pt.cex=0, bty='n', cex=1.8)
mtext('b', side=2, line=2, las=2, adj=3.3, padj=-13, cex=1.3, font=2)

# Pathogenicity
par(las=1, mar=c(4.5,5.5,1,1), mgp=c(3.9, 1, 0))
x_coords <- (barplot(paloc, col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
                     ylab='Relative Transcript Abundance', ylim=c(0,2.1), cex.lab=1.5))
abline(h=c(0.7,1.4), lty=2)
barplot(paloc, col=select_palette, space=c(0,1.5),  beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,2.1), add=TRUE, cex.lab=1.5)
box()
axis(side=2, at=c(0,0.7,1.4,2.1), c('0%','0.7%','1.4%','2.1%'), tick=TRUE, las=1, cex.axis=1.3)
text(x=seq(3.7,33,5.5), y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(c('cdtR', 'tcdA', 'tcdB', 'tcdC', 'tcdE', 'tcdR')), srt=45, adj=1, xpd=TRUE, cex=1.6)
segments(x0=x_coords, y0=paloc+paloc_sd, x1=x_coords, y1=paloc-paloc_sd, lwd=1.5)
segments(x0=x_coords-0.2, y0=paloc+paloc_sd, x1=x_coords+0.2, y1=paloc+paloc_sd, lwd=1.5)
segments(x0=x_coords-0.2, y0=paloc-paloc_sd, x1=x_coords+0.2, y1=paloc-paloc_sd, lwd=1.5)
legend('topright', legend='Pathogenicity', pt.cex=0, bty='n', cex=1.8)
mtext('c', side=2, line=2, las=2, adj=3.3, padj=-13, cex=1.3, font=2)

# Sigma factors
par(las=1, mar=c(4.5,5,1,1), mgp=c(3.9, 1, 0))
x_coords <- (barplot(sigma, col=select_palette, space=c(0,1.5), beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,27), cex.lab=1.5))
abline(h=c(9, 18), lty=2)
barplot(sigma, col=select_palette, space=c(0,1.5), beside=TRUE, xaxt='n', yaxt='n', 
        ylab='Relative Transcript Abundance', ylim=c(0,27), add=TRUE, cex.lab=1.5)
box()
axis(side=2, at=c(0,9,18,27), c('0%','9%','18%','27%'), tick=TRUE, las=1, cex.axis=1.3)
text(x=seq(3.7,55,5.5), y=par()$usr[3]-0.035*(par()$usr[4]-par()$usr[3]),
     labels=make.italic(c('spo0A', 'ccpA', 'codY', 'prdR', 'rex', 'sigE', 'sigF', 'sigH', 'sigG', 'sigK')), 
     srt=45, adj=1, xpd=TRUE, cex=1.6)
segments(x0=x_coords, y0=sigma+sigma_sd, x1=x_coords, y1=sigma-sigma_sd, lwd=1.5)
segments(x0=x_coords-0.2, y0=sigma+sigma_sd, x1=x_coords+0.2, y1=sigma+sigma_sd, lwd=1.5)
segments(x0=x_coords-0.2, y0=sigma-sigma_sd, x1=x_coords+0.2, y1=sigma-sigma_sd, lwd=1.5)
legend('topright', legend='Sigma factors', pt.cex=0, bty='n', cex=1.8)
mtext('d', side=2, line=2, las=2, adj=3.3, padj=-13, cex=1.3, font=2)

dev.off()

#--------------------------------------------------------------------------------------------------------------#

# Clean up
rm(quorum, sigma, sporulation, paloc, quorum_sd, sigma_sd, sporulation_sd, paloc_sd, plot_file, select_palette, make.italic, x_coords)
for (dep in deps){
  pkg <- paste('package:', dep,sep='')
   detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

