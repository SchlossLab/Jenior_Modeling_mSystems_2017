
# Load dependencies
deps <- c('wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Define variables
nmds_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.unique_list.thetayc.0.03.lt.ave.nmds.axes'
summary_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.unique_list.groups.summary'
metadata_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'
shared_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
taxonomy_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.cons.taxonomy'

# Load in data
nmds <- read.delim(nmds_file, sep='\t', header=T, row.names=1)
nmds <- nmds[!rownames(nmds) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
summary <- read.delim(summary_file, sep='\t', header=T, row.names=2)
summary <- summary[!rownames(summary) %in% c('CefC5M2'), ] # Remove contaminated sample
taxonomy <- read.delim(taxonomy_file, sep='\t', header=T, row.names=1)
shared <- read.delim(shared_file, sep='\t', header=T, row.names=2)
shared <- shared[!rownames(shared) %in% c('CefC5M2'), ]  # Remove contaminated sample
shared$numOtus <- NULL
shared$label <- NULL
rm(nmds_file, summary_file, metadata_file, shared_file, taxonomy_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data

# Combine the metadata with axes and format
nmds <- merge(metadata, nmds, by='row.names')
mock_axes <- subset(nmds, infection == 'mock')
mock_axes <- subset(mock_axes, type != 'germfree')
infected_axes <- subset(nmds, infection == '630')
infected_axes <- subset(infected_axes, type != 'germfree')

# Combine the metadata with summary, format, and get median invsimpson diversity
summary <- merge(metadata, summary, by='row.names')
summary <- subset(summary, type != 'germfree')
summary$abx <- factor(summary$abx, levels=c('streptomycin', 'cefoperazone', 'clindamycin', 'none'))
strep_div <- as.numeric(median(summary[summary$abx == 'streptomycin', 13]))
cef_div <- as.numeric(median(summary[summary$abx == 'cefoperazone', 13]))
clinda_div <- as.numeric(median(summary[summary$abx == 'clindamycin', 13]))
conv_div <- as.numeric(median(summary[summary$abx == 'none', 13]))

# Calculate C. diffile reads vs the rest of the community in each condition
shared <- merge(metadata, shared, by='row.names')
rownames(shared) <- shared$Row.names
shared$Row.names <- NULL
shared$cage <- NULL
shared$mouse <- NULL
shared$gender <- NULL
shared$type <- NULL
shared <- subset(shared, infection == '630')
shared$infection <- NULL
strep_shared <- subset(shared, abx == 'streptomycin')
strep_shared$abx <- NULL
cdiff <- rowSums(cbind(strep_shared$Otu0004, strep_shared$Otu0309))
strep_shared$Otu0004 <- NULL
strep_shared$Otu0309 <- NULL
other <- rowSums(strep_shared)
strep_shared <- as.data.frame(cbind(cdiff, other))
colnames(strep_shared) <- c('cdifficile', 'other')
strep_shared$abx <- rep('strep', nrow(strep_shared))
cef_shared <- subset(shared, abx == 'cefoperazone')
cef_shared$abx <- NULL
cdiff <- rowSums(cbind(cef_shared$Otu0004, cef_shared$Otu0309))
cef_shared$Otu0004 <- NULL
cef_shared$Otu0309 <- NULL
other <- rowSums(cef_shared)
cef_shared <- as.data.frame(cbind(cdiff, other))
colnames(cef_shared) <- c('cdifficile', 'other')
cef_shared$abx <- rep('cef', nrow(cef_shared))
cef <- cef[!rownames(cef) %in% c('CefC1M1', 'CefC1M2', 'CefC1M3'), ]
clinda_shared <- subset(shared, abx == 'clindamycin')
clinda_shared$abx <- NULL
cdiff <- rowSums(cbind(clinda_shared$Otu0004, clinda_shared$Otu0309))
clinda_shared$Otu0004 <- NULL
clinda_shared$Otu0309 <- NULL
other <- rowSums(clinda_shared)
clinda_shared <- as.data.frame(cbind(cdiff, other))
colnames(clinda_shared) <- c('cdifficile', 'other')
clinda_shared$abx <- rep('clinda', nrow(clinda_shared))
gf_shared <- subset(shared, abx == 'germfree')
gf_shared$abx <- NULL
cdiff <- rowSums(cbind(gf_shared$Otu0004, gf_shared$Otu0309))
gf_shared$Otu0004 <- NULL
gf_shared$Otu0309 <- NULL
other <- rowSums(gf_shared)
gf_shared <- as.data.frame(cbind(cdiff, other))
colnames(gf_shared) <- c('cdifficile', 'other')
gf_shared$abx <- rep('gf', nrow(gf_shared))
shared <- rbind(strep_shared,
                cef_shared,
                clinda_shared,
                gf_shared)
shared$cdifficile <- as.numeric(as.character(shared$cdifficile))
shared$cdifficile <- log10(shared$cdifficile + 1)
shared$other <- as.numeric(as.character(shared$other))
shared$other <- log10(shared$other + 1)

# Shared file summary stats
strep_diff <- paste('~', as.character(round((as.numeric(mean(strep_shared[,1])) / mean(quantile(strep_shared[,2]))) * 100, digits=2)), '%', sep='')
cef_diff <- paste('~', as.character(round((as.numeric(mean(cef_shared[,1])) / mean(quantile(cef_shared[,2]))) * 100, digits=2)), '%', sep='')
clinda_diff <- paste('~', as.character(round((as.numeric(mean(clinda_shared[,1])) / mean(quantile(clinda_shared[,2]))) * 100, digits=2)), '%', sep='')
gf_diff <- paste('~', as.character(round((as.numeric(mean(gf_shared[,1])) / mean(quantile(gf_shared[,2]))) * 100, digits=2)), '%', sep='')

rm(metadata, cdiff, other, strep_shared, cef_shared, clinda_shared, gf_shared)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S1.pdf'
pdf(file=plot_file, width=17, height=7)
layout(matrix(c(1,1,2,3,3), nrow=1, ncol=5, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# NMDS of treatment groups
par(las=1, mar=c(5,5,1,1))
plot(metadata_axes$axis1, metadata_axes$axis2, pch=21, cex=0,
     xlim=c(-0.8,0.8), ylim=c(-0.8,0.8), cex.lab=2, cex.axis=1.7,
     xlab='NMDS Axis 1', ylab='NMDS Axis 2')

# add mock points
points(x=mock_axes$axis1, y=mock_axes$axis2, 
       col=c(wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black', wes_palette("FantasticFox")[1])[mock_axes$abx], 
       pch=1, lwd=3, cex=3.5)
# add mock infected points
points(x=infected_axes$axis1, y=infected_axes$axis2, 
       col=c(wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black', wes_palette("FantasticFox")[1])[infected_axes$abx], 
       pch=1, lwd=3, cex=3.5)

# Add legends
legend('topleft', legend=c('Streptomycin-treated', 'Cefoperzone-treated', 'Clindamycin-treated', 'No Antibiotics'), 
       col=c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'black'), 
       pch=15, cex=1.9, pt.cex=2.8, bty='n')
#legend('bottomleft', legend=c('Mock Infected', '630 Infected'), 
#       col='black', pt.bg=c('white','black'), pch=21, cex=2, pt.cex=3, bty='n')

mtext('A', side=2, line=2, las=2, adj=1.7, padj=-18.1, cex=1.6)

#-----------------------#

# Stripchart of Inverse simpson diversity
par(las=1, mar=c(1,4,1,1), mgp=c(2.5,0.7,0))
stripchart(invsimpson~abx, data=metadata_summary, vertical=T, pch=2, lwd=2.5,
           ylim=c(0,20), xaxt='n', cex=2, 
           col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5],'black'),
           ylab='Inv. Simpson Diversity', method='jitter', jitter=0.15, cex.axis=1.7, cex.lab=2)

# add medians
segments(x0=0.7, y0=strep_div, x1=1.3, y1=strep_div, lwd=2.5)
segments(x0=1.7, y0=cef_div, x1=2.3, y1=cef_div, lwd=2.5)
segments(x0=2.7, y0=clinda_div, x1=3.3, y1=clinda_div, lwd=2.5)
segments(x0=3.7, y0=conv_div, x1=4.3, y1=conv_div, lwd=2.5)
text(x=4, y=19, '*', font=2, cex=3) # add significance

# add legend
legend('topleft', legend=c('Streptomycin-treated', 'Cefoperzone-treated', 'Clindamycin-treated', 'No Antibiotics'), 
       col=c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'black'), 
       pch=15, cex=1.9, pt.cex=2.8, bty='n')

mtext('B', side=2, line=2, las=2, adj=1.7, padj=-19.5, cex=1.6)

#-----------------------#

# C. difficile read abundance relative to community
par(las=1, mar=c(3,4,1,1), mgp=c(2.5,1,0))
boxplot(cdifficile~abx, data=shared, cex=0, lwd=2.5, xlim=c(0.5,11.5), ylim=c(0,4), at=c(1,4,7,10),
           col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5],'forestgreen'),
           ylab='Read Abundance', staplewex=0.8, boxwex=1, lty=1, medlwd=2.5, xaxt='n', yaxt='n')
boxplot(other~abx, data=shared, cex=0, lwd=2.5, xlim=c(0.5,11.5), ylim=c(0,4), at=c(2,5,8,11),
           col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5],'forestgreen'),
           ylab='Read Abundance', staplewex=0.8, boxwex=1, lty=1, medlwd=2.5, add=TRUE, xaxt='n', yaxt='n')
text(x=c(1,2,4,5,7,8,10,11), y=4, labels=c('630','Other','630','Other','630','Other','630','Other'), cex=0.8)
axis(side=1, at=c(1.5,4.5,7.5,10.5), labels=c('Streptomycin','Cefoperazone','Clindamycin','Gnotobiotic'), cex=0.9, tick=FALSE)
labelsY <- c(0, parse(text=paste(rep(10,4), '^', seq(1,4,1), sep='')))
axis(side=2, at=c(0:4), labelsY, tick=TRUE)
abline(v=c(3,6,9), lty=2)
text(x=c(2,5,8,11), y=0, labels=c(strep_diff,cef_diff,clinda_diff,gf_diff), cex=0.8)
mtext('C', side=2, line=2, las=2, adj=1.7, padj=-19.5, cex=1.6)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps) {
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only=TRUE)
}
rm(list=ls())
gc()
