
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
#shared_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.family.subsample.shared'
#taxonomy_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.family.cons.taxonomy'
metadata_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'

# Load in data
nmds <- read.delim(nmds_file, sep='\t', header=T, row.names=1)
nmds <- nmds[!rownames(nmds) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
summary <- read.delim(summary_file, sep='\t', header=T, row.names=2)
summary <- summary[!rownames(summary) %in% c('CefC5M2'), ] # Remove contaminated sample
#taxonomy <- read.delim(taxonomy_file, sep='\t', header=T, row.names=1)
#shared <- read.delim(shared_file, sep='\t', header=T, row.names=2)
#shared <- shared[!rownames(shared) %in% c('CefC5M2'), ]  # Remove contaminated sample
#shared$numOtus <- NULL
#shared$label <- NULL
rm(nmds_file, summary_file, metadata_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data

# Combine the metadata with axes and format
metadata_axes <- merge(metadata, nmds, by='row.names')
mock_axes <- subset(metadata_axes, infection == 'mock')
mock_axes <- subset(mock_axes, type != 'germfree')
infected_axes <- subset(metadata_axes, infection == '630')
infected_axes <- subset(infected_axes, type != 'germfree')

# Combine the metadata with summary, format, and get median invsimpson diversity
metadata_summary <- merge(metadata, summary, by='row.names')
metadata_summary <- subset(metadata_summary, type != 'germfree')
metadata_summary$abx <- factor(metadata_summary$abx, levels=c('streptomycin', 'cefoperazone', 'clindamycin', 'none'))
strep_div <- as.numeric(median(metadata_summary[metadata_summary$abx == 'streptomycin', 13]))
cef_div <- as.numeric(median(metadata_summary[metadata_summary$abx == 'cefoperazone', 13]))
clinda_div <- as.numeric(median(metadata_summary[metadata_summary$abx == 'clindamycin', 13]))
conv_div <- as.numeric(median(metadata_summary[metadata_summary$abx == 'none', 13]))

rm(nmds, summary, metadata)
#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S1.pdf'
pdf(file=plot_file, width=12, height=8)
layout(matrix(c(1,1,2),
              nrow=1, ncol=3, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# NMDS of treatment groups
par(las=1, mar=c(5,5,1,1))
plot(metadata_axes$axis1, metadata_axes$axis2, pch=21, cex=0,
     xlim=c(-0.8,0.8), ylim=c(-0.8,0.8), cex.lab=2, cex.axis=1.7,
     xlab='NMDS Axis 1', ylab='NMDS Axis 2')

# add mock points
points(x=mock_axes$axis1, y=mock_axes$axis2, 
       col=c(wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black', wes_palette("FantasticFox")[1])[mock_axes$abx], 
       pch=1, lwd=3, cex=3)
# add mock infected points
points(x=infected_axes$axis1, y=infected_axes$axis2, 
       col=c(wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black', wes_palette("FantasticFox")[1])[infected_axes$abx], 
       pch=1, lwd=3, cex=3)

# Add legends
legend('topleft', legend=c('Streptomycin-treated', 'Cefoperzone-treated', 'Clindamycin-treated', 'No Antibiotics'), 
       col=c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'black'), 
       pch=15, cex=1.9, pt.cex=2.8, bty='n')
#legend('bottomleft', legend=c('Mock Infected', '630 Infected'), 
#       col='black', pt.bg=c('white','black'), pch=21, cex=2, pt.cex=3, bty='n')

mtext('a', side=2, line=2, las=2, adj=1.7, padj=-18.1, cex=1.6, font=2)

#-----------------------#

# Stripchart of Inverse simpson diversity
par(las=1, mar=c(1,4,1,1), mgp=c(2.5,0.7,0), yaxs='i')
stripchart(invsimpson~abx, data=metadata_summary, vertical=T, pch=2, lwd=2.5,
           ylim=c(0,20), xaxt='n', cex=2, 
           col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5],'black'),
           ylab='Inv. Simpson Diversity', method='jitter', jitter=0.15, cex.axis=1.7, cex.lab=2)

# add medians
segments(x0=0.7, y0=strep_div, x1=1.3, y1=strep_div, lwd=2.5)
segments(x0=1.7, y0=cef_div, x1=2.3, y1=cef_div, lwd=2.5)
segments(x0=2.7, y0=clinda_div, x1=3.3, y1=clinda_div, lwd=2.5)
segments(x0=3.7, y0=conv_div, x1=4.3, y1=conv_div, lwd=2.5)
text(x=4, y=19, '***', font=2, cex=3) # add significance

# add legend
legend('topleft', legend=c('Streptomycin-treated', 'Cefoperzone-treated', 'Clindamycin-treated', 'No Antibiotics'), 
       col=c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'black'), 
       pch=15, cex=1.9, pt.cex=2.8, bty='n')

mtext('b', side=2, line=2, las=2, adj=1.7, padj=-19.5, cex=1.6, font=2)

#-----------------------#

# Possible family-level phylotype bar chart...

#mtext('c', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.1, font=2)

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
dev.off()


rm(plot_file, metadata_axes, mock_axes, infected_axes, metadata_summary, 
   strep_div, cef_div, clinda_div, conv_div)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()