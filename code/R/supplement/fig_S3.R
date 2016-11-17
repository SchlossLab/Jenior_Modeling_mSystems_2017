
# Load dependencies
deps <- c('wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Select files
acetate_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/cef_acetate_630.txt'

# Read in data
scfa <- read.delim(acetate_file, sep='\t', header=T)
rm(acetate_file)

# Format the data
scfa$group <- factor(scfa$group, levels=c('mock', 'infected'))
scfa$acetate <- as.numeric(as.character(scfa$acetate))

# Subset for stats
mock <- as.numeric(scfa[scfa$group == 'mock', 2])
infected <- as.numeric(scfa[scfa$group == 'infected', 2])

# Calculate significance
wilcox.test(mock, infected, exact=F)
# p-value = 0.01219 *

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S3.pdf'
pdf(file=plot_file, width=5, height=5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Acetate
par(las=1, mar=c(2,4,0.7,1), mgp=c(2.3,0.7,0), xpd=FALSE, yaxs='i')
stripchart(acetate~group, data=scfa, vertical=T, pch=c(1,19), lwd=3,
           ylim=c(0,10), xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1],
           cex=1.8, ylab='nmol Acetate per mg Cecal Content', method='jitter', jitter=0.25)
abline(h=c(2,4,6,8), lty=2)
axis(side=1, at=c(1,2), c('Mock Infected', '630 Infected'), tick=FALSE)
axis(side=2, at=seq(0,10,2), labels=c('0.0','2.0','4.0','6.0','8.0','10.0'))
segments(x0=c(0.7,1.7), y0=c(median(mock),median(infected)), x1=c(1.3,2.3), y1=, lwd=3)
text(2, median(infected) + 1.5, labels='*', font=2, cex=2)
#legend('topright', legend=expression(paste(italic(p),'-value = 0.01219')), bty='n', pt.cex=0)

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
dev.off()
rm(scfa, mock, infected, plot_file)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()
