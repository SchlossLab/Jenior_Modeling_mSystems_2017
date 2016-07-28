
deps <- c('shape', 'wesanderson', 'dplyr', 'NMF', 'gplots');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Select files
metabolites_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/'

# Read in data
metabolites <- read.delim(metabolites_file, sep='\t', header=TRUE)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format the data


  
  
#-------------------------------------------------------------------------------------------------------------------------------------#
  
# Subset groups of interest





#-------------------------------------------------------------------------------------------------------------------------------------#

# Calculate significance





#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_5.pdf'
pdf(file=plot_file, width=12, height=14)
layout(matrix(c(1,1,
                1,1,
                2,2),
              nrow=3, ncol=2, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# A - Heatmap of metabolites

heatmap.2()


mtext('A', side=2, line=2, las=2, adj=1.6, padj=-10, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# B - Bar chart comparison of specific compounds





mtext('B', side=2, line=2, las=2, adj=1.6, padj=-10, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
dev.off()
rm(metabolites)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

