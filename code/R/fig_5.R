
deps <- c('shape', 'wesanderson', 'dplyr', 'NMF');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Select files


# Read in data


#-------------------------------------------------------------------------------------------------------------------------------------#

# Format the data


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


#-------------------------------------------------------------------------------------------------------------------------------------#

# B - Bar chart comparison of specific compounds


#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
dev.off()
rm()
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

