
# Load dependencies
deps <- c('grid', 'gridExtra');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in and format antibiotic table
abx_table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/abx.tsv'
abx_table <- read.delim(abx_table_file, sep='\t', header=TRUE)

abx_table[,2] <- gsub('_', ' ', abx_table[,2])
abx_table[,3] <- gsub('_', ' ', abx_table[,3])
abx_table[,4] <- gsub('_', ' ', abx_table[,4])
abx_table[,5] <- gsub('_', ' ', abx_table[,5])
abx_table[,6] <- gsub('_', ' ', abx_table[,6])

abx_table$Target <- sapply(lapply(abx_table$Target, strwrap, width=40), paste, collapse="\n")
abx_table$Activity <- sapply(lapply(abx_table$Activity, strwrap, width=40), paste, collapse="\n")

#-------------------------------------------------------------------------------------------------------------------------------------#

plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/tables/table_1.pdf'
pdf(file=plot_file, width=14, height=3)

grid.table(abx_table, rows=NULL)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
rm(abx_table_file, abx_table, plot_file)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()
