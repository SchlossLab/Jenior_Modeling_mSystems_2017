
deps <- c('vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Define input file names
cefoperazone_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/cdifficile630/genes/cefoperazone_630.RNA_reads2cdf630.pool.norm.txt'
clindamycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/cdifficile630/genes/clindamycin_630.RNA_reads2cdf630.pool.norm.txt'
streptomycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/cdifficile630/genes/streptomycin_630.RNA_reads2cdf630.pool.norm.txt'
germfree_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/cdifficile630/genes/germfree.RNA_reads2cdf630.pool.norm.txt'

# Load in data
cefoperazone <- read.delim(metagenome_file, sep='\t', header=TRUE, row.names=4)
clindamycin <- read.delim(metatranscriptome1_file, sep='\t', header=TRUE, row.names=4)
streptomycin <- read.delim(metatranscriptome2_file, sep='\t', header=TRUE, row.names=4)
germfree <- read.delim(metatranscriptome2_file, sep='\t', header=TRUE, row.names=4)

#-------------------------------------------------------------------------------------------------------------------------#

