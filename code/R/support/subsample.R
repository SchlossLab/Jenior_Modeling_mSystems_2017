
# Load dependencies
deps <- c('vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Define input file names
cefoperazone_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/all_genes/cefoperazone_630.RNA_reads2cdf630.norm.annotated.txt'
clindamycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/all_genes/clindamycin_630.RNA_reads2cdf630.norm.annotated.txt'
streptomycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/all_genes/streptomycin_630.RNA_reads2cdf630.norm.annotated.txt'
germfree_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/all_genes/germfree.RNA_reads2cdf630.norm.annotated.txt'

# Load in data
cefoperazone <- read.delim(cefoperazone_file, sep='\t', header=FALSE, row.names=1)
colnames(cefoperazone) <- c('Cefoperazone', 'ko', 'gene', 'pathway')
clindamycin <- read.delim(clindamycin_file, sep='\t', header=FALSE, row.names=1)
colnames(clindamycin) <- c('Clindamycin', 'ko', 'gene', 'pathway')
streptomycin <- read.delim(streptomycin_file, sep='\t', header=FALSE, row.names=1)
colnames(streptomycin) <- c('Streptomycin', 'ko', 'gene', 'pathway')
germfree <- read.delim(germfree_file, sep='\t', header=FALSE, row.names=1)
colnames(germfree) <- c('Germfree', 'ko', 'gene', 'pathway')
rm(cefoperazone_file, clindamycin_file, streptomycin_file, germfree_file)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data for merging
cefoperazone$ko <- NULL
cefoperazone$gene <- NULL
cefoperazone$pathway <- NULL
clindamycin$ko <- NULL
clindamycin$gene <- NULL
clindamycin$pathway <- NULL
streptomycin$ko <- NULL
streptomycin$gene <- NULL
streptomycin$pathway <- NULL

# Merge tables
combined_mapping <- merge(cefoperazone, clindamycin, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
combined_mapping <- merge(combined_mapping, streptomycin, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
combined_mapping <- merge(combined_mapping, germfree, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
rm(cefoperazone, clindamycin, streptomycin, germfree)

# Rarefy mappings to be equal within effectiveness of sequencing effort
sub_size <- round(min(colSums(combined_mapping[,1:3])) * 0.9) # 97930
gf_sub_size <- round(min(sum(combined_mapping[,4])) * 0.9) # 34050
cefoperazone_sub <- t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size))
clindamycin_sub <- t(rrarefy(combined_mapping$Clindamycin, sample=sub_size))
streptomycin_sub <- t(rrarefy(combined_mapping$Streptomycin, sample=sub_size))
germfree_sub <- t(rrarefy(combined_mapping$Germfree, sample=gf_sub_size))
for (index in 1:999) {
  cefoperazone_sub <- cbind(cefoperazone_sub, t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size)))
  clindamycin_sub <- cbind(clindamycin_sub, t(rrarefy(combined_mapping$Clindamycin, sample=sub_size)))
  streptomycin_sub <- cbind(streptomycin_sub, t(rrarefy(combined_mapping$Streptomycin, sample=sub_size)))
  germfree_sub <- cbind(germfree_sub, t(rrarefy(combined_mapping$Germfree, sample=gf_sub_size)))
}

# Calculate medians for subsampling distribution
combined_mapping$Cefoperazone <- apply(cefoperazone_sub, 1, median)
combined_mapping$Clindamycin <- apply(clindamycin_sub, 1, median)
combined_mapping$Streptomycin <- apply(streptomycin_sub, 1, median)
combined_mapping$Germfree <- apply(germfree_sub, 1, median)
rm(index, sub_size, gf_sub_size, cefoperazone_sub, clindamycin_sub, streptomycin_sub, germfree_sub)

# Write to an output file
output_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/all_genes/all_infections.ko.RNA.subsample.txt'
write.table(combined_mapping, file=output_file, quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

# Clean up
rm(output_file, combined_mapping)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()
