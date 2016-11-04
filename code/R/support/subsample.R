
# Load dependencies
deps <- c('vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Define input files
cefoperazone_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/cefoperazone_630.RNA_reads2cdf630.norm.ko.txt'
clindamycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/clindamycin_630.RNA_reads2cdf630.norm.ko.txt'
streptomycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/streptomycin_630.RNA_reads2cdf630.norm.ko.txt'
germfree_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/germfree.RNA_reads2cdf630.norm.ko.txt'

# Load in data
cefoperazone <- read.delim(cefoperazone_file, sep='\t', header=FALSE)
clindamycin <- read.delim(clindamycin_file, sep='\t', header=FALSE)
streptomycin <- read.delim(streptomycin_file, sep='\t', header=FALSE)
germfree <- read.delim(germfree_file, sep='\t', header=FALSE)
rm(cefoperazone_file, clindamycin_file, streptomycin_file, germfree_file)

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate subsampling level
sub_size <- round(min(c(sum(cefoperazone[,2]), sum(clindamycin[,2]), sum(streptomycin[,2]), sum(germfree[,2]))) * 0.9)
# = 27664

# Rarefy mappings to be equal within effectiveness of sequencing effort
cefoperazone_sub <- t(rrarefy(cefoperazone[,2], sample=sub_size))
clindamycin_sub <- t(rrarefy(clindamycin[,2], sample=sub_size))
streptomycin_sub <- t(rrarefy(streptomycin[,2], sample=sub_size))
germfree_sub <- t(rrarefy(germfree[,2], sample=sub_size))
for (index in 1:999) {
  cefoperazone_sub <- cbind(cefoperazone_sub, t(rrarefy(cefoperazone[,2], sample=sub_size)))
  clindamycin_sub <- cbind(clindamycin_sub, t(rrarefy(clindamycin[,2], sample=sub_size)))
  streptomycin_sub <- cbind(streptomycin_sub, t(rrarefy(streptomycin[,2], sample=sub_size)))
  germfree_sub <- cbind(germfree_sub, t(rrarefy(germfree[,2], sample=sub_size)))
}

# Calculate medians for subsampling distribution
cefoperazone[,2] <- apply(cefoperazone_sub, 1, median)
clindamycin[,2] <- apply(clindamycin_sub, 1, median)
streptomycin[,2] <- apply(streptomycin_sub, 1, median)
germfree[,2] <- apply(germfree_sub, 1, median)
rm(index, cefoperazone_sub, clindamycin_sub, streptomycin_sub, germfree_sub)

# Final coverage for KOs
sum(cefoperazone[,2]) / length(cefoperazone[,2]) # ~18x

# Write to output files
cefoperazone_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/cefoperazone_630.RNA_reads2cdf630.norm.ko.pick.txt'
clindamycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/clindamycin_630.RNA_reads2cdf630.norm.ko.pick.txt'
streptomycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/streptomycin_630.RNA_reads2cdf630.norm.ko.pick.txt'
germfree_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/germfree.RNA_reads2cdf630.norm.ko.pick.txt'

write.table(cefoperazone, file=cefoperazone_file, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
write.table(clindamycin, file=clindamycin_file, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
write.table(streptomycin, file=streptomycin_file, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
write.table(germfree, file=germfree_file, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

# Clean up
rm(cefoperazone, clindamycin, streptomycin, germfree)
rm(cefoperazone_file, clindamycin_file, streptomycin_file, germfree_file)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()
