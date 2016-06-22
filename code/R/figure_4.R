
deps <- c('vegan', 'VennDiagram');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Import and format transcript mapping data
cefoperazone_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/cefoperazone_630.RNA_reads2cdf630.norm.annotated.txt'
clindamycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/clindamycin_630.RNA_reads2cdf630.norm.annotated.txt'
streptomycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/streptomycin_630.RNA_reads2cdf630.norm.annotated.txt'
germfree_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/germfree.RNA_reads2cdf630.norm.annotated.txt'

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

# Rarefy mappings
sub_size <- round(min(colSums(combined_mapping[,1:4])) * 0.9) # 30645
combined_mapping$Cefoperazone <- t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size))
combined_mapping$Clindamycin <- t(rrarefy(combined_mapping$Clindamycin, sample=sub_size))
combined_mapping$Streptomycin <- t(rrarefy(combined_mapping$Streptomycin, sample=sub_size))
combined_mapping$Germfree <- t(rrarefy(combined_mapping$Germfree, sample=sub_size))
rm(sub_size)

# Remove zeroes before transformation
combined_mapping[combined_mapping == 0] <- 1

# Transform mapping values
combined_mapping$Cefoperazone <- log10(combined_mapping$Cefoperazone)
combined_mapping$Clindamycin <- log10(combined_mapping$Clindamycin)
combined_mapping$Streptomycin <- log10(combined_mapping$Streptomycin)
combined_mapping$Germfree <- log10(combined_mapping$Germfree)

# Seperate into each comparison table
cef_gf_comparison <- combined_mapping[,c(1,4,6,7)]
clinda_gf_comparison <- combined_mapping[,c(2,4,6,7)]
strep_gf_comparison <- combined_mapping[,c(3,4,6,7)]
rm(combined_mapping)

# Eliminate rows with no transcription
cef_gf_comparison <- subset(cef_gf_comparison, (cef_gf_comparison$Cefoperazone+cef_gf_comparison$Germfree) != 0)
clinda_gf_comparison <- subset(clinda_gf_comparison, (clinda_gf_comparison$Clindamycin+clinda_gf_comparison$Germfree) != 0)
strep_gf_comparison <- subset(strep_gf_comparison, (strep_gf_comparison$Streptomycin+strep_gf_comparison$Germfree) != 0)

#-------------------------------------------------------------------------------------------------------------------------#

# Plot
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_4.pdf'

# Generate figure
pdf(file=plot_file, width=12, height=10)




# generate a figure like cody's paper



