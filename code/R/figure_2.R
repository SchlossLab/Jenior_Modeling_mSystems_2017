
deps <- c('wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

#--------------------------------------------------------------------------------------------------------------#

# Define variables
cefoperazone_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/cefoperazone_630.RNA_reads2cdf630.norm.annotated.txt'
clindamycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/clindamycin_630.RNA_reads2cdf630.norm.annotated.txt'
streptomycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/streptomycin_630.RNA_reads2cdf630.norm.annotated.txt'
germfree_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/germfree.RNA_reads2cdf630.norm.annotated.txt'

# Open files
cefoperazone <- read.delim(cefoperazone_file, sep='\t', header=FALSE, row.names=1)
colnames(cefoperazone) <- c('Cefoperazone', 'ko', 'gene', 'pathway')
clindamycin <- read.delim(clindamycin_file, sep='\t', header=FALSE, row.names=1)
colnames(clindamycin) <- c('Clindamycin', 'ko', 'gene', 'pathway')
streptomycin <- read.delim(streptomycin_file, sep='\t', header=FALSE, row.names=1)
colnames(streptomycin) <- c('Streptomycin', 'ko', 'gene', 'pathway')
germfree <- read.delim(germfree_file, sep='\t', header=FALSE, row.names=1)
colnames(germfree) <- c('Germfree', 'ko', 'gene', 'pathway')

# Clean up
rm(cefoperazone_file, clindamycin_file, streptomycin_file, germfree_file)

#--------------------------------------------------------------------------------------------------------------#

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

# Rarefy mappings to be equal within sequencing type
sub_size <- round(min(colSums(combined_mapping[,1:4])) * 0.9) # 34050
combined_mapping$Cefoperazone <- t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size))
combined_mapping$Clindamycin <- t(rrarefy(combined_mapping$Clindamycin, sample=sub_size))
combined_mapping$Streptomycin <- t(rrarefy(combined_mapping$Streptomycin, sample=sub_size))
combined_mapping$Germfree <- t(rrarefy(combined_mapping$Germfree, sample=sub_size))

# Eliminate genes with no transcripts mapping
combined_mapping <- combined_mapping[rowSums(combined_mapping[,1:4]) != 0, ] 

# Clean up
rm(cefoperazone, clindamycin, streptomycin, germfree, sub_size)

# Remove unused catagories
combined_mapping$ko <- NULL
combined_mapping$gene <- NULL

#--------------------------------------------------------------------------------------------------------------#

# Subset by KEGG catagory and pool

# Metabolism - not working
carbohydrate <- subset(combined_mapping, grepl('Carbohydrate_metabolism', combined_mapping$pathway))
carbohydrate <- colSums(carbohydrate[,1:4])
energy <- subset(combined_mapping, grepl('Energy_metabolism', combined_mapping$pathway))
energy <- colSums(energy[,1:4])
lipid <- subset(combined_mapping, grepl('Lipid_metabolism', combined_mapping$pathway))
lipid <- colSums(lipid[,1:4])
nucleotide <- subset(combined_mapping, grepl('Nucleotide_metabolism', combined_mapping$pathway))
nucleotide <- colSums(nucleotide[,1:4])
amino_acid <- subset(combined_mapping, grepl('Amino_acid_metabolism', combined_mapping$pathway))
amino_acid <- colSums(amino_acid[,1:4])
other_amino_acids <- subset(combined_mapping, grepl('Metabolism_of_other_amino_acids', combined_mapping$pathway))
other_amino_acids <- colSums(other_amino_acids[,1:4])
glycan <- subset(combined_mapping, grepl('Glycan_biosynthesis_and_metabolism', combined_mapping$pathway))
glycan <- colSums(glycan[,1:4])
cofactors_and_vitamins <- subset(combined_mapping, grepl('Metabolism_of_cofactors_and_vitamins', combined_mapping$pathway))
cofactors_and_vitamins <- colSums(cofactors_and_vitamins[,1:4])
terpenoids_and_polyketides <- subset(combined_mapping, grepl('Metabolism_of_terpenoids_and_polyketides', combined_mapping$pathway))
terpenoids_and_polyketides <- colSums(terpenoids_and_polyketides[,1:4])
secondary_metabolites <- subset(combined_mapping, grepl('Biosynthesis_of_other_secondary_metabolites', combined_mapping$pathway))
secondary_metabolites <- colSums(secondary_metabolites[,1:4])
xenobiotics <- subset(combined_mapping, grepl('Xenobiotics_biodegradation_and_metabolism', combined_mapping$pathway))
xenobiotics <- colSums(xenobiotics[,1:4])
genetics <- subset(combined_mapping, grepl('Genetic_Information_Processing', combined_mapping$pathway))
genetics <- colSums(genetics[,1:4])
motility <- subset(combined_mapping, grepl('Cell_motility', combined_mapping$pathway))
motility <- colSums(motility[,1:4])
rm(combined_mapping)

# Assemble final table
pooled_mapping <- rbind(carbohydrate, energy, lipid, nucleotide, amino_acid, other_amino_acids, glycan,
                        cofactors_and_vitamins, terpenoids_and_polyketides, secondary_metabolites,
                        xenobiotics, genetics, motility)

#--------------------------------------------------------------------------------------------------------------#

# Pooling pathways
cefoperazone_pooled_pathways <- aggregate(transcripts ~ pathway_name, FUN = sum, data = cefoperazone_pathways)  
rownames(cefoperazone_pooled_pathways) <- cefoperazone_pooled_pathways[,1]
cefoperazone_pooled_pathways[,1] <- NULL
colnames(cefoperazone_pooled_pathways) <- 'Cefoperazone'






#--------------------------------------------------------------------------------------------------------------#

fox <- wes_palette("FantasticFox")
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_2.pdf'

# Open a PDF
pdf(file=plot_file, width=12, height=7)

# Create layout for multi-plot
layout(mat=matrix(c(1,1, 
                    2,2), 
                  nrow=2, ncol=2, byrow=TRUE))

# Lowly abundant
par(mar=c(3,2,1,1))
barchart()




# Highly abundant
par(mar=c(3,2,1,1))
barchart()









