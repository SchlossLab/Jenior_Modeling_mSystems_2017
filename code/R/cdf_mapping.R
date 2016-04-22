
deps <- c('vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

# Define input file names
cefoperazone_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/cdifficile630/genes/cefoperazone_630.RNA_reads2cdf630.pool.norm.annotated.txt'
clindamycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/cdifficile630/genes/clindamycin_630.RNA_reads2cdf630.pool.norm.annotated.txt'
streptomycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/cdifficile630/genes/streptomycin_630.RNA_reads2cdf630.pool.norm.annotated.txt'
germfree_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/cdifficile630/genes/germfree.RNA_reads2cdf630.pool.norm.annotated.txt'

# Load in data
cefoperazone <- read.delim(cefoperazone_file, sep='\t', header=TRUE, row.names=2)
clindamycin <- read.delim(clindamycin_file, sep='\t', header=TRUE, row.names=2)
streptomycin <- read.delim(streptomycin_file, sep='\t', header=TRUE, row.names=2)
germfree <- read.delim(germfree_file, sep='\t', header=TRUE, row.names=2)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data for merging
cefoperazone$gene_annotation <- NULL
cefoperazone$pathway_annotation <- NULL
clindamycin$gene_annotation <- NULL
clindamycin$pathway_annotation <- NULL
streptomycin$gene_annotation <- NULL
streptomycin$pathway_annotation <- NULL

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
colnames(combined_mapping) <- c('cefoperazone', 'clindamycin', 'streptomycin', 'germfree', 'gene_annotation', 'pathway_annotation')

# Remove file names and single tables from memory
rm(cefoperazone_file, clindamycin_file, streptomycin_file, germfree_file)
rm(cefoperazone, clindamycin, streptomycin, germfree)

# Rarefy mappings to be equal within sequencing type
sub_size <- round(min(colSums(combined_mapping[,c(1:4)]))*0.9)
combined_mapping$cefoperazone <- t(rrarefy(combined_mapping$cefoperazone, sample=sub_size))
combined_mapping$clindamycin <- t(rrarefy(combined_mapping$clindamycin, sample=sub_size))
combined_mapping$streptomycin <- t(rrarefy(combined_mapping$streptomycin, sample=sub_size))
combined_mapping$germfree <- t(rrarefy(combined_mapping$germfree, sample=sub_size))

# Eliminate genes with no transcripts mapping
combined_mapping <- combined_mapping[rowSums(combined_mapping[,c(1:4)]) != 0, ] 

# Log10 transform the data
combined_mapping[combined_mapping == 0] <- 1
combined_mapping$cefoperazone <- log10(combined_mapping$cefoperazone)
combined_mapping$clindamycin <- log10(combined_mapping$clindamycin)
combined_mapping$streptomycin <- log10(combined_mapping$streptomycin)
combined_mapping$germfree <- log10(combined_mapping$germfree)
norm_max <- max(combined_mapping[,c(1:4)])

# Remove genes with no pathway annotation
unannotated <- subset(combined_mapping, combined_mapping$pathway_annotation == 'none')
#combined_mapping <- subset(combined_mapping, pathway_annotation != 'none')

# Sort by pathway annotation
combined_mapping <- combined_mapping[order(combined_mapping$pathway_annotation),] 

#-------------------------------------------------------------------------------------------------------------------------#

# Plotting 



# Subset data to color specific groups
# Highest resolution - Known C. difficile 630 carbon sources
glycolysis <- subset(combined_mapping, grepl('*Glycolysis_/_Gluconeogenesis*', combined_mapping$pathway_annotation))
amino_sugar <- subset(combined_mapping, grepl('*Amino_sugar*', combined_mapping$pathway_annotation))
galactose <- subset(combined_mapping, grepl('*Galactose*', combined_mapping$pathway_annotation))
fructose_mannose <- subset(combined_mapping, grepl('*Fructose_and_mannose*', combined_mapping$pathway_annotation))
starch_sucrose <- subset(combined_mapping, grepl('*Starch_and_sucrose*', combined_mapping$pathway_annotation))
glycan_degradation <- subset(combined_mapping, grepl('*glycan_degradation*', combined_mapping$pathway_annotation))
cdf_carbohydates <- rbind(glycolysis, amino_sugar, galactose, fructose_mannose, starch_sucrose, glycan_degradation)

cysteine_methionine <- subset(combined_mapping, grepl('*Cysteine_and_methionine*', combined_mapping$pathway_annotation))
valine_leucine_isoleucine <- subset(combined_mapping, grepl('*Valine,_leucine_and_isoleucine*', combined_mapping$pathway_annotation))
glycine_serine_threonine <- subset(combined_mapping, grepl('*Glycine,_serine_and_threonine*', combined_mapping$pathway_annotation))
alanine_aspartate_glutamate <- subset(combined_mapping, grepl('*Alanine,_aspartate_and_glutamate*', combined_mapping$pathway_annotation))
cdf_amino_acids <- rbind(cysteine_methionine, valine_leucine_isoleucine, glycine_serine_threonine, alanine_aspartate_glutamate)

cdf_carbon_sources <- rbind(cdf_carbohydates, cdf_amino_acids)




butanoate <- subset(combined_mapping, grepl('*Butanoate*', combined_mapping$pathway_annotation))
pentose_glucuronate <- subset(combined_mapping, grepl('*Pentose_and_glucuronate*', combined_mapping$pathway_annotation))
phenylalanine_tyrosine_tryptophan <- subset(combined_mapping, grepl('*Phenylalanine,_tyrosine_and_tryptophan*', combined_mapping$pathway_annotation))
tca_cycle <- subset(combined_mapping, grepl('*TCA_cycle*', combined_mapping$pathway_annotation))
lipopolysaccharide <- subset(combined_mapping, grepl('*Lipopolysaccharide*', combined_mapping$pathway_annotation))
oxidative_phosphorylation <- subset(combined_mapping, grepl('*Oxidative_phosphorylation*', combined_mapping$pathway_annotation))
glyoxylate_dicarboxylate <- subset(combined_mapping, grepl('*Glyoxylate_and_dicarboxylate*', combined_mapping$pathway_annotation))
antibiotics <- subset(combined_mapping, grepl('*antibiotic*', combined_mapping$pathway_annotation))
polyketide_sugar_biosynthesis <- subset(combined_mapping, grepl('*Polyketide_sugar_unit_biosynthesis*', combined_mapping$pathway_annotation))
lysine <- subset(combined_mapping, grepl('*Lysine*', combined_mapping$pathway_annotation))
arginine_proline <- subset(combined_mapping, grepl('*Arginine_and_proline*', combined_mapping$pathway_annotation))
pyruvate <- subset(combined_mapping, grepl('*Pyruvate*', combined_mapping$pathway_annotation))
peptidoglycan <- subset(combined_mapping, grepl('*Peptidoglycan*', combined_mapping$pathway_annotation))
secondary_bile <- subset(combined_mapping, grepl('*Secondary_bile*', combined_mapping$pathway_annotation))
# More general scale
energy_metabolism <- subset(combined_mapping, grepl('*Energy_metabolism*', combined_mapping$pathway_annotation))
carbohydrate_metabolism <- subset(combined_mapping, grepl('*Carbohydrate_metabolism*', combined_mapping$pathway_annotation))
amino_acid_metabolism <- subset(combined_mapping, grepl('*Amino_acid_metabolism*', combined_mapping$pathway_annotation))
nucleotide_metabolism <- subset(combined_mapping, grepl('*Nucleotide_metabolism*', combined_mapping$pathway_annotation))
lipid_metabolism <- subset(combined_mapping, grepl('*Lipid_metabolism*', combined_mapping$pathway_annotation))
cofactors_and_vitamins <- subset(combined_mapping, grepl('*cofactors_and_vitamins*', combined_mapping$pathway_annotation))
diverse_environments <- subset(combined_mapping, grepl('*diverse_environments*', combined_mapping$pathway_annotation))
sulfur_metabolism <- subset(combined_mapping, grepl('*Sulfur_metabolism*', combined_mapping$pathway_annotation))
secondary_metabolites <- subset(combined_mapping, grepl('*secondary_metabolites*', combined_mapping$pathway_annotation))



#-------------------------------------------------------------------------------------------------------------------------#

# Define which pathway to plot and the ouput file name
pathway <- cdf_carbohydates[,c(1,2)]
point_color <- 'firebrick1'
#pathway_name <- 'Carbohydrate Metabolism'
#plot_file <- '~/Desktop/figures/cdf_amino_acids.streptomycin.pdf'

# Plot it!
#pdf(file=plot_file, width=7, height=6)
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0))
plot(x=combined_mapping$cefoperazone, y=combined_mapping$clindamycin, 
     xlim=c(0,4.1), ylim=c(0,4.1), pch=20, col='gray25', xaxt='n', yaxt='n', cex.lab=1.3, 
     xlab='Normalized cDNA Read Abundance', ylab='Normalized cDNA Read Abundance')
#legend('topleft', legend=pathway_name, bty='n', cex=1.5) 
segments(-2, -2, 8, 8, lwd=2, lty=2)
axis_labels <- parse(text=paste(rep(10,4), '^', seq(0,4,1), sep=''))
axis(side=1, at=c(0:4), axis_labels, tick=TRUE, cex.axis=1.2)
axis(side=2, at=c(0:4), axis_labels, tick=TRUE, las=1, cex.axis=1.2)
points(pathway, cex=1.5, pch=21, bg=point_color, col='black')
#dev.off()




