
deps <- c('vegan', 'rgl', 'klaR');
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
#germfree_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/cdifficile630/genes/germfree.RNA_reads2cdf630.pool.norm.annotated.txt'

# Load in data
cefoperazone <- read.delim(cefoperazone_file, sep='\t', header=TRUE, row.names=2)
clindamycin <- read.delim(clindamycin_file, sep='\t', header=TRUE, row.names=2)
streptomycin <- read.delim(streptomycin_file, sep='\t', header=TRUE, row.names=2)
#germfree <- read.delim(germfree_file, sep='\t', header=TRUE, row.names=2)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data for merging
cefoperazone$gene_annotation <- NULL
cefoperazone$pathway_annotation <- NULL
clindamycin$gene_annotation <- NULL
clindamycin$pathway_annotation <- NULL
#streptomycin$gene_annotation <- NULL
#streptomycin$pathway_annotation <- NULL

# Merge tables
combined_mapping <- merge(cefoperazone, clindamycin, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
combined_mapping <- merge(combined_mapping, streptomycin, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
colnames(combined_mapping) <- c('cefoperazone', 'clindamycin', 'streptomycin', 'gene_annotation', 'pathway_annotation')
#combined_mapping <- merge(combined_mapping, germfree, by='row.names')
#rownames(combined_mapping) <- combined_mapping$Row.names
#combined_mapping$Row.names <- NULL
#colnames(combined_mapping) <- c('cefoperazone', 'clindamycin', 'streptomycin', 'germfree', 'gene_annotation', 'pathway_annotation')

# Remove file names and single tables from memory
rm(cefoperazone_file, clindamycin_file, streptomycin_file)
rm(cefoperazone, clindamycin, streptomycin)

# Rarefy mappings to be equal within sequencing type
sub_size <- round(min(colSums(combined_mapping[,c(1:3)]))*0.9)
combined_mapping$cefoperazone <- t(rrarefy(combined_mapping$cefoperazone, sample=sub_size))
combined_mapping$clindamycin <- t(rrarefy(combined_mapping$clindamycin, sample=sub_size))
combined_mapping$streptomycin <- t(rrarefy(combined_mapping$streptomycin, sample=sub_size))
#combined_mapping$germfree <- t(rrarefy(combined_mapping$germfree, sample=sub_size))

# Eliminate genes with no transcripts mapping
combined_mapping <- combined_mapping[rowSums(combined_mapping[,c(1:3)]) != 0, ] 

# Calculate the factors to scale points size by
averages <- log10(rowSums(combined_mapping[,c(1:3)]) / 3) * 2.5

# Convert each gene into the fraction of the transcription for that gene across treatments
combined_mapping[combined_mapping == 0] <- 1
combined_mapping[,c(1:3)] <- combined_mapping[,c(1:3)] / rowSums(combined_mapping[,c(1:3)])

# Remove genes with unannotated pathways
combined_mapping <- subset(combined_mapping, pathway_annotation != 'none')

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

xenobiotics <- subset(combined_mapping, grepl('*Xenobiotics*', combined_mapping$pathway_annotation))



#----------------------------------------------------------------------------------------------------------------------------#

# amino sugars
acetylglucosamine <- subset(combined_mapping, grepl('*N-acetylglucosamine*', combined_mapping$gene_annotation))
acetylmannosamine <- subset(combined_mapping, grepl('*N-acetylmannosamine*', combined_mapping$gene_annotation))
acetylmuramate <- subset(combined_mapping, grepl('*N-acetylmuramate*', combined_mapping$gene_annotation))
amino_sugars <- rbind(acetylglucosamine, acetylmannosamine, acetylmuramate)

# Stickland substrates
proline <- subset(combined_mapping, grepl('*proline*', combined_mapping$gene_annotation))
glycine <- subset(combined_mapping, grepl('*glycine*', combined_mapping$gene_annotation))
arginine <- subset(combined_mapping, grepl('*arginine*', combined_mapping$gene_annotation))
threonine <- subset(combined_mapping, grepl('*threonine*', combined_mapping$gene_annotation))
methionine <- subset(combined_mapping, grepl('*methionine*', combined_mapping$gene_annotation))
serine <- subset(combined_mapping, grepl('*serine*', combined_mapping$gene_annotation))
alanine <- subset(combined_mapping, grepl('*alanine*', combined_mapping$gene_annotation))
stickland <- rbind(proline, glycine, arginine, threonine, methionine, serine, alanine)

# monosaccharides
galactose <- subset(combined_mapping, grepl('*galactose*', combined_mapping$gene_annotation))
tagatose <- subset(combined_mapping, grepl('*tagatose*', combined_mapping$gene_annotation))
trehalose <- subset(combined_mapping, grepl('*trehalose*', combined_mapping$gene_annotation))
mannose <- subset(combined_mapping, grepl('*mannose*', combined_mapping$gene_annotation))
xylose <- subset(combined_mapping, grepl('*xylose*', combined_mapping$gene_annotation))
ribose <- subset(combined_mapping, grepl('*ribose*', combined_mapping$gene_annotation))
fructose <- subset(combined_mapping, grepl('*fructose*', combined_mapping$gene_annotation))
glucose <- subset(combined_mapping, grepl('*glucose*', combined_mapping$gene_annotation))
maltose <- subset(combined_mapping, grepl('*maltose*', combined_mapping$gene_annotation))
lactose <- subset(combined_mapping, grepl('*lactose*', combined_mapping$gene_annotation))
sucrose <- subset(combined_mapping, grepl('*sucrose*', combined_mapping$gene_annotation))
monosaccharides <- rbind(galactose, tagatose, trehalose, mannose, xylose, ribose, fructose, glucose, maltose, lactose, sucrose)

# sugar alcohols
ribitol <- subset(combined_mapping, grepl('*ribitol*', combined_mapping$gene_annotation))
sorbitol <- subset(combined_mapping, grepl('*sorbitol*', combined_mapping$gene_annotation))
mannitol <- subset(combined_mapping, grepl('*mannitol*', combined_mapping$gene_annotation))
sugar_alcohols <- rbind(ribitol, sorbitol, mannitol)

# nucleosides
uridine <- subset(combined_mapping, grepl('*uridine*', combined_mapping$gene_annotation))
inosine <- subset(combined_mapping, grepl('*inosine*', combined_mapping$gene_annotation))
adenosine <- subset(combined_mapping, grepl('*adenosine*', combined_mapping$gene_annotation))
nucleosides <- rbind(uridine, inosine, adenosine)

# short-chain fatty acids
butyrate <- subset(combined_mapping, grepl('*butyrate*', combined_mapping$gene_annotation))
valerate <- subset(combined_mapping, grepl('*valerate*', combined_mapping$gene_annotation))
acetate <- subset(combined_mapping, grepl('*acetate*', combined_mapping$gene_annotation))
scfas <- rbind(butyrate, valerate, acetate)

#-------------------------------------------------------------------------------------------------------------------------#

# Define which pathway to plot and the ouput file name
pathway1 <- cdf_carbohydates
pathway2 <- cdf_amino_acids
point_color1 <- 'blue2'
point_color2 <- 'red2'
pathway_name1 <- 'cdf_carbohydates'
pathway_name2 <- 'cdf_amino_acids'
#plot_file <- '~/Desktop/figures/cdf_amino_acids.streptomycin.pdf'

# Plot it!
#pdf(file=plot_file, width=7, height=6)
triplot(x=combined_mapping$cefoperazone, y=combined_mapping$clindamycin, z=combined_mapping$streptomycin,
        label=c('Cefoperazone', 'Clindamycin', 'Streptomycin'), pch=16, col='gray25', grid=FALSE, center=TRUE)
tripoints(x=pathway1$cefoperazone, y=pathway1$clindamycin, z=pathway1$streptomycin, cex=averages, pch=21, bg=point_color1, col='black')
tripoints(x=pathway2$cefoperazone, y=pathway2$clindamycin, z=pathway2$streptomycin, cex=averages, pch=21, bg=point_color2, col='black')
legend('topleft', legend=c(pathway_name1, pathway_name2), bty='n', cex=1.5, ncol=1, pch=21, pt.cex=2.5, pt.bg=c(point_color1, point_color2), col='black')



#dev.off()








selected <- subset(combined_mapping, grepl('*rginine*', combined_mapping$pathway_annotation))



