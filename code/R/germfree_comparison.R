
deps <- c('vegan', 'wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

# Define input file names
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

# Rarefy mappings to be equal within sequencing type
sub_size <- round(min(colSums(combined_mapping[,1:4])) * 0.9) # 34050
combined_mapping$Cefoperazone <- t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size))
combined_mapping$Clindamycin <- t(rrarefy(combined_mapping$Clindamycin, sample=sub_size))
combined_mapping$Streptomycin <- t(rrarefy(combined_mapping$Streptomycin, sample=sub_size))
combined_mapping$Germfree <- t(rrarefy(combined_mapping$Germfree, sample=sub_size))
rm(sub_size)

# Log transform columns
combined_mapping[combined_mapping == 0] <- 1
combined_mapping$Cefoperazone <- log10(combined_mapping$Cefoperazone)
combined_mapping$Clindamycin <- log10(combined_mapping$Clindamycin)
combined_mapping$Streptomycin <- log10(combined_mapping$Streptomycin)
combined_mapping$Germfree <- log10(combined_mapping$Germfree)

# Eliminate genes with no transcripts mapping
combined_mapping <- combined_mapping[rowSums(combined_mapping[,1:4]) != 0.0, ] 

# Free up memory
rm(cefoperazone_file, clindamycin_file, streptomycin_file, germfree_file)
rm(cefoperazone, clindamycin, streptomycin, germfree)

#----------------------------------------------------------------------------------------------------------------------------#

# Subset by gene annotations

# amino sugars
acetylglucosamine <- subset(combined_mapping, grepl('*N-acetylglucosamine*', combined_mapping$gene))
acetylneuraminate <- subset(combined_mapping, grepl('*N-acetylneuraminate*', combined_mapping$gene))

# Stickland substrates
proline <- subset(combined_mapping, grepl('*proline*', combined_mapping$gene))
glycine <- subset(combined_mapping, grepl('*glycine*', combined_mapping$gene))
threonine <- subset(combined_mapping, grepl('*threonine*', combined_mapping$gene))
methionine <- subset(combined_mapping, grepl('*methionine*', combined_mapping$gene))
serine <- subset(combined_mapping, grepl('*serine*', combined_mapping$gene))
alanine <- subset(combined_mapping, grepl('*alanine*', combined_mapping$gene))

# Saccharides
# hexose
galactose <- subset(combined_mapping, grepl('*galactose*', combined_mapping$gene))
mannose <- subset(combined_mapping, grepl('*mannose*', combined_mapping$gene))
glucose <- subset(combined_mapping, grepl('*glucose*', combined_mapping$gene))
tagatose <- subset(combined_mapping, grepl('*tagatose*', combined_mapping$gene))
fructose <- subset(combined_mapping, grepl('*fructose*', combined_mapping$gene)) #keto-
# pentose
xylose <- subset(combined_mapping, grepl('*xylose*', combined_mapping$gene))
ribose <- subset(combined_mapping, grepl('*ribose*', combined_mapping$gene))
# disaccharides
lactose <- subset(combined_mapping, grepl('*lactose*', combined_mapping$gene))
maltose <- subset(combined_mapping, grepl('*maltose*', combined_mapping$gene))
trehalose <- subset(combined_mapping, grepl('*trehalose*', combined_mapping$gene))

# sugar alcohols
sorbitol <- subset(combined_mapping, grepl('*sorbitol*', combined_mapping$gene))
mannitol <- subset(combined_mapping, grepl('*mannitol*', combined_mapping$gene))

# nucleosides
uridine <- subset(combined_mapping, grepl('*uridine*', combined_mapping$gene))
adenosine <- subset(combined_mapping, grepl('*adenosine*', combined_mapping$gene))

# short-chain fatty acids
butyrate <- subset(combined_mapping, grepl('*butyrate*', combined_mapping$gene))
acetate <- subset(combined_mapping, grepl('*acetate*', combined_mapping$gene))
valerate <- subset(combined_mapping, grepl('*acetate*', combined_mapping$gene))

#-------------------------------------------------------------------------------------------------------------------------#

# or pathways...

# Highest resolution - Known C. difficile 630 carbon sources
glycolysis <- subset(combined_mapping, grepl('*Glycolysis_/_Gluconeogenesis*', combined_mapping$pathway))
amino_sugar <- subset(combined_mapping, grepl('*Amino_sugar*', combined_mapping$pathway))
galactose <- subset(combined_mapping, grepl('*Galactose*', combined_mapping$pathway))
fructose_mannose <- subset(combined_mapping, grepl('*Fructose_and_mannose*', combined_mapping$pathway))
starch_sucrose <- subset(combined_mapping, grepl('*Starch_and_sucrose*', combined_mapping$pathway))
cysteine_methionine <- subset(combined_mapping, grepl('*Cysteine_and_methionine*', combined_mapping$pathway))
valine_leucine_isoleucine <- subset(combined_mapping, grepl('*Valine,_leucine_and_isoleucine*', combined_mapping$pathway))
glycine_serine_threonine <- subset(combined_mapping, grepl('*Glycine,_serine_and_threonine*', combined_mapping$pathway))
alanine_aspartate_glutamate <- subset(combined_mapping, grepl('*Alanine,_aspartate_and_glutamate*', combined_mapping$pathway))



butanoate <- subset(combined_mapping, grepl('*Butanoate*', combined_mapping$pathway))
pentose_glucuronate <- subset(combined_mapping, grepl('*Pentose_and_glucuronate*', combined_mapping$pathway))
phenylalanine_tyrosine_tryptophan <- subset(combined_mapping, grepl('*Phenylalanine,_tyrosine_and_tryptophan*', combined_mapping$pathway))
tca_cycle <- subset(combined_mapping, grepl('*TCA_cycle*', combined_mapping$pathway))
oxidative_phosphorylation <- subset(combined_mapping, grepl('*Oxidative_phosphorylation*', combined_mapping$pathway))
glyoxylate_dicarboxylate <- subset(combined_mapping, grepl('*Glyoxylate_and_dicarboxylate*', combined_mapping$pathway))
lysine <- subset(combined_mapping, grepl('*Lysine*', combined_mapping$pathway))
arginine_proline <- subset(combined_mapping, grepl('*Arginine_and_proline*', combined_mapping$pathway))
pyruvate <- subset(combined_mapping, grepl('*Pyruvate*', combined_mapping$pathway))
peptidoglycan <- subset(combined_mapping, grepl('*Peptidoglycan*', combined_mapping$pathway))



rm(combined_mapping)




#-------------------------------------------------------------------------------------------------------------------------#

# Defineplot details
fox <- wes_palette("FantasticFox")
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/cdf_gf_comparison.pdf'

# Open a PDF
#pdf(file=plot_file, width=9, height=6.5)
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))


# dot plots (stacked)





























# Add the legend
legend('bottomright', legend=c('Antibiotic', 'Gnotobiotic'), cex=1, ncol=1, pch=21, pt.cex=2, col='black', pt.bg=c('firebrick','darkblue'))


legend('topright', legend=c(,,,), 
       cex=1, ncol=1, pch=21, pt.cex=2, col='black', pt.bg=c(fox[1],fox[2],fox[3],fox[4],fox[5]))



# Add figure label
text(x=-0.8, y=0.75, labels='B', font=2, cex=2)

#dev.off()



