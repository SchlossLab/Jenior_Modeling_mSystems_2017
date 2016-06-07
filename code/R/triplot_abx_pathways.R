
deps <- c('vegan', 'plotrix', 'wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

# Define and check color palette
palette_plot <- function(col, border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}

# Define input file names
cefoperazone_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/cefoperazone_630.RNA_reads2cdf630.norm.annotated.txt'
clindamycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/clindamycin_630.RNA_reads2cdf630.norm.annotated.txt'
streptomycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/streptomycin_630.RNA_reads2cdf630.norm.annotated.txt'
germfree_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/germfree.RNA_reads2cdf630.norm.annotated.txt'

# Define which pathway to plot and the ouput file name
plot_file <- '~/Desktop/cdf_nutrients.pdf'

# Load in data
cefoperazone <- read.delim(cefoperazone_file, sep='\t', header=FALSE, row.names=1)
colnames(cefoperazone) <- c('cefoperazone', 'ko', 'gene', 'pathway')
clindamycin <- read.delim(clindamycin_file, sep='\t', header=FALSE, row.names=1)
colnames(clindamycin) <- c('clindamycin', 'ko', 'gene', 'pathway')
streptomycin <- read.delim(streptomycin_file, sep='\t', header=FALSE, row.names=1)
colnames(streptomycin) <- c('streptomycin', 'ko', 'gene', 'pathway')
germfree <- read.delim(germfree_file, sep='\t', header=FALSE, row.names=1)
colnames(germfree) <- c('germfree', 'ko', 'gene', 'pathway')

#-------------------------------------------------------------------------------------------------------------------------#

# Format data for merging
cefoperazone$ko <- NULL
cefoperazone$gene <- NULL
cefoperazone$pathway <- NULL
clindamycin$ko <- NULL
clindamycin$gene <- NULL
clindamycin$pathway <- NULL

# Merge tables
combined_mapping <- merge(cefoperazone, clindamycin, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
combined_mapping <- merge(combined_mapping, streptomycin, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL

# Rarefy mappings to be equal within sequencing type
sub_size <- round(min(colSums(combined_mapping[,c(1:3)]))*0.9) # 34050
combined_mapping$cefoperazone <- t(rrarefy(combined_mapping$cefoperazone, sample=sub_size))
combined_mapping$clindamycin <- t(rrarefy(combined_mapping$clindamycin, sample=sub_size))
combined_mapping$streptomycin <- t(rrarefy(combined_mapping$streptomycin, sample=sub_size))

# Eliminate genes with no transcripts mapping
combined_mapping <- combined_mapping[rowSums(combined_mapping[,c(1:3)]) != 0, ] 

# Calculate the factors to scale points size by
averages <- log10(rowSums(combined_mapping[,c(1:3)]) / 3) * 2.5

# Convert each gene into the fraction of the transcription for that gene across treatments
combined_mapping[combined_mapping == 0] <- 1
combined_mapping[,c(1:3)] <- combined_mapping[,c(1:3)] / rowSums(combined_mapping[,c(1:3)])

# Remove genes with unannotated pathways
combined_mapping <- subset(combined_mapping, pathway != 'none')

# Free up memory
rm(cefoperazone_file, clindamycin_file, streptomycin_file, germfree_file)
rm(cefoperazone, clindamycin, streptomycin, germfree)

#-------------------------------------------------------------------------------------------------------------------------#

# Pathways 

# Highest resolution - Known C. difficile 630 carbon sources
glycolysis <- subset(combined_mapping, grepl('*Glycolysis_/_Gluconeogenesis*', combined_mapping$pathway))
amino_sugar <- subset(combined_mapping, grepl('*Amino_sugar*', combined_mapping$pathway))
galactose <- subset(combined_mapping, grepl('*Galactose*', combined_mapping$pathway))
fructose_mannose <- subset(combined_mapping, grepl('*Fructose_and_mannose*', combined_mapping$pathway))
starch_sucrose <- subset(combined_mapping, grepl('*Starch_and_sucrose*', combined_mapping$pathway))
glycan_degradation <- subset(combined_mapping, grepl('*glycan_degradation*', combined_mapping$pathway))
cdf_carbohydates <- rbind(glycolysis, amino_sugar, galactose, fructose_mannose, starch_sucrose, glycan_degradation)

cysteine_methionine <- subset(combined_mapping, grepl('*Cysteine_and_methionine*', combined_mapping$pathway))
valine_leucine_isoleucine <- subset(combined_mapping, grepl('*Valine,_leucine_and_isoleucine*', combined_mapping$pathway))
glycine_serine_threonine <- subset(combined_mapping, grepl('*Glycine,_serine_and_threonine*', combined_mapping$pathway))
alanine_aspartate_glutamate <- subset(combined_mapping, grepl('*Alanine,_aspartate_and_glutamate*', combined_mapping$pathway))
cdf_amino_acids <- rbind(cysteine_methionine, valine_leucine_isoleucine, glycine_serine_threonine, alanine_aspartate_glutamate)

cdf_carbon_sources <- rbind(cdf_carbohydates, cdf_amino_acids)


butanoate <- subset(combined_mapping, grepl('*Butanoate*', combined_mapping$pathway))
pentose_glucuronate <- subset(combined_mapping, grepl('*Pentose_and_glucuronate*', combined_mapping$pathway))
phenylalanine_tyrosine_tryptophan <- subset(combined_mapping, grepl('*Phenylalanine,_tyrosine_and_tryptophan*', combined_mapping$pathway))
tca_cycle <- subset(combined_mapping, grepl('*TCA_cycle*', combined_mapping$pathway))
lipopolysaccharide <- subset(combined_mapping, grepl('*Lipopolysaccharide*', combined_mapping$pathway))
oxidative_phosphorylation <- subset(combined_mapping, grepl('*Oxidative_phosphorylation*', combined_mapping$pathway))
glyoxylate_dicarboxylate <- subset(combined_mapping, grepl('*Glyoxylate_and_dicarboxylate*', combined_mapping$pathway))
antibiotics <- subset(combined_mapping, grepl('*antibiotic*', combined_mapping$pathway))
polyketide_sugar_biosynthesis <- subset(combined_mapping, grepl('*Polyketide_sugar_unit_biosynthesis*', combined_mapping$pathway))
lysine <- subset(combined_mapping, grepl('*Lysine*', combined_mapping$pathway))
arginine_proline <- subset(combined_mapping, grepl('*Arginine_and_proline*', combined_mapping$pathway))
pyruvate <- subset(combined_mapping, grepl('*Pyruvate*', combined_mapping$pathway))
peptidoglycan <- subset(combined_mapping, grepl('*Peptidoglycan*', combined_mapping$pathway))
secondary_bile <- subset(combined_mapping, grepl('*Secondary_bile*', combined_mapping$pathway))
# More general scale
energy_metabolism <- subset(combined_mapping, grepl('*Energy_metabolism*', combined_mapping$pathway))
carbohydrate_metabolism <- subset(combined_mapping, grepl('*Carbohydrate_metabolism*', combined_mapping$pathway))
amino_acid_metabolism <- subset(combined_mapping, grepl('*Amino_acid_metabolism*', combined_mapping$pathway))
nucleotide_metabolism <- subset(combined_mapping, grepl('*Nucleotide_metabolism*', combined_mapping$pathway))
lipid_metabolism <- subset(combined_mapping, grepl('*Lipid_metabolism*', combined_mapping$pathway))
cofactors_and_vitamins <- subset(combined_mapping, grepl('*cofactors_and_vitamins*', combined_mapping$pathway))
diverse_environments <- subset(combined_mapping, grepl('*diverse_environments*', combined_mapping$pathway))
sulfur_metabolism <- subset(combined_mapping, grepl('*Sulfur_metabolism*', combined_mapping$pathway))
secondary_metabolites <- subset(combined_mapping, grepl('*secondary_metabolites*', combined_mapping$pathway))

xenobiotics <- subset(combined_mapping, grepl('*Xenobiotics*', combined_mapping$pathway))



#----------------------------------------------------------------------------------------------------------------------------#

# Remove non-useful groups
combined_mapping$ko <- NULL

# amino sugars
acetylglucosamine <- subset(combined_mapping, grepl('*N-acetylglucosamine*', combined_mapping$gene))
acetylmannosamine <- subset(combined_mapping, grepl('*N-acetylmannosamine*', combined_mapping$gene))
acetylmuramate <- subset(combined_mapping, grepl('*N-acetylmuramate*', combined_mapping$gene))
amino_sugars <- rbind(acetylglucosamine, acetylmannosamine, acetylmuramate)
rm(acetylglucosamine, acetylmannosamine, acetylmuramate)

# Stickland substrates
proline <- subset(combined_mapping, grepl('*proline*', combined_mapping$gene))
glycine <- subset(combined_mapping, grepl('*glycine*', combined_mapping$gene))
arginine <- subset(combined_mapping, grepl('*arginine*', combined_mapping$gene))
threonine <- subset(combined_mapping, grepl('*threonine*', combined_mapping$gene))
methionine <- subset(combined_mapping, grepl('*methionine*', combined_mapping$gene))
serine <- subset(combined_mapping, grepl('*serine*', combined_mapping$gene))
alanine <- subset(combined_mapping, grepl('*alanine*', combined_mapping$gene))
stickland <- rbind(proline, glycine, arginine, threonine, methionine, serine, alanine)
rm(proline, glycine, arginine, threonine, methionine, serine, alanine)

# Saccharides
# hexose
galactose <- subset(combined_mapping, grepl('*galactose*', combined_mapping$gene))
mannose <- subset(combined_mapping, grepl('*mannose*', combined_mapping$gene))
glucose <- subset(combined_mapping, grepl('*glucose*', combined_mapping$gene))
tagatose <- subset(combined_mapping, grepl('*tagatose*', combined_mapping$gene))
fructose <- subset(combined_mapping, grepl('*fructose*', combined_mapping$gene)) #keto-
hexose <- rbind(galactose, mannose, glucose, tagatose, fructose)
# pentose
xylose <- subset(combined_mapping, grepl('*xylose*', combined_mapping$gene))
ribose <- subset(combined_mapping, grepl('*ribose*', combined_mapping$gene))
pentose <- rbind(xylose, ribose)
# disaccharides
sucrose <- subset(combined_mapping, grepl('*sucrose*', combined_mapping$gene))
lactose <- subset(combined_mapping, grepl('*lactose*', combined_mapping$gene))
maltose <- subset(combined_mapping, grepl('*maltose*', combined_mapping$gene))
trehalose <- subset(combined_mapping, grepl('*trehalose*', combined_mapping$gene))
disaccharides <- rbind(sucrose, lactose, maltose, trehalose)
# all
saccharides <- rbind(galactose, tagatose, trehalose, mannose, xylose, ribose, fructose, glucose, maltose, lactose, sucrose)
rm(galactose, tagatose, trehalose, mannose, xylose, ribose, glucose, maltose, lactose, sucrose)

# sugar alcohols
ribitol <- subset(combined_mapping, grepl('*ribitol*', combined_mapping$gene))
sorbitol <- subset(combined_mapping, grepl('*sorbitol*', combined_mapping$gene))
mannitol <- subset(combined_mapping, grepl('*mannitol*', combined_mapping$gene))
sugar_alcohols <- rbind(ribitol, sorbitol, mannitol)
rm(ribitol, sorbitol, mannitol)

# nucleosides
uridine <- subset(combined_mapping, grepl('*uridine*', combined_mapping$gene))
inosine <- subset(combined_mapping, grepl('*inosine*', combined_mapping$gene))
adenosine <- subset(combined_mapping, grepl('*adenosine*', combined_mapping$gene))
nucleosides <- rbind(uridine, inosine, adenosine)
rm(uridine, inosine, adenosine)

# short-chain fatty acids
butyrate <- subset(combined_mapping, grepl('*butyrate*', combined_mapping$gene))
valerate <- subset(combined_mapping, grepl('*valerate*', combined_mapping$gene))
acetate <- subset(combined_mapping, grepl('*acetate*', combined_mapping$gene))
scfas <- rbind(butyrate, valerate, acetate)
rm(butyrate, valerate, acetate)

#-------------------------------------------------------------------------------------------------------------------------#

# Define a color palette
rainbow <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
fox <- wes_palette("FantasticFox")
tick_labels <- c('10%','','30%','','50%','','70%','','90%')
plotting_map_data <- combined_mapping
plotting_map_data$ko <- NULL
plotting_map_data$gene <- NULL
plotting_map_data$pathway <- NULL
#palette_plot(rainbow)

# Plot it!
#pdf(file=plot_file, width=7, height=6)
triax.plot(x=plotting_map_data, pch=20, lty.grid=2, no.add=FALSE, 
           axis.labels=c('Cefoperazone ', 'Clindamycin', 'Streptomycin'),  
           tick.labels=list(l=tick_labels, r=tick_labels, b=tick_labels), 
           cex.axis=1.5, cex.ticks=1, cc.axes=TRUE, align.labels=TRUE, 
           show.grid=TRUE, mar=c(4.5,0,0,0), at=seq(0.1,0.9, by=0.1))
lines(x=c(0.25,0.75), y=c(0.433,0.433))
lines(x=c(0.25,0.5), y=c(0.433,0))
lines(x=c(0.5,0.75), y=c(0,0.433))

# Plot points with color and size
# Simple saccharides
# hexose, pentose, disaccharides, all = saccharides
selected <- hexose
selected$gene <- NULL
current_color <- fox[1]
triax.points(selected, cex=1.7, col.symbols='black', pch=21, bg.symbols=current_color)
#triax.points(selected, cex=averages, col.symbols='black', pch=21, bg.symbols=current_color)
selected <- pentose 
selected$gene <- NULL
current_color <- fox[2]
triax.points(selected, cex=1.7, col.symbols='black', pch=21, bg.symbols=current_color)
#triax.points(selected, cex=averages, col.symbols='black', pch=21, bg.symbols=current_color)
#selected <- saccharides # all
#selected$gene <- NULL
#current_color <- rainbow[14]
#triax.points(selected, cex=averages, col.symbols='black', pch=21, bg.symbols=current_color)

# peptides
selected <- stickland
selected$gene <- NULL
current_color <- fox[3]
triax.points(selected, cex=1.7, col.symbols='black', pch=21, bg.symbols=current_color)
#triax.points(selected, cex=averages, col.symbols='black', pch=21, bg.symbols=current_color)

# sugar alcohols
selected <- sugar_alcohols
selected$gene <- NULL
current_color <- fox[4]
triax.points(selected, cex=1.7, col.symbols='black', pch=21, bg.symbols=current_color)
#triax.points(selected, cex=averages, col.symbols='black', pch=21, bg.symbols=current_color)

# short-schain fatty acids
selected <- scfas
selected$gene <- NULL
current_color <- fox[5]
triax.points(selected, cex=1.7, col.symbols='black', pch=21, bg.symbols=current_color)
#triax.points(selected, cex=averages, col.symbols='black', pch=21, bg.symbols=current_color)

#dev.off()




selected <- glycolysis
selected$ko <- NULL
selected$gene <- NULL
selected$pathway <- NULL
current_color <- fox[1]
triax.points(selected, cex=1.7, col.symbols='black', pch=21, bg.symbols=current_color)
selected <- amino_sugar
selected$ko <- NULL
selected$gene <- NULL
selected$pathway <- NULL
current_color <- fox[2]
triax.points(selected, cex=1.7, col.symbols='black', pch=21, bg.symbols=current_color)
selected <- galactose
selected$ko <- NULL
selected$gene <- NULL
selected$pathway <- NULL
current_color <- fox[3]
triax.points(selected, cex=1.7, col.symbols='black', pch=21, bg.symbols=current_color)
selected <- fructose_mannose
selected$ko <- NULL
selected$gene <- NULL
selected$pathway <- NULL
current_color <- fox[4]
triax.points(selected, cex=1.7, col.symbols='black', pch=21, bg.symbols=current_color)
selected <- starch_sucrose
selected$ko <- NULL
selected$gene <- NULL
selected$pathway <- NULL
current_color <- fox[5]
triax.points(selected, cex=1.7, col.symbols='black', pch=21, bg.symbols=current_color)


#pdf(file=plot_file, width=6, height=5)
#plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-4,4), ylim=c(-4,4))
legend(x=0, y=1, legend=c('pathway1','pathway2'), 
    cex=1.5, ncol=1, pch=21, pt.cex=2.5, col='black', 
    pt.bg=fox[1:5])
#dev.off()



