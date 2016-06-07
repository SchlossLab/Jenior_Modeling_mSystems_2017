
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

# Load in data
cefoperazone <- read.delim(cefoperazone_file, sep='\t', header=FALSE, row.names=1)
colnames(cefoperazone) <- c('Cefoperazone', 'ko', 'gene', 'pathway')
clindamycin <- read.delim(clindamycin_file, sep='\t', header=FALSE, row.names=1)
colnames(clindamycin) <- c('Clindamycin', 'ko', 'gene', 'pathway')
streptomycin <- read.delim(streptomycin_file, sep='\t', header=FALSE, row.names=1)
colnames(streptomycin) <- c('Streptomycin', 'ko', 'gene', 'pathway')

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
sub_size <- round(min(colSums(combined_mapping[,1:3])) * 0.9) # 97930
combined_mapping$Cefoperazone <- t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size))
combined_mapping$Clindamycin <- t(rrarefy(combined_mapping$Clindamycin, sample=sub_size))
combined_mapping$Streptomycin <- t(rrarefy(combined_mapping$Streptomycin, sample=sub_size))

# Eliminate genes with no transcripts mapping
combined_mapping <- combined_mapping[rowSums(combined_mapping[,1:3]) != 0, ] 

# Calculate the factors to scale points size by
averages <- log10(rowSums(combined_mapping[, 1:3]) / 3) * 2.5

# Convert each gene into the fraction of the transcription for that gene across treatments
combined_mapping[combined_mapping == 0] <- 1
combined_mapping[,1:3] <- combined_mapping[, 1:3] / rowSums(combined_mapping[,1:3])

# Remove genes with unannotated pathways
combined_mapping <- subset(combined_mapping, pathway != 'none')

# Free up memory
rm(cefoperazone_file, clindamycin_file, streptomycin_file)
rm(cefoperazone, clindamycin, streptomycin)

#----------------------------------------------------------------------------------------------------------------------------#

# Pathways 

# Highest resolution - Known C. difficile 630 carbon sources
glycolysis <- subset(combined_mapping[,1:3], grepl('*Glycolysis_/_Gluconeogenesis*', combined_mapping$pathway))
amino_sugar <- subset(combined_mapping[,1:3], grepl('*Amino_sugar*', combined_mapping$pathway))
galactose <- subset(combined_mapping[,1:3], grepl('*Galactose*', combined_mapping$pathway))
fructose_mannose <- subset(combined_mapping[,1:3], grepl('*Fructose_and_mannose*', combined_mapping$pathway))
starch_sucrose <- subset(combined_mapping[,1:3], grepl('*Starch_and_sucrose*', combined_mapping$pathway))
glycan_degradation <- subset(combined_mapping[,1:3], grepl('*glycan_degradation*', combined_mapping$pathway))
cdf_carbohydates <- rbind(glycolysis, amino_sugar, galactose, fructose_mannose, starch_sucrose, glycan_degradation)

cysteine_methionine <- subset(combined_mapping[,1:3], grepl('*Cysteine_and_methionine*', combined_mapping$pathway))
valine_leucine_isoleucine <- subset(combined_mapping[,1:3], grepl('*Valine,_leucine_and_isoleucine*', combined_mapping$pathway))
glycine_serine_threonine <- subset(combined_mapping[,1:3], grepl('*Glycine,_serine_and_threonine*', combined_mapping$pathway))
alanine_aspartate_glutamate <- subset(combined_mapping[,1:3], grepl('*Alanine,_aspartate_and_glutamate*', combined_mapping$pathway))
cdf_amino_acids <- rbind(cysteine_methionine, valine_leucine_isoleucine, glycine_serine_threonine, alanine_aspartate_glutamate)

cdf_carbon_sources <- rbind(cdf_carbohydates, cdf_amino_acids)


butanoate <- subset(combined_mapping[,1:3], grepl('*Butanoate*', combined_mapping$pathway))
pentose_glucuronate <- subset(combined_mapping[,1:3], grepl('*Pentose_and_glucuronate*', combined_mapping$pathway))
phenylalanine_tyrosine_tryptophan <- subset(combined_mapping[,1:3], grepl('*Phenylalanine,_tyrosine_and_tryptophan*', combined_mapping$pathway))
tca_cycle <- subset(combined_mapping[,1:3], grepl('*TCA_cycle*', combined_mapping$pathway))
lipopolysaccharide <- subset(combined_mapping[,1:3], grepl('*Lipopolysaccharide*', combined_mapping$pathway))
oxidative_phosphorylation <- subset(combined_mapping[,1:3], grepl('*Oxidative_phosphorylation*', combined_mapping$pathway))
glyoxylate_dicarboxylate <- subset(combined_mapping[,1:3], grepl('*Glyoxylate_and_dicarboxylate*', combined_mapping$pathway))
antibiotics <- subset(combined_mapping[,1:3], grepl('*antibiotic*', combined_mapping$pathway))
polyketide_sugar_biosynthesis <- subset(combined_mapping[,1:3], grepl('*Polyketide_sugar_unit_biosynthesis*', combined_mapping$pathway))
lysine <- subset(combined_mapping[,1:3], grepl('*Lysine*', combined_mapping$pathway))
arginine_proline <- subset(combined_mapping[,1:3], grepl('*Arginine_and_proline*', combined_mapping$pathway))
pyruvate <- subset(combined_mapping[,1:3], grepl('*Pyruvate*', combined_mapping$pathway))
peptidoglycan <- subset(combined_mapping[,1:3], grepl('*Peptidoglycan*', combined_mapping$pathway))
secondary_bile <- subset(combined_mapping[,1:3], grepl('*Secondary_bile*', combined_mapping$pathway))
# More general scale
energy_metabolism <- subset(combined_mapping[,1:3], grepl('*Energy_metabolism*', combined_mapping$pathway))
carbohydrate_metabolism <- subset(combined_mapping[,1:3], grepl('*Carbohydrate_metabolism*', combined_mapping$pathway))
amino_acid_metabolism <- subset(combined_mapping[,1:3], grepl('*Amino_acid_metabolism*', combined_mapping$pathway))
nucleotide_metabolism <- subset(combined_mapping[,1:3], grepl('*Nucleotide_metabolism*', combined_mapping$pathway))
lipid_metabolism <- subset(combined_mapping[,1:3], grepl('*Lipid_metabolism*', combined_mapping$pathway))
cofactors_and_vitamins <- subset(combined_mapping[,1:3], grepl('*cofactors_and_vitamins*', combined_mapping$pathway))
diverse_environments <- subset(combined_mapping[,1:3], grepl('*diverse_environments*', combined_mapping$pathway))
sulfur_metabolism <- subset(combined_mapping[,1:3], grepl('*Sulfur_metabolism*', combined_mapping$pathway))
secondary_metabolites <- subset(combined_mapping[,1:3], grepl('*secondary_metabolites*', combined_mapping$pathway))

xenobiotics <- subset(combined_mapping[,1:3], grepl('*Xenobiotics*', combined_mapping$pathway))

#-------------------------------------------------------------------------------------------------------------------------#

# Defineplot details
rainbow <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
fox <- wes_palette("FantasticFox")
tick_labels <- c('10%','','30%','','50%','','70%','','90%')
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/cdf_abx_pathways.pdf'
#palette_plot(rainbow)

# Plot it!
pdf(file=plot_file, width=10, height=7)

triax.plot(x=combined_mapping[,1:3], pch=20, lty.grid=2, no.add=FALSE, 
           tick.labels=list(l=tick_labels,r=tick_labels,b=tick_labels),
           cex.axis=1.5, cex.ticks=1, cc.axes=TRUE, align.labels=TRUE, 
           show.grid=TRUE, mar=c(4.5,0,0,0), at=seq(0.1,0.9, by=0.1))
lines(x=c(0.25,0.75), y=c(0.433,0.433))
lines(x=c(0.25,0.5), y=c(0.433,0))
lines(x=c(0.5,0.75), y=c(0,0.433))
lines(x=c(0,0.75), y=c(0,0.433))
lines(x=c(0.5,0.5), y=c(0,0.865))
lines(x=c(0.25,1), y=c(0.433,0))

# Plot points with color +/- relative size
triax.points(hexose[,c(2,1,3)], cex=1.7, col.symbols='black', pch=21, bg.symbols=fox[1])
#triax.points(hexose[,c(2,1,3)], cex=averages, col.symbols='black', pch=21, bg.symbols=fox[1])
triax.points(pentose[,c(2,1,3)], cex=1.7, col.symbols='black', pch=21, bg.symbols=rainbow[7])
#triax.points(pentose[,c(2,1,3)], cex=averages, col.symbols='black', pch=21, bg.symbols=rainbow[7])
triax.points(stickland[,c(2,1,3)], cex=1.7, col.symbols='black', pch=21, bg.symbols=fox[3])
#triax.points(stickland[,c(2,1,3)], cex=averages, col.symbols='black', pch=21, bg.symbols=fox[3])
triax.points(sugar_alcohols[,c(2,1,3)], cex=1.7, col.symbols='black', pch=21, bg.symbols=rainbow[1])
#triax.points(sugar_alcohols[,c(2,1,3)], cex=averages, col.symbols='black', pch=21, bg.symbols=rainbow[1])
triax.points(butyrate[,c(2,1,3)], cex=1.7, col.symbols='black', pch=21, bg.symbols=fox[5])
#triax.points(butyrate[,c(2,1,3)], cex=averages, col.symbols='black', pch=21, bg.symbols=fox[5])

legend(x=-0.3, y=1, legend=c('6-carbon sugars','5-carbon sugars', 'Stickland substrates', 'Sugar alcohols', 'Butyrate production'), 
       cex=1.5, ncol=1, pch=21, pt.cex=2.5, col='black', pt.bg=c(fox[1],rainbow[7],fox[3],rainbow[1],fox[5]))

dev.off()



