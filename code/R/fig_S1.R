
deps <- c('wesanderson','vegan');
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
rm(cefoperazone, clindamycin, streptomycin, germfree)

# Remove genes not in a pathway
combined_mapping <- subset(combined_mapping, pathway != 'none')

# Rarefy mappings to be equal within sequencing type
sub_size <- round(min(colSums(combined_mapping[,1:4])) * 0.9)
combined_mapping$Cefoperazone <- t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size))
combined_mapping$Clindamycin <- t(rrarefy(combined_mapping$Clindamycin, sample=sub_size))
combined_mapping$Streptomycin <- t(rrarefy(combined_mapping$Streptomycin, sample=sub_size))
combined_mapping$Germfree <- t(rrarefy(combined_mapping$Germfree, sample=sub_size))

# Eliminate genes with no transcripts mapping
combined_mapping <- combined_mapping[rowSums(combined_mapping[,1:4]) != 0, ] 

#--------------------------------------------------------------------------------------------------------------#

# Subset by KEGG catagory and pool
all_carbohydrate <- subset(combined_mapping, grepl('Carbohydrate_metabolism', combined_mapping$pathway))
carbohydrate <- t(as.data.frame(colSums(all_carbohydrate[,1:4])))
energy <- subset(combined_mapping, grepl('Energy_metabolism', combined_mapping$pathway))
energy <- t(as.data.frame(colSums(energy[,1:4])))
lipid <- subset(combined_mapping, grepl('Lipid_metabolism', combined_mapping$pathway))
lipid <- t(as.data.frame(colSums(lipid[,1:4])))
nucleotide <- subset(combined_mapping, grepl('Nucleotide_metabolism', combined_mapping$pathway))
nucleotide <- t(as.data.frame(colSums(nucleotide[,1:4])))
amino_acid <- subset(combined_mapping, grepl('Amino_acid_metabolism', combined_mapping$pathway))
other_amino_acids <- subset(combined_mapping, grepl('Metabolism_of_other_amino_acids', combined_mapping$pathway))
all_amino_acid <- rbind(amino_acid, other_amino_acids)
amino_acid <- t(as.data.frame(colSums(amino_acid[,1:4])))
other_amino_acids <- t(as.data.frame(colSums(other_amino_acids[,1:4])))
amino_acids <- amino_acid + other_amino_acids
rm(amino_acid, other_amino_acids)
glycan <- subset(combined_mapping, grepl('Glycan_biosynthesis_and_metabolism', combined_mapping$pathway))
glycan <- t(as.data.frame(colSums(glycan[,1:4])))
cofactors_and_vitamins <- subset(combined_mapping, grepl('Metabolism_of_cofactors_and_vitamins', combined_mapping$pathway))
cofactors_and_vitamins <- t(as.data.frame(colSums(cofactors_and_vitamins[,1:4])))
terpenoids_and_polyketides <- subset(combined_mapping, grepl('Metabolism_of_terpenoids_and_polyketides', combined_mapping$pathway))
terpenoids_and_polyketides <- t(as.data.frame(colSums(terpenoids_and_polyketides[,1:4])))
secondary_metabolites <- subset(combined_mapping, grepl('Biosynthesis_of_other_secondary_metabolites', combined_mapping$pathway))
secondary_metabolites <- t(as.data.frame(colSums(secondary_metabolites[,1:4])))
xenobiotics <- subset(combined_mapping, grepl('Xenobiotics_biodegradation_and_metabolism', combined_mapping$pathway))
xenobiotics <- t(as.data.frame(colSums(xenobiotics[,1:4])))
genetics <- subset(combined_mapping, grepl('Genetic_Information_Processing', combined_mapping$pathway))
genetics <- t(as.data.frame(colSums(genetics[,1:4])))
motility <- subset(combined_mapping, grepl('Cell_motility', combined_mapping$pathway))
motility <- t(as.data.frame(colSums(motility[,1:4])))
rm(combined_mapping)

# Assemble final table
pooled_mapping <- rbind(carbohydrate, energy, lipid, nucleotide, amino_acids, glycan,
                        cofactors_and_vitamins, terpenoids_and_polyketides, secondary_metabolites,
                        xenobiotics, genetics, motility)
rm(carbohydrate, energy, lipid, nucleotide, amino_acids, glycan, cofactors_and_vitamins, terpenoids_and_polyketides, secondary_metabolites, xenobiotics, genetics, motility)
rownames(pooled_mapping) <- c('Carbohydrate metabolism','Energy metabolism','Lipid metabolism','Nucleotide metabolism',
                              'Amino acid metabolism','Glycan metabolism', 'Cofactor/Vitamin metabolism','Terpenoid/Polyketide metabolism',
                              'Secondary metabolite synthesis', 'Xenobiotic metabolism','Genetic Information Processing','Cell motility')
transformed_mapping <- log10(pooled_mapping)
rm(pooled_mapping)

# Agregate carbohydrate and amino acid tables seperately

# Carbohydrates
glycolysis <- subset(all_carbohydrate, grepl('Glycolysis_/_Gluconeogenesis', all_carbohydrate$pathway))
glycolysis <- t(as.data.frame(colSums(glycolysis[,1:4])))
starch <- subset(all_carbohydrate, grepl('Starch_and_sucrose_metabolism', all_carbohydrate$pathway))
starch <- t(as.data.frame(colSums(starch[,1:4])))
citrate_cycle <- subset(all_carbohydrate, grepl('Citrate_cycle', all_carbohydrate$pathway))
citrate_cycle <- t(as.data.frame(colSums(citrate_cycle[,1:4])))
fructose <- subset(all_carbohydrate, grepl('Fructose_and_mannose_metabolism', all_carbohydrate$pathway))
fructose <- t(as.data.frame(colSums(fructose[,1:4])))
butanoate <- subset(all_carbohydrate, grepl('Butanoate_metabolism', all_carbohydrate$pathway))
butanoate <- t(as.data.frame(colSums(butanoate[,1:4])))
pentose_phosphate <- subset(all_carbohydrate, grepl('Pentose_phosphate_pathway', all_carbohydrate$pathway))
pentose_phosphate <- t(as.data.frame(colSums(pentose_phosphate[,1:4])))
amino_sugar <- subset(all_carbohydrate, grepl('Amino_sugar_and_nucleotide_sugar_metabolism', all_carbohydrate$pathway))
amino_sugar <- t(as.data.frame(colSums(amino_sugar[,1:4])))
pyruvate <- subset(all_carbohydrate, grepl('Pyruvate_metabolism', all_carbohydrate$pathway))
pyruvate <- t(as.data.frame(colSums(pyruvate[,1:4])))
glyoxylate <- subset(all_carbohydrate, grepl('Glyoxylate_and_dicarboxylate_metabolism', all_carbohydrate$pathway))
glyoxylate <- t(as.data.frame(colSums(glyoxylate[,1:4])))
fatty_acid <- subset(all_carbohydrate, grepl('Fatty_acid_degradation', all_carbohydrate$pathway))
fatty_acid <- t(as.data.frame(colSums(fatty_acid[,1:4])))
ascorbate <- subset(all_carbohydrate, grepl('Ascorbate_and_aldarate_metabolism', all_carbohydrate$pathway))
ascorbate <- t(as.data.frame(colSums(ascorbate[,1:4])))
galactose <- subset(all_carbohydrate, grepl('Galactose_metabolism', all_carbohydrate$pathway))
galactose <- t(as.data.frame(colSums(galactose[,1:4])))
pooled_carbohydrate <- rbind(glycolysis, starch, citrate_cycle, fructose, butanoate, pentose_phosphate, 
                             amino_sugar, pyruvate, glyoxylate, fatty_acid, ascorbate, galactose)
rm(glycolysis, starch, citrate_cycle, fructose, butanoate, pentose_phosphate, 
   amino_sugar, pyruvate, glyoxylate, fatty_acid, ascorbate, galactose)
rownames(pooled_carbohydrate) <- c('Glycolysis/Gluconeogenesis', 'Starch/Sucrose metabolism', 'Citrate acid cycle', 
                                   'Fructose/Mannose metabolism', 'Butanoate metabolism', 'Pentose phosphate pathway', 
                                   'Amino/Nucleotide sugar metabolism', 'Pyruvate metabolism', 'Glyoxylate/Dicarboxylate metabolism', 
                                   'Fatty acid degradation', 'Ascorbate/Aldarate metabolism', 'Galactose metabolism')

# Log transform mappings
transformed_carbohydrate <- log10(pooled_carbohydrate)
rm(pooled_carbohydrate)

# Amino acids
alanine_aspartate_glutamate <- subset(all_amino_acid, grepl('Alanine,_aspartate_and_glutamate_metabolism', all_amino_acid$pathway))
alanine_aspartate_glutamate <- t(as.data.frame(colSums(alanine_aspartate_glutamate[,1:4])))
alanine <- subset(all_amino_acid, grepl('Alanine_metabolism', all_amino_acid$pathway))
alanine <- t(as.data.frame(colSums(alanine[,1:4])))
glutamine <- subset(all_amino_acid, grepl('D-Glutamine_and_D-glutamate_metabolism', all_amino_acid$pathway))
glutamine <- t(as.data.frame(colSums(glutamine[,1:4])))
alanine_aspartate_glutamate <- alanine_aspartate_glutamate + alanine + glutamine
rm(alanine, glutamine)
cysteine <- subset(all_amino_acid, grepl('Cysteine_and_methionine_metabolism', all_amino_acid$pathway))
cysteine <- t(as.data.frame(colSums(cysteine[,1:4])))
arginine <- subset(all_amino_acid, grepl('Arginine_and_proline_metabolism', all_amino_acid$pathway))
arginine <- t(as.data.frame(colSums(arginine[,1:4])))
histidine <- subset(all_amino_acid, grepl('Histidine_metabolism', all_amino_acid$pathway))
histidine <- t(as.data.frame(colSums(histidine[,1:4])))
valine <- subset(all_amino_acid, grepl('Valine,_leucine_and_isoleucine_biosynthesis', all_amino_acid$pathway))
valine <- t(as.data.frame(colSums(valine[,1:4])))
glycine <- subset(all_amino_acid, grepl('Glycine,_serine_and_threonine_metabolism', all_amino_acid$pathway))
glycine <- t(as.data.frame(colSums(glycine[,1:4])))
lysine <- subset(all_amino_acid, grepl('Lysine_biosynthesis', all_amino_acid$pathway))
lysine <- t(as.data.frame(colSums(lysine[,1:4])))
phenylalanine_tyrosine_tryptophan <- subset(all_amino_acid, grepl('Phenylalanine,_tyrosine_and_tryptophan_biosynthesis', all_amino_acid$pathway))
phenylalanine_tyrosine_tryptophan <- t(as.data.frame(colSums(phenylalanine_tyrosine_tryptophan[,1:4])))
phenylalanine <- subset(all_amino_acid, grepl('Phenylalanine_metabolism', all_amino_acid$pathway))
phenylalanine <- t(as.data.frame(colSums(phenylalanine[,1:4])))
tyrosine <- subset(all_amino_acid, grepl('Tyrosine_metabolism', all_amino_acid$pathway))
tyrosine <- t(as.data.frame(colSums(tyrosine[,1:4])))
taurine <- subset(all_amino_acid, grepl('Taurine_and_hypotaurine_metabolism', all_amino_acid$pathway))
taurine <- t(as.data.frame(colSums(taurine[,1:4])))
selenocompound <- subset(all_amino_acid, grepl('Selenocompound_metabolism', all_amino_acid$pathway))
selenocompound <- t(as.data.frame(colSums(selenocompound[,1:4])))

pooled_amino_acid <- rbind(histidine, cysteine, arginine, alanine_aspartate_glutamate, valine, glycine, 
                           lysine, phenylalanine_tyrosine_tryptophan, phenylalanine, tyrosine, taurine, selenocompound)
rm(alanine_aspartate_glutamate, cysteine, arginine, histidine, valine, glycine, 
   lysine, phenylalanine_tyrosine_tryptophan, phenylalanine, tyrosine, taurine, selenocompound)
rownames(pooled_amino_acid) <- c('Histidine metabolism', 'Cysteine/Methionine metabolism', 
                                 'Arginine/Proline metabolism', 'Alanine/Aspartate/Glutamate metabolism', 'Valine/Leucine/Isoleucine synthesis', 
                                 'Glycine/Serine/Threonine metabolism', 'Lysine synthesis', 'Phenylalanine/Tyrosine/Tryptophan synthesis', 
                                 'Phenylalanine metabolism', 'Tyrosine metabolism', 'Taurine/Hypotaurine metabolism', 'Selenocompound metabolism')

# Log transform mappings
transformed_amino_acid <- log10(pooled_amino_acid)
rm(pooled_amino_acid)

#--------------------------------------------------------------------------------------------------------------#

# Set the color palette and output file name
darjeeling <- wes_palette("Darjeeling")
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figure_S1.pdf'

# Generate figure
pdf(file=plot_file, width=12, height=10)
layout(matrix(c(1,1,
                2,3), 
              nrow=2, ncol=2, byrow = TRUE))

# Global pathway annotations
par(las=1, mar=c(8,4,1,1))
barplot(t(transformed_mapping), col=c(darjeeling[1],darjeeling[2],darjeeling[4],darjeeling[5]), 
        beside=TRUE, xaxt='n', yaxt='n', ylab='Transcript Abundance (Log10)', ylim=c(0,5))
box()
axis(side=2, at=c(1:4), parse(text=paste(rep(10,4), '^', seq(1,4,1), sep='')), tick=TRUE, las=1)
abline(h=c(1:4), lty=2)
legend('topleft', legend=c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Gnotobiotic'), pt.cex=2, bty='n',
       pch=22, col='black', pt.bg=c(darjeeling[1],darjeeling[2],darjeeling[4],darjeeling[5]), ncol=2)
text(x=seq(4,59,5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(transformed_mapping), srt=45, adj=1, xpd=TRUE, cex=0.8)
mtext('A', side=2, line=2, las=2, adj=1, padj=-9, cex=1.5)

# Subpathways in Carbohydrate metabolism
par(las=1, mar=c(11,4,1,1))
barplot(t(transformed_carbohydrate), col=c(darjeeling[1],darjeeling[2],darjeeling[4],darjeeling[5]), 
        beside=TRUE, xaxt='n', yaxt='n', ylab='Transcript Abundance (Log10)', ylim=c(0,4))
box()
axis(side=2, at=c(1:3), parse(text=paste(rep(10,3), '^', seq(1,3,1), sep='')), tick=TRUE, las=1)
abline(h=c(1:3), lty=2)
text(x=seq(4,59,5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(transformed_carbohydrate), srt=45, adj=1, xpd=TRUE, cex=0.8)
mtext('B', side=2, line=2, las=2, adj=1, padj=-8, cex=1.5)
legend('topleft', legend='Carbohydrate Metabolism', bty='n')

# Subpathways in Amino Acid metabolism
par(las=1, mar=c(11,4,1,1))
barplot(t(transformed_amino_acid), col=c(darjeeling[1],darjeeling[2],darjeeling[4],darjeeling[5]), 
        beside=TRUE, xaxt='n', yaxt='n', ylab='Transcript Abundance (Log10)', ylim=c(0,4))
box()
axis(side=2, at=c(1:3), parse(text=paste(rep(10,3), '^', seq(1,3,1), sep='')), tick=TRUE, las=1)
abline(h=c(1:3), lty=2)
text(x=seq(4,59,5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(transformed_amino_acid), srt=45, adj=1, xpd=TRUE, cex=0.8)
mtext('C', side=2, line=2, las=2, adj=1, padj=-8, cex=1.5)
legend('topleft', legend='Amino Acid Metabolism', bty='n')

dev.off()
