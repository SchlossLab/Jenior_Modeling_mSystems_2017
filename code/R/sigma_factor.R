
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



# SUBSET BY SIGMA FACTOR HERE


agr # quorum sensing
ccp

# add virulence gene plot to supplement, tcdE





all_carbohydrate <- subset(combined_mapping, grepl('Carbohydrate_metabolism', combined_mapping$pathway))
carbohydrate <- t(as.data.frame(colSums(all_carbohydrate[,1:4])))


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

#--------------------------------------------------------------------------------------------------------------#

# Set the color palette and output file name
darjeeling <- wes_palette("Darjeeling")
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figure_3.pdf'

# Generate figure
pdf(file=plot_file, width=12, height=10)

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

dev.off()
