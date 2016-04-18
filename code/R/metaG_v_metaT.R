

# Define file variables
metagenome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metagenome/Cefoperazone.DNA_reads2metaG.all.pool.norm.remove.annotated.txt'
metatranscriptome1_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/remove_Cdifficile/cefoperazone_630.RNA_reads2metaG.cdf.all.pool.norm.remove.annotated.txt'
metatranscriptome2_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/remove_Cdifficile/cefoperazone_mock.RNA_reads2metaG.cdf.all.pool.norm.remove.annotated.txt'
plot_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/cefoperazone.metaTvsmetaG.pdf'

# Load in data
metagenome <- read.delim(metagenome_file, sep='\t', header=TRUE, row.names=4)
metatranscriptome1 <- read.delim(metatranscriptome1_file, sep='\t', header=TRUE, row.names=4)
metatranscriptome2 <- read.delim(metatranscriptome2_file, sep='\t', header=TRUE, row.names=4)

# Format data for merging
metagenome$gene_annotation <- NULL
metagenome$evalue <- NULL
metagenome$KEGG_ortholog <- NULL
metagenome$pathway_annotation <- NULL
metatranscriptome1$gene_annotation <- NULL
metatranscriptome1$evalue <- NULL
metatranscriptome1$KEGG_ortholog <- NULL
metatranscriptome1$pathway_annotation <- NULL

# Merge tables
combined_mapping <- merge(metagenome, metatranscriptome1, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
combined_mapping <- merge(combined_mapping, metatranscriptome2, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
colnames(combined_mapping) <- c('metagenome', 'infected_metatranscriptome', 'mock_metatranscriptome', 'gene_annotation', 'evalue', 'KEGG_ortholog', 'pathway_annotation')
combined_mapping$gene_annotation <- NULL
combined_mapping$evalue <- NULL
combined_mapping$KEGG_ortholog <- NULL



# Grep specific strings out of pathway annotation field

# Define points for coloring later
abc_transporters_630 <- subset(combined_data_630, metadata1 == 'ABC transporters', select = c(G.relabund, T.relabund))
starch_sucrose_630 <- subset(combined_data_630, metadata1 == 'Starch and sucrose metabolism', select = c(G.relabund, T.relabund))
carbon_630 <- subset(combined_data_630, metadata1 == 'Carbon metabolism', select = c(G.relabund, T.relabund))
aromatic_compounds_630 <- subset(combined_data_630, metadata1 == 'Degradation of aromatic compounds', select = c(G.relabund, T.relabund))
amino_sugar_630 <- subset(combined_data_630, metadata1 == 'Amino sugar and nucleotide sugar metabolism', select = c(G.relabund, T.relabund))
secondary_metabolites_630 <- subset(combined_data_630, metadata1 == 'Biosynthesis of secondary metabolites', select = c(G.relabund, T.relabund))
diverse_environments_630 <- subset(combined_data_630, metadata1 == 'Microbial metabolism in diverse environments', select = c(G.relabund, T.relabund))
two_component_630 <- subset(combined_data_630, metadata1 == 'Two-component system', select = c(G.relabund, T.relabund))
arginine_proline_630 <- subset(combined_data_630, metadata1 == 'Arginine and proline metabolism', select = c(G.relabund, T.relabund))
pyruvate_630 <- subset(combined_data_630, metadata1 == 'Pyruvate metabolism', select = c(G.relabund, T.relabund))
fatty_acid_630 <- subset(combined_data_630, metadata1 == 'Fatty acid metabolism', select = c(G.relabund, T.relabund))
glycerophospholipid_630 <- subset(combined_data_630, metadata1 == 'Glycerophospholipid metabolism', select = c(G.relabund, T.relabund))
peptidoglycan_630 <- subset(combined_data_630, metadata1 == 'Peptidoglycan biosynthesis', select = c(G.relabund, T.relabund))

abc_transporters_mock <- subset(combined_data_mock, metadata1 == 'ABC transporters', select = c(G.relabund, T.relabund))
starch_sucrose_mock <- subset(combined_data_mock, metadata1 == 'Starch and sucrose metabolism', select = c(G.relabund, T.relabund))
carbon_mock <- subset(combined_data_mock, metadata1 == 'Carbon metabolism', select = c(G.relabund, T.relabund))
aromatic_compounds_mock <- subset(combined_data_mock, metadata1 == 'Degradation of aromatic compounds', select = c(G.relabund, T.relabund))
amino_sugar_mock <- subset(combined_data_mock, metadata1 == 'Amino sugar and nucleotide sugar metabolism', select = c(G.relabund, T.relabund))
secondary_metabolites_mock <- subset(combined_data_mock, metadata1 == 'Biosynthesis of secondary metabolites', select = c(G.relabund, T.relabund))
diverse_environments_mock <- subset(combined_data_mock, metadata1 == 'Microbial metabolism in diverse environments', select = c(G.relabund, T.relabund))
two_component_mock <- subset(combined_data_mock, metadata1 == 'Two-component system', select = c(G.relabund, T.relabund))
arginine_proline_mock <- subset(combined_data_mock, metadata1 == 'Arginine and proline metabolism', select = c(G.relabund, T.relabund))
pyruvate_mock <- subset(combined_data_mock, metadata1 == 'Pyruvate metabolism', select = c(G.relabund, T.relabund))
fatty_acid_mock <- subset(combined_data_mock, metadata1 == 'Fatty acid metabolism', select = c(G.relabund, T.relabund))
glycerophospholipid_mock <- subset(combined_data_mock, metadata1 == 'Glycerophospholipid metabolism', select = c(G.relabund, T.relabund))
peptidoglycan_mock <- subset(combined_data_mock, metadata1 == 'Peptidoglycan biosynthesis', select = c(G.relabund, T.relabund))

# Plot the most interesting pathways
pdf(file=plot_file, width=25, height=12)
layout(matrix(1:12, nrow=2, byrow=T))
axis_labels <- parse(text=paste(rep(10,7), '^', seq(-6,0,1), sep=''))
par(mar=c(5, 4, 4, 4) + 0.1)
# 630 Infected
plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', main='ABC transporters', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(abc_transporters_630), cex=1.5, pch=21, bg='steelblue3', col='black')
plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', main='Starch and sucrose metabolism', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(starch_sucrose_630), cex=1.5, pch=21, bg='yellowgreen', col='black')
plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', main='Carbon metabolism', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(carbon_630), cex=1.5, pch=21, bg='yellow1', col='black')
plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', main='Amino sugar and nucleotide sugar metabolism', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(amino_sugar_630), cex=1.5, pch=21, bg='orange', col='black')
plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', main='Biosynthesis of secondary metabolites', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(secondary_metabolites_630), cex=1.5, pch=21, bg='blueviolet', col='black')
plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', main='Microbial metabolism in diverse environments', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(diverse_environments_630), cex=1.5, pch=21, bg='blue4', col='black')
mtext('630 Infected', side=4, cex=1.5, padj=1)
#plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(two_component_630), cex=1.5, pch=21, bg='cyan3', col='black')
#plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(arginine_proline_630), cex=1.5, pch=21, bg='darkgoldenrod1', col='black')
#plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(pyruvate_630), cex=1.5, pch=21, bg='orangered3', col='black')
#plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(fatty_acid_630), cex=1.5, pch=21, bg='firebrick2', col='black')
#plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(glycerophospholipid_630), cex=1.5, pch=21, bg='forestgreen', col='black')
#plot(log10(combined_data_630), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(peptidoglycan_630), cex=1.5, pch=21, bg='chocolate4', col='black')

# Mock
plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(abc_transporters_mock), cex=1.5, pch=21, bg='steelblue3', col='black')
plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(starch_sucrose_mock), cex=1.5, pch=21, bg='yellowgreen', col='black')
plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(carbon_mock), cex=1.5, pch=21, bg='yellow1', col='black')
plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(amino_sugar_mock), cex=1.5, pch=21, bg='orange', col='black')
plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(secondary_metabolites_mock), cex=1.5, pch=21, bg='blueviolet', col='black')
plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25', xaxt='n', yaxt='n')
segments(-7, -7, 1, 1, lwd=2)
axis(side=1, at=c(-6:0), axis_labels, tick=TRUE)
axis(side=2, at=c(-6:0), axis_labels, tick=TRUE, las=1)
points(log10(diverse_environments_mock), cex=1.5, pch=21, bg='blue4', col='black')
mtext('Mock Treated', side=4, cex=1.5, padj=1)
#plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(two_component_mock), cex=1.5, pch=21, bg='cyan3', col='black')
#plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(arginine_proline_mock), cex=1.5, pch=21, bg='darkgoldenrod1', col='black')
#plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(pyruvate_mock), cex=1.5, pch=21, bg='orangered3', col='black')
#plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(fatty_acid_mock), cex=1.5, pch=21, bg='firebrick2', col='black')
#plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(glycerophospholipid_mock), cex=1.5, pch=21, bg='forestgreen', col='black')
#plot(log10(combined_data_mock), xlim=c(-6,0), ylim=c(-6,0), pch=20, col='gray25')
#segments(-7, -7, 1, 1, lwd=2)
#points(log10(peptidoglycan_mock), cex=1.5, pch=21, bg='chocolate4', col='black')
dev.off()







