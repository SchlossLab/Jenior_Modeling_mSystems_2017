
# Load packages
source("http://bioconductor.org/biocLite.R")
biocLite(c("Rgraphviz", "png", "KEGGgraph", "org.Hs.eg.db"))

install.packages("pathview", repos="http://R-Forge.R-project.org")
library(pathview)

#---------------------------------------------------------------------------------------------------#

# Set working directory where R will generate plots
setwd('~/Desktop/Matt/pathview_test/')

#---------------------------------------------------------------------------------------------------#

# Define expression file variables
transcripts_file <- 'test.lst'
inputs_file <- 'cefoperazone_630.input_score.txt'


#---------------------------------------------------------------------------------------------------#

# Read in data
transcripts <- read.delim(transcripts_file, sep='\t', header=FALSE, row.names=2)
inputs <- read.delim(inputs_file, sep='\t', header=FALSE, row.names=1)

#---------------------------------------------------------------------------------------------------#

# Format expression data
expression <- as.vector(as.numeric(transcripts[,1]))
names(expression) <- rownames(transcripts)
expression <- log2(expression)

compounds <- inputs
compounds <- as.vector(as.numeric(compounds[,2]))
names(compounds) <- rownames(inputs)
compounds <- log2(compounds)

#---------------------------------------------------------------------------------------------------#

# Generate plot - cpd.data = compounds, cpd.idtype = "kegg", 

# Carbohydrate metabolism
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00010", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # Glycolysis
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00020", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # tca_cycle
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00030", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # pentose_phosphate
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00040", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # pentose_to_glucuronate
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00051", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # fructose_mannose
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00052", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # galactose
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00053", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # ascorbate_aldarate
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00500", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # starch_sucrose 
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00520", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # amino_sugar
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00620", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # pyruvate
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00630", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # glyoxylate_dicarboxylate 
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00640", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # propanoate
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00650", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # butanoate
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00660", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # c5_dibasic_acid

# Energy metabolism
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00190", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # oxidative_phosphorylation
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00680", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # methane
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00910", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # nitrogen
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00920", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # sulfur

# Lipid metabolism
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00061", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # fatty_acid_biosynthesis
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00071", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # fatty_acid_degradation
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00072", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # ketone
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00561", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # glycerolipid
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00564", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # glycerophospholipid 
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "01040", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # unsaturated_fatty_acids

# Nucleotide metabolism
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00230", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # purine
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00240", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # pyrimidine

# Amino acid metabolism
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00250", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # alanine_aspartate_glutamate
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00260", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # glycine_serine_threonine
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00270", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # cysteine_methionine
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00280", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # valine_leucine_isoleucine1
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00290", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # valine_leucine_isoleucine2 
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00300", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # lysine1
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00310", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # lysine2
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00330", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # arginine_proline
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00340", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # histidine
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00350", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # tyrosine
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00360", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # phenylalanine
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00380", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # tryptophan
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00400", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # phenylalanine_tyrosine_tryptophan

# Metabolism of other amino acids
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00410", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # beta_alanine
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00430", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # taurine_hypotaurine
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00440", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # phosphonate_phosphinate
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00450", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # selenocompound 
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00460", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # cyanoamino_acid 
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00471", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # d_glutamine_d_glutamate
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00473", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # d_alanine
pathview(gene.data = expression, gene.idtype = "KEGG", pathway.id = "00480", cpd.data = compounds, cpd.idtype = "kegg", species = "cdf", kegg.dir = ".", file.type = 'png', same.layer=TRUE) # glutathione

#---------------------------------------------------------------------------------------------------#























































#---------------------------------------------------------------------------------------------------#
# Download KEGG pathway information without mapping

# Carbohydrate metabolism
glycolysis <- download.kegg(pathway.id = "00010", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
tca_cycle <- download.kegg(pathway.id = "00020", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
pentose_phosphate <- download.kegg(pathway.id = "00030", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
pentose_to_glucuronate <- download.kegg(pathway.id = "00040", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
fructose_mannose<- download.kegg(pathway.id = "00051", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
galactose <- download.kegg(pathway.id = "00052", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
ascorbate_aldarate <- download.kegg(pathway.id = "00053", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
starch_sucrose <- download.kegg(pathway.id = "00500", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
amino_sugar <- download.kegg(pathway.id = "00520", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
pyruvate <- download.kegg(pathway.id = "00620", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
glyoxylate_dicarboxylate <- download.kegg(pathway.id = "00630", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
propanoate <- download.kegg(pathway.id = "00640", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
butanoate <- download.kegg(pathway.id = "00650", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
c5_dibasic_acid <- download.kegg(pathway.id = "00660", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))

# Energy metabolism
oxidative_phosphorylation <- download.kegg(pathway.id = "00190", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
methane <- download.kegg(pathway.id = "00680", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
nitrogen <- download.kegg(pathway.id = "00910", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
sulfur <- download.kegg(pathway.id = "00920", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))

# Lipid metabolism
fatty_acid_biosynthesis <- download.kegg(pathway.id = "00061", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
fatty_acid_degradation <- download.kegg(pathway.id = "00071", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
ketone <- download.kegg(pathway.id = "00072", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
glycerolipid <- download.kegg(pathway.id = "00561", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
glycerophospholipid <- download.kegg(pathway.id = "00564", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
unsaturated_fatty_acids <- download.kegg(pathway.id = "01040", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))

# Nucleotide metabolism
purine <- download.kegg(pathway.id = "00230", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
pyrimidine <- download.kegg(pathway.id = "00240", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))

# Amino acid metabolism
alanine_aspartate_glutamate <- download.kegg(pathway.id = "00250", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
glycine_serine_threonine <- download.kegg(pathway.id = "00260", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
cysteine_methionine <- download.kegg(pathway.id = "00270", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
valine_leucine_isoleucine1 <- download.kegg(pathway.id = "00280", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
valine_leucine_isoleucine2 <- download.kegg(pathway.id = "00290", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
lysine1 <- download.kegg(pathway.id = "00300", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
lysine2 <- download.kegg(pathway.id = "00310", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
arginine_proline <- download.kegg(pathway.id = "00330", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
histidine <- download.kegg(pathway.id = "00340", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
tyrosine <- download.kegg(pathway.id = "00350", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
phenylalanine <- download.kegg(pathway.id = "00360", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
tryptophan <- download.kegg(pathway.id = "00380", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
phenylalanine_tyrosine_tryptophan <- download.kegg(pathway.id = "00400", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))

# Metabolism of other amino acids
beta_alanine <- download.kegg(pathway.id = "00410", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
taurine_hypotaurine <- download.kegg(pathway.id = "00430", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
phosphonate_phosphinate <- download.kegg(pathway.id = "00440", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
selenocompound <- download.kegg(pathway.id = "00450", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
cyanoamino_acid <- download.kegg(pathway.id = "00460", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
d_glutamine_d_glutamate <- download.kegg(pathway.id = "00471", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
d_alanine <- download.kegg(pathway.id = "00473", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))
glutathione <- download.kegg(pathway.id = "00480", species = "cdf", kegg.dir = ".", file.type=c("xml", "png"))

