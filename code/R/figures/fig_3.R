
# Load dependencies
deps <- c('vegan', 'klaR', 'wesanderson', 'scatterplot3d');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Define input file names
cefoperazone_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/all_genes/cefoperazone_630.RNA_reads2cdf630.norm.annotated.txt'
clindamycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/all_genes/clindamycin_630.RNA_reads2cdf630.norm.annotated.txt'
streptomycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/all_genes/streptomycin_630.RNA_reads2cdf630.norm.annotated.txt'

# Load in data
cefoperazone <- read.delim(cefoperazone_file, sep='\t', header=FALSE, row.names=1)
colnames(cefoperazone) <- c('Cefoperazone', 'ko', 'gene', 'pathway')
clindamycin <- read.delim(clindamycin_file, sep='\t', header=FALSE, row.names=1)
colnames(clindamycin) <- c('Clindamycin', 'ko', 'gene', 'pathway')
streptomycin <- read.delim(streptomycin_file, sep='\t', header=FALSE, row.names=1)
colnames(streptomycin) <- c('Streptomycin', 'ko', 'gene', 'pathway')
rm(cefoperazone_file, clindamycin_file, streptomycin_file)

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
rm(cefoperazone, clindamycin, streptomycin)

# Rarefy mappings to be equal within sequencing type
sub_size <- round(min(colSums(combined_mapping[,1:3])) * 0.9) # 97930
cefoperazone_sub <- t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size))
clindamycin_sub <- t(rrarefy(combined_mapping$Clindamycin, sample=sub_size))
streptomycin_sub <- t(rrarefy(combined_mapping$Streptomycin, sample=sub_size))
for (index in 1:999) {
  cefoperazone_sub <- cbind(cefoperazone_sub, t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size)))
  clindamycin_sub <- cbind(clindamycin_sub, t(rrarefy(combined_mapping$Clindamycin, sample=sub_size)))
  streptomycin_sub <- cbind(streptomycin_sub, t(rrarefy(combined_mapping$Streptomycin, sample=sub_size)))
}
combined_mapping$Cefoperazone <- round(rowMeans(cefoperazone_sub), digits=0)
combined_mapping$Clindamycin <- round(rowMeans(clindamycin_sub), digits=0)
combined_mapping$Streptomycin <- round(rowMeans(streptomycin_sub), digits=0)
rm(index, sub_size, cefoperazone_sub, clindamycin_sub, streptomycin_sub)

# Eliminate genes with no transcripts mapping
combined_mapping <- combined_mapping[rowSums(combined_mapping[,1:3]) != 0.0, ] 

# Transform the data for comparison
combined_mapping[,1:3] <- log10(combined_mapping[,1:3] + 1)

#----------------------------------------------------------------------------------------------------------------------------#

# Subset by gene annotations and calculate relative abundances

# Amino sugar catabolism
mur <- rbind(subset(combined_mapping, grepl('mur.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('mur;', combined_mapping$gene))) # Muramidase
nag <- rbind(subset(combined_mapping, grepl('nag.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('nag;', combined_mapping$gene))) # Acetylglucosaminidase
acd <- rbind(subset(combined_mapping, grepl('acd.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('acd;', combined_mapping$gene))) # Peptidoglycan hydrolase 
ldt <- rbind(subset(combined_mapping, grepl('ldt.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('ldt;', combined_mapping$gene))) # l,d-transpeptidase 
gne <- rbind(subset(combined_mapping, grepl('gne.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('gne;', combined_mapping$gene))) # Glucosamine (UDP-N-Acetyl)-2-Epimerase/N-Acetylmannosamine Kinase
nan <- rbind(subset(combined_mapping, grepl('nan.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('nan;', combined_mapping$gene))) # Sialidase
glm <- rbind(subset(combined_mapping, grepl('glm.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('glm;', combined_mapping$gene))) # Glutamate mutase
amino_sugars <- rbind(mur, nag, acd, ldt, gne, nan, glm)
rm(mur, nag, acd, ldt, gne, nan, glm)
amino_sugars$grouping <- rep('Amino sugar catabolism', nrow(amino_sugars))
amino_sugars_relabund <- amino_sugars[,1:3] / rowSums(amino_sugars[,1:3])
gene_table <- amino_sugars
amino_sugars$ko <- NULL
amino_sugars$gene <- NULL
amino_sugars$pathway <- NULL

# Amino acid catabolism (including Stickland fermentation)
pep <- rbind(subset(combined_mapping, grepl('pep.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('pep;', combined_mapping$gene))) # general peptidases
prd <- rbind(subset(combined_mapping, grepl('prd.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('prd;', combined_mapping$gene))) # proline fermentation
grd <- rbind(subset(combined_mapping, grepl('grd.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('grd;', combined_mapping$gene))) # glycine fermentation
kamA <- subset(combined_mapping, grepl('kam.;', combined_mapping$gene)) # lycine fermentation
had <- subset(combined_mapping, grepl('had.;', combined_mapping$gene)) # Leucine fermentation
tdcB <- subset(combined_mapping, grepl('tdcB;', combined_mapping$gene)) # threonine dehydratase
fdh <- subset(combined_mapping, grepl('fdh.;', combined_mapping$gene)) # Formate dehydrogenase
panB  <- subset(combined_mapping, grepl('panB;', combined_mapping$gene)) # Butanoate formation
arg <- subset(combined_mapping, grepl('arg.;', combined_mapping$gene)) # arginine deaminase operon
sdaB <- subset(combined_mapping, grepl('sdaB;', combined_mapping$gene)) # L-serine dehydratase (serine catabolism)
stickland <- rbind(pep, prd, grd, had, tdcB, fdh, panB, arg, sdaB, kamA)
rm(pep, prd, grd, had, tdcB, fdh, panB, arg, sdaB, kamA)
stickland$grouping <- rep('Stickland reactions', nrow(stickland))
stickland_relabund <- stickland[,1:3] / rowSums(stickland[,1:3])
gene_table <- rbind(gene_table, stickland)
stickland$ko <- NULL
stickland$gene <- NULL
stickland$pathway <- NULL

# Monosaccharide catabolism
gap <- subset(combined_mapping, grepl('gap.;', combined_mapping$gene)) # Glyceraldehyde 3-phosphate dehydrogenase (Glycolysis)
gpmI <- subset(combined_mapping, grepl('gpmI;', combined_mapping$gene)) # Phosphoglyceromutase (Glycolysis)
pfk <- rbind(subset(combined_mapping, grepl('pfkA;', combined_mapping$gene)),
              subset(combined_mapping, grepl('phosphofructokinase', combined_mapping$gene)))# Phosphofructokinase (Glycolysis)
tpi <- subset(combined_mapping, grepl('tpi;', combined_mapping$gene)) # Triosephosphate isomerase (Glycolysis)
pyk <- subset(combined_mapping, grepl('pyk;', combined_mapping$gene)) # Pyruvate kinase (Glycolysis)
eno <- subset(combined_mapping, grepl('eno;', combined_mapping$gene)) # Enolase (Glycolysis)
pgm <- rbind(subset(combined_mapping, grepl('pgm;', combined_mapping$gene)),
             subset(combined_mapping, grepl('phosphoglycerate_mutase', combined_mapping$gene))) # Phosphoglycerate mutase (Glycolysis)
galactose <- rbind(subset(combined_mapping, grepl('galE;', combined_mapping$gene)), 
                   subset(combined_mapping, grepl('galactose-1-phosphate_uridylyltransferase', combined_mapping$gene))) # hexose
mannose <- rbind(subset(combined_mapping, grepl('mng.;', combined_mapping$gene)),
                 subset(combined_mapping, grepl('man.;', combined_mapping$gene)),
                 subset(combined_mapping, grepl('pmi;', combined_mapping$gene)))# hexose
tagatose <- subset(combined_mapping, grepl('tagatose', combined_mapping$gene)) # hexose 
fructose <- rbind(subset(combined_mapping, grepl('fbp;', combined_mapping$gene)),
                  subset(combined_mapping, grepl('fru;', combined_mapping$gene)),
                  subset(combined_mapping, grepl('fru.;', combined_mapping$gene)),
                  subset(combined_mapping, grepl('fru...;', combined_mapping$gene)),
                  subset(combined_mapping, grepl('fba;', combined_mapping$gene))) # hexose
xylose <- subset(combined_mapping, grepl('xylose', combined_mapping$gene)) # pentose
monosaccharides <- rbind(gap, gpmI, pfk, tpi, pyk, eno, pgm, galactose, mannose, tagatose, fructose, xylose)
rm(gap, gpmI, pfk, tpi, pyk, eno, pgm, galactose, mannose, tagatose, fructose, xylose)
monosaccharides$grouping <- rep('Monosaccharide catabolism', nrow(monosaccharides))
monosaccharides_relabund <- monosaccharides[,1:3] / rowSums(monosaccharides[,1:3])
gene_table <- rbind(gene_table, monosaccharides)
monosaccharides$ko <- NULL
monosaccharides$gene <- NULL
monosaccharides$pathway <- NULL

# Polysaccharide catabolism
sucrose <- subset(combined_mapping, grepl('scr.;', combined_mapping$gene))
maltose <- rbind(subset(combined_mapping, grepl('maltose-6\'-phosphate_glucosidase', combined_mapping$gene)),
                 subset(combined_mapping, grepl('maa;', combined_mapping$gene)),
                 subset(combined_mapping, grepl('map.;', combined_mapping$gene)),
                 subset(combined_mapping, grepl('maltose_O-acetyltransferase', combined_mapping$gene)))
tre <- subset(combined_mapping, grepl('tre.;', combined_mapping$gene)) # Trehalose utilization operon
glucosidase <- subset(combined_mapping, grepl('glucosidase', combined_mapping$gene))
cel <- rbind(subset(combined_mapping, grepl('celG;', combined_mapping$gene)),
             subset(combined_mapping, grepl('celC;', combined_mapping$gene)))
polysaccharides <- rbind(sucrose, maltose, tre, glucosidase, cel)
rm(sucrose, maltose, tre, glucosidase, cel)
polysaccharides$grouping <- rep('Polysaccharide catabolism', nrow(polysaccharides))
polysaccharides_relabund <- polysaccharides[,1:3] / rowSums(polysaccharides[,1:3])
gene_table <- rbind(gene_table, polysaccharides)
polysaccharides$ko <- NULL
polysaccharides$gene <- NULL
polysaccharides$pathway <- NULL

# PTS systems
PTS <- rbind(subset(combined_mapping, grepl('PTS_system', combined_mapping$gene)),
             subset(combined_mapping, grepl('pyridoxal_phosphate-dependent_transferase', combined_mapping$gene)))
PTS$grouping <- rep('PEP group translocators', nrow(PTS))
PTS_relabund <- PTS[,1:3] / rowSums(PTS[,1:3])
gene_table <- rbind(gene_table, PTS)
PTS$ko <- NULL
PTS$gene <- NULL
PTS$pathway <- NULL

# ABC transporters
ABC <- subset(combined_mapping, grepl('ABC_transporter_sugar', combined_mapping$gene))
ABC$grouping <- rep('ABC sugar transporters', nrow(ABC))
ABC[,1:3] <- ABC[,1:3] / rowSums(ABC[,1:3])
gene_table <- rbind(gene_table, ABC)
ABC$ko <- NULL
ABC$gene <- NULL
ABC$pathway <- NULL

# Sugar alcohols
srl <- subset(combined_mapping, grepl('srl.;', combined_mapping$gene)) # sorbitol utilization locus
srlE <- subset(combined_mapping, grepl('srlE.;', combined_mapping$gene)) # sorbitol import
mtl <- subset(combined_mapping, grepl('mtl.;', combined_mapping$gene)) # mannitol utilization locus
sugar_alcohols <- rbind(srl, srlE, mtl)
rm(srl, srlE, mtl)
sugar_alcohols$grouping <- rep('Sugar alcohol catabolism', nrow(sugar_alcohols))
sugar_alcohols_relabund <- sugar_alcohols[,1:3] / rowSums(sugar_alcohols[,1:3])
gene_table <- rbind(gene_table, sugar_alcohols)
sugar_alcohols$ko <- NULL
sugar_alcohols$gene <- NULL
sugar_alcohols$pathway <- NULL

# Fermentation genes
buk <- rbind(subset(combined_mapping, grepl('buk;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('buk.;', combined_mapping$gene))) # Butyrate kinase
ptb <- rbind(subset(combined_mapping, grepl('ptb;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('ptb.;', combined_mapping$gene))) # phosphate butyryltransferase
acetate <- subset(combined_mapping, grepl('acetate', combined_mapping$gene))
valerate <- subset(combined_mapping, grepl('valerate', combined_mapping$gene))
sucD <- subset(combined_mapping, grepl('sucD;', combined_mapping$gene)) # succinate-semialdehyde dehydrogenase
adh <- subset(combined_mapping, grepl('adh.;', combined_mapping$gene)) # Alcohol dehydrogenase
cat <- subset(combined_mapping, grepl('cat.;', combined_mapping$gene)) # Acetate to butyrate conversion
abfD <- subset(combined_mapping, grepl('abfD;', combined_mapping$gene)) # gamma-aminobutyrate metabolism dehydratase
hbd <- subset(combined_mapping, grepl('hbd;', combined_mapping$gene)) # 3-hydroxybutyryl-CoA dehydrogenase
fermentation <- rbind(buk, ptb, acetate, valerate, sucD, adh, cat, abfD, hbd)
rm(buk, ptb, acetate, valerate, sucD, adh, cat, abfD, hbd)
fermentation$grouping <- rep('Fermentation end steps', nrow(fermentation))
fermentation_relabund <- fermentation[,1:3] / rowSums(fermentation[,1:3])
gene_table <- rbind(gene_table, fermentation)
fermentation$ko <- NULL
fermentation$gene <- NULL
fermentation$pathway <- NULL

# Prep the data from and write it to a file
gene_table$KEGG_code <- rownames(gene_table)
table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S2.tsv'
write.table(gene_table, file=table_file, sep='\t', row.names=FALSE, quote=FALSE)
rm(table_file, gene_table)

# Calculate rleative abundance for the rest of the genes
combined_mapping[,1:3] <- combined_mapping[,1:3] / rowSums(combined_mapping[,1:3])

#----------------------------------------------------------------------------------------------------------------------------#

# Create supplementary 3D scatterplots



scatterplot3d(x=monosaccharides[,1], y=monosaccharides[,2], z=monosaccharides[,3], 
              xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), 
              highlight.3d=TRUE, col.grid="gray68", pch=19)








#----------------------------------------------------------------------------------------------------------------------------#

# Calculate Bray Curtis dissimilarities and significant differences in similarity matricies

# Combine correct groups for comparison
transport <- rbind(ABC, PTS)
transport_dis <- vegdist(x=transport[,1:3], method='bray', binary=FALSE, diag=FALSE, upper=FALSE)
anosim(transport_dis, grouping=transport$grouping, permutations=1000, distance='bray')
# R = 0.1582
# p-value = 0.00999
# adjusted p-value = 0.039960 *
rm(transport, transport_dis)

sugar_alcohols$association <- ifelse(sugar_alcohols$Streptomycin >= 0.5, 'Streptomycin', 'Cefoperazone')
sugar_alcohols_dis <- vegdist(x=sugar_alcohols[,1:3], method='bray', binary=FALSE, diag=FALSE, upper=FALSE)
anosim(sugar_alcohols_dis, grouping=sugar_alcohols$association, permutations=1000, distance='bray')
# R = 1
# p-value = 0.004995
# adjusted p-value = 0.024975 *
rm(sugar_alcohols_dis)
sugar_alcohols$association <- NULL

polysaccharides$association <- ifelse(polysaccharides$Clindamycin >= 0.5, 'Clindamycin', 'Cefoperazone')
polysaccharides_dis <- vegdist(x=polysaccharides[,1:3], method='bray', binary=FALSE, diag=FALSE, upper=FALSE)
anosim(polysaccharides_dis, grouping=polysaccharides$association, permutations=1000, distance='bray')
# R = 0.1494 
# p-value = 0.12488 
# adjusted p-value = 0.249760 n.s.
rm(polysaccharides_dis)
polysaccharides$association <- NULL

carb_type <- rbind(monosaccharides, fermentation)
carb_type_dis <- vegdist(x=carb_type[,1:3], method='bray', binary=FALSE, diag=FALSE, upper=FALSE)
anosim(carb_type_dis, grouping=carb_type$grouping, permutations=1000, distance='bray')
# R = 0.1132 
# p-value = 0.025974
# adjusted p-value = 0.077922 n.s.
rm(carb_type)

carb_type <- rbind(polysaccharides, fermentation)
carb_type_dis <- vegdist(x=carb_type[,1:3], method='bray', binary=FALSE, diag=FALSE, upper=FALSE)
anosim(carb_type_dis, grouping=carb_type$grouping, permutations=1000, distance='bray')
# R = 0.3028 
# p-value = 0.000999
# adjusted p-value = 0.005994 **
rm(carb_type, carb_type_dis)

ferm_type <- rbind(stickland, fermentation)
ferm_type_dis <- vegdist(x=ferm_type[,1:3], method='bray', binary=FALSE, diag=FALSE, upper=FALSE)
anosim(ferm_type_dis, grouping=ferm_type$grouping, permutations=1000, distance='bray')
# R = 0.02785 
# p-value = 0.2008
# adjusted p-value = 0.249760 n.s.
rm(ferm_type, ferm_type_dis)



stickland$association <- ifelse(stickland$Clindamycin >= 0.5, 'Clindamycin', 'Cefoperazone')
stickland_dis <- vegdist(x=stickland[,1:3], method='bray', binary=FALSE, diag=FALSE, upper=FALSE)
anosim(stickland_dis, grouping=stickland$association, permutations=1000, distance='bray')
# R =  
# p-value = 
# adjusted p-value = 
rm(stickland_dis)
stickland$association <- NULL
stickland$association <- ifelse(stickland$Clindamycin >= 0.5, 'Clindamycin', 'Streptomycin')
stickland_dis <- vegdist(x=stickland[,1:3], method='bray', binary=FALSE, diag=FALSE, upper=FALSE)
anosim(stickland_dis, grouping=stickland$association, permutations=1000, distance='bray')
# R =  
# p-value = 
# adjusted p-value = 
rm(stickland_dis)
stickland$association <- NULL
stickland$association <- ifelse(stickland$Streptomycin >= 0.5, 'Streptomycin', 'Cefoperazone')
stickland_dis <- vegdist(x=stickland[,1:3], method='bray', binary=FALSE, diag=FALSE, upper=FALSE)
anosim(stickland_dis, grouping=stickland$association, permutations=1000, distance='bray')
# R =  
# p-value = 
# adjusted p-value = 
rm(stickland_dis)
stickland$association <- NULL




# Delete abundances
rm(amino_sugars, stickland, monosaccharides, polysaccharides, ABC, PTS, sugar_alcohols, fermentation)

# Correct p-values
p_values <- c(0.00999, 0.004995, 0.12488, 0.025974, 0.000999, 0.2008,)
corrected_p_values <- p.adjust(p_values, method='holm')
rm(p_values)

# Create supplementary table for corrected p-values
p_table <- cbind(c('ABC_transporters','Sugar_alcohol_Strep','Polysaccharide_Clinda','Monosaccharide','Polysaccharide','Amino_acid'), 
                 c('PTS_transporters','Sugar_alcohol_Cef','Polysaccharide_Cef','SCFA_production','SCFA_production','SCFA_production'), 
                 corrected_p_values)
colnames(p_table) <- c('Group_1', 'Group_2', 'Corrected_p_value_HOLM')
rm(corrected_p_values)

# Write table
table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/pvalues.tsv'
write.table(p_table, file=table_file, sep='\t', row.names=FALSE, quote=FALSE)
rm(table_file, p_table)

#-------------------------------------------------------------------------------------------------------------------------#

# Define plot details
rainbow <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
fox <- wes_palette("FantasticFox")
tick_labels <- c('10%','','30%','','50%','','70%','','90%')
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_3.pdf'

# Open a PDF
pdf(file=plot_file, width=9, height=9)

# Create layout for multi-plot
layout(mat=matrix(c(1,1,1,2, 
                    1,1,1,3, 
                    1,1,1,4, 
                    5,6,7,8), nrow=4, ncol=4, byrow=TRUE))

# Generate raw plot
par(mar=c(0,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('','',''), grid=seq(0.1,0.9,by=0.1), cex=0)

# Center triangle
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')

# 50% lines
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))

# Add points
tripoints(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], cex=0.8)

# Axis labels
text(x=-0.35, y=-0.41, labels='Cefoperazone', cex=1.4)
text(x=-0.21, y=0.48, labels='Clindamycin', cex=1.4, srt=60)
text(x=0.55, y=-0.12, labels='Streptomycin', cex=1.4, srt=-60)

# Left axis - Clindmycin
lines(x=c(-0.52,-0.54), y=c(-0.233,-0.233))
lines(x=c(-0.462,-0.492), y=c(-0.133,-0.133))
lines(x=c(-0.404,-0.424), y=c(-0.033,-0.033))
lines(x=c(-0.346,-0.366), y=c(0.067,0.067))
lines(x=c(-0.289,-0.309), y=c(0.167,0.167))
lines(x=c(-0.23,-0.25), y=c(0.267,0.267))
lines(x=c(-0.173,-0.193), y=c(0.367,0.367))
lines(x=c(-0.115,-0.135), y=c(0.467,0.467))
lines(x=c(-0.057,-0.077), y=c(0.567,0.567))
text(x=c(-0.57,-0.522,-0.454,-0.396,-0.339,-0.28,-0.223,-0.165,-0.107), 
     y=c(-0.233,-0.133,-0.033,0.067,0.167,0.267,0.367,0.467,0.567), 
     labels=tick_labels, cex=0.9)

# Right axis - Streptomycin
lines(x=c(0.52,0.54), y=c(-0.233,-0.233))
lines(x=c(0.462,0.482), y=c(-0.133,-0.133))
lines(x=c(0.404,0.424), y=c(-0.033,-0.033))
lines(x=c(0.346,0.366), y=c(0.067,0.067))
lines(x=c(0.289,0.309), y=c(0.167,0.167))
lines(x=c(0.23,0.25), y=c(0.267,0.267))
lines(x=c(0.173,0.193), y=c(0.367,0.367))
lines(x=c(0.115,0.135), y=c(0.467,0.467))
lines(x=c(0.057,0.077), y=c(0.567,0.567))
text(x=c(0.57,0.522,0.454,0.396,0.339,0.28,0.223,0.165,0.107), 
     y=c(-0.233,-0.133,-0.033,0.067,0.167,0.267,0.367,0.467,0.567), 
     labels=rev(tick_labels), cex=0.9)

# Bottom axis - Cefoperzone
lines(x=c(-0.462,-0.462), y=c(-0.333,-0.353))
lines(x=c(-0.346,-0.346), y=c(-0.333,-0.353))
lines(x=c(-0.231,-0.231), y=c(-0.333,-0.353))
lines(x=c(-0.115,-0.115), y=c(-0.333,-0.353))
lines(x=c(0,0), y=c(-0.333,-0.353))
lines(x=c(0.116,0.116), y=c(-0.333,-0.353))
lines(x=c(0.232,0.232), y=c(-0.333,-0.353))
lines(x=c(0.347,0.347), y=c(-0.333,-0.353))
lines(x=c(0.463,0.463), y=c(-0.333,-0.353))
text(x=c(-0.462,-0.346,-0.231,-0.115,0,0.116,0.232,0.347,0.463), 
     y=c(-0.373,-0.373,-0.373,-0.373,-0.373,-0.373,-0.373,-0.373,-0.373), 
     labels=rev(tick_labels), cex=0.9)

# Color points by substrate
tripoints(x=PTS[,1], y=PTS[,2], z=PTS[,3], pch=21, cex=2.3, bg=fox[3])
tripoints(x=ABC[,1], y=ABC[,2], z=ABC[,3], pch=21, cex=2.3, bg=rainbow[7])
tripoints(x=monosaccharides[,1], y=monosaccharides[,2], z=monosaccharides[,3], pch=21, cex=2.3, bg=fox[1])
tripoints(x=stickland[,1], y=stickland[,2], z=stickland[,3], pch=21, cex=2.3, bg=fox[2])
tripoints(x=sugar_alcohols[,1], y=sugar_alcohols[,2], z=sugar_alcohols[,3], pch=21, cex=2.3, bg='darkorchid3')
tripoints(x=fermentation[,1], y=fermentation[,2], z=fermentation[,3], pch=21, cex=2.3, bg=fox[5])
tripoints(x=polysaccharides[,1], y=polysaccharides[,2], z=polysaccharides[,3], pch=21, cex=2.3, bg='blue3')

# Add the legend
legend(x=0.26, y=0.63, legend=c('Monosaccharide catabolism', 'Polysaccharide catabolism', 'Sugar alcohol catabolism', 'Amino acid catabolism', 'SCFA production', 'PEP group translocators', 'ABC sugar transporters', 'All other genes'), 
    ncol=1, pch=21, cex=1.4, pt.cex=c(2.5,2.5,2.5,2.5,2.5,2.5,2.5,1.4), col='black', pt.bg=c(fox[1],'blue3','darkorchid3',fox[2],fox[5],fox[3],rainbow[7], 'white'), bty='n')
# Add figure label
legend(x=-0.6, y=0.6, legend='A', cex=2, bty='n')

# monosaccharides alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=seq(0.1,0.9,by=0.1), cex=0.3, col='gray75')
legend('topleft', legend='B', cex=2, bty='n')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=monosaccharides[,1], y=monosaccharides[,2], z=monosaccharides[,3], pch=21, cex=2, bg=fox[1])
text(x=0, y=-0.48, labels='Monosaccharide catabolism', cex=1.3)

# polysaccharides alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=seq(0.1,0.9,by=0.1), cex=0.3, col='gray75')
legend('topleft', legend='C', cex=2, bty='n')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=polysaccharides[,1], y=polysaccharides[,2], z=polysaccharides[,3], pch=21, cex=2, bg='blue3')
text(x=0, y=-0.48, labels='Polysaccharide catabolism', cex=1.3)

# sugar_alcohols alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=seq(0.1,0.9,by=0.1), cex=0.3, col='gray75')
legend('topleft', legend='D', cex=2, bty='n')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=sugar_alcohols[,1], y=sugar_alcohols[,2], z=sugar_alcohols[,3], pch=21, cex=2, bg='darkorchid3')
text(x=0, y=-0.48, labels='Sugar alcohol catabolism', cex=1.3)

# stickland alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=seq(0.1,0.9,by=0.1), cex=0.3, col='gray75')
legend('topleft', legend='E', cex=2, bty='n')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=stickland[,1], y=stickland[,2], z=stickland[,3], pch=21, cex=2, bg=fox[2])
text(x=0, y=-0.48, labels='Amino acid catabolism', cex=1.3)

# fermentation
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=seq(0.1,0.9,by=0.1), cex=0.3, col='gray75')
legend('topleft', legend='F', cex=2, bty='n')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=fermentation[,1], y=fermentation[,2], z=fermentation[,3], pch=21, cex=2, bg=fox[5])
text(x=0, y=-0.48, labels='SCFA production', cex=1.3)

# PTS alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=seq(0.1,0.9,by=0.1), cex=0.3, col='gray75')
legend('topleft', legend='G', cex=2, bty='n')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=PTS[,1], y=PTS[,2], z=PTS[,3], pch=21, cex=2, bg=fox[3])
text(x=0, y=-0.48, labels='PEP group translocators', cex=1.3)

# ABC alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=seq(0.1,0.9,by=0.1), cex=0.3, col='gray75')
legend('topleft', legend='H', cex=2, bty='n')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=ABC[,1], y=ABC[,2], z=ABC[,3], pch=21, cex=2, bg=rainbow[7])
text(x=0, y=-0.48, labels='ABC sugar transporters', cex=1.3)

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
rm(ABC, amino_sugars, combined_mapping, fermentation, monosaccharides, 
   polysaccharides, PTS, stickland, sugar_alcohols, plot_file, 
   rainbow, tick_labels, fox)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

