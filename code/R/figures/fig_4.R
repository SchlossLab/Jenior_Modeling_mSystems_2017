
# Load dependencies
deps <- c('vegan', 'klaR', 'wesanderson', 'scatterplot3d', 'scales');
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

# Calculate coverage
sum(combined_mapping$Cefoperazone) / length(combined_mapping$Cefoperazone) # ~24x

# Eliminate genes with no transcripts mapping
combined_mapping <- combined_mapping[rowSums(combined_mapping[,1:3]) > 5, ] 

#----------------------------------------------------------------------------------------------------------------------------#

# Subset by gene annotations and calculate relative abundances

# Amino sugar catabolism
mur <- rbind(subset(combined_mapping, grepl('murA;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('murB;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('murG;', combined_mapping$gene))) # Muramidase
nag <- rbind(subset(combined_mapping, grepl('nag.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('nag;', combined_mapping$gene))) # Acetylglucosaminidase
acd <- rbind(subset(combined_mapping, grepl('acd.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('acd;', combined_mapping$gene))) # Peptidoglycan hydrolase 
ldt <- rbind(subset(combined_mapping, grepl('ldt.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('ldt;', combined_mapping$gene))) # l,d-transpeptidase 
gne <- rbind(subset(combined_mapping, grepl('GNE;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('gne;', combined_mapping$gene))) # Glucosamine (UDP-N-Acetyl)-2-Epimerase/N-Acetylmannosamine Kinase
nan <- rbind(subset(combined_mapping, grepl('nan.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('nan;', combined_mapping$gene))) # Sialidase
glm <- rbind(subset(combined_mapping, grepl('glm.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('glm;', combined_mapping$gene))) # Glutamate mutase
glucosamine <- rbind(subset(combined_mapping, grepl('^N-acetylglucosamine-6-phosphate_deacetylase', combined_mapping$gene)))
amino_sugars <- rbind(mur, nag, acd, ldt, gne, nan, glm, glucosamine)
rm(mur, nag, acd, ldt, gne, nan, glm, glucosamine)
amino_sugars$grouping <- rep('Amino sugar catabolism', nrow(amino_sugars))
amino_sugars_relabund <- amino_sugars[,1:3] / rowSums(amino_sugars[,1:3])
gene_table <- amino_sugars
amino_sugars$ko <- NULL
amino_sugars$gene <- NULL
amino_sugars$pathway <- NULL
amino_sugars$grouping <- NULL
amino_sugars[,1:3] <- log10(amino_sugars[,1:3] + 1)

# Amino acid catabolism (including Stickland fermentation)
pep <- rbind(subset(combined_mapping, grepl('pep.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('pep;', combined_mapping$gene))) # general peptidases
prd <- rbind(subset(combined_mapping, grepl('prd.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('prd;', combined_mapping$gene))) # proline fermentation
grd <- rbind(subset(combined_mapping, grepl('grd.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('grd;', combined_mapping$gene))) # glycine fermentation
kamA <- subset(combined_mapping, grepl('kam.;', combined_mapping$gene)) # lycine fermentation
map <- rbind(subset(combined_mapping, grepl('map1;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('map2;', combined_mapping$gene))) # methionine-directed peptidases
had <- subset(combined_mapping, grepl('had.;', combined_mapping$gene)) # Leucine fermentation
tdcB <- subset(combined_mapping, grepl('tdcB;', combined_mapping$gene)) # threonine dehydratase
fdh <- subset(combined_mapping, grepl('fdh.;', combined_mapping$gene)) # Formate dehydrogenase
panB  <- subset(combined_mapping, grepl('panB;', combined_mapping$gene)) # Butanoate formation
arg <- subset(combined_mapping, grepl('arg.;', combined_mapping$gene)) # arginine deaminase operon
sdaB <- subset(combined_mapping, grepl('sdaB;', combined_mapping$gene)) # L-serine dehydratase (serine catabolism)
stickland <- rbind(pep, prd, grd, had, tdcB, fdh, panB, arg, sdaB, kamA, map)
rm(pep, prd, grd, had, tdcB, fdh, panB, arg, sdaB, kamA, map)
stickland$grouping <- rep('Stickland reactions', nrow(stickland))
stickland_relabund <- stickland[,1:3] / rowSums(stickland[,1:3])
gene_table <- rbind(gene_table, stickland)
stickland$ko <- NULL
stickland$gene <- NULL
stickland$pathway <- NULL
stickland$grouping <- NULL
stickland[,1:3] <- log10(stickland[,1:3] + 1)

# Monosaccharide catabolism (glycolysis)
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
mannose <- rbind(subset(combined_mapping, grepl('man.;', combined_mapping$gene)),
                 subset(combined_mapping, grepl('pmi;', combined_mapping$gene)))# hexose
tagatose <- subset(combined_mapping, grepl('tagatose', combined_mapping$gene)) # hexose 
fructose <- rbind(subset(combined_mapping, grepl('fbp;', combined_mapping$gene)),
                  subset(combined_mapping, grepl('fru;', combined_mapping$gene)),
                  subset(combined_mapping, grepl('fru.;', combined_mapping$gene)),
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
monosaccharides$grouping <- NULL
monosaccharides[,1:3] <- log10(monosaccharides[,1:3] + 1)

# Polysaccharide catabolism
sucrose <- subset(combined_mapping, grepl('scr.;', combined_mapping$gene))
maltose <- rbind(subset(combined_mapping, grepl('maltose-6\'-phosphate_glucosidase', combined_mapping$gene)),
                 subset(combined_mapping, grepl('maa;', combined_mapping$gene)),
                 subset(combined_mapping, grepl('mapA;', combined_mapping$gene)),
                 subset(combined_mapping, grepl('maltose_O-acetyltransferase', combined_mapping$gene)))
tre <- subset(combined_mapping, grepl('tre.;', combined_mapping$gene)) # Trehalose utilization operon
glucosidase <- subset(combined_mapping, grepl('glucosidase', combined_mapping$gene))
cel <- rbind(subset(combined_mapping, grepl('celG;', combined_mapping$gene)))
polysaccharides <- rbind(sucrose, maltose, tre, glucosidase, cel)
rm(sucrose, maltose, tre, glucosidase, cel)
polysaccharides$grouping <- rep('Polysaccharide catabolism', nrow(polysaccharides))
polysaccharides_relabund <- polysaccharides[,1:3] / rowSums(polysaccharides[,1:3])
gene_table <- rbind(gene_table, polysaccharides)
polysaccharides$ko <- NULL
polysaccharides$gene <- NULL
polysaccharides$pathway <- NULL
polysaccharides$grouping <- NULL
polysaccharides[,1:3] <- log10(polysaccharides[,1:3] + 1)

# PTS systems
PTS <- rbind(subset(combined_mapping, grepl('PTS_system', combined_mapping$gene)),
             subset(combined_mapping, grepl('pyridoxal_phosphate-dependent_transferase', combined_mapping$gene)))
PTS$grouping <- rep('PEP group translocators', nrow(PTS))
PTS_relabund <- PTS[,1:3] / rowSums(PTS[,1:3])
gene_table <- rbind(gene_table, PTS)
PTS$ko <- NULL
PTS$gene <- NULL
PTS$pathway <- NULL
PTS$grouping <- NULL
PTS[,1:3] <- log10(PTS[,1:3] + 1)

# ABC transporters
ABC <- subset(combined_mapping, grepl('ABC_transporter_sugar', combined_mapping$gene))
ABC$grouping <- rep('ABC sugar transporters', nrow(ABC))
ABC_relabund <- ABC[,1:3] / rowSums(ABC[,1:3])
gene_table <- rbind(gene_table, ABC)
ABC$ko <- NULL
ABC$gene <- NULL
ABC$pathway <- NULL
ABC$grouping <- NULL
ABC[,1:3] <- log10(ABC[,1:3] + 1)

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
sugar_alcohols$grouping <- NULL
sugar_alcohols[,1:3] <- log10(sugar_alcohols[,1:3] + 1)

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
fermentation$grouping <- NULL
fermentation[,1:3] <- log10(fermentation[,1:3] + 1)

# Prep the data from and write it to a file
gene_table$KEGG_code <- rownames(gene_table)
table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S1.tsv'
write.table(gene_table, file=table_file, sep='\t', row.names=FALSE, quote=FALSE)
rm(table_file, gene_table)

# Calculate rleative abundance for the rest of the genes
combined_mapping[,1:3] <- combined_mapping[,1:3] / rowSums(combined_mapping[,1:3])

#-------------------------------------------------------------------------------------------------------------------------#

# Remove gene sets from combined mapping
combined_mapping <- as.data.frame(subset(combined_mapping, !(rownames(combined_mapping) %in% rownames(amino_sugars))))
combined_mapping <- as.data.frame(subset(combined_mapping, !(rownames(combined_mapping) %in% rownames(stickland))))
combined_mapping <- as.data.frame(subset(combined_mapping, !(rownames(combined_mapping) %in% rownames(monosaccharides))))
combined_mapping <- as.data.frame(subset(combined_mapping, !(rownames(combined_mapping) %in% rownames(polysaccharides))))
combined_mapping <- as.data.frame(subset(combined_mapping, !(rownames(combined_mapping) %in% rownames(PTS))))
combined_mapping <- as.data.frame(subset(combined_mapping, !(rownames(combined_mapping) %in% rownames(ABC))))
combined_mapping <- as.data.frame(subset(combined_mapping, !(rownames(combined_mapping) %in% rownames(sugar_alcohols))))
combined_mapping <- as.data.frame(subset(combined_mapping, !(rownames(combined_mapping) %in% rownames(fermentation))))

#-------------------------------------------------------------------------------------------------------------------------#

# Define plot details
rainbow <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", 
             "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
fox <- wes_palette("FantasticFox")
tick_labels <- c('10%','30%','50%','70%','90%')
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_4.pdf'

# Open a PDF
pdf(file=plot_file, width=9, height=13)

# Create layout for multi-plot
layout(mat=matrix(c(1,1,1,1, 
                    1,1,1,1, 
                    1,1,1,1,
                    1,1,1,1,
                    2,3,4,5,
                    6,7,8,9), nrow=6, ncol=4, byrow=TRUE))

# Generate raw plot
par(mar=c(0,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('','',''), grid=FALSE, cex=0)

# Center triangle
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')

# 50% lines
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))

# Axis labels
text(x=-0.35, y=-0.43, labels='Cefoperazone (SPF)', cex=1.4)
text(x=-0.22, y=0.49, labels='Clindamycin (SPF)', cex=1.4, srt=60)
text(x=0.56, y=-0.12, labels='Streptomycin (SPF)', cex=1.4, srt=-60)

# Left axis - Clindmycin
lines(x=c(-0.52,-0.54), y=c(-0.233,-0.22))
lines(x=c(-0.462,-0.482), y=c(-0.133,-0.12))
lines(x=c(-0.404,-0.424), y=c(-0.033,-0.02))
lines(x=c(-0.346,-0.366), y=c(0.067,0.08))
lines(x=c(-0.289,-0.309), y=c(0.167,0.18))
lines(x=c(-0.23,-0.25), y=c(0.267,0.28))
lines(x=c(-0.173,-0.193), y=c(0.367,0.38))
lines(x=c(-0.115,-0.135), y=c(0.467,0.48))
lines(x=c(-0.057,-0.077), y=c(0.567,0.58))
text(x=c(-0.552, -0.436, -0.321, -0.205, -0.089),  
     y=c(-0.215, -0.015, 0.195, 0.395, 0.595),  
     labels=tick_labels, srt=60)

# Right axis - Streptomycin
lines(x=c(0.52,0.54), y=c(-0.233,-0.22))
lines(x=c(0.462,0.482), y=c(-0.133,-0.12))
lines(x=c(0.404,0.424), y=c(-0.033,-0.02))
lines(x=c(0.346,0.366), y=c(0.067,0.08))
lines(x=c(0.289,0.309), y=c(0.167,0.18))
lines(x=c(0.23,0.25), y=c(0.267,0.28))
lines(x=c(0.173,0.193), y=c(0.367,0.38))
lines(x=c(0.115,0.135), y=c(0.467,0.48))
lines(x=c(0.057,0.077), y=c(0.567,0.58))
text(x=c(0.553, 0.44, 0.324, 0.209, 0.092), 
     y=c(-0.212, -0.02, 0.185, 0.379, 0.584),
     labels=rev(tick_labels), srt=-60) 

# Bottom axis - Cefoperzone
lines(x=c(-0.462,-0.462), y=c(-0.333,-0.356))
lines(x=c(-0.346,-0.346), y=c(-0.333,-0.356))
lines(x=c(-0.231,-0.231), y=c(-0.333,-0.356))
lines(x=c(-0.115,-0.115), y=c(-0.333,-0.356))
lines(x=c(0,0), y=c(-0.333,-0.356))
lines(x=c(0.116,0.116), y=c(-0.333,-0.356))
lines(x=c(0.232,0.232), y=c(-0.333,-0.356))
lines(x=c(0.347,0.347), y=c(-0.333,-0.356))
lines(x=c(0.463,0.463), y=c(-0.333,-0.356))
text(x=c(-0.462, -0.231, 0, 0.232, 0.463),  
     y=c(-0.373, -0.373, -0.373, -0.373, -0.373),  
     labels=rev(tick_labels))

# Add points
tripoints(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], cex=0.8, col='gray65', pch=19)

# Color points by substrate
tripoints(x=PTS_relabund[,1], y=PTS_relabund[,2], z=PTS_relabund[,3], pch=21, cex=apply(PTS, 1, max)*3.5, bg=alpha(fox[3],0.7))
tripoints(x=ABC_relabund[,1], y=ABC_relabund[,2], z=ABC_relabund[,3], pch=21, cex=apply(ABC, 1, max)*3.5, bg=alpha(rainbow[7],0.7))
tripoints(x=monosaccharides_relabund[,1], y=monosaccharides_relabund[,2], z=monosaccharides_relabund[,3], pch=21, cex=apply(monosaccharides, 1, max)*3.5, bg=alpha(fox[1],0.7))
tripoints(x=stickland_relabund[,1], y=stickland_relabund[,2], z=stickland_relabund[,3], pch=21, cex=apply(stickland, 1, max)*3.5, bg=alpha(fox[2],0.7))
tripoints(x=sugar_alcohols_relabund[,1], y=sugar_alcohols_relabund[,2], z=sugar_alcohols_relabund[,3], pch=21, cex=apply(sugar_alcohols, 1, max)*3.5, bg=alpha('darkorchid3',0.7))
tripoints(x=fermentation_relabund[,1], y=fermentation_relabund[,2], z=fermentation_relabund[,3], pch=21, cex=apply(fermentation, 1, max)*3.5, bg=alpha(fox[5],0.7))
tripoints(x=polysaccharides_relabund[,1], y=polysaccharides_relabund[,2], z=polysaccharides_relabund[,3], pch=21, cex=apply(polysaccharides, 1, max)*3.5, bg=alpha('blue3',0.7))
tripoints(x=amino_sugars_relabund[,1], y=amino_sugars_relabund[,2], z=amino_sugars_relabund[,3], pch=21, cex=apply(amino_sugars, 1, max)*3.5, bg=alpha('firebrick1',0.7))

# Add the legend
legend(x=0.3, y=0.51, legend=c('Amino acid catabolism',
                               'Amino sugar metabolism',
                               'Glycolysis-associated',
                               'Polysaccharide metabolism',
                               'PTS transporters', 
                               'ABC sugar transporters', 
                               'Sugar alcohol metabolism', 
                               'Fermentation product synthesis', 
                               'All genes'), 
       ncol=1, pch=21, cex=1.4, pt.cex=c(2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,1.4), col=c('black','black','black','black','black','black','black','black','gray65'), 
       pt.bg=c(fox[2],
               'firebrick1',
               fox[1],
               'blue3', 
               fox[3], 
               rainbow[7], 
               'darkorchid3', 
               fox[5], 
               'gray65'), bty='n')
# Size legend
legend(x=-0.6, y=0.41, legend=c('     500 transcripts','','','   50 transcripts','',' 5 transcripts'), pch=21, col='black', pt.bg='gray87', pt.cex=c(9.446395,0,0,5.946395,0,2.446395), cex=1.4, bty='n')

# Add figure label
#text(x=-0.5, y=0.5, labels='a', font=2, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------#

# Individual plots

# amino acid (stickland) alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=FALSE, cex=0.3, col='gray75')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=stickland_relabund[,1], y=stickland_relabund[,2], z=stickland_relabund[,3], pch=21, cex=2, bg=fox[2])
text(x=0, y=-0.48, labels='Amino acid catabolism', cex=1.3)
#text(x=-0.5, y=0.5, labels='b', font=2, cex=1.6)

# amino sugars alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=FALSE, cex=0.3, col='gray75')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=amino_sugars_relabund[,1], y=amino_sugars_relabund[,2], z=amino_sugars_relabund[,3], 
          pch=21, cex=2, bg='firebrick1')
text(x=0, y=-0.48, labels='Amino sugar metabolism', cex=1.3)
#text(x=-0.5, y=0.5, labels='c', font=2, cex=1.6)

# glycolysis alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=FALSE, cex=0.3, col='gray75')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=monosaccharides_relabund[,1], y=monosaccharides_relabund[,2], z=monosaccharides_relabund[,3], 
          pch=21, cex=2, bg=fox[1])
text(x=0, y=-0.48, labels='Glycolysis-associated', cex=1.3)
#text(x=-0.5, y=0.5, labels='d', font=2, cex=1.6)

# polysaccharides alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=FALSE, cex=0.3, col='gray75')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=polysaccharides_relabund[,1], y=polysaccharides_relabund[,2], z=polysaccharides_relabund[,3], 
          pch=21, cex=2, bg='blue3')
text(x=0, y=-0.48, labels='Polysaccharide metabolism', cex=1.3)
#text(x=-0.5, y=0.5, labels='e', font=2, cex=1.6)

# PTS alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=FALSE, cex=0.3, col='gray75')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=PTS_relabund[,1], y=PTS_relabund[,2], z=PTS_relabund[,3], 
          pch=21, cex=2, bg=fox[3])
text(x=0, y=-0.48, labels='PTS transporters', cex=1.3)
#text(x=-0.5, y=0.5, labels='f', font=2, cex=1.6)

# ABC alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=FALSE, cex=0.3, col='gray75')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=ABC_relabund[,1], y=ABC_relabund[,2], z=ABC_relabund[,3], 
          pch=21, cex=2, bg=rainbow[7])
text(x=0, y=-0.48, labels='ABC sugar transporters', cex=1.3)
#text(x=-0.5, y=0.5, labels='g', font=2, cex=1.6)

# sugar alcohols alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=FALSE, cex=0.3, col='gray75')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=sugar_alcohols_relabund[,1], y=sugar_alcohols_relabund[,2], z=sugar_alcohols_relabund[,3], 
          pch=21, cex=2, bg='darkorchid3')
text(x=0, y=-0.48, labels='Sugar alcohol metabolism', cex=1.3)
#text(x=-0.5, y=0.5, labels='h', font=2, cex=1.6)

# fermentation alone
par(mar=c(1,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cef','Clinda','Strep'), grid=FALSE, cex=0.3, col='gray75')
lines(x=c(-0.288,0), y=c(0.1665,-0.333), col='gray68')
lines(x=c(-0.288,0.288), y=c(0.1665,0.1665), col='gray68')
lines(x=c(0,0.288), y=c(-0.333,0.1665), col='gray68')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=fermentation_relabund[,1], y=fermentation_relabund[,2], z=fermentation_relabund[,3], 
          pch=21, cex=2, bg=fox[5])
text(x=0, y=-0.48, labels='Fermentation product synthesis', cex=1.3)
#text(x=-0.5, y=0.5, labels='i', font=2, cex=1.6)

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
rm(amino_sugars, stickland, monosaccharides, polysaccharides, ABC, PTS, sugar_alcohols, fermentation, 
   ABC_relabund, amino_sugars_relabund, combined_mapping, fermentation_relabund, monosaccharides_relabund, 
   polysaccharides_relabund, PTS_relabund, stickland_relabund, sugar_alcohols_relabund, plot_file, 
   rainbow, tick_labels, fox)

for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()
