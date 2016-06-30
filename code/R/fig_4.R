
deps <- c('vegan', 'klaR', 'wesanderson');
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
#cefoperazone_file <- '/media/mjenior/Data/mapping/cdifficile630/cefoperazone_630.RNA_reads2cdf630.norm.annotated.txt'
#clindamycin_file <- '/media/mjenior/Data/mapping/cdifficile630/clindamycin_630.RNA_reads2cdf630.norm.annotated.txt'
#streptomycin_file <- '/media/mjenior/Data/mapping/cdifficile630/streptomycin_630.RNA_reads2cdf630.norm.annotated.txt'

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

# Calculate the factors to scale points size by
averages <- log10(rowSums(combined_mapping[, 1:3]) / 3) * 2.5

# Convert each gene into the fraction of the transcription for that gene across treatments
combined_mapping[combined_mapping == 0] <- 1
combined_mapping[,1:3] <- combined_mapping[, 1:3] / rowSums(combined_mapping[,1:3])

# Eliminate genes with no transcripts mapping
combined_mapping <- combined_mapping[rowSums(combined_mapping[,1:3]) != 0, ] 

# Free up memory
rm(cefoperazone_file, clindamycin_file, streptomycin_file)
rm(cefoperazone, clindamycin, streptomycin)

#----------------------------------------------------------------------------------------------------------------------------#

# Subset by gene annotations

# Amino sugar catabolism
mur <- rbind(subset(combined_mapping, grepl('mur.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('mur;', combined_mapping$gene))) # 
nag <- rbind(subset(combined_mapping, grepl('nag.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('nag;', combined_mapping$gene))) # 
acd <- rbind(subset(combined_mapping, grepl('acd.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('acd;', combined_mapping$gene))) # 
ldt <- rbind(subset(combined_mapping, grepl('ldt.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('ldt;', combined_mapping$gene))) # 
gne <- rbind(subset(combined_mapping, grepl('gne.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('gne;', combined_mapping$gene))) # 
anm <- rbind(subset(combined_mapping, grepl('anm.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('anm;', combined_mapping$gene))) # 
nan <- rbind(subset(combined_mapping, grepl('nan.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('nan;', combined_mapping$gene))) #
glm <- rbind(subset(combined_mapping, grepl('glm.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('glm;', combined_mapping$gene)))
amino_sugars <- rbind(mur, nag, acd, ldt, gne, anm, nan, glm)

# Stickland fermentation
pep <- rbind(subset(combined_mapping, grepl('pep.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('pep;', combined_mapping$gene))) # general peptidases
prd <- rbind(subset(combined_mapping, grepl('prd.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('prd;', combined_mapping$gene))) # proline fermentation
grd <- rbind(subset(combined_mapping, grepl('grd.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('grd;', combined_mapping$gene))) # glycine fermentation
kamA <- subset(combined_mapping, grepl('kam.;', combined_mapping$gene)) # lycine fermentation
ldhA <- subset(combined_mapping, grepl('ldhA;', combined_mapping$gene)) # (R)-2-hydroxyisocaproate dehydrogenase (leucine degradation)
had <- subset(combined_mapping, grepl('had.;', combined_mapping$gene)) # Leucine fermentation
tdcB <- subset(combined_mapping, grepl('tdcB;', combined_mapping$gene)) # threonine dehydratase
fdh <- subset(combined_mapping, grepl('fdh.;', combined_mapping$gene)) # Formate dehydrogenase
panB  <- subset(combined_mapping, grepl('panB;', combined_mapping$gene)) # Butanoate formation
arg <- subset(combined_mapping, grepl('arg.;', combined_mapping$gene)) # arginine deaminase operon
serA <- subset(combined_mapping, grepl('serA;', combined_mapping$gene)) # D-3-phosphoglycerate dehydrogenase (serine catabolism)
sdaB <- subset(combined_mapping, grepl('sdaB;', combined_mapping$gene)) # L-serine dehydratase (serine catabolism)
stickland <- rbind(pep, prd, grd, ldhA, had, tdcB, fdh, panB, arg, serA, sdaB, kamA)

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
glycolysis <- rbind(gap, gpmI, pfk, tpi, pyk, eno, pgm)
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

# PTS systems
PTS <- rbind(subset(combined_mapping, grepl('PTS_system', combined_mapping$gene)),
             subset(combined_mapping, grepl('pyridoxal_phosphate-dependent_transferase', combined_mapping$gene)))

# ABC transporters
ABC <- subset(combined_mapping, grepl('ABC_transporter_sugar', combined_mapping$gene))

# sugar alcohols
srl <- subset(combined_mapping, grepl('srl.;', combined_mapping$gene)) # sorbitol utilization locus
srlE <- subset(combined_mapping, grepl('srlE.;', combined_mapping$gene)) # sorbitol import
mtl <- subset(combined_mapping, grepl('mtl.;', combined_mapping$gene)) # mannitol utilization locus
sugar_alcohols <- rbind(srl, srlE, mtl)

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

#-------------------------------------------------------------------------------------------------------------------------#

# Define plot details
rainbow <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
#palette_plot(rainbow)
fox <- wes_palette("FantasticFox")
tick_labels <- c('10%','','30%','','50%','','70%','','90%')
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_4.pdf'

# Open a PDF
pdf(file=plot_file, width=9, height=9)

# Create layout for multi-plot
layout(mat=matrix(c(1,1,1,2, 
                    1,1,1,3, 
                    1,1,1,4, 
                    5,6,7,8), 
                  nrow=4, ncol=4, byrow=TRUE))

# Generate raw plot
par(mar=c(0,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('','',''), grid=seq(0.1,0.9,by=0.1), cex=0.8)

# 50% lines
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))

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
legend(x=0.26, y=0.63, legend=c('Monosaccharide catabolism', 'Polysaccharide catabolism', 'Sugar alcohol catabolism', 'Stickland reactions', 'Fermentation end steps', 'PEP group translocators', 'ABC sugar transporters', 'Other'), 
    ncol=1, pch=21, cex=1.4, pt.cex=c(2.5,2.5,2.5,2.5,2.5,2.5,2.5,1.4), col='black', pt.bg=c(fox[1],'blue3','darkorchid3',fox[2],fox[5],fox[3],rainbow[7], 'white'), bty='n')
# Add figure label
legend(x=-0.6, y=0.6, legend='A', cex=2, bty='n')

# monosaccharides alone
par(mar=c(0,0,0,0))
triplot(x=monosaccharides[,1], y=monosaccharides[,2], z=monosaccharides[,3], 
        frame=TRUE, label=c('','',''), grid=seq(0.1,0.9,by=0.1), cex=0.8)
legend('topleft', legend='B', cex=2, bty='n')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=monosaccharides[,1], y=monosaccharides[,2], z=monosaccharides[,3], pch=21, cex=2, bg=fox[1])
text(x=0, y=-0.41, labels='Monosaccharide catabolism')

# polysaccharides alone
par(mar=c(0,0,0,0))
triplot(x=polysaccharides[,1], y=polysaccharides[,2], z=polysaccharides[,3], 
        frame=TRUE, label=c('','',''), grid=seq(0.1,0.9,by=0.1), cex=0.8)
legend('topleft', legend='C', cex=2, bty='n')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=polysaccharides[,1], y=polysaccharides[,2], z=polysaccharides[,3], pch=21, cex=2, bg='blue3')
text(x=0, y=-0.41, labels='Polysaccharide catabolism')

# sugar_alcohols alone
par(mar=c(0,0,0,0))
triplot(x=sugar_alcohols[,1], y=sugar_alcohols[,2], z=sugar_alcohols[,3], 
        frame=TRUE, label=c('','',''), grid=seq(0.1,0.9,by=0.1), cex=0.8)
legend('topleft', legend='D', cex=2, bty='n')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=sugar_alcohols[,1], y=sugar_alcohols[,2], z=sugar_alcohols[,3], pch=21, cex=2, bg='darkorchid3')
text(x=0, y=-0.41, labels='Sugar alcohol catabolism')

# stickland alone
par(mar=c(0,0,0,0))
triplot(x=stickland[,1], y=stickland[,2], z=stickland[,3], 
        frame=TRUE, label=c('','',''), grid=seq(0.1,0.9,by=0.1), cex=0.8)
legend('topleft', legend='E', cex=2, bty='n')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=stickland[,1], y=stickland[,2], z=stickland[,3], pch=21, cex=2, bg=fox[2])
text(x=0, y=-0.41, labels='Stickland reactions')

# fermentation
par(mar=c(0,0,0,0))
triplot(x=fermentation[,1], y=fermentation[,2], z=fermentation[,3], 
        frame=TRUE, label=c('','',''), grid=seq(0.1,0.9,by=0.1), cex=0.8)
legend('topleft', legend='F', cex=2, bty='n')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=fermentation[,1], y=fermentation[,2], z=fermentation[,3], pch=21, cex=2, bg=fox[5])
text(x=0, y=-0.41, labels='Fermentation end steps')

# PTS alone
par(mar=c(0,0,0,0))
triplot(x=PTS[,1], y=PTS[,2], z=PTS[,3], 
        frame=TRUE, label=c('','',''), grid=seq(0.1,0.9,by=0.1), cex=0.8)
legend('topleft', legend='G', cex=2, bty='n')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=PTS[,1], y=PTS[,2], z=PTS[,3], pch=21, cex=2, bg=fox[3])
text(x=0, y=-0.41, labels='PEP group translocators')

# ABC alone
par(mar=c(0,0,0,0))
triplot(x=ABC[,1], y=ABC[,2], z=ABC[,3], 
        frame=TRUE, label=c('','',''), grid=seq(0.1,0.9,by=0.1), cex=0.8)
legend('topleft', legend='H', cex=2, bty='n')
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))
tripoints(x=ABC[,1], y=ABC[,2], z=ABC[,3], pch=21, cex=2, bg=rainbow[7])
text(x=0, y=-0.41, labels='ABC sugar transporters')

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
rm(abfD, ABC, acd, acetate, adh, amino_sugars, anm, arg, buk, cat, cel, combined_mapping, eno, fdh, fermentation, fructose,
   galactose, glycolysis, gne, gpmI, grd, had, hbd, kamA, ldhA, ldt, maltose, mannose, monosaccharides, mtl, mur, nag, nan,
   panB, pep, pfk, pgm, polysaccharides, prd, ptb, PTS, pyk, sdaB, serA, srl, srlE, stickland, sucD, sucrose, sugar_alcohols,
   tagatose, tdcB, tpi, tre, valerate, xylose, averages, fox, plot_file, rainbow, sub_size, tick_labels, palette_plot, gap,
   glm, glucosidase)
