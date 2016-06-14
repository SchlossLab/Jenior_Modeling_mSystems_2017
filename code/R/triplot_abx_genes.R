
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
cefoperazone_file <- '/media/mjenior/Data/mapping/cdifficile630/cefoperazone_630.RNA_reads2cdf630.norm.annotated.txt'
clindamycin_file <- '/media/mjenior/Data/mapping/cdifficile630/clindamycin_630.RNA_reads2cdf630.norm.annotated.txt'
streptomycin_file <- '/media/mjenior/Data/mapping/cdifficile630/streptomycin_630.RNA_reads2cdf630.norm.annotated.txt'

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

# amino sugars
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
rm(mur, nag, acd, ldt, gne, anm, nan, glm)

# Stickland fermentation
pep <- rbind(subset(combined_mapping, grepl('pep.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('pep;', combined_mapping$gene))) # general peptidases
prd <- rbind(subset(combined_mapping, grepl('prd.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('prd;', combined_mapping$gene))) # proline fermentation
grd <- rbind(subset(combined_mapping, grepl('grd.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('grd;', combined_mapping$gene))) # glycine fermentation
ldhA <- subset(combined_mapping, grepl('ldhA;', combined_mapping$gene)) # (R)-2-hydroxyisocaproate dehydrogenase (leucine degradation)
had <- subset(combined_mapping, grepl('had.;', combined_mapping$gene)) # Leucine fermentation
tdcB <- subset(combined_mapping, grepl('tdcB;', combined_mapping$gene)) # threonine dehydratase
fdh <- subset(combined_mapping, grepl('fdh.;', combined_mapping$gene)) # Formate dehydrogenase
panB  <- subset(combined_mapping, grepl('panB;', combined_mapping$gene)) # Butanoate formation
arg <- subset(combined_mapping, grepl('arg.;', combined_mapping$gene)) # arginine deaminase operon
serA <- subset(combined_mapping, grepl('serA;', combined_mapping$gene)) # D-3-phosphoglycerate dehydrogenase (serine catabolism)
sdaB <- subset(combined_mapping, grepl('sdaB;', combined_mapping$gene)) # L-serine dehydratase (serine catabolism)
stickland <- rbind(pep, prd, grd, ldhA, had, tdcB, fdh, panB, arg, serA, sdaB)
rm(pep, prd, grd, ldhA, had, tdcB, fdh, panB, arg, serA, sdaB)

# Monosaccharides
glycolysis <- subset(combined_mapping, grepl('Glycolysis_/_Gluconeogenesis', combined_mapping$pathway))
galactose <- rbind(subset(combined_mapping, grepl('galE;', combined_mapping$gene)), 
                   subset(combined_mapping, grepl('galactose-1-phosphate_uridylyltransferase', combined_mapping$gene))) # hexose
mannose <- rbind(subset(combined_mapping, grepl('mng.;', combined_mapping$gene)),
                 subset(combined_mapping, grepl('man.;', combined_mapping$gene)),
                 subset(combined_mapping, grepl('pmi;', combined_mapping$gene)))# mannose utilization 
tagatose <- subset(combined_mapping, grepl('tagatose', combined_mapping$gene)) # hexose 
fructose <- rbind(subset(combined_mapping, grepl('fbp;', combined_mapping$gene)),
                  subset(combined_mapping, grepl('fru;', combined_mapping$gene)),
                  subset(combined_mapping, grepl('fru.;', combined_mapping$gene)),
                  subset(combined_mapping, grepl('fru...;', combined_mapping$gene))) # hexose
xylose <- subset(combined_mapping, grepl('xylose', combined_mapping$gene)) # pentose
monosaccharides <- rbind(glycolysis, galactose, mng, tagatose, fructose, xylose)
rm(glycolysis, galactose, mng, tagatose, fructose, xylose)

# Disaccharides
sucrose <- subset(combined_mapping, grepl('scr.;', combined_mapping$gene))
maltose <- rbind(subset(combined_mapping, grepl('maltose-6\'-phosphate_glucosidase', combined_mapping$gene)),
                 subset(combined_mapping, grepl('maa;', combined_mapping$gene)),
                 subset(combined_mapping, grepl('map.;', combined_mapping$gene)),
                 subset(combined_mapping, grepl('maltose_O-acetyltransferase', combined_mapping$gene)))
tre <- subset(combined_mapping, grepl('tre.;', combined_mapping$gene)) # Trehalose utilization operon
disaccharides <- rbind(sucrose,  maltose, tre)
rm(sucrose, maltose, tre)

# PTS systems
PTS <- subset(combined_mapping, grepl('PTS_system', combined_mapping$gene))

# ABC transporters
ABC <- subset(combined_mapping, grepl('ABC_transporter', combined_mapping$gene))

# sugar alcohols
srl <- subset(combined_mapping, grepl('srl.;', combined_mapping$gene)) # sorbitol utilization locus
srlE <- subset(combined_mapping, grepl('srlE.;', combined_mapping$gene)) # sorbitol import
mtl <- subset(combined_mapping, grepl('mtl.;', combined_mapping$gene)) # mannitol utilization locus
sugar_alcohols <- rbind(srl, srlE, mtl)
rm(srl, srlE, mtl)

# Butyrate production
buk <- rbind(subset(combined_mapping, grepl('buk;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('buk.;', combined_mapping$gene))) # Butyrate kinase
ptb <- rbind(subset(combined_mapping, grepl('ptb;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('ptb.;', combined_mapping$gene))) # phosphate butyryltransferase
butyrate <- rbind(buk, ptb)
rm(buk, ptb)

# GABA 
#gabT <- subset(combined_mapping, grepl('gabT;', combined_mapping$gene))

#-------------------------------------------------------------------------------------------------------------------------#

# Defineplot details
rainbow <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
#palette_plot(rainbow)
fox <- wes_palette("FantasticFox")
tick_labels <- c('10%','','30%','','50%','','70%','','90%')
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/cdf_abx_genes.pdf'

# Open a PDF
#pdf(file=plot_file, width=9, height=6.5)

# Generate raw plot
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('','',''), grid=seq(0.1,0.9,by=0.1), pch=20)

# 50% lines
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))

# Left axis - Clindmycin
#lines(x=c(-0.52,-0.54), y=c(-0.233,-0.233))
#lines(x=c(-0.462,-0.492), y=c(-0.133,-0.133))
#lines(x=c(-0.404,-0.424), y=c(-0.033,-0.033))
#lines(x=c(-0.346,-0.366), y=c(0.067,0.067))
#lines(x=c(-0.289,-0.309), y=c(0.167,0.167))
#lines(x=c(-0.23,-0.25), y=c(0.267,0.267))
#lines(x=c(-0.173,-0.193), y=c(0.367,0.367))
#lines(x=c(-0.115,-0.135), y=c(0.467,0.467))
#lines(x=c(-0.057,-0.077), y=c(0.567,0.567))
#text(x=c(-0.57,-0.522,-0.454,-0.396,-0.339,-0.28,-0.223,-0.165,-0.107), 
#     y=c(-0.233,-0.133,-0.033,0.067,0.167,0.267,0.367,0.467,0.567), 
#     labels=tick_labels, cex=0.9)
text(x=0, y=0.7, labels='Clindamycin', cex=1.4)

# Right axis - Streptomycin
#lines(x=c(0.52,0.54), y=c(-0.233,-0.233))
#lines(x=c(0.462,0.482), y=c(-0.133,-0.133))
#lines(x=c(0.404,0.424), y=c(-0.033,-0.033))
#lines(x=c(0.346,0.366), y=c(0.067,0.067))
#lines(x=c(0.289,0.309), y=c(0.167,0.167))
#lines(x=c(0.23,0.25), y=c(0.267,0.267))
#lines(x=c(0.173,0.193), y=c(0.367,0.367))
#lines(x=c(0.115,0.135), y=c(0.467,0.467))
#lines(x=c(0.057,0.077), y=c(0.567,0.567))
#text(x=c(0.57,0.522,0.454,0.396,0.339,0.28,0.223,0.165,0.107), 
#     y=c(-0.233,-0.133,-0.033,0.067,0.167,0.267,0.367,0.467,0.567), 
#     labels=rev(tick_labels), cex=0.9)
text(x=0.65, y=-0.38, labels='Streptomycin', cex=1.4)

# Bottom axis - Cefoperzone
#lines(x=c(-0.462,-0.462), y=c(-0.333,-0.353))
#lines(x=c(-0.346,-0.346), y=c(-0.333,-0.353))
#lines(x=c(-0.231,-0.231), y=c(-0.333,-0.353))
#lines(x=c(-0.115,-0.115), y=c(-0.333,-0.353))
#lines(x=c(0,0), y=c(-0.333,-0.353))
#lines(x=c(0.116,0.116), y=c(-0.333,-0.353))
#lines(x=c(0.232,0.232), y=c(-0.333,-0.353))
#lines(x=c(0.347,0.347), y=c(-0.333,-0.353))
#lines(x=c(0.463,0.463), y=c(-0.333,-0.353))
#text(x=c(-0.462,-0.346,-0.231,-0.115,0,0.116,0.232,0.347,0.463), 
#     y=c(-0.373,-0.373,-0.373,-0.373,-0.373,-0.373,-0.373,-0.373,-0.373), 
#     labels=rev(tick_labels), cex=0.9)
text(x=-0.65, y=-0.38, labels='Cefoperzone', cex=1.4)

# Color points by substrate
tripoints(x=PTS[,1], y=PTS[,2], z=PTS[,3], pch=21, cex=2, bg=fox[3])
tripoints(x=ABC[,1], y=ABC[,2], z=ABC[,3], pch=21, cex=2, bg=rainbow[7])
tripoints(x=monosaccharides[,1], y=monosaccharides[,2], z=monosaccharides[,3], pch=21, cex=2, bg=fox[1])
#tripoints(x=disaccharides[,1], y=disaccharides[,2], z=disaccharides[,3], pch=21, cex=2, bg=rainbow[7])
tripoints(x=stickland[,1], y=stickland[,2], z=stickland[,3], pch=21, cex=2, bg=fox[2])
tripoints(x=sugar_alcohols[,1], y=sugar_alcohols[,2], z=sugar_alcohols[,3], pch=21, cex=2, bg=rainbow[1])
tripoints(x=butyrate[,1], y=butyrate[,2], z=butyrate[,3], pch=21, cex=2, bg=fox[5])

# Add the legend
legend('topright', legend=c('Monosaccharide catabolism', 'Sugar alcohol catabolism', 'Stickland fermentation', 'Butyrate production', 'PTS genes', 'ABC transporters', 'Other'), 
    ncol=1, pch=c(21,21,21,21,21,21,20), pt.cex=2, col='black', pt.bg=c(fox[1],rainbow[1],fox[2],fox[5],fox[3],rainbow[7], NA))

# Add figure label
text(x=-0.8, y=0.75, labels='A', font=2, cex=2)

#dev.off()



