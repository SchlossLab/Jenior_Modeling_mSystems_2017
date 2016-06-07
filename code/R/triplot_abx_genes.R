
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

# Subset by gene annotations

# amino sugars
acetylglucosamine <- subset(combined_mapping[,1:3], grepl('*N-acetylglucosamine*', combined_mapping$gene))
acetylmannosamine <- subset(combined_mapping[,1:3], grepl('*N-acetylmannosamine*', combined_mapping$gene))
acetylmuramate <- subset(combined_mapping[,1:3], grepl('*N-acetylmuramate*', combined_mapping$gene))
amino_sugars <- rbind(acetylglucosamine, acetylmannosamine, acetylmuramate)
rm(acetylglucosamine, acetylmannosamine, acetylmuramate)

# Stickland substrates
proline <- subset(combined_mapping[,1:3], grepl('*proline*', combined_mapping$gene))
glycine <- subset(combined_mapping[,1:3], grepl('*glycine*', combined_mapping$gene))
arginine <- subset(combined_mapping[,1:3], grepl('*arginine*', combined_mapping$gene))
threonine <- subset(combined_mapping[,1:3], grepl('*threonine*', combined_mapping$gene))
methionine <- subset(combined_mapping[,1:3], grepl('*methionine*', combined_mapping$gene))
serine <- subset(combined_mapping[,1:3], grepl('*serine*', combined_mapping$gene))
alanine <- subset(combined_mapping[,1:3], grepl('*alanine*', combined_mapping$gene))
stickland <- rbind(proline, glycine, arginine, threonine, methionine, serine, alanine)
rm(proline, glycine, arginine, threonine, methionine, serine, alanine)

# Saccharides
# hexose
galactose <- subset(combined_mapping[,1:3], grepl('*galactose*', combined_mapping$gene))
mannose <- subset(combined_mapping[,1:3], grepl('*mannose*', combined_mapping$gene))
glucose <- subset(combined_mapping[,1:3], grepl('*glucose*', combined_mapping$gene))
tagatose <- subset(combined_mapping[,1:3], grepl('*tagatose*', combined_mapping$gene))
fructose <- subset(combined_mapping[,1:3], grepl('*fructose*', combined_mapping$gene)) #keto-
hexose <- rbind(galactose, mannose, glucose, tagatose, fructose)
# pentose
xylose <- subset(combined_mapping[,1:3], grepl('*xylose*', combined_mapping$gene))
ribose <- subset(combined_mapping[,1:3], grepl('*ribose*', combined_mapping$gene))
pentose <- rbind(xylose, ribose)
# disaccharides
sucrose <- subset(combined_mapping[,1:3], grepl('*sucrose*', combined_mapping$gene))
lactose <- subset(combined_mapping[,1:3], grepl('*lactose*', combined_mapping$gene))
maltose <- subset(combined_mapping[,1:3], grepl('*maltose*', combined_mapping$gene))
trehalose <- subset(combined_mapping[,1:3], grepl('*trehalose*', combined_mapping$gene))
disaccharides <- rbind(sucrose, lactose, maltose, trehalose)
# all
saccharides <- rbind(galactose, tagatose, trehalose, mannose, xylose, ribose, fructose, glucose, maltose, lactose, sucrose)
rm(galactose, tagatose, trehalose, mannose, xylose, ribose, glucose, maltose, lactose, sucrose)

# sugar alcohols
ribitol <- subset(combined_mapping[,1:3], grepl('*ribitol*', combined_mapping$gene))
sorbitol <- subset(combined_mapping[,1:3], grepl('*sorbitol*', combined_mapping$gene))
mannitol <- subset(combined_mapping[,1:3], grepl('*mannitol*', combined_mapping$gene))
sugar_alcohols <- rbind(ribitol, sorbitol, mannitol)
rm(ribitol, sorbitol, mannitol)

# nucleosides
uridine <- subset(combined_mapping[,1:3], grepl('*uridine*', combined_mapping$gene))
inosine <- subset(combined_mapping[,1:3], grepl('*inosine*', combined_mapping$gene))
adenosine <- subset(combined_mapping[,1:3], grepl('*adenosine*', combined_mapping$gene))
nucleosides <- rbind(uridine, inosine, adenosine)
rm(uridine, inosine, adenosine)

# short-chain fatty acids
butyrate <- subset(combined_mapping[,1:3], grepl('*butyrate*', combined_mapping$gene))

#-------------------------------------------------------------------------------------------------------------------------#

# Defineplot details
rainbow <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
fox <- wes_palette("FantasticFox")
tick_labels <- c('10%','','30%','','50%','','70%','','90%')
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/cdf_abx_genes.pdf'
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



