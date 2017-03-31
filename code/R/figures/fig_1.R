
# Start with clean environment
rm(list=ls())
gc()

# Load dependencies
deps <- c('shape', 'wesanderson', 'plotrix');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Select files
cfu_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/cfu.dat'
toxin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/toxin_titer.dat'

# Read in data
cfu <- read.delim(cfu_file, sep='\t', header=T)
toxin <- read.delim(toxin_file, sep='\t', header=T)
rm(cfu_file, toxin_file)

# Format CFU data and collect summary statistics
cfu[cfu == 0] <- 50
cfu$cfu_vegetative <- log10(cfu$cfu_vegetative)
cfu$cfu_spore <- log10(cfu$cfu_spore)
cfu$mouse <- NULL
cfu <- subset(cfu, cage < 4 ) # Remove uninfected controls
cfu$cage <- NULL
cfu$treatment <- factor(cfu$treatment, levels=c('conventional', 'streptomycin', 'cefoperazone', 'clindamycin', 'germfree'))
vegetative_cfu <- cfu
vegetative_cfu$cfu_spore <- NULL
spore_cfu <- cfu
spore_cfu$cfu_vegetative <- NULL
cfu_veg_calc <- vegetative_cfu
cfu_spore_calc <- spore_cfu
cef <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'cefoperazone', 2]))
strep <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'streptomycin', 2]))
clinda <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'clindamycin', 2]))
gf <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'germfree', 2]))
conv <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'conventional', 2]))
vege_medians <- c(conv, strep, cef, clinda, gf)
cef <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'cefoperazone', 2]))
strep <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'streptomycin', 2]))
clinda <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'clindamycin', 2]))
gf <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'germfree', 2]))
conv <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'conventional', 2]))
spore_medians <- c(conv, strep, cef, clinda, gf)
rm(cfu, cef, strep, clinda, gf, conv)

# Format toxin data and find summary statistics
toxin$mouse <- NULL
toxin$cage <- NULL
toxin$treatment <- factor(toxin$treatment, levels=c('Conventional', 'Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree'))
toxin$titer[toxin$titer < 2.3] <- 2.15
toxin_calc <- toxin
cef <- as.numeric(median(toxin[toxin$treatment == 'Cefoperazone', 2]))
strep <- as.numeric(median(toxin[toxin$treatment == 'Streptomycin', 2]))
clinda <- as.numeric(median(toxin[toxin$treatment == 'Clindamycin', 2]))
gf <- as.numeric(median(toxin[toxin$treatment == 'Germfree', 2]))
conv <- as.numeric(median(toxin[toxin$treatment == 'Conventional', 2]))
toxin_medians <- c(conv, strep, cef, clinda, gf)
rm(cef, strep, clinda, gf, conv)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Calculate significant differences

# Vegetative CFU
cefoperazone <- subset(cfu_veg_calc, treatment == 'cefoperazone')$cfu_vegetative
clindamycin <- subset(cfu_veg_calc, treatment == 'clindamycin')$cfu_vegetative
streptomycin <- subset(cfu_veg_calc, treatment == 'streptomycin')$cfu_vegetative
germfree <- subset(cfu_veg_calc, treatment == 'germfree')$cfu_vegetative
p_values <- p.adjust(c(wilcox.test(germfree, cefoperazone, exact=F)$p.value,
           wilcox.test(germfree, clindamycin, exact=F)$p.value,
           wilcox.test(germfree, streptomycin, exact=F)$p.value,
           wilcox.test(clindamycin, cefoperazone, exact=F)$p.value,
           wilcox.test(streptomycin, cefoperazone, exact=F)$p.value,
           wilcox.test(clindamycin, streptomycin, exact=F)$p.value), method='BH')
cat("Vegetative: Germfree vs Cefoperazone: p =", p_values[1],'\n',
    "Vegetative: Germfree vs Clindamycin: p =", p_values[2],'\n',
    "Vegetative: Germfree vs Streptomycin: p =", p_values[3],'\n',
    "Vegetative: Cefoperazone vs Clindamycin: p =", p_values[4],'\n',
    "Vegetative: Cefoperazone vs Streptomycin: p =", p_values[5],'\n',
    "Vegetative: Clindamycin vs Streptomycin: p =", p_values[6])

# Spore CFU
cefoperazone <- subset(cfu_spore_calc, treatment == 'cefoperazone')$cfu_spore
clindamycin <- subset(cfu_spore_calc, treatment == 'clindamycin')$cfu_spore
streptomycin <- subset(cfu_spore_calc, treatment == 'streptomycin')$cfu_spore
germfree <- subset(cfu_spore_calc, treatment == 'germfree')$cfu_spore
p_values <- p.adjust(c(wilcox.test(germfree, cefoperazone, exact=F)$p.value,
           wilcox.test(germfree, clindamycin, exact=F)$p.value,
           wilcox.test(germfree, streptomycin, exact=F)$p.value,
           wilcox.test(clindamycin, cefoperazone, exact=F)$p.value,
           wilcox.test(streptomycin, cefoperazone, exact=F)$p.value,
           wilcox.test(clindamycin, streptomycin, exact=F)$p.value), method='BH')
cat("Spores: Germfree vs Cefoperazone: p =", p_values[1],'\n',
    "Spores: Germfree vs Clindamycin: p =", p_values[2],'\n',
    "Spores: Germfree vs Streptomycin: p =", p_values[3],'\n',
    "Spores: Cefoperazone vs Clindamycin: p =", p_values[4],'\n',
    "Spores: Cefoperazone vs Streptomycin: p =", p_values[5],'\n',
    "Spores: Clindamycin vs Streptomycin: p =", p_values[6])

# Toxin titer
cefoperazone <- subset(toxin_calc, treatment == 'Cefoperazone')$titer
clindamycin <- subset(toxin_calc, treatment == 'Clindamycin')$titer
streptomycin <- subset(toxin_calc, treatment == 'Streptomycin')$titer
germfree <- subset(toxin_calc, treatment == 'Germfree')$titer
p_values <- p.adjust(c(wilcox.test(germfree, cefoperazone, exact=F)$p.value,
           wilcox.test(germfree, clindamycin, exact=F)$p.value,
           wilcox.test(germfree, streptomycin, exact=F)$p.value,
           wilcox.test(cefoperazone, clindamycin, exact=F)$p.value,
           wilcox.test(cefoperazone, streptomycin, exact=F)$p.value,
           wilcox.test(clindamycin, streptomycin, exact=F)$p.value), method='BH')
cat("Toxin: Germfree vs Cefoperazone: p =", p_values[1],'\n',
    "Toxin: Germfree vs Clindamycin: p =", p_values[2],'\n',
    "Toxin: Germfree vs Streptomycin: p =", p_values[3],'\n',
    "Toxin: Cefoperazone vs Clindamycin: p =", p_values[4],'\n',
    "Toxin: Cefoperazone vs Streptomycin: p =", p_values[5],'\n',
    "Toxin: Clindamycin vs Streptomycin: p =", p_values[6])

rm(germfree, cefoperazone, clindamycin, streptomycin)

# Change undetectable points for to points actually on the plot for clarity
toxin$titer[toxin$titer <= 2.3] <- 2.15

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_1.pdf'
select_palette <- c('gray40', wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen')
pdf(file=plot_file, width=6, height=8)
layout(matrix(c(1,
                2,
                3),
              nrow=3, ncol=1, byrow = TRUE))
par(las=1, mgp=c(3,0.7,0), yaxs='i')

#-------------------------------------------------------------------------------------------------------------------------------------#

# A.  Vegetative cell CFU
par(mar=c(0.6,5,1,1))
stripchart(cfu_vegetative~treatment, data=vegetative_cfu, vertical=T, pch=19, lwd=2.2,
           ylim=c(1,9), xaxt='n', yaxt='n', cex=0, col=select_palette,
           ylab='Vegetative cfu/g content', method='jitter', jitter=0.15, cex.lab=1.2)
labelsY <- c(0, parse(text=paste(rep(10,8), '^', seq(2,9,1), sep='')))
axis(side=2, at=c(1:9), labelsY, tick=TRUE, cex.axis=1.2)
abline(h=2, col="black", lty=2, lwd=1.5) # LOD
stripchart(cfu_vegetative~treatment, data=vegetative_cfu, vertical=T, pch=19, lwd=2.5,
           ylim=c(1,9), xaxt='n', yaxt='n', cex=2, col=select_palette,
           ylab='Vegetative cfu/g content', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)

# Draw axis break
axis.break(2, 1.5, style='slash')

# Draw median
segments(0.7, vege_medians[1], 1.3, vege_medians[1], lwd=3)
segments(1.7, vege_medians[2], 2.3, vege_medians[2], lwd=3)
segments(2.7, vege_medians[3], 3.3, vege_medians[3], lwd=3)
segments(3.7, vege_medians[4], 4.3, vege_medians[4], lwd=3)
segments(4.7, vege_medians[5], 5.3, vege_medians[5], lwd=3)

# Adding significance to plot
text(1, 2.4, labels='*', cex=3, font=2, col='gray40') # Resistant mice

mtext('A', side=2, line=2, las=2, adj=2, padj=-6.5, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# B.  Spore CFU
par(mar=c(0.6,5,0.6,1))
stripchart(cfu_spore~treatment, data=spore_cfu, vertical=T, pch=19, lwd=2.2, 
           ylim=c(1,9), xaxt='n', yaxt='n', cex=0, col=select_palette,
           ylab='Spore cfu/g content', method='jitter', jitter=0.15, cex.lab=1.2)
axis(side=2, at=c(1:9), labelsY, tick=TRUE, cex.axis=1.2)
abline(h=2, col="black", lty=2, lwd=1.5)
stripchart(cfu_spore~treatment, data=spore_cfu, vertical=T, pch=19, lwd=2.5, 
           ylim=c(1,9), xaxt='n', yaxt='n', cex=2, col=select_palette,
           ylab='Spore cfu/g content', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)

# Draw axis break
axis.break(2, 1.5, style='slash') 

# Draw median
segments(0.7, spore_medians[1], 1.3, spore_medians[1], lwd=3) # cefoperazone
segments(1.7, spore_medians[2], 2.3, spore_medians[2], lwd=3) # streptomycin
segments(2.7, spore_medians[3], 3.3, spore_medians[3], lwd=3) # clindamycin
segments(3.7, spore_medians[4], 4.3, spore_medians[4], lwd=3) # germfree
segments(4.7, spore_medians[5], 5.3, spore_medians[5], lwd=3) # conventional

# Adding significance to plot
segments(x0=c(2,3,4), y0=c(7,7.5,8), x1=c(5,5,5), y1=c(7,7.5,8), lwd=2)
text(c(3.5,4,4.5), c(7.2,7.7,8.2), labels=c('*','*','*'), cex=2.2)
text(1, 2.4, labels='*', cex=3, font=2, col='gray40') # Resistant mice

mtext('B', side=2, line=2, las=2, adj=2, padj=-6.5, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# C.  Toxin data
par(mar=c(3,5,0.6,1))
stripchart(titer~treatment, data=toxin, vertical=T, pch=19, lwd=2.2,
           ylim=c(1.9,3.8), xlim=c(0.5,5.5), xaxt='n', yaxt='n', col=select_palette,
           ylab='Toxin titer/g content', xlab='',
           method='jitter', jitter=0.15, cex=0, cex.lab=1.2)
abline(h=2.3, lty=2, lwd=1.5) # LOD
stripchart(titer~treatment, data=toxin, vertical=T, pch=19, lwd=2.5,
           ylim=c(1.5,3.8), xlim=c(0.5,5.5), xaxt='n', yaxt='n', col=select_palette,
           ylab='Toxin titer/g content', xlab='',
           method='jitter', jitter=0.15, cex=2, cex.lab=1.2, add=TRUE)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-Germfree'), side=1, at=c(1:5), padj=1, cex=0.85)
axis(side=2, at=c(1.9,2.3,2.6,2.9,3.2,3.5,3.8), labels=c('0','200','398','794','1585','3162','6310'), cex.axis=1.2)

# Draw axis break
axis.break(2, 2.05, style='slash') 

# Draw median
segments(0.7, toxin_medians[1], 1.3, toxin_medians[1], lwd=3) # cefoperazone
segments(1.7, toxin_medians[2], 2.3, toxin_medians[2], lwd=3) # streptomycin
segments(2.7, toxin_medians[3], 3.3, toxin_medians[3], lwd=3) # clindamycin
segments(3.7, toxin_medians[4], 4.3, toxin_medians[4], lwd=3) # germfree
segments(4.7, toxin_medians[5], 5.3, toxin_medians[5], lwd=3) # conventional

# Adding significance to plot
segments(x0=c(2,3,4), y0=c(3.1,3.25,3.4), x1=c(5,5,5), y1=c(3.1,3.25,3.4), lwd=2)
text(c(3.5,4,4.5), c(3.16,3.31,3.46), labels=c('*','*','*'), cex=2.2)
text(1, 2.4, labels='*', cex=3, font=2, col='gray40') # Resistant mice

# Plot label
mtext('C', side=2, line=2, las=2, adj=2, padj=-5.7, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
dev.off()
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

