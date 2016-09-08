
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
cfu_raw <- cfu
toxin <- read.delim(toxin_file, sep='\t', header=T)
toxin_raw <- toxin
rm(cfu_file, toxin_file)

# Format CFU data and collect summary statistics
cfu[cfu == 0] <- 100
cfu$cfu_vegetative <- log10(cfu$cfu_vegetative)
cfu$cfu_spore <- log10(cfu$cfu_spore)
cfu$mouse <- NULL
cfu <- subset(cfu, cage < 4 ) # Remove uninfected controls
cfu$cage <- NULL
cfu$treatment <- factor(cfu$treatment, levels=c('streptomycin', 'cefoperazone', 'clindamycin', 'germfree', 'conventional'))
vegetative_cfu <- cfu
vegetative_cfu$cfu_spore <- NULL
spore_cfu <- cfu
spore_cfu$cfu_vegetative <- NULL
cef <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'cefoperazone', 2]))
strep <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'streptomycin', 2]))
clinda <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'clindamycin', 2]))
gf <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'germfree', 2]))
conv <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'conventional', 2]))
vege_medians <- c(strep, cef, clinda, gf, conv)
cef <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'cefoperazone', 2]))
strep <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'streptomycin', 2]))
clinda <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'clindamycin', 2]))
gf <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'germfree', 2]))
conv <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'conventional', 2]))
spore_medians <- c(strep, cef, clinda, gf, conv)
rm(cfu, cef, strep, clinda, gf, conv)
vegetative_cfu$color <- ifelse(vegetative_cfu$cfu_vegetative == 2.0, 'gray50', 'black')
vegetative_cfu$cfu_vegetative[vegetative_cfu$cfu_vegetative == 2.0] <- 1.6
spore_cfu$cfu_spore[spore_cfu$cfu_spore == 2.0] <- 1.6

# Format toxin data and find summary statistics
toxin$mouse <- NULL
toxin$cage <- NULL
toxin$treatment <- factor(toxin$treatment, levels=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree', 'Conventional'))
cef <- as.numeric(median(toxin[toxin$treatment == 'Cefoperazone', 2]))
strep <- as.numeric(median(toxin[toxin$treatment == 'Streptomycin', 2]))
clinda <- as.numeric(median(toxin[toxin$treatment == 'Clindamycin', 2]))
gf <- as.numeric(median(toxin[toxin$treatment == 'Germfree', 2]))
conv <- as.numeric(median(toxin[toxin$treatment == 'Conventional', 2]))
toxin_medians <- c(strep, cef, clinda, gf, conv)
rm(cef, strep, clinda, gf, conv)
toxin$titer[toxin$titer <= 2.0] <- 1.9

#-------------------------------------------------------------------------------------------------------------------------------------#

# Calculate significant differences

toxin_raw$titer[toxin_raw$titer < 2.0] <- NA
cefoperazone <- subset(toxin_raw, treatment == 'Cefoperazone')$titer
clindamycin <- subset(toxin_raw, treatment == 'Clindamycin')$titer
streptomycin <- subset(toxin_raw, treatment == 'Streptomycin')$titer
germfree <- subset(toxin_raw, treatment == 'Germfree')$titer
rm(toxin_raw)

wilcox.test(germfree, cefoperazone, exact=F) # p-value = 0.001533 ***
wilcox.test(germfree, clindamycin, exact=F) # p-value = 0.0001881 ***
wilcox.test(germfree, streptomycin, exact=F) # p-value = 0.05 *
wilcox.test(cefoperazone, clindamycin, exact=F) # n.s.
wilcox.test(cefoperazone, streptomycin, exact=F) # n.s.
wilcox.test(clindamycin, streptomycin, exact=F) # n.s.

p_values <- c(0.001533, 0.0001881, 0.05, 0.1, 0.1, 0.1)
corrected_p_values <- p.adjust(p_values, method='holm')

cfu_raw[cfu_raw == 0.0] <- NA
cefoperazone <- subset(cfu_raw, treatment == 'cefoperazone')$cfu_spore
clindamycin <- subset(cfu_raw, treatment == 'clindamycin')$cfu_spore
streptomycin <- subset(cfu_raw, treatment == 'streptomycin')$cfu_spore
germfree <- subset(cfu_raw, treatment == 'germfree')$cfu_spore
rm(cfu_raw)

wilcox.test(germfree, cefoperazone, exact=F) # p-value = 0.02341 *
wilcox.test(germfree, clindamycin, exact=F) # p-value = 0.0005699 ***
wilcox.test(germfree, streptomycin, exact=F) # p-value = 0.001764 **
wilcox.test(clindamycin, cefoperazone, exact=F) # n.s.
wilcox.test(streptomycin, cefoperazone, exact=F) # n.s.
wilcox.test(clindamycin, cefoperazone, exact=F) # n.s.
rm(germfree, cefoperazone, clindamycin, streptomycin)

p_values <- c(0.02341, 0.0005699, 0.001764, 0.1, 0.1, 0.1)
p.adjust(p_values, method='holm')
rm(p_values)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_1.pdf'
select_palette <- c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black')
pdf(file=plot_file, width=12, height=6)
layout(matrix(c(1,2,
                3,4), 
              nrow=2, ncol=2, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# A.  Timeline of mouse experiments

# Create an empty plot
par(mar=c(1,2,1,1))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-4.8,4), ylim=c(-2,5))

# Abx in drinking water timeline
rect(xleft=-4, ybottom=2.8, xright=0, ytop=3.2, col='darkorchid2', border='black')
Arrows(x0=-4, y0=3, x1=3.5, y1=3, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,0,2,2.75), y0=c(3.5,3.5,3.5,3.5), x1=c(-4,0,2,2.75), y1=c(2.5,2.5,2.5,2.5), lwd=4)
segments(x0=c(-4,-3,-2,-1,1), y0=c(3.25,3.25,3.25,3.25,3.25), x1=c(-4,-3,-2,-1,1), y1=c(2.75,2.75,2.75,2.75,2.75), lwd=2)
points(x=c(2,2.75), y=c(4,4), pch=25, bg=c('white','black'), col='black', cex=2.5)
text(x=c(-4,0,2,2.75), y=c(2.2,2.2,2.2,2.2), c('Day -7', 'Day -2', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.6, y=3.2, 'Cefoperazone', cex=0.7)
text(x=-4.6, y=2.95, 'or', font=2, cex=0.8)
text(x=-4.6, y=2.7, 'Streptomycin', cex=0.7)

# IP injection abx timeline
Arrows(x0=-4, y0=0, x1=-1.5, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,-3,-2.25), y0=c(-0.5,-0.5,-0.5), x1=c(-4,-3,-2.25), y1=c(0.5,0.5,0.5), lwd=4)
points(x=c(-4,-3,-2.25), y=c(1,1,1), pch=c(25,25,25), bg=c('gray60','white','black'), col='black', cex=2.5)
text(x=c(-4,-3,-2.25), y=c(-0.8,-0.8,-0.8), c('Day -1', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.6, y=0, 'Clindamycin', cex=0.7)

# Legend
legend(x=0, y=1.3, legend=expression('Antibiotic in Drinking Water', 'IP Injection of Antibiotic', paste(italic('C. difficile'), ' Spore Gavage'), 'Euthanize & Necropsy'), 
       pt.bg=c('darkorchid2','gray60','white','black'), pch=c(22,25,25,25), pt.cex=c(2.5,2,2,2), bty='n')

# Plot label
mtext('A', side=2, line=2, las=2, adj=-0.4, padj=-7.1, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# B.  Vegetative cell CFU
par(las=1, mar=c(2,4,1,1), mgp=c(2.5,0.7,0), yaxs='i')
stripchart(cfu_vegetative~treatment, data=vegetative_cfu, vertical=T, pch=1, lwd=2,
           ylim=c(1,9), xaxt='n', yaxt='n', cex=1.5, col=select_palette,
           ylab='Vegetative CFU/g Cecal Content', method='jitter', jitter=0.25, cex.lab=0.9)
axis(side=1, at=c(1:5), c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Gnotobiotic', 'No Antibiotics'), tick = FALSE)
labelsY <- c(0, parse(text=paste(rep(10,8), '^', seq(2,9,1), sep='')))
axis(side=2, at=c(1:9), labelsY, tick=TRUE)
abline(h=2, col="black", lty=2, lwd=1.5)

# Draw axis break
axis.break(2, 1.5, style='slash')

# Draw median
segments(0.6, vege_medians[1], 1.4, vege_medians[1], lwd=3) # cefoperazone
segments(1.6, vege_medians[2], 2.4, vege_medians[2], lwd=3) # streptomycin
segments(2.6, vege_medians[3], 3.4, vege_medians[3], lwd=3) # clindamycin
segments(3.6, vege_medians[4], 4.4, vege_medians[4], lwd=3) # germfree
segments(4.6, vege_medians[5], 5.4, vege_medians[5], lwd=3) # conventional

mtext('B', side=2, line=2, las=2, adj=1.5, padj=-6.7, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# C.  Spore CFU
par(las=1, mar=c(2,4,1,1), mgp=c(2.5,0.7,0), yaxs='i')
stripchart(cfu_spore~treatment, data=spore_cfu, vertical=T, pch=1, lwd=2, 
           ylim=c(1,9), xaxt='n', yaxt='n', cex=1.5, col=select_palette,
           ylab='Spore CFU/g Cecal Content', method='jitter', jitter=0.25)
axis(side=1, at=c(1:5), c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Gnotobiotic', 'No Antibiotics'), tick = FALSE)
axis(side=2, at=c(1:9), labelsY, tick=TRUE)
abline(h=2, col="black", lty=2, lwd=1.5)

# Draw axis break
axis.break(2, 1.5, style='slash') 

# Draw median
segments(0.6, spore_medians[1], 1.4, spore_medians[1], lwd=3) # cefoperazone
segments(1.6, spore_medians[2], 2.4, spore_medians[2], lwd=3) # streptomycin
segments(2.6, spore_medians[3], 3.4, spore_medians[3], lwd=3) # clindamycin
segments(3.6, spore_medians[4], 4.4, spore_medians[4], lwd=3) # germfree
segments(4.6, spore_medians[5], 5.4, spore_medians[5], lwd=3) # conventional

# Adding significance to plot
segments(x0=c(1,2,3), y0=c(7,7.5,8), x1=c(4,4,4), y1=c(7,7.5,8), lwd=2)
text(c(2.5,3,3.5), c(7.2,7.7,8.2), labels=c('*','*','*'), cex=1.5)

mtext('C', side=2, line=2, las=2, adj=1.5, padj=-6.7, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# D.  Toxin data
par(las=1, mar=c(2,4,1,1), mgp=c(2.5,0.7,0), xpd=FALSE, yaxs='i')
stripchart(titer~treatment, data=toxin, vertical=T, pch=1, lwd=2,
           ylim=c(1.5,3.5), xlim=c(0.5,5.5), xaxt='n', yaxt='n', col=select_palette,
           cex=1.5, ylab=expression(paste('Toxin Titer/g Cecal Content (',Log[10],')')), method='jitter', jitter=0.25)
axis(side=1, at=c(1:5), c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Gnotobiotic', 'No Antibiotics'), tick=FALSE)
axis(side=2, at=c(1.5,2.0,2.5,3.0,3.5), labels=c('0','2.0','2.5','3.0','3.5'))

# Draw axis break
axis.break(2, 1.75, style='slash') 

# Draw limit of detection
abline(h=2, lty=2, lwd=1.5)

# Draw median
segments(0.6, toxin_medians[1], 1.4, toxin_medians[1], lwd=3) # cefoperazone
segments(1.6, toxin_medians[2], 2.4, toxin_medians[2], lwd=3) # streptomycin
segments(2.6, toxin_medians[3], 3.4, toxin_medians[3], lwd=3) # clindamycin
segments(3.6, toxin_medians[4], 4.4, toxin_medians[4], lwd=3) # germfree
segments(4.6, toxin_medians[5], 5.4, toxin_medians[5], lwd=3) # conventional

# Adding significance to plot
segments(x0=c(1,2,3), y0=c(3.2,3.3,3.4), x1=c(4,4,4), y1=c(3.2,3.3,3.4), lwd=2)
text(c(2.5,3,3.5), c(3.25,3.35,3.45), labels=c('*','*','*'), cex=1.5)

# Plot label
mtext('D', side=2, line=2, las=2, adj=2, padj=-6.7, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
dev.off()
rm(labelsY, plot_file, toxin_medians, spore_medians, vege_medians, spore_cfu, toxin, vegetative_cfu, select_palette)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

