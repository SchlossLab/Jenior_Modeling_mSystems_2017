


deps <- c('wesanderson', 'shape');
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

fox <- wes_palette("FantasticFox")
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_1.pdf'
cfu_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Transcriptomics_2015/data/bench_assays/cfu.dat'


# Read, format, and separate the data
cfu_raw <- read.delim(cfu_file, sep='\t', header=T)

cfu <- cfu_raw
cfu[cfu == 0] <- 100
cfu$cfu_vegetative <- log10(cfu$cfu_vegetative)
cfu$cfu_spore <- log10(cfu$cfu_spore)
cfu$mouse <- NULL
cfu <- subset(cfu, cage < 4)
cfu$cage <- NULL
cfu <- subset(cfu, treatment != 'conventional')

vegetative_cfu <- cfu
vegetative_cfu$cfu_spore <- NULL
vegetative_cfu <- droplevels(vegetative_cfu)
spore_cfu <- cfu
spore_cfu$cfu_vegetative <- NULL
spore_cfu <- droplevels(spore_cfu)

vegetative_spore_cfu <- cbind(vegetative_cfu, spore_cfu)
colnames(vegetative_spore_cfu) <- c('group', 'vegetative', 'treatment', 'spores')
vegetative_spore_cfu$treatment <- NULL
vegetative_spore_cfu_median <- ddply(vegetative_spore_cfu,~group,summarise,vegetative=median(vegetative),spores=median(spores))
rownames(vegetative_spore_cfu_median) <- vegetative_spore_cfu_median$group
vegetative_spore_cfu_median$group <- NULL
vegetative_spore_cfu_median <- t(as.matrix(vegetative_spore_cfu_median))

cef_iqr_vege <- quantile(vegetative_cfu[vegetative_cfu$treatment == 'cefoperazone', 2])
clinda_iqr_vege <- quantile(vegetative_cfu[vegetative_cfu$treatment == 'clindamycin', 2])
strep_iqr_vege <- quantile(vegetative_cfu[vegetative_cfu$treatment == 'streptomycin', 2])
gf_iqr_vege <- quantile(vegetative_cfu[vegetative_cfu$treatment == 'germfree', 2])

cef_iqr_spore <- quantile(spore_cfu[spore_cfu$treatment == 'cefoperazone', 2])
clinda_iqr_spore <- quantile(spore_cfu[spore_cfu$treatment == 'clindamycin', 2])
strep_iqr_spore <- quantile(spore_cfu[spore_cfu$treatment == 'streptomycin', 2])
gf_iqr_spore <- quantile(spore_cfu[spore_cfu$treatment == 'germfree', 2])








toxin_file <- '/Users/pschloss/Desktop/Repositories/Jenior_Transcriptomics_2015/data/raw/toxin_titer.dat'

# Read in toxin data
raw_toxin <- read.delim(toxin_file, sep='\t', header=T)
raw_toxin$mouse <- NULL
raw_toxin$cage <- NULL
raw_toxin$treatment <- factor(raw_toxin$treatment, levels=c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree'))
raw_toxin[,2] <- raw_toxin[,2] - 1.5

# Calculate quantiles
cef_iqr <- quantile(raw_toxin[raw_toxin$treatment == 'Cefoperazone', 2])
strep_iqr <- quantile(raw_toxin[raw_toxin$treatment == 'Streptomycin', 2])
clinda_iqr <- quantile(raw_toxin[raw_toxin$treatment == 'Clindamycin', 2])
gf_iqr <- quantile(raw_toxin[raw_toxin$treatment == 'Germfree', 2])

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
pdf(file=plot_file, width=15, height=9)
layout(matrix(c(1,1,
                2,3), 
              nrow=2, ncol=2, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 1A.  Timeline of mouse experiments

# Create an empty plot
par(mar=c(0,0,0,0))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-4.75,4), ylim=c(-2,5))

# Abx in drinking water timeline
rect(xleft=-4, ybottom=2.8, xright=0, ytop=3.2, col=fox[3], border='black')
Arrows(x0=-4, y0=3, x1=3.5, y1=3, lwd=4, arr.type='triangle', arr.length=0.75, arr.width=0.4)
segments(x0=c(-4,0,2,2.75), y0=c(3.5,3.5,3.5,3.5), x1=c(-4,0,2,2.75), y1=c(2.5,2.5,2.5,2.5), lwd=4)
segments(x0=c(-4,-3,-2,-1,1), y0=c(3.25,3.25,3.25,3.25,3.25), x1=c(-4,-3,-2,-1,1), y1=c(2.75,2.75,2.75,2.75,2.75), lwd=2)
points(x=c(2,2.75), y=c(4,4), pch=25, bg=c('white','black'), col='black', cex=2.5)
text(x=c(-4,0,2,2.75), y=c(2.2,2.2,2.2,2.2), c('Day -7', 'Day -2', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.5, y=3.2, 'Cefoperazone', cex=0.8)
text(x=-4.5, y=2.95, 'or', font=2)
text(x=-4.5, y=2.7, 'Streptomycin', cex=0.8)

# IP injection abx timeline
Arrows(x0=-4, y0=0, x1=-1.5, y1=0, lwd=4, arr.type='triangle', arr.length=0.75, arr.width=0.4)
segments(x0=c(-4,-3,-2.25), y0=c(-0.5,-0.5,-0.5), x1=c(-4,-3,-2.25), y1=c(0.5,0.5,0.5), lwd=4)
points(x=c(-4,-3,-2.25), y=c(1,1,1), pch=c(25,25,25), bg=c(fox[5],'white','black'), col='black', cex=2.5)
text(x=c(-4,-3,-2.25), y=c(-0.8,-0.8,-0.8), c('Day -1', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.5, y=0, 'Clindamycin', cex=0.8)

# Legend
legend(x=0, y=1.3, legend=expression('Antibiotic in Drinking Water', 'IP Injection of Antibiotic',paste(italic('C. difficile'), ' Spore Gavage'), 'Sacrifice & Necropsy'), 
       pt.bg=c(fox[3],fox[5],'white','black'), pch=c(22,25,25,25), cex=2, pt.cex=c(4.2,3.2,3.2,3.2), bty='n')

# Plot label
legend('topleft', legend='A', cex=2, bty='n')

#-------------------------------------------------------------------------------------------------------------------------------------#

# Grouped bar plot of both vegetative cells and spores

par(las=1, mar=c(3.5,4,1,1), mgp=c(2.5,0.7,0))
barplot(vegetative_spore_cfu_median, col=c("blue3","firebrick"), beside=TRUE, ylim=c(0,9), xaxt='n', yaxt='n', ylab='CFU per gram cecal content', cex.lab=1.3)
box()
axis(side=1, at=c(2,5,8,11), c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree'), tick = FALSE, font=2, cex.axis=1.4)
mtext(c('0.5 mg/ml DW', '5 mg/ml DW', '10 mg/kg IP', 'NA'), side=1, at=c(2,5,8,11), cex=0.9, padj=3.8)
labelsY <- parse(text=paste(rep(10,7), '^', seq(1,9,1), sep=''))
axis(side=2, at=c(1:9), labelsY, tick=TRUE, cex.axis=1.2)
abline(h=2, col="black", lty=2, lwd=2)
legend('topleft', legend=c('Vegetative cells', 'Spores'), col=c("blue3","firebrick"), pch=15, cex=1, pt.cex=2)
text(11.5, gf_iqr_spore[4]+0.3, '**', font=2, cex=2)

segments(1.2, cef_iqr_vege[4], 1.8, cef_iqr_vege[4], lwd=4)
segments(1.2, cef_iqr_vege[2], 1.8, cef_iqr_vege[2], lwd=4)
segments(1.5, cef_iqr_vege[4], 1.5, cef_iqr_vege[2], lwd=4)
segments(2.2, cef_iqr_spore[4], 2.8, cef_iqr_spore[4], lwd=4)
segments(2.2, cef_iqr_spore[2], 2.8, cef_iqr_spore[2], lwd=4) 
segments(2.5, cef_iqr_spore[4], 2.5, cef_iqr_spore[2], lwd=4)

segments(4.2, strep_iqr_vege[4], 4.8, strep_iqr_vege[4], lwd=4)
segments(4.2, strep_iqr_vege[2], 4.8, strep_iqr_vege[2], lwd=4)
segments(4.5, strep_iqr_vege[4], 4.5, strep_iqr_vege[2], lwd=4)
segments(5.2, strep_iqr_spore[4], 5.8, strep_iqr_spore[4], lwd=4)
segments(5.2, strep_iqr_spore[2], 5.8, strep_iqr_spore[2], lwd=4) 
segments(5.5, strep_iqr_spore[4], 5.5, strep_iqr_spore[2], lwd=4)

segments(7.2, clinda_iqr_vege[4], 7.8, clinda_iqr_vege[4], lwd=4)
segments(7.2, clinda_iqr_vege[2], 7.8, clinda_iqr_vege[2], lwd=4)
segments(7.5, clinda_iqr_vege[4], 7.5, clinda_iqr_vege[2], lwd=4)
segments(8.2, clinda_iqr_spore[4], 8.8, clinda_iqr_spore[4], lwd=4)
segments(8.2, clinda_iqr_spore[2], 8.8, clinda_iqr_spore[2], lwd=4) 
segments(8.5, clinda_iqr_spore[4], 8.5, clinda_iqr_spore[2], lwd=4)

segments(10.2, gf_iqr_vege[4], 10.8, gf_iqr_vege[4], lwd=4)
segments(10.2, gf_iqr_vege[2], 10.8, gf_iqr_vege[2], lwd=4)
segments(10.5, gf_iqr_vege[4], 10.5, gf_iqr_vege[2], lwd=4)
segments(11.2, gf_iqr_spore[4], 11.8, gf_iqr_spore[4], lwd=4)
segments(11.2, gf_iqr_spore[2], 11.8, gf_iqr_spore[2], lwd=4) 
segments(11.5, gf_iqr_spore[4], 11.5, gf_iqr_spore[2], lwd=4)

# Plot label
legend('topleft', legend='B', cex=2, bty='n')

#-------------------------------------------------------------------------------------------------------------------------------------#


# Plot the data
par(las=1, mar=c(4,4,1,1), mgp=c(2.5,0.7,0), xpd=FALSE)
stripchart(titer~treatment, data=raw_toxin, vertical=T, pch=21, ylim=c(0.5,1.75), xlim=c(0.5,4.5), yaxt='n', xaxt='n', cex=3, bg='forestgreen', ylab='', method='jitter', jitter=0.25)
mtext('Toxin Titer (log10)', side=2, at=1.125, cex=1.5, padj=-3, las=0)
axis(side=1, at=c(1:4), c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree'), tick = FALSE, font=2, cex.axis=1.5)
mtext(c('0.5 mg/ml DW', '5 mg/ml DW', '10 mg/kg IP', ''), side=1, at=c(1:4), cex=1, padj=3.5)
axis(side=2, at=c(0.5, 0.75, 1, 1.25, 1.5, 1.75), c(2, 2.25, 2.5, 2.75, 3, 3.25), tick=TRUE, cex=1.2)

# Draw limit of detection
abline(h=0.5, lty=2, lwd=3)

# Draw median
segments(0.6, cef_iqr[3], 1.4, cef_iqr[3], lwd=8) # cefoperazone
segments(1.6, strep_iqr[3], 2.4, strep_iqr[3], lwd=8) # streptomycin
segments(2.6, clinda_iqr[3], 3.4, clinda_iqr[3], lwd=8) # clindamycin
segments(3.6, gf_iqr[3], 4.4, gf_iqr[3], lwd=8) # germfree

# Adding significance to plot
text(4, gf_iqr[4] + 0.1, labels='***', cex=3, font=2)


# Plot label
legend('topleft', legend='C', cex=2, bty='n')


dev.off()

