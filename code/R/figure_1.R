deps <- c('wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

cfu_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/cfu.dat'
toxin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/toxin_titer.dat'

# Read in data
cfu <- read.delim(cfu_file, sep='\t', header=T)
toxin <- read.delim(toxin_file, sep='\t', header=T)

# Format CFU data and collect summary statistics
cfu[cfu == 0] <- 100
cfu$cfu_vegetative <- log10(cfu$cfu_vegetative)
cfu$cfu_spore <- log10(cfu$cfu_spore)
cfu$mouse <- NULL
cfu <- subset(cfu, cage < 4 ) # Remove uninfected controls
cfu$cage <- NULL
cfu$treatment <- factor(cfu$treatment, levels=c('cefoperazone', 'streptomycin', 'clindamycin', 'germfree', 'conventional'))
vegetative_cfu <- cfu
vegetative_cfu$cfu_spore <- NULL
spore_cfu <- cfu
spore_cfu$cfu_vegetative <- NULL
rm(cfu_file, cfu)

# Format toxin data and find summary statistics
toxin$mouse <- NULL
toxin$cage <- NULL
toxin$treatment <- factor(toxin$treatment, levels=c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree', 'Conventional'))
cef <- as.numeric(median(toxin[toxin$treatment == 'Cefoperazone', 2]))
strep <- as.numeric(median(toxin[toxin$treatment == 'Streptomycin', 2]))
clinda <- as.numeric(median(toxin[toxin$treatment == 'Clindamycin', 2]))
gf <- as.numeric(median(toxin[toxin$treatment == 'Germfree', 2]))
conv <- as.numeric(median(toxin[toxin$treatment == 'Conventional', 2]))
toxin_medians <- c(cef, strep, clinda, gf, conv)
rm(cef, strep, clinda, gf, conv, toxin_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_1.pdf'
fox <- wes_palette('FantasticFox')
pdf(file=plot_file, width=12, height=6)
layout(matrix(c(1,2,
                3,4), 
              nrow=2, ncol=2, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 1A.  Timeline of mouse experiments

# Create an empty plot
par(mar=c(1,2,1,1))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-4.75,4), ylim=c(-2,5))

# Abx in drinking water timeline
rect(xleft=-4, ybottom=2.8, xright=0, ytop=3.2, col=fox[3], border='black')
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
points(x=c(-4,-3,-2.25), y=c(1,1,1), pch=c(25,25,25), bg=c(fox[5],'white','black'), col='black', cex=2.5)
text(x=c(-4,-3,-2.25), y=c(-0.8,-0.8,-0.8), c('Day -1', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.6, y=0, 'Clindamycin', cex=0.7)

# Legend
legend(x=0, y=1.3, legend=expression('Antibiotic in Drinking Water', 'IP Injection of Antibiotic',paste(italic('C. difficile'), ' Spore Gavage'), 'Sacrifice & Necropsy'), 
       pt.bg=c(fox[3],fox[5],'white','black'), pch=c(22,25,25,25), pt.cex=c(2.5,2,2,2), bty='n')

# Plot label
mtext('A', side=2, line=2, las=2, adj=-0.3, padj=-6.5, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Plot toxin data
par(las=1, mar=c(2,4,1,1), mgp=c(2.5,0.7,0), xpd=FALSE)
stripchart(titer~treatment, data=toxin, vertical=T, pch=20, 
           ylim=c(0,3.5), xlim=c(0.5,5.5), xaxt='n', 
           cex=2, col='black', ylab='Toxin Titer (log10)', method='jitter', jitter=0.25)
axis(side=1, at=c(1:5), c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Gnotobiotic', 'Conventional'), tick=FALSE)

# Draw limit of detection
abline(h=2, lty=2, lwd=1.5)

# Draw median
segments(0.6, toxin_medians[1], 1.4, toxin_medians[1], lwd=3) # cefoperazone
segments(1.6, toxin_medians[2], 2.4, toxin_medians[2], lwd=3) # streptomycin
segments(2.6, toxin_medians[3], 3.4, toxin_medians[3], lwd=3) # clindamycin
segments(3.6, toxin_medians[4], 4.4, toxin_medians[4], lwd=3) # germfree
segments(4.6, toxin_medians[5], 5.4, toxin_medians[5], lwd=3) # conventional

# Adding significance to plot
text(4, toxin_medians[4] + 0.1, labels='***', cex=2, font=2)

# Plot label
mtext('B', side=2, line=2, las=2, adj=2, padj=-6.2, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Plot formatted data - vegetative
par(las=1, mar=c(2,4,1,1), mgp=c(2.5,0.7,0))
boxplot(cfu_vegetative~treatment, data=vegetative_cfu, 
        ylim=c(1,9), xaxt='n', yaxt='n', boxlwd=3, outwex=2, whisklwd=3,
        staplelwd=3, outline=FALSE, range=2, ylab='Vegetative CFU/g Cecal Content')
axis(side=1, at=c(1:5), c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Gnotobiotic', 'Conventional'), tick = FALSE)
labelsY <- parse(text=paste(rep(10,9), '^', seq(1,9,1), sep=''))
axis(side=2, at=c(1:9), labelsY, tick=TRUE)
abline(h=2, col="black", lty=2, lwd=1.5)
mtext('C', side=2, line=2, las=2, adj=1.5, padj=-6.7, cex=1.5)

# Plot formatted data - spores
par(las=1, mar=c(2,4,1,1), mgp=c(2.5,0.7,0))
boxplot(cfu_spore~treatment, data=spore_cfu, col='gray64', 
        ylim=c(1,9), xaxt='n', yaxt='n', boxlwd=3, outwex=2, whisklwd=3, 
        staplelwd=3, outline=FALSE, range=2, ylab='Spore CFU/g Cecal Content')
axis(side=1, at=c(1:5), c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Gnotobiotic', 'Conventional'), tick = FALSE)
labelsY <- parse(text=paste(rep(10,9), '^', seq(1,9,1), sep=''))
axis(side=2, at=c(1:9), labelsY, tick=TRUE)
abline(h=2, col="black", lty=2, lwd=1.5)
text(4, 6.8, '***', cex=2, font=2)
mtext('D', side=2, line=2, las=2, adj=1.2, padj=-6.7, cex=1.5)


dev.off()



