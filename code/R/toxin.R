
# Read in toxin data
raw_toxin <- read.delim('/Users/mljenior/Desktop/repositories/Jenior_Transcriptomics_2015.docx/data/raw/toxin_data_2015.txt', sep='\t', header=T)
raw_toxin$mouse <- NULL
raw_toxin$cage <- NULL
colnames(raw_toxin) <- c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Germfree')


# Plot the data
stripchart(raw_toxin, vertical=T, method='jitter', jitter=.1, pch=16, ylim=c(2,3.3), xaxt='n', cex=1.4, ylab=substitute(paste(italic('C. diffile'), ' 630 Toxin Titer')))
axis(side = 1, at = 1:4, colnames(raw_toxin), tick = FALSE, font=2)
abline(h=2, col='black', lty=2) # Limit of detection


# Plot medians
segments(x0=0.7, x1=1.3, y0=median(raw_toxin$Cefoperazone), y1=median(raw_toxin$Cefoperazone), col="black", lwd=3) # cefoperazone
segments(x0=1.7, x1=2.3, y0=median(raw_toxin$Clindamycin), y1=median(raw_toxin$Clindamycin), col="black", lwd=3) # clindamycin
segments(x0=2.7, x1=3.3, y0=median(raw_toxin$Streptomycin), y1=median(raw_toxin$Streptomycin), col="black", lwd=3) # streptomycin
segments(x0=3.7, x1=4.3, y0=median(raw_toxin$Gnotobiotic), y1=median(raw_toxin$Gnotobiotic), col="black", lwd=3) # gnotobiotic


# Compare groups
wilcox.test(raw_toxin$Cefoperazone, raw_toxin$Clindamycin, exact=F) # W = 42, p-value = 0.9189  NS
wilcox.test(raw_toxin$Cefoperazone, raw_toxin$Streptomycin, exact=F) # W = 56.5, p-value = 0.1242  NS
wilcox.test(raw_toxin$Cefoperazone, raw_toxin$Gnotobiotic, exact=F) # W = 4.5, p-value = 0.000476  ***
wilcox.test(raw_toxin$Clindamycin, raw_toxin$Streptomycin, exact=F) # W = 59, p-value = 0.07892  NS
wilcox.test(raw_toxin$Clindamycin, raw_toxin$Gnotobiotic, exact=F) # W = 0, p-value = 0.0001152  ***
wilcox.test(raw_toxin$Streptomycin, raw_toxin$Gnotobiotic, exact=F) # W = 4.5, p-value = 0.0003599  ***


# Adding significance to plot - only significant
segments(x0=1, x1=1, y0=3.04, y1=3.08, col="black", lwd=2) # Cefoperazone vs Gnotobiotic
segments(x0=1, x1=4, y0=3.08, y1=3.08, col="black", lwd=2)
segments(x0=4, x1=4, y0=3.04, y1=3.08, col="black", lwd=2)
text(2.5, 3.1, labels=substitute(paste(italic('p'), ' < 0.0005')), cex=0.8)

segments(x0=2, x1=2, y0=3.14, y1=3.18, col="black", lwd=2) # Clindamycin vs Gnotobiotic
segments(x0=2, x1=4, y0=3.18, y1=3.18, col="black", lwd=2)
segments(x0=4, x1=4, y0=3.14, y1=3.18, col="black", lwd=2)
text(3, 3.2, labels=substitute(paste(italic('p'), ' < 0.0005')), cex=0.8)

segments(x0=3, x1=3, y0=3.24, y1=3.28, col="black", lwd=2) # Streptomycin vs Gnotobiotic
segments(x0=3, x1=4, y0=3.28, y1=3.28, col="black", lwd=2)
segments(x0=4, x1=4, y0=3.24, y1=3.28, col="black", lwd=2)
text(3.5, 3.3, labels=substitute(paste(italic('p'), ' < 0.0005')), cex=0.8)





