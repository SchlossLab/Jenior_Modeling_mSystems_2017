
# Read in and format data - spores vs vegetative %
cfu_spore_only <- read.delim('/Users/mljenior/Desktop/repositories/Jenior_Transcriptomics_2015.docx/data/raw/composite_%_sporeVvege.tsv', sep='\t', header=T)
cfu_spore_only$mouse <- NULL
cfu_spore_only$cage <- NULL
colnames(cfu_spore_only) <- c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Gnotobiotic')


# Plot the cfus 
par(las=1)
stripchart(cfu_spore_only, vertical=T, method='jitter', jitter=.1, pch=16, ylim=c(0,10), yaxt='n', xaxt='n', cex=1.4, ylab=substitute(paste('Percent Spores of Cultureable ', italic('C. diffile'), ' 630')))
axis(side = 1, at = 1:4, colnames(cfu_spore_only), tick = FALSE, font=2)
labelsY=c('0%', '', '2%', '', '4%', '', '6%', '', '8%', '', '10%')
axis(side = 2, at = 0:10, labelsY, tick = TRUE)
abline(h=0, col='black', lty=2)
abline(h=2, col='black', lty=2)
abline(h=4, col='black', lty=2)
abline(h=6, col='black', lty=2)
abline(h=8, col='black', lty=2)
abline(h=10, col='black', lty=2)

# Spore vs spore
wilcox.test(cfu_sub_spore$Cefoperazone, cfu_sub_spore$Clindamycin, exact=F) # W = 31, p-value = 0.4241  NS
wilcox.test(cfu_sub_spore$Cefoperazone, cfu_sub_spore$Streptomycin, exact=F) # W = 33.5, p-value = 0.5617  NS
wilcox.test(cfu_sub_spore$Cefoperazone, cfu_sub_spore$Gnotobiotic, exact=F) # W = 5, p-value = 0.001892  **
wilcox.test(cfu_sub_spore$Clindamycin, cfu_sub_spore$Streptomycin, exact=F) # W = 41, p-value = 1  NS
wilcox.test(cfu_sub_spore$Clindamycin, cfu_sub_spore$Gnotobiotic, exact=F) # W = 1, p-value = 0.0005699  ***
wilcox.test(cfu_sub_spore$Streptomycin, cfu_sub_spore$Gnotobiotic, exact=F) # W = 3, p-value = 0.001086  **


# Plot medians
segments(x0=0.7, x1=1.3, y0=median(cfu_spore_only$Cefoperazone), y1=median(cfu_spore_only$Cefoperazone), col="black", lwd=3) # cefoperazone median
segments(x0=1.7, x1=2.3, y0=median(cfu_spore_only$Clindamycin), y1=median(cfu_spore_only$Clindamycin), col="black", lwd=3) # clindamycin median
segments(x0=2.7, x1=3.3, y0=median(cfu_spore_only$Streptomycin), y1=median(cfu_spore_only$Streptomycin), col="black", lwd=3) # streptomycin median
segments(x0=3.7, x1=4.3, y0=median(cfu_spore_only$Gnotobiotic), y1=median(cfu_spore_only$Gnotobiotic), col="black", lwd=3) # gnotobiotic median

# Lines denoting significance
segments(x0=1, x1=4, y0=7, y1=7, col="black", lwd=2) # Gnoto vs Cef
segments(x0=1, x1=1, y0=7, y1=6.6, col="black", lwd=2)
segments(x0=4, x1=4, y0=7, y1=6.6, col="black", lwd=2)
text(2.5, 7.25, labels=substitute(paste(italic('p'), ' = 0.002')), cex=0.8)

segments(x0=2, x1=4, y0=8, y1=8, col="black", lwd=2) # Gnoto vs Clinda
segments(x0=2, x1=2, y0=8, y1=7.6, col="black", lwd=2)
segments(x0=4, x1=4, y0=8, y1=7.6, col="black", lwd=2)
text(3, 8.25, labels=substitute(paste(italic('p'), ' < 0.001')), cex=0.8)

segments(x0=3, x1=4, y0=9, y1=9, col="black", lwd=2) # Gnoto vs Strep
segments(x0=3, x1=3, y0=9, y1=8.6, col="black", lwd=2)
segments(x0=4, x1=4, y0=9, y1=8.6, col="black", lwd=2)
text(3.5, 9.25, labels=substitute(paste(italic('p'), ' = 0.001')), cex=0.8)





