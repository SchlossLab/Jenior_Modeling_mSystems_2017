
# This script plots spores only quantifications across all conditions

# Read in and format data - spores
cfu_spore <- read.delim('/Users/mljenior/Desktop/repositories/Jenior_Transcriptomics_2015.docx/data/raw/composite_cfu_spore.tsv', sep='\t', header=T)
cfu_sub_spore <- subset(cfu_spore, cage < 4)
cfu_sub_spore$mouse <- NULL
cfu_sub_spore$cage <- NULL
colnames(cfu_sub_spore) <- c('Conventional', 'Cefoperazone', 'Clindamycin', 'Streptomycin', 'Germfree')
cfu_sub_spore_log <- log10(cfu_sub_spore)

# Plot the cfus 
par(las=1)
stripchart(cfu_sub_spore_log, vertical=T, method='jitter', jitter=.1, pch=16, ylim=c(1,8), xaxt='n', yaxt='n', cex=1.4, ylab='Spores per gram cecal content')
axis(side = 1, at = 1:5, colnames(cfu_sub_spore), tick = FALSE, font=2)
labelsY=parse(text=paste(rep(10,8), "^", seq(1,8,1), sep=""))
axis(side = 2, at = 1:8, labelsY, tick = TRUE)
abline(h=2, col='black', lty=2) # Limit of detection

# Plot medians and IQRs
segments(x0=0.7, x1=1.3, y0=median(cfu_sub_spore_log$Conventional), y1=median(cfu_sub_spore_log$Conventional), col="black", lwd=2) 

segments(x0=1.7, x1=2.3, y0=median(cfu_sub_spore_log$Cefoperazone), y1=median(cfu_sub_spore_log$Cefoperazone), col="black", lwd=2) 
segments(x0=1.8, x1=2.2, y0=quantile(cfu_sub_spore_log$Cefoperazone)[4], y1=quantile(cfu_sub_spore_log$Cefoperazone)[4], col="black", lwd=2) 
segments(x0=1.8, x1=2.2, y0=quantile(cfu_sub_spore_log$Cefoperazone)[2], y1=quantile(cfu_sub_spore_log$Cefoperazone)[2], col="black", lwd=2) 
segments(x0=2, x1=2, y0=quantile(cfu_sub_spore_log$Cefoperazone)[4], y1=quantile(cfu_sub_spore_log$Cefoperazone)[2], col="black", lwd=2)

segments(x0=2.7, x1=3.3, y0=median(cfu_sub_spore_log$Clindamycin), y1=median(cfu_sub_spore_log$Clindamycin), col="black", lwd=2) 
segments(x0=2.8, x1=3.2, y0=quantile(cfu_sub_spore_log$Clindamycin)[4], y1=quantile(cfu_sub_spore_log$Clindamycin)[4], col="black", lwd=2) 
segments(x0=2.8, x1=3.2, y0=quantile(cfu_sub_spore_log$Clindamycin)[2], y1=quantile(cfu_sub_spore_log$Clindamycin)[2], col="black", lwd=2)
segments(x0=3, x1=3, y0=quantile(cfu_sub_spore_log$Clindamycin)[4], y1=quantile(cfu_sub_spore_log$Clindamycin)[2], col="black", lwd=2)

segments(x0=3.7, x1=4.3, y0=median(cfu_sub_spore_log$Streptomycin), y1=median(cfu_sub_spore_log$Streptomycin), col="black", lwd=2) 
segments(x0=3.8, x1=4.2, y0=quantile(cfu_sub_spore_log$Streptomycin)[4], y1=quantile(cfu_sub_spore_log$Streptomycin)[4], col="black", lwd=2) 
segments(x0=3.8, x1=4.2, y0=quantile(cfu_sub_spore_log$Streptomycin)[2], y1=quantile(cfu_sub_spore_log$Streptomycin)[2], col="black", lwd=2) 
segments(x0=4, x1=4, y0=quantile(cfu_sub_spore_log$Streptomycin)[4], y1=quantile(cfu_sub_spore_log$Streptomycin)[2], col="black", lwd=2)

segments(x0=4.7, x1=5.3, y0=median(cfu_sub_spore_log$Germfree), y1=median(cfu_sub_spore_log$Germfree), col="black", lwd=2)
segments(x0=4.8, x1=5.2, y0=quantile(cfu_sub_spore_log$Germfree)[4], y1=quantile(cfu_sub_spore_log$Germfree)[4], col="black", lwd=2) 
segments(x0=4.8, x1=5.2, y0=quantile(cfu_sub_spore_log$Germfree)[2], y1=quantile(cfu_sub_spore_log$Germfree)[2], col="black", lwd=2) 
segments(x0=5, x1=5, y0=quantile(cfu_sub_spore_log$Germfree)[4], y1=quantile(cfu_sub_spore_log$Germfree)[2], col="black", lwd=2)


# Calculate differences
wilcox.test(cfu_sub_spore$Cefoperazone, cfu_sub_spore$Clindamycin, exact=F) # W = 31, p-value = 0.4241
wilcox.test(cfu_sub_spore$Cefoperazone, cfu_sub_spore$Streptomycin, exact=F) # W = 33.5, p-value = 0.5617
wilcox.test(cfu_sub_spore$Cefoperazone, cfu_sub_spore$Germfree, exact=F) # V = 45, p-value = 0.008551
wilcox.test(cfu_sub_spore$Clindamycin, cfu_sub_spore$Streptomycin, exact=F) # W = 41, p-value = 1
wilcox.test(cfu_sub_spore$Clindamycin, cfu_sub_spore$Germfree, exact=F) # V = 45, p-value = 0.009091
wilcox.test(cfu_sub_spore$Streptomycin, cfu_sub_spore$CGermfree, exact=F) # V = 45, p-value = 0.009152



