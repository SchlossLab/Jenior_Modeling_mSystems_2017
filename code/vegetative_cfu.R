

# Read in and format data - vegetative cells
cfu_vege <- read.delim('~/Desktop/composite_cfu_vegetative.tsv', sep='\t', header=T)
cfu_sub_vege <- subset(cfu_vege, cage < 4)
cfu_sub_vege$mouse <- NULL
cfu_sub_vege$cage <- NULL
colnames(cfu_sub_vege) <- c('Conventional', 'Cefoperazone', 'Clindamycin', 'Streptomycin', 'Gnotobiotic')
cfu_sub_vege_log <- log10(cfu_sub_vege)

# Plot the cfus 
par(las=1)
stripchart(cfu_sub_vege_log, vertical=T, method='jitter', jitter=.1, pch=16, ylim=c(1,10), xaxt='n', yaxt='n', cex=1.4)
axis(side = 1, at = 1:5, colnames(cfu_sub_vege), tick = FALSE, font=2)
labelsY=parse(text=paste(rep(10,10), "^", seq(1,10,1), sep=""))
axis(side = 2, at = 1:10, labelsY, tick = TRUE)
abline(h=2, col='black', lty=2) # Limit of detection


# Plot median and iqr
segments(x0=.7, x1=1.3, y0=median(cfu_sub_vege_log$Conventional), y1=median(cfu_sub_vege_log$Conventional), col="black", lwd=2) # conventional median
segments(x0=.8, x1=1.2, y0=quantile(cfu_sub_vege_log$Conventional)[4], y1=quantile(cfu_sub_vege_log$Conventional)[4], col="black", lwd=2) # conventional iqr 1
segments(x0=.8, x1=1.2, y0=quantile(cfu_sub_vege_log$Conventional)[2], y1=quantile(cfu_sub_vege_log$Conventional)[2], col="black", lwd=2) # conventional iqr 2
segments(x0=1, x1=1, y0=quantile(cfu_sub_vege_log$Conventional)[4], y1=quantile(cfu_sub_vege_log$Conventional)[2], col="black", lwd=2)

segments(x0=1.7, x1=2.3, y0=median(cfu_sub_vege_log$Cefoperazone), y1=median(cfu_sub_vege_log$Cefoperazone), col="black", lwd=2) # cefoperazone median
segments(x0=1.8, x1=2.2, y0=quantile(cfu_sub_vege_log$Cefoperazone)[4], y1=quantile(cfu_sub_vege_log$Cefoperazone)[4], col="black", lwd=2) # cefoperazone iqr 1
segments(x0=1.8, x1=2.2, y0=quantile(cfu_sub_vege_log$Cefoperazone)[2], y1=quantile(cfu_sub_vege_log$Cefoperazone)[2], col="black", lwd=2) # cefoperazone iqr 2
segments(x0=2, x1=2, y0=quantile(cfu_sub_vege_log$Cefoperazone)[4], y1=quantile(cfu_sub_vege_log$Cefoperazone)[2], col="black", lwd=2)

segments(x0=2.7, x1=3.3, y0=median(cfu_sub_vege_log$Clindamycin), y1=median(cfu_sub_vege_log$Clindamycin), col="black", lwd=2) # clindamycin median
segments(x0=2.8, x1=3.2, y0=quantile(cfu_sub_vege_log$Clindamycin)[4], y1=quantile(cfu_sub_vege_log$Clindamycin)[4], col="black", lwd=2) # clindamycin iqr 1
segments(x0=2.8, x1=3.2, y0=quantile(cfu_sub_vege_log$Clindamycin)[2], y1=quantile(cfu_sub_vege_log$Clindamycin)[2], col="black", lwd=2) # clindamycin iqr 2
segments(x0=3, x1=3, y0=quantile(cfu_sub_vege_log$Clindamycin)[4], y1=quantile(cfu_sub_vege_log$Clindamycin)[2], col="black", lwd=2)

segments(x0=3.7, x1=4.3, y0=median(cfu_sub_vege_log$Streptomycin), y1=median(cfu_sub_vege_log$Streptomycin), col="black", lwd=2) # streptomycin median
segments(x0=3.8, x1=4.2, y0=quantile(cfu_sub_vege_log$Streptomycin)[4], y1=quantile(cfu_sub_vege_log$Streptomycin)[4], col="black", lwd=2) # streptomycin iqr 1
segments(x0=3.8, x1=4.2, y0=quantile(cfu_sub_vege_log$Streptomycin)[2], y1=quantile(cfu_sub_vege_log$Streptomycin)[2], col="black", lwd=2) # streptomycin iqr 2
segments(x0=4, x1=4, y0=quantile(cfu_sub_vege_log$Streptomycin)[4], y1=quantile(cfu_sub_vege_log$Streptomycin)[2], col="black", lwd=2)

segments(x0=4.7, x1=5.3, y0=median(cfu_sub_vege_log$Gnotobiotic), y1=median(cfu_sub_vege_log$Gnotobiotic), col="black", lwd=2) # germfree median
segments(x0=4.8, x1=5.2, y0=quantile(cfu_sub_vege_log$Gnotobiotic)[4], y1=quantile(cfu_sub_vege_log$Gnotobiotic)[4], col="black", lwd=2) # germfree iqr 1
segments(x0=4.8, x1=5.2, y0=quantile(cfu_sub_vege_log$Gnotobiotic)[2], y1=quantile(cfu_sub_vege_log$Gnotobiotic)[2], col="black", lwd=2) # germfree iqr 2
segments(x0=5, x1=5, y0=quantile(cfu_sub_vege_log$Gnotobiotic)[4], y1=quantile(cfu_sub_vege_log$Gnotobiotic)[2], col="black", lwd=2)


# Calculate significant differences
wilcox.test(cfu_sub_vege$Cefoperazone, cfu_sub_vege$Clindamycin, exact=F) # W = 53, p-value = 0.2868  NS
wilcox.test(cfu_sub_vege$Cefoperazone, cfu_sub_vege$Streptomycin, exact=F) # W = 33.5, p-value = 0.563  NS
wilcox.test(cfu_sub_vege$Cefoperazone, cfu_sub_vege$Gnotobiotic, exact=F) # W = 12, p-value = 0.01197  *
wilcox.test(cfu_sub_vege$Clindamycin, cfu_sub_vege$Streptomycin, exact=F) # W = 35.5, p-value = 0.6908  NS
wilcox.test(cfu_sub_vege$Clindamycin, cfu_sub_vege$Gnotobiotic, exact=F) # W = 7, p-value = 0.0035  **
wilcox.test(cfu_sub_vege$Streptomycin, cfu_sub_vege$Gnotobiotic, exact=F) # W = 33.5, p-value = 0.5644  NS


# Adding significance to plot
segments(x0=2, x1=2, y0=quantile(cfu_sub_vege_log$Cefoperazone)[4]+1.3, y1=quantile(cfu_sub_vege_log$Cefoperazone)[4]+1.5, col="black", lwd=2) # Cef vs gnoto
segments(x0=2, x1=5, y0=quantile(cfu_sub_vege_log$Cefoperazone)[4]+1.5, y1=quantile(cfu_sub_vege_log$Cefoperazone)[4]+1.5, col="black", lwd=2)
segments(x0=5, x1=5, y0=quantile(cfu_sub_vege_log$Cefoperazone)[4]+1.3, y1=quantile(cfu_sub_vege_log$Cefoperazone)[4]+1.5, col="black", lwd=2)
text(3.5, quantile(cfu_sub_vege_log$Cefoperazone)[4]+1.65, labels=substitute(paste(italic('p'), ' = 0.01197')), cex=0.8)

segments(x0=3, x1=3, y0=quantile(cfu_sub_vege_log$Cefoperazone)[4]+2.0, y1=quantile(cfu_sub_vege_log$Cefoperazone)[4]+2.2, col="black", lwd=2) # Clinda vs gnoto
segments(x0=3, x1=5, y0=quantile(cfu_sub_vege_log$Cefoperazone)[4]+2.2, y1=quantile(cfu_sub_vege_log$Cefoperazone)[4]+2.2, col="black", lwd=2)
segments(x0=5, x1=5, y0=quantile(cfu_sub_vege_log$Cefoperazone)[4]+2.0, y1=quantile(cfu_sub_vege_log$Cefoperazone)[4]+2.2, col="black", lwd=2)
text(4, quantile(cfu_sub_vege_log$Cefoperazone)[4]+2.35, labels=substitute(paste(italic('p'), ' = 0.0035')), cex=0.8)



