
# Define variables
toxin_file <- '/Users/schloss/Desktop/Jenior_812/data/toxin_titer.dat'
figure_file <- '/Users/schloss/Desktop/toxin.pdf'

# Read in toxin data
raw_toxin <- read.delim(toxin_file, sep='\t', header=T)
raw_toxin$mouse <- NULL
raw_toxin$cage <- NULL
colnames(raw_toxin) <- c('Cefoperazone','Clindamycin','Streptomycin','Germfree')

# Calculate quantiles
cef_iqr <- quantile(raw_toxin$Cefoperazone) - 1.5
strep_iqr <- quantile(raw_toxin$Streptomycin) - 1.5
clinda_iqr <- quantile(raw_toxin$Clindamycin) - 1.5
gf_iqr <- quantile(raw_toxin$Germfree) - 1.5
median_vector <- c(cef_iqr[3],strep_iqr[3],clinda_iqr[3],gf_iqr[3])


# Plot the data
pdf(file=figure_file, width=10, height=7)
par(las=1, mar=c(4,3.5,1,1), mgp=c(2.5,0.7,0), xpd=FALSE)
barplot(median_vector, space=0.5, col='forestgreen', ylim=c(0,1.75), yaxt='n', xaxt='n', ylab='Toxin Titer (log10)')
axis(side=1, at=c(1,2.5,4,5.5), c('Cefoperazone','Streptomycin','Clindamycin','Germfree'), tick = FALSE, font=2, cex.axis=1.4)
mtext(c('0.5 mg/ml DW', '0.5 mg/ml DW', '10 mg/kg IP', ''), 
      side=1, at=c(1,2.5,4,5.5), cex=0.9, padj=3.5)
axis(side=2, at=seq(0,1.75,0.25), c('1.50','1.75','2.00','2.25','2.50','2.75','3.00','3.25'), tick=TRUE, cex.axis=0.9)
abline(h=0.5, col='black', lty=2)


# Draw IQR
segments(0.7, cef_iqr[2], 1.3, cef_iqr[2], lwd=3) # cefoperazone
segments(0.7, cef_iqr[4], 1.3, cef_iqr[4], lwd=3)
segments(1, cef_iqr[2], 1, cef_iqr[4], lwd=3)
segments(2.2, strep_iqr[2], 2.8, strep_iqr[2], lwd=3) # streptomycin
segments(2.2, strep_iqr[4], 2.8, strep_iqr[4], lwd=3)
segments(2.5, strep_iqr[2], 2.5, strep_iqr[4], lwd=3)
segments(3.7, clinda_iqr[2], 4.3, clinda_iqr[2], lwd=3) # clindamycin
segments(3.7, clinda_iqr[4], 4.3, clinda_iqr[4], lwd=3)
segments(4, clinda_iqr[2], 4, clinda_iqr[4], lwd=3)
segments(5.2, gf_iqr[2], 5.8, gf_iqr[2], lwd=3) # germfree
segments(5.2, gf_iqr[4], 5.8, gf_iqr[4], lwd=3)
segments(5.5, gf_iqr[2], 5.5, gf_iqr[4], lwd=3)


# Adding significance to plot
text(5.5, 1.6, labels='***', cex=2.5, font=2)
dev.off()

# Test for normal distribution
shapiro.test(raw_toxin$Cefoperazone) # p-value = 0.003751
shapiro.test(raw_toxin$Streptomycin) # p-value = 1.155e-05
shapiro.test(raw_toxin$Clindamycin) # p-value = 0.0001526
# Germfre not possible because all values are identical


# Compare groups - use rank test for non-normallly distributed data
wilcox.test(raw_toxin$Cefoperazone, raw_toxin$Clindamycin, exact=F) # W = 42, p-value = 0.9189  NS
wilcox.test(raw_toxin$Cefoperazone, raw_toxin$Streptomycin, exact=F) # W = 56.5, p-value = 0.1242  NS
wilcox.test(raw_toxin$Cefoperazone, raw_toxin$Gnotobiotic, exact=F) # W = 4.5, p-value = 0.000476  ***
wilcox.test(raw_toxin$Clindamycin, raw_toxin$Streptomycin, exact=F) # W = 59, p-value = 0.07892  NS
wilcox.test(raw_toxin$Clindamycin, raw_toxin$Gnotobiotic, exact=F) # W = 0, p-value = 0.0001152  ***
wilcox.test(raw_toxin$Streptomycin, raw_toxin$Gnotobiotic, exact=F) # W = 4.5, p-value = 0.0003599  ***

