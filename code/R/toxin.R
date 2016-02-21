
# Define variables
toxin_file <- '/Users/pschloss/Desktop/Repositories/Jenior_Transcriptomics_2015/data/raw/toxin_titer.dat'
figure_file <- '/Users/pschloss/Desktop/toxin.pdf'

# Read in toxin data
raw_toxin <- read.delim(toxin_file, sep='\t', header=T)
raw_toxin$mouse <- NULL
raw_toxin$cage <- NULL
raw_toxin$treatment <- factor(raw_toxin$treatment, levels=c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Germfree'))
raw_toxin[,2] <- raw_toxin[,2] - 1.5

# Calculate quantiles
cef_iqr <- quantile(raw_toxin[raw_toxin$treatment == 'Cefoperazone', 2])
clinda_iqr <- quantile(raw_toxin[raw_toxin$treatment == 'Clindamycin', 2])
strep_iqr <- quantile(raw_toxin[raw_toxin$treatment == 'Streptomycin', 2])
gf_iqr <- quantile(raw_toxin[raw_toxin$treatment == 'Germfree', 2])

# Initialize PDF
pdf(file=figure_file, width=10, height=7)

# Plot the data
par(las=1, mar=c(4,4,1,1), mgp=c(2.5,0.7,0), xpd=FALSE)
stripchart(titer~treatment, data=raw_toxin, vertical=T, method='jitter', jitter=.1, pch=21, ylim=c(0.5,1.75), yaxt='n', xaxt='n', cex=2, bg='forestgreen', ylab='')
mtext('Toxin Titer (log10)', side=2, at=1.125, cex=1.5, padj=-3, las=0)
axis(side=1, at=c(1:4), c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Germfree'), tick = FALSE, font=2, cex.axis=1.5)
mtext(c('0.5 mg/ml DW', '10 mg/kg IP', '5 mg/ml DW', ''), side=1, at=c(1:4), cex=1, padj=3.5)
axis(side=2, at=c(0.5, 0.75, 1, 1.25, 1.5, 1.75), c(2, 2.25, 2.5, 2.75, 3, 3.25), tick=TRUE, cex=1.2)

# Draw limit of detection
abline(h=0.5, lty=3, lwd=3)

# Draw median
segments(0.8, cef_iqr[3], 1.2, cef_iqr[3], lwd=5) # cefoperazone
segments(1.8, clinda_iqr[3], 2.2, clinda_iqr[3], lwd=5) # clindamycin
segments(2.8, strep_iqr[3], 3.2, strep_iqr[3], lwd=5) # streptomycin
segments(3.8, gf_iqr[3], 4.2, gf_iqr[3], lwd=5) # germfree

# Adding significance to plot
text(4, gf_iqr[4] + 0.1, labels='***', cex=2.5, font=2)

# Close PDF
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

