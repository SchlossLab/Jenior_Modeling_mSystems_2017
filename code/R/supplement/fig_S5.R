
pdf(file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S5.pdf', width=14, height=5)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow = TRUE))


# Read in score distribution file
all_scores <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/test_distribution.txt', header=FALSE)
metabolite_1 <- unique(all_scores[,1])
metabolite_1_score <- 6.329
metabolite_2 <- unique(all_scores[,2])
metabolite_2_score <- -4.297

# Create a curve for the sample distribution
par(mar=c(4,4,1,1))
d <- density(metabolite_1)
plot(d, type='n', main='', xlab='Metabolite Importance Score', ylab='Score Frequency (Log10)', las=1, ylim=c(0,0.12))
polygon(d, col='gray')

# Calculate stats
all_quantile = quantile(metabolite_1)
score_median = all_quantile[3]
iqr = as.numeric(all_quantile[4] - all_quantile[2])
numerator = 1.25 * iqr
denominator = 1.35 * sqrt(length(metabolite_1))
range_factor = numerator / denominator
range_95 = 1.6 * range_factor
lower_95 = score_median - range_95
upper_95 = score_median + range_95

# Summary lines
abline(v=score_median, lwd=2, col='black') # Median
abline(v=c(lower_95,upper_95), lty=2, lwd=2, col='red')

# Actual score for sorbitol
arrows(x0=metabolite_1_score, y0=0.11, x1=metabolite_1_score, y1=0.089, col='blue', length=0.2, angle=20, lwd=3)
mtext('a', side=2, line=2, las=2, adj=3, padj=-15, cex=1.1, font=2)

#-----------------------------#

# Create a curve for the sample distribution
par(mar=c(4,4,1,1))
d <- density(metabolite_2)
plot(d, type='n', main='', xlab='Metabolite Importance Score', ylab='Score Frequency (Log10)', las=1, ylim=c(0,0.12))
polygon(d, col='gray')

# Calculate stats
all_quantile = quantile(metabolite_2)
score_median = all_quantile[3]
iqr = all_quantile[4] - all_quantile[2]
numerator = 1.25 * iqr
denominator = 1.35 * sqrt(length(metabolite_2))
range_factor = numerator / denominator
range_95 = 1.6 * range_factor
lower_95 = score_median - range_95
upper_95 = score_median + range_95

# Summary lines
abline(v=score_median, lwd=2, col='black') # Median
abline(v=c(lower_95,upper_95), lty=2, lwd=2, col='red')

# Actual score for aspartate
arrows(x0=metabolite_2_score, y0=0.095, x1=metabolite_2_score, y1=0.075, col='blue', length=0.2, angle=20, lwd=3)
mtext('b', side=2, line=2, las=2, adj=3, padj=-15, cex=1.1, font=2)

dev.off()
