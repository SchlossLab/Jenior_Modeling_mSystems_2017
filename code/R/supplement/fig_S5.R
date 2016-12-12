
pdf(file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S5.pdf', width=14, height=5)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow = TRUE))


# Read in score distribution file
all_scores <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/test_distribution.txt', header=FALSE)
fructose <- unique(all_scores[,1])
propanoate <- unique(all_scores[,2])

# Create a curve for the sample distribution
par(mar=c(4,4,1,1))
scores <- fructose
score_hist <- hist(scores, breaks=10, main='', xlab='Metabolite Importance Score', ylab='Score Frequency', 
                   xlim=c(-9,9), ylim=c(0,500), las=1)
box()

# Calculate stats
all_quantile = quantile(scores)
score_median = all_quantile[3]
iqr = all_quantile[4] - all_quantile[2]
numerator = 1.25 * iqr
denominator = 1.35 * sqrt(length(scores))
range_factor = numerator / denominator
range_95 = 1.6 * range_factor
lower_95 = score_median - range_95
upper_95 = score_median + range_95

# Summary lines
abline(v=score_median, lwd=2, col='black') # Median
abline(v=c(lower_95,upper_95), lty=2, lwd=2, col='red')

# Actual score for fructose 1-phosphate
arrows(x0=5.946, y0=75, x1=5.946, y1=50, col='blue', length=0.2, angle=20, lwd=4)
mtext('a', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.1, font=2)


# Create a curve for the sample distribution
par(mar=c(4,4,1,1))
scores <- propanoate
score_hist <- hist(scores, breaks=10, main='', xlab='Metabolite Importance Score', ylab='Score Frequency', 
                   xlim=c(-9,9), ylim=c(0,500), las=1)
box()

# Calculate stats
all_quantile = quantile(scores)
score_median = all_quantile[3]
iqr = all_quantile[4] - all_quantile[2]
numerator = 1.25 * iqr
denominator = 1.35 * sqrt(length(scores))
range_factor = numerator / denominator
range_95 = 1.6 * range_factor
lower_95 = score_median - range_95
upper_95 = score_median + range_95

# Summary lines
abline(v=score_median, lwd=2, col='black') # Median
abline(v=c(lower_95,upper_95), lty=2, lwd=2, col='red')

# Actual score for propanoate
arrows(x0=5.946, y0=75, x1=5.946, y1=50, col='blue', length=0.2, angle=20, lwd=4)
mtext('b', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.1, font=2)

dev.off()
