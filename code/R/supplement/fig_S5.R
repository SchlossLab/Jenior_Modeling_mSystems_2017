
# Read in score distribution file
scores <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/C01094_distribution.txt', header=FALSE)[,1]

pdf(file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S5.pdf', width=7, height=5)

# Create a curve for the sample distribution
par(mar=c(4,4,1,1))
score_hist <- hist(scores, breaks=50, main='', xlab='Metabolite Importance Score', ylab='Score Frequency', xlim=c(-7,7), ylim=c(0,500), axes=FALSE, border='white')
box()
axis(1, at=seq(-9,9,3), labels=seq(-9,9,3))
axis(2, at=seq(0,500,100), labels=seq(0,500,100), las=1)

xfit <- seq(min(scores), max(scores), length=10000)
yfit <- dnorm(xfit, mean=mean(scores), sd=sd(scores)) * diff(score_hist$mids[1:2]) * length(scores) 
lines(xfit, yfit, lwd=2)

# Calulacte condfidence
score_median <- quantile(scores)[3]
lower_iqr <- quantile(scores)[2]
upper_iqr <- quantile(scores)[4]
lower_95 <- score_median - abs(1.7 * (lower_iqr / sqrt(length(scores))))
upper_95 <- score_median + abs(1.7 * (upper_iqr / sqrt(length(scores))))
lower_99 <- score_median - abs(1.95 * (lower_iqr / sqrt(length(scores))))
upper_99 <- score_median + abs(1.95 * (upper_iqr / sqrt(length(scores))))

# Summary lines
abline(v=score_median, lwd=2, col='black') # Median
abline(v=c(lower_99,upper_99), lty=2, lwd=2, col='red') # 0.99

# Actual score for fructose 1-phosphate
arrows(x0=5.946, y0=250, x1=5.946, y1=175, col='blue', length=0.2, angle=20, lwd=4)
text(x=5.946, y=270, 'Measured Score')

dev.off()
