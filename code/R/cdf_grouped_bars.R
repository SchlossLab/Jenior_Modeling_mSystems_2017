
# File variables
cef_expression_file <- '/Users/pschloss/Desktop/Jenior_812/cefoperazone.cdf_genes.txt'
strep_expression_file <- '/Users/pschloss/Desktop/Jenior_812/streptomycin.cdf_genes.txt'
clinda_expression_file <- '/Users/pschloss/Desktop/Jenior_812/clindamycin.cdf_genes.txt'
gf_expression_file <- '/Users/pschloss/Desktop/Jenior_812/germfree.cdf_genes.txt'
figure_file <- '/Users/pschloss/Desktop/cdf_expression.pdf'
legend_file <- '/Users/pschloss/Desktop/cdf_expression_legend.pdf'

# Read in expression data
cef_expression <- read.delim(cef_expression_file, header = T, sep='\t', row.names=1)
strep_expression <- read.delim(strep_expression_file, header = T, sep='\t', row.names=1)
clinda_expression <- read.delim(clinda_expression_file, header = T, sep='\t', row.names=1)
gf_expression <- read.delim(gf_expression_file, header = T, sep='\t', row.names=1)

# Pool within a gene catagory
pool_cef_expression <- aggregate(cefoperazone ~ family, FUN = sum, data = cef_expression)
pool_strep_expression <- aggregate(streptomycin ~ family, FUN = sum, data = strep_expression)
pool_clinda_expression <- aggregate(clindamycin ~ family, FUN = sum, data = clinda_expression)
pool_gf_expression <- aggregate(germfree ~ family, FUN = sum, data = gf_expression)

# Merge datasets
all_pool_expression <- merge(pool_cef_expression, pool_strep_expression, by='family')
all_pool_expression <- merge(all_pool_expression, pool_clinda_expression, by='family')
all_pool_expression <- merge(all_pool_expression, pool_gf_expression, by='family')
rownames(all_pool_expression) <- all_pool_expression$family
all_pool_expression$family <- NULL

# Transform merged data
log_all_pool_expression <- all_pool_expression
log_all_pool_expression[log_all_pool_expression  == 0] <- 1
log_all_pool_expression <- log10(log_all_pool_expression)
rownames(log_all_pool_expression) <- c('Adhesion', 'Antibiotic resitance', 'Glycan metabolism', 'Sulfur metabolism', 'Flagella', 'Global regulators', 'Iron uptake', 'PaLoc genes', 'Sporulation', 'AA metabolism')

# Plot the transformed data
pdf(file=figure_file, width=14, height=10)
par(las=3, mar=c(9,5,1,1))
barplot(t(log_all_pool_expression), col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), 
        beside=TRUE, ylim=c(0.001, 4.5), ylab='Log10 Transcript Abundance', 
        cex.lab=1.8, yaxt='n', font=2)
box()
abline(h=seq(1:4), lty=3, lwd=2)
axis(2, c(1:4), las=1, cex.axis=1.5)
dev.off()

pdf(file=legend_file, width=5, height=5)
plot(10, xlim=c(-10,10), ylim=c(-10,10))
legend('center', legend=c('Cefoperzone', 'Streptomycin', 'Clindamycin', 'Germfree'), 
       col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), pch=15, cex=2, pt.cex=3.5)
dev.off()