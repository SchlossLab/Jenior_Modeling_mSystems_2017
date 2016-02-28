
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
plot(0, xlim=c(-5,5), ylim=c(-5,5), axes=FALSE, frame.plot=FALSE, xaxt='n', yaxt='n', xlab='', ylab='', type='n')
legend('center', legend=c('Cefoperzone', 'Streptomycin', 'Clindamycin', 'Germfree'), 
       col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), pch=15, cex=2, pt.cex=3.5)
dev.off()






#--------------------------------------------------------------------------------------------------------------------------#

# Each expression group seperately

# File variables
cef_expression_file <- '/Users/pschloss/Desktop/Jenior_812/cefoperazone.cdf_genes.txt'
strep_expression_file <- '/Users/pschloss/Desktop/Jenior_812/streptomycin.cdf_genes.txt'
clinda_expression_file <- '/Users/pschloss/Desktop/Jenior_812/clindamycin.cdf_genes.txt'
gf_expression_file <- '/Users/pschloss/Desktop/Jenior_812/germfree.cdf_genes.txt'
#figure_file <- '/Users/pschloss/Desktop/cdf_expression.pdf'
#legend_file <- '/Users/pschloss/Desktop/cdf_expression_legend.pdf'


# Read in expression data
cef_expression <- read.delim(cef_expression_file, header = T, sep='\t', row.names=1)
strep_expression <- read.delim(strep_expression_file, header = T, sep='\t', row.names=1)
clinda_expression <- read.delim(clinda_expression_file, header = T, sep='\t', row.names=1)
gf_expression <- read.delim(gf_expression_file, header = T, sep='\t', row.names=1)


# Merge and format the data
all_expression <- cbind(cef_expression, strep_expression$streptomycin, clinda_expression$clindamycin, gf_expression$germfree)
colnames(all_expression) <- c('annotation', 'family', 'cefoperazone', 'streptomycin', 'clindamycin', 'germfree')


# Subset by gene family
adhesion_all_expression <- subset(all_expression, family == 'Adhesion')
flagella_all_expression <- subset(all_expression, family == 'Flagella')
sporulation_all_expression <- subset(all_expression, family == 'Sporulation')
antibiotic_all_expression <- subset(all_expression, family == 'Antibiotic_resitance')
global_all_expression <- subset(all_expression, family == 'Global_regulators')
aminoacid_all_expression <- subset(all_expression, family == 'Stickland_fermentation_/_Amino_acid_metabolism')
carbohydrate_all_expression <- subset(all_expression, family == 'Carbohydrate_metabolism_/_Fermentation')
iron_all_expression <- subset(all_expression, family == 'Iron_uptake')
cysteine_all_expression <- subset(all_expression, family == 'Cysteine/Sulfur_metabolism')
paloc_all_expression <- subset(all_expression, family == 'PaLoc_genes')


# Format and transform each new table
adhesion_all_expression$annotation <- NULL
adhesion_all_expression$family <- NULL
adhesion_all_expression[adhesion_all_expression  == 0] <- 1
adhesion_all_expression <- log10(adhesion_all_expression)

flagella_all_expression$annotation <- NULL
flagella_all_expression$family <- NULL
flagella_all_expression[flagella_all_expression  == 0] <- 1
flagella_all_expression <- log10(flagella_all_expression)

sporulation_all_expression$annotation <- NULL
sporulation_all_expression$family <- NULL
sporulation_all_expression[sporulation_all_expression  == 0] <- 1
sporulation_all_expression <- log10(sporulation_all_expression)

antibiotic_all_expression$annotation <- NULL
antibiotic_all_expression$family <- NULL
antibiotic_all_expression[antibiotic_all_expression  == 0] <- 1
antibiotic_all_expression <- log10(antibiotic_all_expression)

global_all_expression$annotation <- NULL
global_all_expression$family <- NULL
global_all_expression[global_all_expression  == 0] <- 1
global_all_expression <- log10(global_all_expression)

aminoacid_all_expression$annotation <- NULL
aminoacid_all_expression$family <- NULL
aminoacid_all_expression[aminoacid_all_expression  == 0] <- 1
aminoacid_all_expression <- log10(aminoacid_all_expression)

carbohydrate_all_expression$annotation <- NULL
carbohydrate_all_expression$family <- NULL
carbohydrate_all_expression[carbohydrate_all_expression  == 0] <- 1
carbohydrate_all_expression <- log10(carbohydrate_all_expression)

iron_all_expression$annotation <- NULL
iron_all_expression$family <- NULL
iron_all_expression[iron_all_expression  == 0] <- 1
iron_all_expression <- log10(iron_all_expression)

cysteine_all_expression$annotation <- NULL
cysteine_all_expression$family <- NULL
cysteine_all_expression[cysteine_all_expression  == 0] <- 1
cysteine_all_expression <- log10(cysteine_all_expression)

paloc_all_expression$annotation <- NULL
paloc_all_expression$family <- NULL
paloc_all_expression[paloc_all_expression  == 0] <- 1
paloc_all_expression <- log10(paloc_all_expression)


# Plot each new group seperately
pdf(file='/Users/pschloss/Desktop/adhesion_expr.pdf', width=14, height=10)
par(las=3, mar=c(5,5,3,1))
barplot(t(adhesion_all_expression), col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), 
        beside=TRUE, ylim=c(0.001, 4.5), ylab='Log10 Transcript Abundance', 
        cex.lab=1.8, yaxt='n', font=2, main='Adhesion')
box()
abline(h=seq(1:4), lty=3, lwd=2)
axis(2, c(1:4), las=1, cex.axis=1.5)
dev.off()
#-------------------------#
pdf(file='/Users/pschloss/Desktop/flagella_expr.pdf', width=14, height=10)
par(las=3, mar=c(5,5,3,1))
barplot(t(flagella_all_expression), col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), 
        beside=TRUE, ylim=c(0.001, 4.5), ylab='Log10 Transcript Abundance', 
        cex.lab=1.8, yaxt='n', font=2, main='Flagella')
box()
abline(h=seq(1:4), lty=3, lwd=2)
axis(2, c(1:4), las=1, cex.axis=1.5)
dev.off()
#-------------------------#
pdf(file='/Users/pschloss/Desktop/sporulation_expr.pdf', width=14, height=10)
par(las=3, mar=c(5,5,3,1))
barplot(t(sporulation_all_expression), col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), 
        beside=TRUE, ylim=c(0.001, 4.5), ylab='Log10 Transcript Abundance', 
        cex.lab=1.8, yaxt='n', font=2, main='Sporulation')
box()
abline(h=seq(1:4), lty=3, lwd=2)
axis(2, c(1:4), las=1, cex.axis=1.5)
dev.off()
#-------------------------#
pdf(file='/Users/pschloss/Desktop/antibiotic_expr.pdf', width=14, height=10)
par(las=3, mar=c(5,5,3,1))
barplot(t(antibiotic_all_expression), col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), 
        beside=TRUE, ylim=c(0.001, 4.5), ylab='Log10 Transcript Abundance', 
        cex.lab=1.8, yaxt='n', font=2, main='Antibiotic Resistance')
box()
abline(h=seq(1:4), lty=3, lwd=2)
axis(2, c(1:4), las=1, cex.axis=1.5)
dev.off()
#-------------------------#
pdf(file='/Users/pschloss/Desktop/global_expr.pdf', width=14, height=10)
par(las=3, mar=c(5,5,3,1))
barplot(t(global_all_expression), col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), 
        beside=TRUE, ylim=c(0.001, 4.5), ylab='Log10 Transcript Abundance', 
        cex.lab=1.8, yaxt='n', font=2, main='Global Regulators')
box()
abline(h=seq(1:4), lty=3, lwd=2)
axis(2, c(1:4), las=1, cex.axis=1.5)
dev.off()
#-------------------------#
pdf(file='/Users/pschloss/Desktop/aminoacid_expr.pdf', width=14, height=10)
par(las=3, mar=c(5,5,3,1))
barplot(t(aminoacid_all_expression), col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), 
        beside=TRUE, ylim=c(0.001, 4.5), ylab='Log10 Transcript Abundance', 
        cex.lab=1.8, yaxt='n', font=2, main='Stickland Fermentation/Amino Acid Metabolism')
box()
abline(h=seq(1:4), lty=3, lwd=2)
axis(2, c(1:4), las=1, cex.axis=1.5)
dev.off()
#-------------------------#
pdf(file='/Users/pschloss/Desktop/carbohydrate_expr.pdf', width=14, height=10)
par(las=3, mar=c(5,5,3,1))
barplot(t(carbohydrate_all_expression), col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), 
        beside=TRUE, ylim=c(0.001, 4.5), ylab='Log10 Transcript Abundance', 
        cex.lab=1.8, yaxt='n', font=2, main='Glycolysis / Carbohydrate Fermentation')
box()
abline(h=seq(1:4), lty=3, lwd=2)
axis(2, c(1:4), las=1, cex.axis=1.5)
dev.off()
#-------------------------#
pdf(file='/Users/pschloss/Desktop/iron_expr.pdf', width=14, height=10)
par(las=3, mar=c(5,5,3,1))
barplot(t(iron_all_expression), col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), 
        beside=TRUE, ylim=c(0.001, 4.5), ylab='Log10 Transcript Abundance', 
        cex.lab=1.8, yaxt='n', font=2, main='Iron Uptake')
box()
abline(h=seq(1:4), lty=3, lwd=2)
axis(2, c(1:4), las=1, cex.axis=1.5)
dev.off()
#-------------------------#
pdf(file='/Users/pschloss/Desktop/cysteine_expr.pdf', width=14, height=10)
par(las=3, mar=c(5,5,3,1))
barplot(t(cysteine_all_expression), col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), 
        beside=TRUE, ylim=c(0.001, 4.5), ylab='Log10 Transcript Abundance', 
        cex.lab=1.8, yaxt='n', font=2, main='Cysteine/Sulfur Metabolism')
box()
abline(h=seq(1:4), lty=3, lwd=2)
axis(2, c(1:4), las=1, cex.axis=1.5)
dev.off()
#-------------------------#
pdf(file='/Users/pschloss/Desktop/paloc_expr.pdf', width=14, height=10)
par(las=3, mar=c(5,5,3,1))
barplot(t(paloc_all_expression), col=c('firebrick2', 'blue2', 'chartreuse4', 'goldenrod2'), 
        beside=TRUE, ylim=c(0.001, 4.5), ylab='Log10 Transcript Abundance', 
        cex.lab=1.8, yaxt='n', font=2, main='PaLoc Genes')
box()
abline(h=seq(1:4), lty=3, lwd=2)
axis(2, c(1:4), las=1, cex.axis=1.5)
dev.off()


