
# Load dependencies
deps <- c('wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Function for population variance of columns in matrix
pop_var <- function(data) {
  vars <- c()
  for (x in 1:ncol(data)){
    vars[x] <- sum((data[,x] - mean(data[,x]))^2) / length(data[,x])
  }
  return(vars)
}

# Define files
metadata <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'
metabolome <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/metabolomics.tsv'
shared <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'

#----------------------------------------#

# Read data
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination
metabolome <- read.delim(metabolome, sep='\t', header=T, row.names=1)
metabolome <- metabolome[, !colnames(metabolome) %in% c('CefC5M2','StrepC4M1')] # Remove possible contamination
shared <- read.delim(shared, sep='\t', header=T, row.names=2)
shared <- shared[!rownames(shared) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination

# Format data
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$type <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome <- as.data.frame(t(metabolome))
metabolites <- colnames(metabolome)
shared$label <- NULL
shared$numOtus <- NULL
shared <- log10(shared + 10)
otus <- colnames(shared)

# Merge datasets
shared <- merge(metadata, shared, by='row.names')
rownames(shared) <- shared$Row.names
shared$Row.names <- NULL
metabolome <- merge(metadata, metabolome, by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
rm(metadata)

# Subset datasets
strep <- subset(metabolome, abx == 'streptomycin')
strep$abx <- NULL
strep_metabolome_mock <- subset(strep, infection == 'mock')
strep_metabolome_mock$infection <- NULL
strep_metabolome_630 <- subset(strep, infection == '630')
strep_metabolome_630$infection <- NULL
cef <- subset(metabolome, abx == 'cefoperazone')
cef$abx <- NULL
cef_metabolome_mock <- subset(cef, infection == 'mock')
cef_metabolome_mock$infection <- NULL
cef_metabolome_630 <- subset(cef, infection == '630')
cef_metabolome_630$infection <- NULL
clinda <- subset(metabolome, abx == 'clindamycin')
clinda$abx <- NULL
clinda_metabolome_mock <- subset(clinda, infection == 'mock')
clinda_metabolome_mock$infection <- NULL
clinda_metabolome_630 <- subset(clinda, infection == '630')
clinda_metabolome_630$infection <- NULL
conv <- subset(metabolome, abx == 'none')
conv$abx <- NULL
conv_metabolome_mock <- subset(conv, infection == 'mock')
conv_metabolome_mock$infection <- NULL
rm(metabolome)
strep <- subset(shared, abx == 'streptomycin')
strep$abx <- NULL
strep_shared_mock <- subset(strep, infection == 'mock')
strep_shared_mock$infection <- NULL
strep_shared_630 <- subset(strep, infection == '630')
strep_shared_630$infection <- NULL
cef <- subset(shared, abx == 'cefoperazone')
cef$abx <- NULL
cef_shared_mock <- subset(cef, infection == 'mock')
cef_shared_mock$infection <- NULL
cef_shared_630 <- subset(cef, infection == '630')
cef_shared_630$infection <- NULL
clinda <- subset(shared, abx == 'clindamycin')
clinda$abx <- NULL
clinda_shared_mock <- subset(clinda, infection == 'mock')
clinda_shared_mock$infection <- NULL
clinda_shared_630 <- subset(clinda, infection == '630')
clinda_shared_630$infection <- NULL
conv <- subset(shared, abx == 'none')
conv$abx <- NULL
conv_shared_mock <- subset(conv, infection == 'mock')
conv_shared_mock$infection <- NULL
rm(shared)
rm(strep, cef, clinda, conv)

#----------------------------------------#

# Calculate population variance
strep_metabolome_mock <- cbind(pop_var(strep_metabolome_mock), rep('strep_metabolome_mock', ncol(strep_metabolome_mock)))
strep_metabolome_630 <- cbind(pop_var(strep_metabolome_630), rep('strep_metabolome_630', ncol(strep_metabolome_630)))
cef_metabolome_mock <- cbind(pop_var(cef_metabolome_mock), rep('cef_metabolome_mock', ncol(cef_metabolome_mock)))
cef_metabolome_630 <- cbind(pop_var(cef_metabolome_630), rep('cef_metabolome_630', ncol(cef_metabolome_630)))
clinda_metabolome_mock <- cbind(pop_var(clinda_metabolome_mock), rep('clinda_metabolome_mock', ncol(clinda_metabolome_mock)))
clinda_metabolome_630 <- cbind(pop_var(clinda_metabolome_630), rep('clinda_metabolome_630', ncol(clinda_metabolome_630)))
conv_metabolome_mock <- cbind(pop_var(conv_metabolome_mock), rep('conv_metabolome_mock', ncol(conv_metabolome_mock)))
strep_shared_mock <- cbind(pop_var(strep_shared_mock), rep('strep_shared_mock', ncol(strep_shared_mock)))
strep_shared_630 <- cbind(pop_var(strep_shared_630), rep('strep_shared_630', ncol(strep_shared_630)))
cef_shared_mock <- cbind(pop_var(cef_shared_mock), rep('cef_shared_mock', ncol(cef_shared_mock)))
cef_shared_630 <- cbind(pop_var(cef_shared_630), rep('cef_shared_630', ncol(cef_shared_630)))
clinda_shared_mock <- cbind(pop_var(clinda_shared_mock), rep('clinda_shared_mock', ncol(clinda_shared_mock)))
clinda_shared_630 <- cbind(pop_var(clinda_shared_630), rep('clinda_shared_630', ncol(clinda_shared_630)))
conv_shared_mock <- cbind(pop_var(conv_shared_mock), rep('conv_shared_mock', ncol(conv_shared_mock)))

# Assemble tables of plotting
metabolome_mock <- rbind(strep_metabolome_mock, cef_metabolome_mock, clinda_metabolome_mock, conv_metabolome_mock)
colnames(metabolome_mock) <-  c('pop_variance','treatment')
metabolome_630 <- rbind(strep_metabolome_630, cef_metabolome_630, clinda_metabolome_630, conv_metabolome_630)
colnames(metabolome_630) <-  c('pop_variance','treatment')
shared_mock <- rbind(strep_shared_mock, cef_shared_mock, clinda_shared_mock, conv_shared_mock)
colnames(shared_mock) <-  c('pop_variance','treatment')
shared_630 <- rbind(strep_shared_630, cef_shared_630, clinda_shared_630, conv_metabolome_630)
colnames(shared_630) <-  c('pop_variance','treatment')

# Calculate differences
#pvalues_mock <- p.adjust(c(wilcox.test(strep_metabolome_mock, cef_metabolome_mock, exact=F)$p.value,
#                           wilcox.test(strep_metabolome_mock, clinda_metabolome_mock, exact=F)$p.value,
#                           wilcox.test(strep_metabolome_mock, conv_metabolome_mock, exact=F)$p.value,
#                           wilcox.test(cef_metabolome_mock, clinda_metabolome_mock, exact=F)$p.value,
#                           wilcox.test(cef_metabolome_mock, conv_metabolome_mock, exact=F)$p.value,
#                           wilcox.test(clinda_metabolome_mock, conv_metabolome_mock, exact=F)$p.value), method='BH')
#pvalues_630 <- p.adjust(c(wilcox.test(strep_metabolome_630, cef_metabolome_630, exact=F)$p.value,
#                          wilcox.test(strep_metabolome_630, clinda_metabolome_630, exact=F)$p.value,
#                          wilcox.test(strep_metabolome_630, conv_metabolome_630, exact=F)$p.value,
#                          wilcox.test(cef_metabolome_630, clinda_metabolome_630, exact=F)$p.value,
#                          wilcox.test(cef_metabolome_630, conv_metabolome_630, exact=F)$p.value,
#                          wilcox.test(clinda_metabolome_630, conv_metabolome_630, exact=F)$p.value), method='BH')

#----------------------------------------#

# Generate plot
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S6.pdf'
pdf(file=plot_file, width=7, height=4)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))
par(las=1, mar=c(4,4,1,1), mgp=c(2.5,0.7,0))


# still need to define ylim for both plots

# Metabolome
boxplot(pop_variance~treatment, data=shared_mock, cex=0, lwd=2.5, xlim=c(0.5,11.5), at=c(1,4,7,10),
        col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5], 'gray50'),
        ylab='16S rRNA Gene Amplicon Read Abundance', staplewex=0.8, boxwex=1, lty=1, medlwd=2.5, xaxt='n', yaxt='n')
boxplot(pop_variance~treatment, data=shared_630, cex=0, lwd=2.5, xlim=c(0.5,11.5), at=c(2,5,8,11),
        col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5], 'gray50'),
        ylab='16S rRNA Gene Amplicon Read Abundance', staplewex=0.8, boxwex=1, lty=1, medlwd=2.5, xaxt='n', yaxt='n', add=TRUE)
axis(side=1, labels=c('Streptomycin','Cefoperazone','Clindamycin','No Antibiotics'), at=c(1.5,4.5,7.5,10.5), cex=1.2)
legend('topleft', legend='16S Abundance', bty='n', cex=1.2)
mtext('A', side=2, line=2, las=2, adj=1.3, padj=-7.5, cex=1.5)

# 16S
boxplot(pop_variance~treatment, data=metabolome_mock, cex=0, lwd=2.5, xlim=c(0.5,11.5), at=c(1,4,7,10),
        col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5], 'gray50'),
        ylab='16S rRNA Gene Amplicon Read Abundance', staplewex=0.8, boxwex=1, lty=1, medlwd=2.5, xaxt='n', yaxt='n')
boxplot(pop_variance~treatment, data=metabolome_630, cex=0, lwd=2.5, xlim=c(0.5,11.5), at=c(2,5,8,11),
        col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5], 'gray50'),
        ylab='16S rRNA Gene Amplicon Read Abundance', staplewex=0.8, boxwex=1, lty=1, medlwd=2.5, xaxt='n', yaxt='n', add=TRUE)
axis(side=1, labels=c('Streptomycin','Cefoperazone','Clindamycin','No Antibiotics'), at=c(1.5,4.5,7.5,10.5), cex=1.2)
legend('topleft', legend='Metabolome', bty='n', cex=1.2)
mtext('B', side=2, line=2, las=2, adj=1.3, padj=-7.5, cex=1.5)

dev.off()

#----------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
#rm(list=ls())
#gc()
