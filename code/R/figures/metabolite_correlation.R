
# Load dependencies
deps <- c('wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Define files
metabolome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolome/metabolomics.tsv'
cef_importances <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/importances.tsv'
strep_importances <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/streptomycin_630.bipartite.files/importances.tsv'
clinda_importances <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/clindamycin_630.bipartite.files/importances.tsv'
gf_importances <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/germfree_630.bipartite.files/importances.tsv'
metadata <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'

# Metabolomic data
metabolome <- read.delim(metabolome, sep='\t', header=T)
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome <- subset(metabolome, KEGG != 'NA')
metabolome <- metabolome[match(unique(metabolome$KEGG), metabolome$KEGG),]
rownames(metabolome) <- metabolome$KEGG
metabolome$KEGG <- NULL
metabolome$BIOCHEMICAL <- NULL

# Combine with metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metabolome <- merge(t(metabolome), metadata, by='row.names')
rm(metadata)
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
metabolome <- subset(metabolome, infection == '630')
metabolome <- subset(metabolome, abx != 'none')
metabolome$cage <- NULL
metabolome$mouse <- NULL
metabolome$gender <- NULL
metabolome$type <- NULL
metabolome$infection <- NULL
metabolome <- aggregate(metabolome[, 1:398], list(metabolome$abx), median)
rownames(metabolome) <- metabolome$Group.1
metabolome$Group.1 <- NULL
metabolome <- as.data.frame(t(metabolome))
colnames(metabolome) <- c('cefoperazone_conc', 'clindamycin_conc', 'germfree_conc', 'streptomycin_conc')
metabolome <- log10(metabolome + 1)

# Merge importances to one table
cef_importances <- read.delim(cef_importances, sep='\t', header=T, row.names=1)
cef_importances$p_value <- NULL
clinda_importances <- read.delim(clinda_importances, sep='\t', header=T, row.names=1)
clinda_importances$p_value <- NULL
clinda_importances$Compound_name <- NULL
strep_importances <- read.delim(strep_importances, sep='\t', header=T, row.names=1)
strep_importances$p_value <- NULL
strep_importances$Compound_name <- NULL
gf_importances <- read.delim(gf_importances, sep='\t', header=T, row.names=1)
gf_importances$p_value <- NULL
gf_importances$Compound_name <- NULL
importances <- merge(cef_importances, clinda_importances, by='row.names')
rownames(importances) <- importances$Row.names
importances$Row.names <- NULL
importances <- merge(importances, strep_importances, by='row.names')
rownames(importances) <- importances$Row.names
importances$Row.names <- NULL
importances <- merge(importances, gf_importances, by='row.names')
rownames(importances) <- importances$Row.names
importances$Row.names <- NULL
colnames(importances) <- c('Compound_name', 'cefoperazone_score', 'clindamycin_score', 'streptomycin_score', 'germfree_score')
importances$Compound_name <- gsub('_',' ', importances$Compound_name)
rm(cef_importances, clinda_importances, strep_importances, gf_importances)

# Merge metabolome medians and importance values
combined <- merge(importances, metabolome, by='row.names')
rownames(combined) <- combined$Row.names
combined$Row.names <- NULL
rm(importances, metabolome)

# Separate treatment groups
cef <- as.data.frame(cbind(combined$cefoperazone_score, combined$cefoperazone_conc))
rownames(cef) <- rownames(combined)
colnames(cef) <- c('score', 'conc')
cef$name <- combined$Compound_name
cef <- subset(cef, cef[,1] != 0)
clinda <- as.data.frame(cbind(combined$clindamycin_score, combined$clindamycin_conc))
rownames(clinda) <- rownames(combined)
colnames(clinda) <- c('score', 'conc')
clinda$name <- combined$Compound_name
clinda <- subset(clinda, clinda[,1] != 0)
strep <- as.data.frame(cbind(combined$streptomycin_score, combined$streptomycin_conc))
rownames(strep) <- rownames(combined)
colnames(strep) <- c('score', 'conc')
strep$name <- combined$Compound_name
strep <- subset(strep, strep[,1] != 0)
germfree <- as.data.frame(cbind(combined$germfree_score, combined$germfree_conc))
rownames(germfree) <- rownames(combined)
colnames(germfree) <- c('score', 'conc')
germfree$name <- combined$Compound_name
germfree <- subset(germfree, germfree[,1] != 0)
rm(combined)

# Calculate Spearman correlation
p.adjust(c(cor.test(cef[,1], cef[,2], method='spearman', exact=FALSE)$p.value, 
           cor.test(clinda[,1], clinda[,2], method='spearman', exact=FALSE)$p.value, 
           cor.test(strep[,1], strep[,2], method='spearman', exact=FALSE)$p.value, 
           cor.test(germfree[,1], germfree[,2], method='spearman', exact=FALSE)$p.value), method='BH')
c(cor.test(cef[,1], cef[,2], method='spearman', exact=FALSE)$estimate, 
  cor.test(clinda[,1], clinda[,2], method='spearman', exact=FALSE)$estimate, 
  cor.test(strep[,1], strep[,2], method='spearman', exact=FALSE)$estimate, 
  cor.test(germfree[,1], germfree[,2], method='spearman', exact=FALSE)$estimate)

# Identify outliers
cef_outliers <- subset(cef, conc >= (as.numeric(quantile(cef$conc)[4]) + (1.5 * IQR(cef$conc))) & score > 0)
clinda_outliers <- subset(clinda, conc >= (as.numeric(quantile(clinda$conc)[4]) + (1.5 * IQR(clinda$conc))) & score > 0)
strep_outliers <- subset(strep, conc >= (as.numeric(quantile(strep$conc)[4]) + (1.5 * IQR(strep$conc))) & score > 0)
germfree_outliers <- subset(germfree, conc >= (as.numeric(quantile(germfree$conc)[4]) + (1.5 * IQR(germfree$conc))) & score > 0)

#----------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_6.pdf'
select_palette <- c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black')
pdf(file=plot_file, width=8, height=7)
layout(matrix(c(1,2,
                3,4), 
              nrow=2, ncol=2, byrow=TRUE))
par(las=1, mar=c(4,4,1,1), mgp=c(2.2,1,0))

# Plot the data and correlations
plot(strep[,1], strep[,2], xlab=expression(paste('Importance Score (',log[2],')')), 
     ylab=expression(paste('Scaled Intensity (',log[10],')')), pch=19, cex=0.9, 
     xlim=c(-6,8), ylim=c(-0.1,1.1), col=wes_palette("FantasticFox")[1]) 
abline(lm(conc ~ score, data=strep), col='gray30', lwd=2)
points(strep_outliers[,1], strep_outliers[,2], pch=21, bg=wes_palette("FantasticFox")[1], cex=1.5, lwd=2) # color outlier
text(x=-0.8, y=0.9743735, strep_outliers$name, cex=0.9) # label outlier
mtext('A', side=2, line=2, las=2, adj=1.7, padj=-9, cex=1.2)
legend('topright', legend=c(expression(paste(italic('r'),' = -0.249')), expression(paste(italic('P'),' < 0.05'))), pt.cex=0, bty='n', cex=1.1)

plot(cef[,1], cef[,2], xlab=expression(paste('Importance Score (',log[2],')')), 
     ylab=expression(paste('Scaled Intensity (',log[10],')')), pch=19, cex=0.9, 
     xlim=c(-6,8), ylim=c(0,2.1), col=wes_palette("FantasticFox")[3]) 
abline(lm(conc ~ score, data=cef), col='gray30', lwd=2)
points(cef_outliers[,1], cef_outliers[,2], pch=21, bg=wes_palette("FantasticFox")[3], cex=1.5, lwd=2) # color outliers
text(x=c(0.913,1.287,1,-1.5), y=c(0.5920879,0.8211674,1.5861168,2.0454313), cef_outliers$name, cex=0.9) # label outliers
mtext('B', side=2, line=2, las=2, adj=1.7, padj=-9, cex=1.2)
legend('topright', legend=c(expression(paste(italic('r'),' = -0.121')), expression(paste(italic('P'),' = n.s.'))), pt.cex=0, bty='n', cex=1.1)

plot(clinda[,1], clinda[,2], xlab=expression(paste('Importance Score (',log[2],')')), 
     ylab=expression(paste('Scaled Intensity (',log[10],')')), pch=19, cex=0.9, 
     xlim=c(-7,7), ylim=c(0,0.8), col=wes_palette("FantasticFox")[5]) 
abline(lm(conc ~ score, data=clinda), col='gray30', lwd=2)
mtext('C', side=2, line=2, las=2, adj=1.7, padj=-9, cex=1.2)
legend('topright', legend=c(expression(paste(italic('r'),' = -0.341')), expression(paste(italic('P'),' < 0.05'))), pt.cex=0, bty='n', cex=1.1)

plot(germfree[,1], germfree[,2], xlab=expression(paste('Importance Score (',log[2],')')), 
     ylab=expression(paste('Scaled Intensity (',log[10],')')), pch=19, cex=0.9, 
     xlim=c(-6.5,6.5), ylim=c(0,1.6), col='forestgreen') 
abline(lm(conc ~ score, data=germfree), col='gray30', lwd=2)
points(germfree_outliers[,1], germfree_outliers[,2], pch=21, bg='forestgreen', cex=1.5, lwd=2) # color outliers
text(x=c(0.185,-0.2), y=c(1.164246,1.523506), germfree_outliers$name, cex=0.9) # label outliers
mtext('D', side=2, line=2, las=2, adj=1.7, padj=-9, cex=1.2)
legend('topright', legend=c(expression(paste(italic('r'),' = -0.117')), expression(paste(italic('P'),' = n.s.'))), pt.cex=0, bty='n', cex=1.1)

dev.off()

#----------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

