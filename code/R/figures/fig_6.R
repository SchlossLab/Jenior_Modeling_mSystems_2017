
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
metabolome <- subset(metabolome, KEGG != 'C00337') ###
metabolome <- subset(metabolome, KEGG != 'C00438') ###
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
metabolome <- subset(metabolome, abx != 'none')

metabolome_630 <- subset(metabolome, infection == '630')
metabolome_630$cage <- NULL
metabolome_630$mouse <- NULL
metabolome_630$gender <- NULL
metabolome_630$type <- NULL
metabolome_630$infection <- NULL
metabolome_630 <- aggregate(metabolome_630[, 1:396], list(metabolome_630$abx), median)
rownames(metabolome_630) <- metabolome_630$Group.1
metabolome_630$Group.1 <- NULL
metabolome_630 <- as.data.frame(t(metabolome_630))
colnames(metabolome_630) <- c('cef_630','clinda_630','gf_630','strep_630')
metabolome_mock <- subset(metabolome, infection == 'mock')
metabolome_mock$cage <- NULL
metabolome_mock$mouse <- NULL
metabolome_mock$gender <- NULL
metabolome_mock$type <- NULL
metabolome_mock$infection <- NULL
metabolome_mock <- aggregate(metabolome_mock[, 1:396], list(metabolome_mock$abx), median)
rownames(metabolome_mock) <- metabolome_mock$Group.1
metabolome_mock$Group.1 <- NULL
metabolome_mock <- as.data.frame(t(metabolome_mock))
colnames(metabolome_mock) <- c('cef_mock','clinda_mock','gf_mock','strep_mock')
metabolome <- metabolome_mock - metabolome_630
colnames(metabolome) <- c('cefoperazone_conc', 'clindamycin_conc', 'germfree_conc', 'streptomycin_conc')
rm(metabolome_mock, metabolome_630)

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

# Identify outliers with a positive importance score
cef_outliers <- rbind(subset(cef, conc >= (as.numeric(quantile(cef$conc)[4]) + (1.5 * IQR(cef$conc))) & score > 0), 
                      subset(cef, conc <= (as.numeric(quantile(cef$conc)[2]) - (1.5 * IQR(cef$conc))) & score > 0))
clinda_outliers <- rbind(subset(clinda, conc >= (as.numeric(quantile(clinda$conc)[4]) + (1.5 * IQR(clinda$conc))) & score > 0), 
                         subset(clinda, conc <= (as.numeric(quantile(clinda$conc)[2]) - (1.5 * IQR(clinda$conc))) & score > 0))
strep_outliers <- rbind(subset(strep, conc >= (as.numeric(quantile(strep$conc)[4]) + (1.5 * IQR(strep$conc))) & score > 0), 
                        subset(strep, conc <= (as.numeric(quantile(strep$conc)[2]) - (1.5 * IQR(strep$conc))) & score > 0))
germfree_outliers <- rbind(subset(germfree, conc >= (as.numeric(quantile(germfree$conc)[4]) + (1.5 * IQR(germfree$conc))) & score > 0), 
                           subset(germfree, conc <= (as.numeric(quantile(germfree$conc)[2]) - (1.5 * IQR(germfree$conc))) & score > 0))

#----------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_6.pdf'
select_palette <- c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black')
pdf(file=plot_file, width=3, height=10.5)
layout(matrix(c(1,
                2,
                3,
                4), 
              nrow=4, ncol=1, byrow=TRUE))
par(las=1, mar=c(3,4,1,1), mgp=c(2,0.7,0), xaxs='i', yaxs='i')

# Plot the data and correlations
plot(strep[,1], strep[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Scaled Intensity')), 
     pch=19, cex=0.9, xlim=c(-9,9), ylim=c(-8,6), col=wes_palette("FantasticFox")[1], xaxt='n')
axis(side=1, at=seq(-9,9,3), labels=seq(-9,9,3))
abline(v=0, lty=2, col='gray30')
abline(h=0, lty=2, col='gray30')
abline(lm(conc ~ score, data=strep), col='black', lwd=2)
mtext('A', side=2, line=2, las=2, adj=1.7, padj=-7, cex=1.2)
legend('topleft', legend=c(expression(paste(italic('r'),' = 0.108')), expression(paste(italic('P'),' = 0.406'))), pt.cex=0, bty='n', cex=1.1)
# Label outliers
points(strep_outliers[,1], strep_outliers[,2], pch=21, bg=wes_palette("FantasticFox")[1], cex=1.5, lwd=2)
text(x=c(4.5,5.2,5.8,3.5,7,4,-1.5,6.2), 
     y=c(1.9862,4.6330,-2.5686,-7.1469,-3.5686,-5.8,-2.3,-4.8), 
     strep_outliers$name, cex=0.8)
segments(x0=c(2.6,2.47,3.6), y0=c(-3.5686,-4.57,-4.3), 
         x1=c(5.2,3,4.2), y1=c(-3.5686,-5.6,-4.6))


plot(cef[,1], cef[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Scaled Intensity')), 
     pch=19, cex=0.9, xlim=c(-10,8), ylim=c(-4,4), col=wes_palette("FantasticFox")[3], xaxt='n') 
axis(side=1, at=seq(-10,8,2), labels=seq(-10,8,2))
abline(v=0, lty=2, col='gray30')
abline(h=0, lty=2, col='gray30')
abline(lm(conc ~ score, data=cef), col='black', lwd=2)
mtext('B', side=2, line=2, las=2, adj=1.7, padj=-7, cex=1.2)
legend('topleft', legend=c(expression(paste(italic('r'),' = 0.051')), expression(paste(italic('P'),' = 0.643'))), pt.cex=0, bty='n', cex=1.1)
# Label outliers
points(cef_outliers[,1], cef_outliers[,2], pch=21, bg=wes_palette("FantasticFox")[3], cex=1.5, lwd=2)
#text(x=c(6.4,4,5.2,2,4), 
#     y=c(6,14,-4,9,49), 
#     cef_outliers$name, cex=0.8)
#segments(x0=c(2,4.1,4.9), y0=c(7,12,-0.5), 
#         x1=c(3.4,4.1,5.1), y1=c(4,5,-2.5))


plot(clinda[,1], clinda[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Scaled Intensity')), 
     pch=19, cex=0.9, xlim=c(-10,8), ylim=c(-4,2), col=wes_palette("FantasticFox")[5], xaxt='n')
axis(side=1, at=seq(-10,8,2), labels=seq(-10,8,2))
abline(v=0, lty=2, col='gray30')
abline(h=0, lty=2, col='gray30')
abline(lm(conc ~ score, data=clinda), col='black', lwd=2)
mtext('C', side=2, line=2, las=2, adj=1.7, padj=-7, cex=1.2)
legend('topleft', legend=c(expression(paste(italic('r'),' = 0.268')), expression(paste(italic('P'),' = 0.022*'))), pt.cex=0, bty='n', cex=1.1)
# Label outliers
points(clinda_outliers[,1], clinda_outliers[,2], pch=21, bg=wes_palette("FantasticFox")[5], cex=1.5, lwd=2)
text(x=c(6.3,-4,1,3.9,1.585,4.4,-3.1), 
     y=c(1.5,0.75,1.3,-1.45,-2,-0.7,-0.7368), 
     clinda_outliers$name, cex=0.8)
segments(x0=c(-0.4,3,1.6), y0=c(0.75,1.15,-1.1), 
         x1=c(0.7,3.9,1.6), y1=c(0.65,0.75,-1.8))

plot(germfree[,1], germfree[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Scaled Intensity')), 
     pch=19, cex=0.9, xlim=c(-6,8), ylim=c(-6,15), col='forestgreen', yaxt='n')
axis(side=2, at=seq(-6,15,3), labels=seq(-6,15,3))

abline(v=0, lty=2, col='gray30')
abline(h=0, lty=2, col='gray30')
abline(lm(conc ~ score, data=germfree), col='black', lwd=2)
mtext('D', side=2, line=2, las=2, adj=1.7, padj=-7, cex=1.2)
legend('topleft', legend=c(expression(paste(italic('r'),' = 0.371')), expression(paste(italic('P'),' = 0.022*'))), pt.cex=0, bty='n', cex=1.1)
# Label outliers
points(germfree_outliers[,1], germfree_outliers[,2], pch=21, bg='forestgreen', cex=1.5, lwd=2) # color outliers
text(x=c(1.950,4.2,3.358,3,3.4), 
     y=c(4,3.2,6.7,9.7564,12.0140), 
     germfree_outliers$name, cex=0.8)

dev.off()

#----------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

