
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
annotation <- metabolome[,1:5]
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome <- subset(metabolome, KEGG != 'NA')
metabolome <- metabolome[match(unique(metabolome$KEGG), metabolome$KEGG),]
rownames(metabolome) <- metabolome$KEGG
metabolome$KEGG <- NULL
metabolome$BIOCHEMICAL <- NULL
annotation <- subset(annotation, KEGG != 'NA')
annotation <- annotation[match(unique(annotation$KEGG), annotation$KEGG),]
rownames(annotation) <- annotation$KEGG
annotation$KEGG <- NULL
annotation$SUB_PATHWAY <- NULL
annotation$PUBCHEM <- NULL
annotation$BIOCHEMICAL <- NULL
annotation$SUPER_PATHWAY <- gsub('_',' ', annotation$SUPER_PATHWAY)

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
metabolome_630 <- aggregate(metabolome_630[, 1:398], list(metabolome_630$abx), median)
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
metabolome_mock <- aggregate(metabolome_mock[, 1:398], list(metabolome_mock$abx), median)
rownames(metabolome_mock) <- metabolome_mock$Group.1
metabolome_mock$Group.1 <- NULL
metabolome_mock <- as.data.frame(t(metabolome_mock))
colnames(metabolome_mock) <- c('cef_mock','clinda_mock','gf_mock','strep_mock')
metabolome <- metabolome_mock / metabolome_630
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
cef <- subset(cef, cef[,1] != 0) # remove metabolites with 0 importance
clinda <- as.data.frame(cbind(combined$clindamycin_score, combined$clindamycin_conc))
rownames(clinda) <- rownames(combined)
colnames(clinda) <- c('score', 'conc')
clinda$name <- combined$Compound_name
clinda <- subset(clinda, clinda[,1] != 0) # remove metabolites with 0 importance
strep <- as.data.frame(cbind(combined$streptomycin_score, combined$streptomycin_conc))
rownames(strep) <- rownames(combined)
colnames(strep) <- c('score', 'conc')
strep$name <- combined$Compound_name
strep <- subset(strep, strep[,1] != 0) # remove metabolites with 0 importance
germfree <- as.data.frame(cbind(combined$germfree_score, combined$germfree_conc))
rownames(germfree) <- rownames(combined)
colnames(germfree) <- c('score', 'conc')
germfree$name <- combined$Compound_name
germfree <- subset(germfree, germfree[,1] != 0)
combined <- rbind(cef, clinda, strep, germfree)

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
cef_outliers <- rbind(subset(cef, conc >= (as.numeric(quantile(cef$conc)[4]) + (1.5 * IQR(cef$conc)))), 
                      subset(cef, conc <= (as.numeric(quantile(cef$conc)[2]) - (1.5 * IQR(cef$conc)))))
cef_outliers <- subset(cef_outliers, score > 0)
clinda_outliers <- rbind(subset(clinda, conc >= (as.numeric(quantile(clinda$conc)[4]) + (1.5 * IQR(clinda$conc)))), 
                         subset(clinda, conc <= (as.numeric(quantile(clinda$conc)[2]) - (1.5 * IQR(clinda$conc)))))
clinda_outliers <- subset(clinda_outliers, score > 0)
strep_outliers <- rbind(subset(strep, conc >= (as.numeric(quantile(strep$conc)[4]) + (1.5 * IQR(strep$conc)))), 
                        subset(strep, conc <= (as.numeric(quantile(strep$conc)[2]) - (1.5 * IQR(strep$conc)))))
strep_outliers <- subset(strep_outliers, score > 0)
germfree_outliers <- rbind(subset(germfree, conc >= (as.numeric(quantile(germfree$conc)[4]) + (1.5 * IQR(germfree$conc)))), 
                           subset(germfree, conc <= (as.numeric(quantile(germfree$conc)[2]) - (1.5 * IQR(germfree$conc)))))
germfree_outliers <- subset(germfree_outliers, score > 0)
germfree_outliers <- subset(germfree_outliers, conc >= 3.910228)
combined_outliers <- rbind(subset(combined, conc >= (as.numeric(quantile(combined$conc)[4]) + (1.5 * IQR(combined$conc)))), 
                           subset(combined, conc <= (as.numeric(quantile(combined$conc)[2]) - (1.5 * IQR(combined$conc)))))

# Add pathways to combined analysis
combined_annotation <- merge(annotation, combined_outliers, by='row.names')
rownames(combined_annotation) <- combined_annotation$Row.names
combined_annotation$Row.names <- NULL

# not merging properly...





# Calculate stats
test <- as.data.frame(cbind(round(c(lm(conc ~ score, data=strep)$coefficients[[2]],
                                    lm(conc ~ score, data=cef)$coefficients[[2]],
                                    lm(conc ~ score, data=clinda)$coefficients[[2]],
                                    lm(conc ~ score, data=germfree)$coefficients[[2]],
                                    lm(conc ~ score, data=combined)$coefficients[[2]]), digits=3),
                            round(c(cor.test(strep[,1], strep[,2], method='spearman', exact=FALSE)$estimate,
                                    cor.test(cef[,1], cef[,2], method='spearman', exact=FALSE)$estimate,
                                    cor.test(clinda[,1], clinda[,2], method='spearman', exact=FALSE)$estimate,
                                    cor.test(germfree[,1], germfree[,2], method='spearman', exact=FALSE)$estimate,
                                    cor.test(combined[,1], combined[,2], method='spearman', exact=FALSE)$estimate), digits=3),
                            round(c(cor.test(strep[,1], strep[,2], method='spearman', exact=FALSE)$p.value,
                                             cor.test(cef[,1], cef[,2], method='spearman', exact=FALSE)$p.value,
                                             cor.test(clinda[,1], clinda[,2], method='spearman', exact=FALSE)$p.value,
                                             cor.test(germfree[,1], germfree[,2], method='spearman', exact=FALSE)$p.value,
                                             cor.test(combined[,1], combined[,2], method='spearman', exact=FALSE)$p.value), digits=3)))
rownames(test) <- c('streptomycin','cefoperazone','clindamycin','germfree','combined')
colnames(test) <- c('m','r', 'p')

#----------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_6.pdf'
select_palette <- c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black')
pdf(file=plot_file, width=6.3, height=9)
layout(matrix(c(1,1,
                2,3,
                4,5), 
              nrow=3, ncol=2, byrow=TRUE))
par(las=1, mar=c(3,4,1,1), mgp=c(1.8,0.7,0))

# Plot the data and correlations
plot(combined[,1], combined[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Scaled Intensity')), 
     pch=19, cex=1.1, xlim=c(-10,10), ylim=c(-1,15), col='gray30')
abline(lm(conc ~ score, data=strep), col='black', lwd=2)
legend('topright', legend='All Treatment Groups', bty='n', cex=1.2)
mtext('A', side=2, line=2, las=2, adj=1.2, padj=-8, cex=1.2)


legend('topleft', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(test$r[5])))), 
                           as.expression(bquote(paste(italic('P'),' = ',.(test$p[5]),'*')))), pt.cex=0, bty='n', cex=1.1)


points(combined_outliers[,1], combined_outliers[,2], pch=21, bg='white', col='gray30', cex=1.7, lwd=2)

# lEGEND COLORING POINTS BASED ON TYPE OF MOLECULE


# Add dotted line for outlier cutoff and frequently appearing outliers






# streptomycin alone
plot(strep[,1], strep[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Scaled Intensity')), 
     pch=19, xlim=c(-10,10), ylim=c(-1,9), col=wes_palette("FantasticFox")[1])
abline(v=0, lty=2, col='gray30')
abline(lm(conc ~ score, data=strep), col='black', lwd=2)
mtext('B', side=2, line=2, las=2, adj=1.2, padj=-8, cex=1.2)
legend('topleft', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(test$r[1])))), 
                           as.expression(bquote(paste(italic('P'),' = ',.(test$p[1]))))), pt.cex=0, bty='n', cex=1.1)
legend('topright', legend='Streptomycin', bty='n', cex=1.1)
points(strep_outliers[,1], strep_outliers[,2], pch=21, bg=wes_palette("FantasticFox")[1], cex=1.7, lwd=2)
text(x=c(8,7,5.5,5,-5.6), 
     y=c(1.8,7.3,3,4.5,2.6), 
     strep_outliers$name, cex=0.9)
segments(x0=c(2.5,-0.8,2.3), y0=c(4.3,2.6,2), 
         x1=c(1.3,1,5.6), y1=c(3.5,2.3,1.9))

# cefoperazone alone
plot(cef[,1], cef[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Scaled Intensity')), 
     pch=19, xlim=c(-10,10), ylim=c(-1,11), col=wes_palette("FantasticFox")[3]) 
abline(v=0, lty=2, col='gray30')
abline(lm(conc ~ score, data=cef), col='black', lwd=2)
mtext('C', side=2, line=2, las=2, adj=1.2, padj=-8, cex=1.2)
legend('topleft', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(test$r[2])))), 
                           as.expression(bquote(paste(italic('P'),' = ',.(test$p[2]))))), pt.cex=0, bty='n', cex=1.1)
legend('topright', legend='Cefoperazone', bty='n', cex=1.1)
points(cef_outliers[,1], cef_outliers[,2], pch=21, bg=wes_palette("FantasticFox")[3], cex=1.7, lwd=2)
text(x=c(5,8,7.3), 
     y=c(5,7.9,3.409411), 
     cef_outliers$name, cex=0.9)
segments(x0=1, y0=3.5, x1=2.5, y1=4.6)

# clindamycin alone
plot(clinda[,1], clinda[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Scaled Intensity')), 
     pch=19, xlim=c(-10,10), ylim=c(-1,5), col=wes_palette("FantasticFox")[5])
abline(v=0, lty=2, col='gray30')
abline(lm(conc ~ score, data=clinda), col='black', lwd=2)
mtext('D', side=2, line=2, las=2, adj=1.2, padj=-8, cex=1.2)
legend('topleft', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(test$r[3])))), 
                           as.expression(bquote(paste(italic('P'),' = ',.(test$p[3]),'*')))), pt.cex=0, bty='n', cex=1.1)
legend('topright', legend='Clindamycin', bty='n', cex=1.1)
points(clinda_outliers[,1], clinda_outliers[,2], pch=21, bg=wes_palette("FantasticFox")[5], cex=1.7, lwd=2)
text(x=c(4.7,7.2,-5,5.7,5.5), 
     y=c(2.2279556,4,2.2,3.2,-0.1677397), 
     clinda_outliers$name, cex=0.9)
segments(x0=-0.5, y0=2.2, x1=0.6, y1=2)

# germfree alone
plot(germfree[,1], germfree[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Scaled Intensity')), 
     pch=19, cex=0.9, xlim=c(-8,8), ylim=c(-1,17), col='forestgreen', xaxt='n')
axis(side=1, at=c(-8,-4,0,4,8), labels=c(-8,-4,0,4,8))
abline(v=0, lty=2, col='gray30')
abline(lm(conc ~ score, data=germfree), col='black', lwd=2)
mtext('E', side=2, line=2, las=2, adj=1.2, padj=-8, cex=1.2)
legend('topleft', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(test$r[4])))), 
                           as.expression(bquote(paste(italic('P'),' = ',.(test$p[4]),'*')))), pt.cex=0, bty='n', cex=1.1)
legend('topright', legend='Gnotobiotic', bty='n', cex=1.1)
points(germfree_outliers[,1], germfree_outliers[,2], pch=21, bg='forestgreen', cex=1.5, lwd=2) # color outliers
text(x=c(4,3.2,5.7,5.9), 
     y=c(4.523075,5.924171,14.439394,11.890565), 
     germfree_outliers$name, cex=0.9)

dev.off()

#----------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
test
#rm(list=ls())
#gc()

