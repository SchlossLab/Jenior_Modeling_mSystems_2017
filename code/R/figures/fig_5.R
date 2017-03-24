
# Start with blank slate
#rm(list=ls())
#gc()

# Load dependency
if ('wesanderson' %in% installed.packages()[,"Package"] == FALSE){
  install.packages(as.character('wesanderson'), quiet=TRUE);
}
library('wesanderson', verbose=FALSE, character.only=TRUE)

# Select files
metabolome <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/metabolomics.scaled_intensities.tsv'
metadata <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'

# Read in data
metabolome <- read.delim(metabolome, sep='\t', header=T, row.names=1)
metabolome <- metabolome[, !colnames(metabolome) %in% c('CefC5M2','StrepC4M1')] # Remove possible contamination
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination

# Merge metabolomics with metadata
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$type <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome <- as.data.frame(t(metabolome))
metabolome <- merge(metadata, metabolome, by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
rm(metadata)

# Subset metabolomics - untreated vs mock infected comparison
mock_metabolome <- subset(metabolome, infection == 'mock')
mock_metabolome$infection <- NULL
acetylglucosamine_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('N-acetylglucosamine/N-acetylgalactosamine')))]
colnames(acetylglucosamine_mock) <- c('abx', 'substrate')
acetylglucosamine_mock$abx <- factor(acetylglucosamine_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
proline_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('proline')))]
colnames(proline_mock) <- c('abx', 'substrate')
proline_mock$abx <- factor(proline_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
mannitolsorbitol_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('mannitol/sorbitol')))]
colnames(mannitolsorbitol_mock) <- c('abx', 'substrate')
mannitolsorbitol_mock$abx <- factor(mannitolsorbitol_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
acetylneuraminate_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('N-acetylneuraminate')))]
colnames(acetylneuraminate_mock) <- c('abx', 'substrate')
acetylneuraminate_mock$abx <- factor(acetylneuraminate_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
glucose_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('glucose')))]
colnames(glucose_mock) <- c('abx', 'substrate')
glucose_mock$abx <- factor(glucose_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
fructose_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('fructose')))]
colnames(fructose_mock) <- c('abx', 'substrate')
fructose_mock$abx <- factor(fructose_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
ribose_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('ribose')))]
colnames(ribose_mock) <- c('abx', 'substrate')
ribose_mock$abx <- factor(ribose_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
rm(mock_metabolome)

# Subset metabolites - mock vs infected comparison
acetylglucosamine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-acetylglucosamine/N-acetylgalactosamine')))]
acetylglucosamine_cef <- subset(acetylglucosamine, abx == 'cefoperazone')
acetylglucosamine_cef$abx <- NULL
colnames(acetylglucosamine_cef) <- c('infection', 'substrate')
acetylglucosamine_cef$infection <- factor(acetylglucosamine_cef$infection, levels=c('mock','630'))
acetylglucosamine_strep <- subset(acetylglucosamine, abx == 'streptomycin')
acetylglucosamine_strep$abx <- NULL
colnames(acetylglucosamine_strep) <- c('infection', 'substrate')
acetylglucosamine_strep$infection <- factor(acetylglucosamine_strep$infection, levels=c('mock','630'))
acetylglucosamine_clinda <- subset(acetylglucosamine, abx == 'clindamycin')
acetylglucosamine_clinda$abx <- NULL
colnames(acetylglucosamine_clinda) <- c('infection', 'substrate')
acetylglucosamine_clinda$infection <- factor(acetylglucosamine_clinda$infection, levels=c('mock','630'))
acetylglucosamine_gf <- subset(acetylglucosamine, abx == 'germfree')
acetylglucosamine_gf$abx <- NULL
colnames(acetylglucosamine_gf) <- c('infection', 'substrate')
acetylglucosamine_gf$infection <- factor(acetylglucosamine_gf$infection, levels=c('mock','630'))
rm(acetylglucosamine)
proline <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('proline')))]
proline_cef <- subset(proline, abx == 'cefoperazone')
proline_cef$abx <- NULL
colnames(proline_cef) <- c('infection', 'substrate')
proline_cef$infection <- factor(proline_cef$infection, levels=c('mock','630'))
proline_strep <- subset(proline, abx == 'streptomycin')
proline_strep$abx <- NULL
colnames(proline_strep) <- c('infection', 'substrate')
proline_strep$infection <- factor(proline_strep$infection, levels=c('mock','630'))
proline_clinda <- subset(proline, abx == 'clindamycin')
proline_clinda$abx <- NULL
colnames(proline_clinda) <- c('infection', 'substrate')
proline_clinda$infection <- factor(proline_clinda$infection, levels=c('mock','630'))
proline_gf <- subset(proline, abx == 'germfree')
proline_gf$abx <- NULL
colnames(proline_gf) <- c('infection', 'substrate')
proline_gf$infection <- factor(proline_gf$infection, levels=c('mock','630'))
rm(proline)
mannitolsorbitol <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('mannitol/sorbitol')))]
mannitolsorbitol_cef <- subset(mannitolsorbitol, abx == 'cefoperazone')
mannitolsorbitol_cef$abx <- NULL
colnames(mannitolsorbitol_cef) <- c('infection', 'substrate')
mannitolsorbitol_cef$infection <- factor(mannitolsorbitol_cef$infection, levels=c('mock','630'))
mannitolsorbitol_strep <- subset(mannitolsorbitol, abx == 'streptomycin')
mannitolsorbitol_strep$abx <- NULL
colnames(mannitolsorbitol_strep) <- c('infection', 'substrate')
mannitolsorbitol_strep$infection <- factor(mannitolsorbitol_strep$infection, levels=c('mock','630'))
mannitolsorbitol_clinda <- subset(mannitolsorbitol, abx == 'clindamycin')
mannitolsorbitol_clinda$abx <- NULL
colnames(mannitolsorbitol_clinda) <- c('infection', 'substrate')
mannitolsorbitol_clinda$infection <- factor(mannitolsorbitol_clinda$infection, levels=c('mock','630'))
mannitolsorbitol_gf <- subset(mannitolsorbitol, abx == 'germfree')
mannitolsorbitol_gf$abx <- NULL
colnames(mannitolsorbitol_gf) <- c('infection', 'substrate')
mannitolsorbitol_gf$infection <- factor(mannitolsorbitol_gf$infection, levels=c('mock','630'))
rm(mannitolsorbitol)
acetylneuraminate <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-acetylneuraminate')))]
acetylneuraminate_cef <- subset(acetylneuraminate, abx == 'cefoperazone')
acetylneuraminate_cef$abx <- NULL
colnames(acetylneuraminate_cef) <- c('infection', 'substrate')
acetylneuraminate_cef$infection <- factor(acetylneuraminate_cef$infection, levels=c('mock','630'))
acetylneuraminate_strep <- subset(acetylneuraminate, abx == 'streptomycin')
acetylneuraminate_strep$abx <- NULL
colnames(acetylneuraminate_strep) <- c('infection', 'substrate')
acetylneuraminate_strep$infection <- factor(acetylneuraminate_strep$infection, levels=c('mock','630'))
acetylneuraminate_clinda <- subset(acetylneuraminate, abx == 'clindamycin')
acetylneuraminate_clinda$abx <- NULL
colnames(acetylneuraminate_clinda) <- c('infection', 'substrate')
acetylneuraminate_clinda$infection <- factor(acetylneuraminate_clinda$infection, levels=c('mock','630'))
acetylneuraminate_gf <- subset(acetylneuraminate, abx == 'germfree')
acetylneuraminate_gf$abx <- NULL
colnames(acetylneuraminate_gf) <- c('infection', 'substrate')
acetylneuraminate_gf$infection <- factor(acetylneuraminate_gf$infection, levels=c('mock','630'))
rm(acetylneuraminate)
glucose <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('glucose')))]
glucose_strep <- subset(glucose, abx == 'streptomycin')
glucose_strep$abx <- NULL
colnames(glucose_strep) <- c('infection', 'substrate')
glucose_strep$infection <- factor(glucose_strep$infection, levels=c('mock','630'))
glucose_cef <- subset(glucose, abx == 'cefoperazone')
glucose_cef$abx <- NULL
colnames(glucose_cef) <- c('infection', 'substrate')
glucose_cef$infection <- factor(glucose_cef$infection, levels=c('mock','630'))
glucose_clinda <- subset(glucose, abx == 'clindamycin')
glucose_clinda$abx <- NULL
colnames(glucose_clinda) <- c('infection', 'substrate')
glucose_clinda$infection <- factor(glucose_clinda$infection, levels=c('mock','630'))
glucose_gf <- subset(glucose, abx == 'germfree')
glucose_gf$abx <- NULL
colnames(glucose_gf) <- c('infection', 'substrate')
glucose_gf$infection <- factor(glucose_gf$infection, levels=c('mock','630'))
rm(glucose)
fructose <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('fructose')))]
fructose_strep <- subset(fructose, abx == 'streptomycin')
fructose_strep$abx <- NULL
colnames(fructose_strep) <- c('infection', 'substrate')
fructose_strep$infection <- factor(fructose_strep$infection, levels=c('mock','630'))
fructose_cef <- subset(fructose, abx == 'cefoperazone')
fructose_cef$abx <- NULL
colnames(fructose_cef) <- c('infection', 'substrate')
fructose_cef$infection <- factor(fructose_cef$infection, levels=c('mock','630'))
fructose_clinda <- subset(fructose, abx == 'clindamycin')
fructose_clinda$abx <- NULL
colnames(fructose_clinda) <- c('infection', 'substrate')
fructose_clinda$infection <- factor(fructose_clinda$infection, levels=c('mock','630'))
fructose_gf <- subset(fructose, abx == 'germfree')
fructose_gf$abx <- NULL
colnames(fructose_gf) <- c('infection', 'substrate')
fructose_gf$infection <- factor(fructose_gf$infection, levels=c('mock','630'))
rm(fructose)
ribose <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('ribose')))]
ribose_strep <- subset(ribose, abx == 'streptomycin')
ribose_strep$abx <- NULL
colnames(ribose_strep) <- c('infection', 'substrate')
ribose_strep$infection <- factor(ribose_strep$infection, levels=c('mock','630'))
ribose_cef <- subset(ribose, abx == 'cefoperazone')
ribose_cef$abx <- NULL
colnames(ribose_cef) <- c('infection', 'substrate')
ribose_cef$infection <- factor(ribose_cef$infection, levels=c('mock','630'))
ribose_clinda <- subset(ribose, abx == 'clindamycin')
ribose_clinda$abx <- NULL
colnames(ribose_clinda) <- c('infection', 'substrate')
ribose_clinda$infection <- factor(ribose_clinda$infection, levels=c('mock','630'))
ribose_gf <- subset(ribose, abx == 'germfree')
ribose_gf$abx <- NULL
colnames(ribose_gf) <- c('infection', 'substrate')
ribose_gf$infection <- factor(ribose_gf$infection, levels=c('mock','630'))
rm(ribose)
rm(metabolome)


#-------------------------------------------------------------------------------------------------------------------------------------#

# Calculate significant differences

# N-acetylglucosamine - Untreated vs Mock
p.adjust(c(wilcox.test(subset(acetylglucosamine_mock, abx=='none')[,2], subset(acetylglucosamine_mock, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_mock, abx=='none')[,2], subset(acetylglucosamine_mock, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_mock, abx=='none')[,2], subset(acetylglucosamine_mock, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_mock, abx=='none')[,2], subset(acetylglucosamine_mock, abx=='germfree')[,2], exact=F)$p.value), method='BH')
# Proline - Untreated vs Mock
p.adjust(c(wilcox.test(subset(proline_mock, abx=='none')[,2], subset(proline_mock, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_mock, abx=='none')[,2], subset(proline_mock, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_mock, abx=='none')[,2], subset(proline_mock, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_mock, abx=='none')[,2], subset(proline_mock, abx=='germfree')[,2], exact=F)$p.value), method='BH')
# Mannitol + Sorbitol - Untreated vs Mock
p.adjust(c(wilcox.test(subset(mannitolsorbitol_mock, abx=='none')[,2], subset(mannitolsorbitol_mock, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_mock, abx=='none')[,2], subset(mannitolsorbitol_mock, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_mock, abx=='none')[,2], subset(mannitolsorbitol_mock, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_mock, abx=='none')[,2], subset(mannitolsorbitol_mock, abx=='germfree')[,2], exact=F)$p.value), method='BH')
# N-acetylneuraminate - Untreated vs Mock
p.adjust(c(wilcox.test(subset(acetylneuraminate_mock, abx=='none')[,2], subset(acetylneuraminate_mock, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_mock, abx=='none')[,2], subset(acetylneuraminate_mock, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_mock, abx=='none')[,2], subset(acetylneuraminate_mock, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_mock, abx=='none')[,2], subset(acetylneuraminate_mock, abx=='germfree')[,2], exact=F)$p.value), method='BH')
# Glucose - Untreated vs Mock
p.adjust(c(wilcox.test(subset(glucose_mock, abx=='none')[,2], subset(glucose_mock, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(glucose_mock, abx=='none')[,2], subset(glucose_mock, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(glucose_mock, abx=='none')[,2], subset(glucose_mock, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(glucose_mock, abx=='none')[,2], subset(glucose_mock, abx=='germfree')[,2], exact=F)$p.value), method='BH')
# Fructose - Untreated vs Mock
p.adjust(c(wilcox.test(subset(fructose_mock, abx=='none')[,2], subset(fructose_mock, abx=='streptomycin')[,2], exact=F)$p.value,
           wilcox.test(subset(fructose_mock, abx=='none')[,2], subset(fructose_mock, abx=='cefoperazone')[,2], exact=F)$p.value,
           wilcox.test(subset(fructose_mock, abx=='none')[,2], subset(fructose_mock, abx=='clindamycin')[,2], exact=F)$p.value,
           wilcox.test(subset(fructose_mock, abx=='none')[,2], subset(fructose_mock, abx=='germfree')[,2], exact=F)$p.value), method='BH')

# N-acetylglucosamine - Mock vs Infected
p.adjust(c(wilcox.test(subset(acetylglucosamine_strep, infection=='630')[,2], subset(acetylglucosamine_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_cef, infection=='630')[,2], subset(acetylglucosamine_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_clinda, infection=='630')[,2], subset(acetylglucosamine_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylglucosamine_gf, infection=='630')[,2], subset(acetylglucosamine_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
# Proline - Mock vs Infected
p.adjust(c(wilcox.test(subset(proline_strep, infection=='630')[,2], subset(proline_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_cef, infection=='630')[,2], subset(proline_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_clinda, infection=='630')[,2], subset(proline_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(proline_gf, infection=='630')[,2], subset(proline_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
# Mannitol + Sorbitol - Mock vs Infected
p.adjust(c(wilcox.test(subset(mannitolsorbitol_strep, infection=='630')[,2], subset(mannitolsorbitol_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_cef, infection=='630')[,2], subset(mannitolsorbitol_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_clinda, infection=='630')[,2], subset(mannitolsorbitol_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(mannitolsorbitol_gf, infection=='630')[,2], subset(mannitolsorbitol_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
# N-acetylneuraminate - Mock vs Infected
p.adjust(c(wilcox.test(subset(acetylneuraminate_strep, infection=='630')[,2], subset(acetylneuraminate_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_cef, infection=='630')[,2], subset(acetylneuraminate_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_clinda, infection=='630')[,2], subset(acetylneuraminate_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(acetylneuraminate_gf, infection=='630')[,2], subset(acetylneuraminate_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
# Glucose - Mock vs Infected
p.adjust(c(wilcox.test(subset(glucose_strep, infection=='630')[,2], subset(glucose_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glucose_cef, infection=='630')[,2], subset(glucose_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glucose_clinda, infection=='630')[,2], subset(glucose_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(glucose_gf, infection=='630')[,2], subset(glucose_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
# Fructose - Mock vs Infected
p.adjust(c(wilcox.test(subset(fructose_strep, infection=='630')[,2], subset(fructose_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(fructose_cef, infection=='630')[,2], subset(fructose_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(fructose_clinda, infection=='630')[,2], subset(fructose_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(fructose_gf, infection=='630')[,2], subset(fructose_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')

#-------------------------------------------------------------------------------------------------------------------------------------#

# Extract just untreated points from mock infection subsets
acetylglucosamine_mock <- subset(acetylglucosamine_mock, abx == 'none')
proline_mock <- subset(proline_mock, abx == 'none')
mannitolsorbitol_mock <- subset(mannitolsorbitol_mock, abx == 'none')
acetylneuraminate_mock <- subset(acetylneuraminate_mock, abx == 'none')
glucose_mock <- subset(glucose_mock, abx == 'none')
fructose_mock <- subset(fructose_mock, abx == 'none')
ribose_mock <- subset(ribose_mock, abx == 'none')

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_5.pdf'
select_palette <- c('gray40', wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen')
pdf(file=plot_file, width=7, height=11)
layout(matrix(c(1,
                2,
                3,
                4,
                5,
                6), nrow=6, ncol=1, byrow=TRUE))

par(mar=c(3,5,1,1), xpd=FALSE, las=1, mgp=c(3,0.7,0))

# N-acetylglucosamine
stripchart(substrate~abx, data=acetylglucosamine_mock, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,14), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
stripchart(substrate~infection, data=acetylglucosamine_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,14), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylglucosamine_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,14), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylglucosamine_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,14), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylglucosamine_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,14), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
axis(side=2, at=seq(0,14,3.5), labels=c('0.0','3.5','7.0','10.5','14.0'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.5, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.9)
mtext('A', side=2, line=2, las=2, adj=2, padj=-4, cex=1.3)
legend('topright', legend='GlcNAc / GalNAc', pt.cex=0, cex=1.2, bty='n')
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(acetylglucosamine_mock[,2]),
              median(subset(acetylglucosamine_strep, infection=='mock')[,2]), median(subset(acetylglucosamine_strep, infection=='630')[,2]),
              median(subset(acetylglucosamine_cef, infection=='mock')[,2]), median(subset(acetylglucosamine_cef, infection=='630')[,2]),
              median(subset(acetylglucosamine_clinda, infection=='mock')[,2]), median(subset(acetylglucosamine_clinda, infection=='630')[,2]),
              median(subset(acetylglucosamine_gf, infection=='mock')[,2]), median(subset(acetylglucosamine_gf, infection=='630')[,2])), 
         y1=c(median(acetylglucosamine_mock[,2]),
              median(subset(acetylglucosamine_strep, infection=='mock')[,2]), median(subset(acetylglucosamine_strep, infection=='630')[,2]),
              median(subset(acetylglucosamine_cef, infection=='mock')[,2]), median(subset(acetylglucosamine_cef, infection=='630')[,2]),
              median(subset(acetylglucosamine_clinda, infection=='mock')[,2]), median(subset(acetylglucosamine_clinda, infection=='630')[,2]),
              median(subset(acetylglucosamine_gf, infection=='mock')[,2]), median(subset(acetylglucosamine_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=12, y0=2.5, x1=13, y1=2.5, lwd=2)
text(x=12.5, y=3.5, '*', font=2, cex=2.5)

#------------------#

# Proline
stripchart(substrate~abx, data=proline_mock, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
stripchart(substrate~infection, data=proline_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=proline_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=proline_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=proline_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
axis(side=2, at=c(0:4), labels=c('0.0','1.0','2.0','3.0', '4.0'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.5, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.9)
mtext('B', side=2, line=2, las=2, adj=2, padj=-4, cex=1.3)
legend('topright', legend='Proline', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(proline_mock[,2]),
              median(subset(proline_strep, infection=='mock')[,2]), median(subset(proline_strep, infection=='630')[,2]),
              median(subset(proline_cef, infection=='mock')[,2]), median(subset(proline_cef, infection=='630')[,2]),
              median(subset(proline_clinda, infection=='mock')[,2]), median(subset(proline_clinda, infection=='630')[,2]),
              median(subset(proline_gf, infection=='mock')[,2]), median(subset(proline_gf, infection=='630')[,2])), 
         y1=c(median(proline_mock[,2]),
              median(subset(proline_strep, infection=='mock')[,2]), median(subset(proline_strep, infection=='630')[,2]),
              median(subset(proline_cef, infection=='mock')[,2]), median(subset(proline_cef, infection=='630')[,2]),
              median(subset(proline_clinda, infection=='mock')[,2]), median(subset(proline_clinda, infection=='630')[,2]),
              median(subset(proline_gf, infection=='mock')[,2]), median(subset(proline_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=c(3,6,9,12), y0=c(2.9,2.8,2.5,2.9), x1=c(4,7,10,13), y1=c(2.9,2.8,2.5,2.9), lwd=2)
text(x=c(3.5,6.5,9.5,12.5), y=c(3.2,3.1,2.8,3.2), '*', font=2, cex=2.5)

#------------------#

# Mannitol / Sorbitol
stripchart(substrate~abx, data=mannitolsorbitol_mock, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,70), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
stripchart(substrate~infection, data=mannitolsorbitol_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,70), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=mannitolsorbitol_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,70), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=mannitolsorbitol_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,70), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=mannitolsorbitol_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,70), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
axis(side=2, at=seq(0,70,17.5), labels=c('0.0','17.5','35','52.5','70.0'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.5, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.9)
mtext('C', side=2, line=2, las=2, adj=2, padj=-4, cex=1.3)
legend('topright', legend='Mannitol / Sorbitol', pt.cex=0, cex=1.2, bty='n')
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(mannitolsorbitol_mock[,2]),
              median(subset(mannitolsorbitol_strep, infection=='mock')[,2]), median(subset(mannitolsorbitol_strep, infection=='630')[,2]),
              median(subset(mannitolsorbitol_cef, infection=='mock')[,2]), median(subset(mannitolsorbitol_cef, infection=='630')[,2]),
              median(subset(mannitolsorbitol_clinda, infection=='mock')[,2]), median(subset(mannitolsorbitol_clinda, infection=='630')[,2]),
              median(subset(mannitolsorbitol_gf, infection=='mock')[,2]), median(subset(mannitolsorbitol_gf, infection=='630')[,2])), 
         y1=c(median(mannitolsorbitol_mock[,2]),
              median(subset(mannitolsorbitol_strep, infection=='mock')[,2]), median(subset(mannitolsorbitol_strep, infection=='630')[,2]),
              median(subset(mannitolsorbitol_cef, infection=='mock')[,2]), median(subset(mannitolsorbitol_cef, infection=='630')[,2]),
              median(subset(mannitolsorbitol_clinda, infection=='mock')[,2]), median(subset(mannitolsorbitol_clinda, infection=='630')[,2]),
              median(subset(mannitolsorbitol_gf, infection=='mock')[,2]), median(subset(mannitolsorbitol_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=12, y0=45, x1=13, y1=45, lwd=2)
text(x=12.5, y=50, '*', font=2, cex=2.5)

#------------------#

# Glucose
stripchart(substrate~abx, data=glucose_mock, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,26), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
stripchart(substrate~infection, data=glucose_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,26), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=glucose_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,26), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=glucose_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,26), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=glucose_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,26), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
axis(side=2, at=seq(0,70,17.5), labels=c('0.0','6.5','13.0','19.5','26.0'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.5, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.9)
mtext('C', side=2, line=2, las=2, adj=2, padj=-4, cex=1.3)
legend('topright', legend='Glucose', pt.cex=0, cex=1.2, bty='n')
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(glucose_mock[,2]),
              median(subset(glucose_strep, infection=='mock')[,2]), median(subset(glucose_strep, infection=='630')[,2]),
              median(subset(glucose_cef, infection=='mock')[,2]), median(subset(glucose_cef, infection=='630')[,2]),
              median(subset(glucose_clinda, infection=='mock')[,2]), median(subset(glucose_clinda, infection=='630')[,2]),
              median(subset(glucose_gf, infection=='mock')[,2]), median(subset(glucose_gf, infection=='630')[,2])), 
         y1=c(median(glucose_mock[,2]),
              median(subset(glucose_strep, infection=='mock')[,2]), median(subset(glucose_strep, infection=='630')[,2]),
              median(subset(glucose_cef, infection=='mock')[,2]), median(subset(glucose_cef, infection=='630')[,2]),
              median(subset(glucose_clinda, infection=='mock')[,2]), median(subset(glucose_clinda, infection=='630')[,2]),
              median(subset(glucose_gf, infection=='mock')[,2]), median(subset(glucose_gf, infection=='630')[,2])),
         lwd=3)

#------------------#

# N-acetylneuraminate
stripchart(substrate~abx, data=acetylneuraminate_mock, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
stripchart(substrate~infection, data=acetylneuraminate_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylneuraminate_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylneuraminate_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=acetylneuraminate_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,4), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
axis(side=2, at=c(0:4), labels=c('0.0','1.0','2.0','3.0', '4.0'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.5, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.9)
mtext('F', side=2, line=2, las=2, adj=2, padj=-4, cex=1.3)
legend('topright', legend='Neu5Ac', pt.cex=0, bty='n', cex=1.2)
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(acetylneuraminate_mock[,2]),
              median(subset(acetylneuraminate_strep, infection=='mock')[,2]), median(subset(acetylneuraminate_strep, infection=='630')[,2]),
              median(subset(acetylneuraminate_cef, infection=='mock')[,2]), median(subset(acetylneuraminate_cef, infection=='630')[,2]),
              median(subset(acetylneuraminate_clinda, infection=='mock')[,2]), median(subset(acetylneuraminate_clinda, infection=='630')[,2]),
              median(subset(acetylneuraminate_gf, infection=='mock')[,2]), median(subset(acetylneuraminate_gf, infection=='630')[,2])), 
         y1=c(median(acetylneuraminate_mock[,2]),
              median(subset(acetylneuraminate_strep, infection=='mock')[,2]), median(subset(acetylneuraminate_strep, infection=='630')[,2]),
              median(subset(acetylneuraminate_cef, infection=='mock')[,2]), median(subset(acetylneuraminate_cef, infection=='630')[,2]),
              median(subset(acetylneuraminate_clinda, infection=='mock')[,2]), median(subset(acetylneuraminate_clinda, infection=='630')[,2]),
              median(subset(acetylneuraminate_gf, infection=='mock')[,2]), median(subset(acetylneuraminate_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=c(6,12), y0=c(2.3,2.3), x1=c(7,13), y1=c(2.3,2.3), lwd=2)
text(x=c(6.5,12.5), y=c(2.6,2.6), '*', font=2, cex=2.5)

#------------------#

# Fructose
stripchart(substrate~abx, data=fructose_mock, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,10), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.25, cex.lab=1.2)
stripchart(substrate~infection, data=fructose_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,10), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=fructose_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,10), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=fructose_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,10), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=fructose_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,10), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.25, cex.lab=1.2, add=TRUE)
axis(side=2, at=seq(0,70,17.5), labels=c('0.0','2.5','5.0','7.5','10.0'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.8)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.5, cex=1.2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-GF'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.9)
mtext('C', side=2, line=2, las=2, adj=2, padj=-4, cex=1.3)
legend('topright', legend='Fructose', pt.cex=0, cex=1.2, bty='n')
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(fructose_mock[,2]),
              median(subset(fructose_strep, infection=='mock')[,2]), median(subset(fructose_strep, infection=='630')[,2]),
              median(subset(fructose_cef, infection=='mock')[,2]), median(subset(fructose_cef, infection=='630')[,2]),
              median(subset(fructose_clinda, infection=='mock')[,2]), median(subset(fructose_clinda, infection=='630')[,2]),
              median(subset(fructose_gf, infection=='mock')[,2]), median(subset(fructose_gf, infection=='630')[,2])), 
         y1=c(median(fructose_mock[,2]),
              median(subset(fructose_strep, infection=='mock')[,2]), median(subset(fructose_strep, infection=='630')[,2]),
              median(subset(fructose_cef, infection=='mock')[,2]), median(subset(fructose_cef, infection=='630')[,2]),
              median(subset(fructose_clinda, infection=='mock')[,2]), median(subset(fructose_clinda, infection=='630')[,2]),
              median(subset(fructose_gf, infection=='mock')[,2]), median(subset(fructose_gf, infection=='630')[,2])),
         lwd=3)
segments(x0=c(3,6,9,12), y0=c(2.9,2.8,2.5,2.9), x1=c(4,7,10,13), y1=c(2.9,2.8,2.5,2.9), lwd=2)
text(x=c(3.5,6.5,9.5,12.5), y=c(3.2,3.1,2.8,3.2), '*', font=2, cex=2.5)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
#detach('package:wesanderson', character.only = TRUE)
#rm(list=ls())
#gc()

