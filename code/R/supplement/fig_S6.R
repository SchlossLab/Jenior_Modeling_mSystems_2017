
# Load dependencies
deps <- c('wesanderson', 'plyr');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Select files
concentrations <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/ms_substrates.tsv'
metadata <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'
importances <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/select_importances.tsv'

# Read in data
concentrations <- read.delim(concentrations, sep='\t', header=T, row.names=1)
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
importances <- read.delim(importances, sep='\t', header=T, row.names=1)

# Format and merge tables
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$type <- NULL
concentrations$SUPER_PATHWAY <- NULL
concentrations$KEGG <- NULL
concentrations <- t(concentrations)
concentrations <- merge(concentrations, metadata, by='row.names')
rownames(concentrations) <- concentrations$Row.names
concentrations$Row.names <- NULL
colnames(concentrations) <- c('acetylglucosamine', 'mannitol_sorbitol', 'proline', 'glycine', 
                              'salicylate', 'acetylneuraminate', 'galactitol', 'glucose', 'abx', 'infection')

# Sort by concentrations antibiotic treatment group
concentrations$abx <- factor(concentrations$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))

# Subset concentrations by infection group
concentrations_630 <- subset(concentrations, infection == '630')
concentrations_630$infection <- NULL
concentrations_mock <- subset(concentrations, infection == 'mock')
concentrations_mock$infection <- NULL
rm(concentrations)

# Calculate concentration medians
acetylglucosamine_630_median <- ddply(concentrations_630, .(abx), summarize, 'median'=median(acetylglucosamine))
rownames(acetylglucosamine_630_median) <- acetylglucosamine_630_median[,1]
acetylglucosamine_630_median[,1] <- NULL
colnames(acetylglucosamine_630_median) <- 'acetylglucosamine'
mannitolsorbitol_630_median <- ddply(concentrations_630, .(abx), summarize, 'median'=median(mannitol_sorbitol))
rownames(mannitolsorbitol_630_median) <- mannitolsorbitol_630_median[,1]
mannitolsorbitol_630_median[,1] <- NULL
colnames(mannitolsorbitol_630_median) <- 'mannitol_sorbitol'
proline_630_median <- ddply(concentrations_630, .(abx), summarize, 'median'=median(proline))
rownames(proline_630_median) <- proline_630_median[,1]
proline_630_median[,1] <- NULL
colnames(proline_630_median) <- 'proline'
glycine_630_median <- ddply(concentrations_630, .(abx), summarize, 'median'=median(glycine))
rownames(glycine_630_median) <- glycine_630_median[,1]
glycine_630_median[,1] <- NULL
colnames(glycine_630_median) <- 'glycine'
salicylate_630_median <- ddply(concentrations_630, .(abx), summarize, 'median'=median(salicylate))
rownames(salicylate_630_median) <- salicylate_630_median[,1]
salicylate_630_median[,1] <- NULL
colnames(salicylate_630_median) <- 'salicylate'
acetylneuraminate_630_median <- ddply(concentrations_630, .(abx), summarize, 'median'=median(acetylneuraminate))
rownames(acetylneuraminate_630_median) <- acetylneuraminate_630_median[,1]
acetylneuraminate_630_median[,1] <- NULL
colnames(acetylneuraminate_630_median) <- 'acetylneuraminate'
galactitol_630_median <- ddply(concentrations_630, .(abx), summarize, 'median'=median(galactitol))
rownames(galactitol_630_median) <- galactitol_630_median[,1]
galactitol_630_median[,1] <- NULL
colnames(galactitol_630_median) <- 'galactitol'
glucose_630_median <- ddply(concentrations_630, .(abx), summarize, 'median'=median(glucose))
rownames(glucose_630_median) <- glucose_630_median[,1]
glucose_630_median[,1] <- NULL
colnames(glucose_630_median) <- 'glucose'
concentrations_630_medians <- rbind(t(acetylglucosamine_630_median), t(mannitolsorbitol_630_median), t(proline_630_median), t(glycine_630_median), 
                                    t(salicylate_630_median), t(acetylneuraminate_630_median), t(galactitol_630_median), t(glucose_630_median))
rm(acetylglucosamine_630_median, proline_630_median, glycine_630_median, mannitolsorbitol_630_median, 
   salicylate_630_median, acetylneuraminate_630_median, galactitol_630_median, glucose_630_median)
acetylglucosamine_mock_median <- ddply(concentrations_mock, .(abx), summarize, 'median'=median(acetylglucosamine))
rownames(acetylglucosamine_mock_median) <- acetylglucosamine_mock_median[,1]
acetylglucosamine_mock_median[,1] <- NULL
colnames(acetylglucosamine_mock_median) <- 'acetylglucosamine'
mannitolsorbitol_mock_median <- ddply(concentrations_mock, .(abx), summarize, 'median'=median(mannitol_sorbitol))
rownames(mannitolsorbitol_mock_median) <- mannitolsorbitol_mock_median[,1]
mannitolsorbitol_mock_median[,1] <- NULL
colnames(mannitolsorbitol_mock_median) <- 'mannitol_sorbitol'
proline_mock_median <- ddply(concentrations_mock, .(abx), summarize, 'median'=median(proline))
rownames(proline_mock_median) <- proline_mock_median[,1]
proline_mock_median[,1] <- NULL
colnames(proline_mock_median) <- 'proline'
glycine_mock_median <- ddply(concentrations_mock, .(abx), summarize, 'median'=median(glycine))
rownames(glycine_mock_median) <- glycine_mock_median[,1]
glycine_mock_median[,1] <- NULL
colnames(glycine_mock_median) <- 'glycine'
salicylate_mock_median <- ddply(concentrations_mock, .(abx), summarize, 'median'=median(salicylate))
rownames(salicylate_mock_median) <- salicylate_mock_median[,1]
salicylate_mock_median[,1] <- NULL
colnames(salicylate_mock_median) <- 'salicylate'
acetylneuraminate_mock_median <- ddply(concentrations_mock, .(abx), summarize, 'median'=median(acetylneuraminate))
rownames(acetylneuraminate_mock_median) <- acetylneuraminate_mock_median[,1]
acetylneuraminate_mock_median[,1] <- NULL
colnames(acetylneuraminate_mock_median) <- 'acetylneuraminate'
galactitol_mock_median <- ddply(concentrations_mock, .(abx), summarize, 'median'=median(galactitol))
rownames(galactitol_mock_median) <- galactitol_mock_median[,1]
galactitol_mock_median[,1] <- NULL
colnames(galactitol_mock_median) <- 'galactitol'
glucose_mock_median <- ddply(concentrations_mock, .(abx), summarize, 'median'=median(glucose))
rownames(glucose_mock_median) <- glucose_mock_median[,1]
glucose_mock_median[,1] <- NULL
colnames(glucose_mock_median) <- 'glucose'
concentrations_mock_medians <- rbind(t(acetylglucosamine_mock_median), t(mannitolsorbitol_mock_median), t(proline_mock_median), t(glycine_mock_median), 
                                     t(salicylate_mock_median), t(acetylneuraminate_mock_median), t(galactitol_mock_median), t(glucose_mock_median))
rm(acetylglucosamine_mock_median, proline_mock_median, glycine_mock_median, mannitolsorbitol_mock_median, 
   salicylate_mock_median, acetylneuraminate_mock_median, galactitol_mock_median, glucose_mock_median)

# Combine importance with concentration median














# Subset each metabolites for boxplots
acetylglucosamine_630 <- concentrations_630[,c(1,9)]
mannitolsorbitol_630 <- concentrations_630[,c(2,9)]
galactitol_630 <- concentrations_630[,c(7,9)]
salicylate_630 <- concentrations_630[,c(5,9)]
acetylneuraminate_630 <- concentrations_630[,c(6,9)]
proline_630 <- concentrations_630[,c(3,9)]
glycine_630 <- concentrations_630[,c(4,9)]
glucose_630 <- concentrations_630[,c(8,9)]
acetylglucosamine_mock <- concentrations_mock[,c(1,9)]
mannitolsorbitol_mock <- concentrations_mock[,c(2,9)]
galactitol_mock <- concentrations_mock[,c(7,9)]
salicylate_mock <- concentrations_mock[,c(5,9)]
acetylneuraminate_mock <- concentrations_mock[,c(6,9)]
proline_mock <- concentrations_mock[,c(3,9)]
glycine_mock <- concentrations_mock[,c(4,9)]
glucose_mock <- concentrations_mock[,c(8,9)]
rm(concentrations_630, concentrations_mock, metadata)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Calculate significant differences for concentrations and correct p-values
acetylglucosamine_mockV630_p <- p.adjust(c(wilcox.test(subset(acetylglucosamine_mock, abx=='none')[,1], subset(acetylglucosamine_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                           wilcox.test(subset(acetylglucosamine_mock, abx=='none')[,1], subset(acetylglucosamine_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                           wilcox.test(subset(acetylglucosamine_mock, abx=='none')[,1], subset(acetylglucosamine_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                           wilcox.test(subset(acetylglucosamine_mock, abx=='none')[,1], subset(acetylglucosamine_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                                         method='BH')
names(acetylglucosamine_mockV630_p) <- c('conv_v_strep', 'conv_v_cef', 'conv_v_clinda', 'conv_v_gf')
acetylglucosamine_abx_p <- p.adjust(c(wilcox.test(subset(acetylglucosamine_mock, abx=='cefoperazone')[,1], subset(acetylglucosamine_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                           wilcox.test(subset(acetylglucosamine_mock, abx=='streptomycin')[,1], subset(acetylglucosamine_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                           wilcox.test(subset(acetylglucosamine_mock, abx=='clindamycin')[,1], subset(acetylglucosamine_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                           wilcox.test(subset(acetylglucosamine_mock, abx=='germfree')[,1], subset(acetylglucosamine_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                                         method='BH')
names(acetylglucosamine_abx_p) <- c('cef_mockV630', 'strep_mockV630', 'clinda_mockV630', 'gf_mockV630')
mannitolsorbitol_mockV630_p <- p.adjust(c(wilcox.test(subset(mannitolsorbitol_mock, abx=='none')[,1], subset(mannitolsorbitol_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                          wilcox.test(subset(mannitolsorbitol_mock, abx=='none')[,1], subset(mannitolsorbitol_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                          wilcox.test(subset(mannitolsorbitol_mock, abx=='none')[,1], subset(mannitolsorbitol_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                          wilcox.test(subset(mannitolsorbitol_mock, abx=='none')[,1], subset(mannitolsorbitol_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                                        method='BH')
names(mannitolsorbitol_mockV630_p) <- c('conv_v_strep', 'conv_v_cef', 'conv_v_clinda', 'conv_v_gf')
mannitolsorbitol_abx_p <- p.adjust(c(wilcox.test(subset(mannitolsorbitol_mock, abx=='cefoperazone')[,1], subset(mannitolsorbitol_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                     wilcox.test(subset(mannitolsorbitol_mock, abx=='streptomycin')[,1], subset(mannitolsorbitol_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                     wilcox.test(subset(mannitolsorbitol_mock, abx=='clindamycin')[,1], subset(mannitolsorbitol_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                     wilcox.test(subset(mannitolsorbitol_mock, abx=='germfree')[,1], subset(mannitolsorbitol_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                                   method='BH')
names(mannitolsorbitol_abx_p) <- c('cef_mockV630', 'strep_mockV630', 'clinda_mockV630', 'gf_mockV630')
galactitol_mockV630_p <- p.adjust(c(wilcox.test(subset(galactitol_mock, abx=='none')[,1], subset(galactitol_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                    wilcox.test(subset(galactitol_mock, abx=='none')[,1], subset(galactitol_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                    wilcox.test(subset(galactitol_mock, abx=='none')[,1], subset(galactitol_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                    wilcox.test(subset(galactitol_mock, abx=='none')[,1], subset(galactitol_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                                  method='BH')
names(galactitol_mockV630_p) <- c('conv_v_strep', 'conv_v_cef', 'conv_v_clinda', 'conv_v_gf')
galactitol_abx_p <- p.adjust(c(wilcox.test(subset(galactitol_mock, abx=='cefoperazone')[,1], subset(galactitol_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                               wilcox.test(subset(galactitol_mock, abx=='streptomycin')[,1], subset(galactitol_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                               wilcox.test(subset(galactitol_mock, abx=='clindamycin')[,1], subset(galactitol_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                               wilcox.test(subset(galactitol_mock, abx=='germfree')[,1], subset(galactitol_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                             method='BH')
names(galactitol_abx_p) <- c('cef_mockV630', 'strep_mockV630', 'clinda_mockV630', 'gf_mockV630')
salicylate_mockV630_p <- p.adjust(c(wilcox.test(subset(salicylate_mock, abx=='none')[,1], subset(salicylate_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                    wilcox.test(subset(salicylate_mock, abx=='none')[,1], subset(salicylate_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                    wilcox.test(subset(salicylate_mock, abx=='none')[,1], subset(salicylate_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                    wilcox.test(subset(salicylate_mock, abx=='none')[,1], subset(salicylate_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                                  method='BH')
names(salicylate_mockV630_p) <- c('conv_v_strep', 'conv_v_cef', 'conv_v_clinda', 'conv_v_gf')
salicylate_abx_p <- p.adjust(c(wilcox.test(subset(salicylate_mock, abx=='cefoperazone')[,1], subset(salicylate_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                               wilcox.test(subset(salicylate_mock, abx=='streptomycin')[,1], subset(salicylate_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                               wilcox.test(subset(salicylate_mock, abx=='clindamycin')[,1], subset(salicylate_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                               wilcox.test(subset(salicylate_mock, abx=='germfree')[,1], subset(salicylate_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                             method='BH')
names(salicylate_abx_p) <- c('cef_mockV630', 'strep_mockV630', 'clinda_mockV630', 'gf_mockV630')
acetylneuraminate_mockV630_p <- p.adjust(c(wilcox.test(subset(acetylneuraminate_mock, abx=='none')[,1], subset(acetylneuraminate_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                           wilcox.test(subset(acetylneuraminate_mock, abx=='none')[,1], subset(acetylneuraminate_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                           wilcox.test(subset(acetylneuraminate_mock, abx=='none')[,1], subset(acetylneuraminate_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                           wilcox.test(subset(acetylneuraminate_mock, abx=='none')[,1], subset(acetylneuraminate_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                                         method='BH')
names(acetylneuraminate_mockV630_p) <- c('conv_v_strep', 'conv_v_cef', 'conv_v_clinda', 'conv_v_gf')
acetylneuraminate_abx_p <- p.adjust(c(wilcox.test(subset(acetylneuraminate_mock, abx=='cefoperazone')[,1], subset(acetylneuraminate_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                      wilcox.test(subset(acetylneuraminate_mock, abx=='streptomycin')[,1], subset(acetylneuraminate_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                      wilcox.test(subset(acetylneuraminate_mock, abx=='clindamycin')[,1], subset(acetylneuraminate_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                      wilcox.test(subset(acetylneuraminate_mock, abx=='germfree')[,1], subset(acetylneuraminate_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                                    method='BH')
names(acetylneuraminate_abx_p) <- c('cef_mockV630', 'strep_mockV630', 'clinda_mockV630', 'gf_mockV630')
proline_mockV630_p <- p.adjust(c(wilcox.test(subset(proline_mock, abx=='none')[,1], subset(proline_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                 wilcox.test(subset(proline_mock, abx=='none')[,1], subset(proline_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                 wilcox.test(subset(proline_mock, abx=='none')[,1], subset(proline_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                 wilcox.test(subset(proline_mock, abx=='none')[,1], subset(proline_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                               method='BH')
names(proline_mockV630_p) <- c('conv_v_strep', 'conv_v_cef', 'conv_v_clinda', 'conv_v_gf')
proline_abx_p <- p.adjust(c(wilcox.test(subset(proline_mock, abx=='cefoperazone')[,1], subset(proline_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                            wilcox.test(subset(proline_mock, abx=='streptomycin')[,1], subset(proline_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                            wilcox.test(subset(proline_mock, abx=='clindamycin')[,1], subset(proline_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                            wilcox.test(subset(proline_mock, abx=='germfree')[,1], subset(proline_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                          method='BH')
names(proline_abx_p) <- c('cef_mockV630', 'strep_mockV630', 'clinda_mockV630', 'gf_mockV630')
glycine_mockV630_p <- p.adjust(c(wilcox.test(subset(glycine_mock, abx=='none')[,1], subset(glycine_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                 wilcox.test(subset(glycine_mock, abx=='none')[,1], subset(glycine_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                 wilcox.test(subset(glycine_mock, abx=='none')[,1], subset(glycine_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                 wilcox.test(subset(glycine_mock, abx=='none')[,1], subset(glycine_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                               method='BH')
names(glycine_mockV630_p) <- c('conv_v_strep', 'conv_v_cef', 'conv_v_clinda', 'conv_v_gf')
glycine_abx_p <- p.adjust(c(wilcox.test(subset(glycine_mock, abx=='cefoperazone')[,1], subset(glycine_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                            wilcox.test(subset(glycine_mock, abx=='streptomycin')[,1], subset(glycine_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                            wilcox.test(subset(glycine_mock, abx=='clindamycin')[,1], subset(glycine_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                            wilcox.test(subset(glycine_mock, abx=='germfree')[,1], subset(glycine_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                          method='BH')
names(glycine_abx_p) <- c('cef_mockV630', 'strep_mockV630', 'clinda_mockV630', 'gf_mockV630')
glucose_mockV630_p <- p.adjust(c(wilcox.test(subset(glucose_mock, abx=='none')[,1], subset(glucose_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                                 wilcox.test(subset(glucose_mock, abx=='none')[,1], subset(glucose_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                                 wilcox.test(subset(glucose_mock, abx=='none')[,1], subset(glucose_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                                 wilcox.test(subset(glucose_mock, abx=='none')[,1], subset(glucose_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                               method='BH')
names(glucose_mockV630_p) <- c('conv_v_strep', 'conv_v_cef', 'conv_v_clinda', 'conv_v_gf')
glucose_abx_p <- p.adjust(c(wilcox.test(subset(glucose_mock, abx=='cefoperazone')[,1], subset(glucose_630, abx=='cefoperazone')[,1], exact=FALSE)$p.value,
                            wilcox.test(subset(glucose_mock, abx=='streptomycin')[,1], subset(glucose_630, abx=='streptomycin')[,1], exact=FALSE)$p.value,
                            wilcox.test(subset(glucose_mock, abx=='clindamycin')[,1], subset(glucose_630, abx=='clindamycin')[,1], exact=FALSE)$p.value,
                            wilcox.test(subset(glucose_mock, abx=='germfree')[,1], subset(glucose_630, abx=='germfree')[,1], exact=FALSE)$p.value), 
                          method='BH')
names(glucose_abx_p) <- c('cef_mockV630', 'strep_mockV630', 'clinda_mockV630', 'gf_mockV630')

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_6.pdf'
select_palette <- c('gray', wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen')
pdf(file=plot_file, width=9, height=10)
layout(matrix(c(1,2,3,
                4,5,6,
                7,8,8),
              nrow=3, ncol=3, byrow = TRUE))

#--------------------------------#

# Acetylglucosamine
par(las=1, mar=c(0.2,4,1,1), mgp=c(2.5,0.7,0))
boxplot(Nacetylglucosamine_Nacetylgalactosamine~abx, data=acetylglucosamine, col=select_palette, ylim=c(0,14), whisklty=1, cex.lab=1.1,
           xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('A', side=2, line=2, las=2, adj=1.7, padj=-5)
legend('topright', 'N-Acetylglucosamine + N-Acetylgalactosamine', bty='n', cex=0.9)
axis(2, at=seq(0,14,3.5), labels=c('0.0','3.5','7.0','10.5','14.0'))

abline(v=2, lty=2) # Line separating resistant from susceptible



text(c(2,3,4,5), c(7,2.3,1.6,2.3), labels=c('*','*','*','*'), cex=2.5)

#--------------------------------#

# Salicylate
par(las=1, mar=c(0.2,4,0.2,1), mgp=c(2.5,0.7,0))
boxplot(salicylate~abx, data=salicylate, col=select_palette, ylim=c(0,5.2), whisklty=1, cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('B', side=2, line=2, las=2, adj=1.7, padj=-5)
legend('topright', 'Salicylate', bty='n', cex=0.9)
axis(2, at=seq(0,5.2,1.3), labels=c('0.0','1.3','2.6','3.9','5.2'))

abline(v=2, lty=2) # Line separating resistant from susceptible

text(c(2,3,4,5), c(2.8,1.8,1.9,0.6), labels=c('*','*','*','*'), cex=2.5)

# B-Glucose
par(las=1, mar=c(0.2,4,0.2,1), mgp=c(2.5,0.7,0))
boxplot(salicylate~abx, data=salicylate, col=select_palette, ylim=c(0,5.2), whisklty=1, cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('B', side=2, line=2, las=2, adj=1.7, padj=-5)
legend('topright', 'Salicylate', bty='n', cex=0.9)
axis(2, at=seq(0,5.2,1.3), labels=c('0.0','1.3','2.6','3.9','5.2'))

abline(v=2, lty=2) # Line separating resistant from susceptible

text(c(2,3,4,5), c(2.8,1.8,1.9,0.6), labels=c('*','*','*','*'), cex=2.5)







#--------------------------------#

# Galactitol
par(las=1, mar=c(0.2,4,0.2,1), mgp=c(2.5,0.7,0))
boxplot(galactitol~abx, data=galactitol, col=select_palette, ylim=c(0,3), whisklty=1, cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('C', side=2, line=2, las=2, adj=1.4, padj=-5)
legend('topright', 'Galactitol', bty='n', cex=0.9)
axis(2, at=seq(0,3,0.75), labels=c('0.0','0.75','1.5','2.25','3.0'))
text(c(3,5), c(2.7,2.35), labels=c('*','*'), cex=2.5)

abline(v=2, lty=2) # Line separating resistant from susceptible



#--------------------------------#

# Mannitol - Sorbitol
par(las=1, mar=c(0.2,4,0.2,1), mgp=c(2.5,0.7,0))
boxplot(mannitol_sorbitol~abx, data=mannitol_sorbitol, col=select_palette, ylim=c(0,50), whisklty=1, cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('D', side=2, line=2, las=2, adj=1.5, padj=-5)
legend('topright', 'Mannitol + Sorbitol', bty='n', cex=0.9)
axis(2, at=seq(0,50,12.5), labels=c('0.0','12.5','25.0','37.5','50.0'))
text(c(3,5), c(47,43), labels=c('*','*'), cex=2.5)

abline(v=2, lty=2) # Line separating resistant from susceptible



#--------------------------------#

# Acetylneuraminate
par(las=1, mar=c(0.2,4,0.2,1), mgp=c(2.5,0.7,0))
boxplot(Nacetylneuraminate~abx, data=acetylneuraminate, col=select_palette, ylim=c(0,3), whisklty=1, cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
mtext('E', side=2, line=2, las=2, adj=1.7, padj=-5)
legend('topright', 'N-Acetylneuraminate', bty='n', cex=0.9)
axis(2, at=seq(0,3,0.75), labels=c('0.0','0.75','1.5','2.25','3.0'))
text(c(3,5), c(2.2,2.3), labels=c('*','*'), cex=2.5)

abline(v=2, lty=2) # Line separating resistant from susceptible



#--------------------------------#

# Proline
par(las=1, mar=c(0.2,4,0.2,1), mgp=c(2.5,0.7,0))
boxplot(proline~abx, data=proline, col=select_palette, whisklty=1, ylim=c(0,3), cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
legend('topright', 'Proline', bty='n', cex=0.9)
mtext('F', side=2, line=2, las=2, adj=1.3, padj=-5)
axis(2, at=seq(0,3,0.75), labels=c('0.0','0.75','1.5','2.25','3.0'))
text(c(2,3,4,5), c(3.0,2.8,2.3,2.85), labels=c('*','*','*','*'), cex=2.5)

abline(v=2, lty=2) # Line separating resistant from susceptible


# Glycine
par(las=1, mar=c(0.2,4,0.2,1), mgp=c(2.5,0.7,0))
boxplot(glycine~abx, data=glycine, col=select_palette, whisklty=1, ylim=c(0,3), cex.lab=1.1,
        xaxt='n', yaxt='n', ylab='Scaled Intensity', boxlwd=2, whisklwd=2, staplelwd=2, outline=FALSE, range=0, medlwd=2)
legend('topright', 'Glycine', bty='n', cex=0.9)
mtext('G', side=2, line=2, las=2, adj=1.3, padj=-5)
axis(2, at=seq(0,3,0.75), labels=c('0.0','0.75','1.5','2.25','3.0'))
text(c(3,5), c(2.5,2.1), labels=c('*','*'), cex=2.5)

abline(v=2, lty=2) # Line separating resistant from susceptible











mtext(c('No Antibiotics\nSPF','Streptomycin\nSPF','Cefoperazone\nSPF','Clindamycin\nSPF','No Antibiotics\nGF'), side=1, at=c(1:5), cex=0.77, padj=1)
mtext('Treatment:', side=1, at=0.14, padj=1.3, cex=0.7)
mtext('Mice:', side=1, at=0.12, padj=3.1, cex=0.7)

#--------------------------------#



# NEED A PANEL ADDRESSING CORRELATION WITH IMPORTANCE SCORES
par(las=1, mar=c(0.2,4,0.2,1), mgp=c(2.5,0.7,0))
plot(1, xlab='Metabolite Concentration', ylab='Importance Score')

# Each treatment group a separate colored line
legend('topright')


mtext('H', side=2, line=2, las=2, adj=1.3, padj=-5)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

