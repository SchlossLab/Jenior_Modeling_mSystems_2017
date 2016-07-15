
deps <- c('vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in substrate importance data
cef_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/cefoperazone_630.monte_carlo.score.txt'
clinda_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/clindamycin_630.bipartite.files/clindamycin_630.monte_carlo.score.txt'
strep_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/streptomycin_630.bipartite.files/streptomycin_630.monte_carlo.score.txt'
gf_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/germfree.bipartite.files/germfree.monte_carlo.score.txt'

cef_importance <- read.table(cef_importance_file, header=TRUE, sep='\t', row.names=1)
clinda_importance <- read.table(clinda_importance_file, header=TRUE, sep='\t', row.names=1)
strep_importance <- read.table(strep_importance_file, header=TRUE, sep='\t', row.names=1)
gf_importance <- read.table(gf_importance_file, header=TRUE, sep='\t', row.names=1)
rm(cef_importance_file, clinda_importance_file, strep_importance_file, gf_importance_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Subset substrates outside of randomized confidence interval
cef_importance <- subset(cef_importance, cef_importance$Std_from_Sim_Mean > 1)
clinda_importance <- subset(clinda_importance, clinda_importance$Std_from_Sim_Mean > 1)
strep_importance <- subset(strep_importance, strep_importance$Std_from_Sim_Mean > 1)
gf_importance <- subset(gf_importance, gf_importance$Std_from_Sim_Mean > 1)

# Rank metabolite importance scores
cef_importance <- cef_importance[order(-cef_importance$Metabolite_score),]
clinda_importance <- clinda_importance[order(-clinda_importance$Metabolite_score),]
strep_importance <- strep_importance[order(-strep_importance$Metabolite_score),]
gf_importance <- gf_importance[order(-gf_importance$Metabolite_score),]

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_6.pdf'
pdf(file=plot_file, width=15, height=7)
layout(matrix(c(1,2,3), nrow=2, ncol=3, byrow=FALSE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Plot it!

# A - Cefoperazone
par(mar=c(4,3,1,1), yaxs='i', xaxs='i')
dotchart(cef_importance$Metabolite_score, labels=cef_importance$Compound_name, 
         xlab='Metabolite Importance Score', xlim=c(-12,12), pch=19, cex=1.5)
segments(x0=rep(-12,14), y0=c(1:14), x1=rep(12,14), y1=c(1:14), lty=2)
abline(v=0)

# Add simulated means and standard deviation
points(x=cef_importance$Sim_Mean, y=c(1:14), cex=1.8, pch='|') # mean
points(x=cef_importance$Sim_Mean - cef_importance$Sim_StD, y=c(1:14), cex=1.6, pch='|') # lower std
points(x=cef_importance$Sim_Mean + cef_importance$Sim_StD, y=c(1:14), cex=1.6, pch='|') # upper std

mtext('A', side=2, line=2, las=2, adj=0.5, padj=-14.5, cex=1.8)

#------------------#

# B - Clindamycin
par(mar=c(4,3,1,1), yaxs='i', xaxs='i')
dotchart(clinda_importance$Metabolite_score, labels=clinda_importance$Compound_name, 
         xlab='Metabolite Importance Score', xlim=c(-12,12), pch=19, cex=1.5)
segments(x0=rep(-12, 5), y0=c(1:5), x1=rep(12, 5), y1=c(1:5), lty=2)
abline(v=0)

# Add simulated means and standard deviation
points(x=clinda_importance$Sim_Mean, y=c(1:5), cex=1.8, pch='|') # mean
points(x=clinda_importance$Sim_Mean - clinda_importance$Sim_StD, y=c(1:5), cex=1.6, pch='|') # lower std
points(x=clinda_importance$Sim_Mean + clinda_importance$Sim_StD, y=c(1:5), cex=1.6, pch='|') # upper std

mtext('B', side=2, line=2, las=2, adj=0.5, padj=-14.5, cex=1.8)

#------------------#

# C - Streptomycin
par(mar=c(4,3,1,1), yaxs='i', xaxs='i')
dotchart(strep_importance$Metabolite_score, labels=strep_importance$Compound_name, 
         xlab='Metabolite Importance Score', xlim=c(-12,12), pch=19, cex=1.5)
segments(x0=rep(-12, 17), y0=c(1:17), x1=rep(12, 17), y1=c(1:17), lty=2)
abline(v=0)

# Add simulated means and standard deviation
points(x=strep_importance$Sim_Mean, y=c(1:17), cex=1.8, pch='|') # mean
points(x=strep_importance$Sim_Mean - strep_importance$Sim_StD, y=c(1:17), cex=1.6, pch='|') # lower std
points(x=strep_importance$Sim_Mean + strep_importance$Sim_StD, y=c(1:17), cex=1.6, pch='|') # upper std

mtext('C', side=2, line=2, las=2, adj=0.5, padj=-14.5, cex=1.8)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
rm(cef_importance, clinda_importance, strep_importance, gf_importance, plot_file)


