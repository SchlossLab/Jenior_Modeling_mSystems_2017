
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

# Format columns
cef_importance$Compound_name <- as.character(cef_importance$Compound_name)
clinda_importance$Compound_name <- as.character(clinda_importance$Compound_name)
strep_importance$Compound_name <- as.character(strep_importance$Compound_name)
gf_importance$Compound_name <- as.character(gf_importance$Compound_name)

# Subset substrates outside of randomized confidence interval
cef_importance <- subset(cef_importance, cef_importance$Std_from_Sim_Mean > 1)
clinda_importance <- subset(clinda_importance, clinda_importance$Std_from_Sim_Mean > 1)
strep_importance <- subset(strep_importance, strep_importance$Std_from_Sim_Mean > 1)
gf_importance <- subset(gf_importance, gf_importance$Std_from_Sim_Mean > 1)

# Rank metabolite importance scores
cef_importance <- cef_importance[order(cef_importance$Metabolite_score),]
clinda_importance <- clinda_importance[order(clinda_importance$Metabolite_score),]
strep_importance <- strep_importance[order(strep_importance$Metabolite_score),]
gf_importance <- gf_importance[order(gf_importance$Metabolite_score),]

# Replace long compound names with shorter versions
cef_importance$Compound_name <- replace(cef_importance$Compound_name, cef_importance$Compound_name=='[Dihydrolipoyllysine-residue_acetyltransferase]_S-acetyldihydrolipoyllysine', 'Acetyldihydrolipoyllysine')
strep_importance$Compound_name <- replace(strep_importance$Compound_name, strep_importance$Compound_name=='3-Hydroxy-3-(4-methylpent-3-en-1-yl)glutaryl-CoA', '3-Hydroxy-3-isohexeneylglutaryl-CoA')

# Replace '_' with ' '
cef_importance$Compound_name <- gsub('_', ' ', cef_importance$Compound_name)
clinda_importance$Compound_name <- gsub('_', ' ', clinda_importance$Compound_name)
strep_importance$Compound_name <- gsub('_', ' ',strep_importance$Compound_name)
gf_importance$Compound_name <- gsub('_', ' ', gf_importance$Compound_name)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figure_S4.pdf'
pdf(file=plot_file, width=10, height=12)
layout(matrix(c(1,2,3), nrow=3, ncol=1, byrow=FALSE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Plot it!

# A - Cefoperazone
par(mar=c(4,3,1,1), yaxs='i', xaxs='i')
dotchart(cef_importance$Metabolite_score, labels=cef_importance$Compound_name, 
         xlab='Metabolite Importance Score', xlim=c(-14,14), pch=19,
         col=c('black','black','black','black','red','black','red',
               'red','black','black','black','black','black','black'))
segments(x0=rep(-14,14), y0=c(1:14), x1=rep(14,14), y1=c(1:14), lty=2)
abline(v=0)

# Add simulated means and standard deviation
#points(x=cef_importance$Sim_Mean, y=c(1:14), cex=2, pch=4) # mean
points(x=cef_importance$Sim_Mean - (2*cef_importance$Sim_StD), y=c(1:14), cex=1.6, pch='|')
points(x=cef_importance$Sim_Mean + (2*cef_importance$Sim_StD), y=c(1:14), cex=1.6, pch='|')
points(x=cef_importance$Sim_Mean - (3*cef_importance$Sim_StD), y=c(1:14), cex=1.6, pch='|')
points(x=cef_importance$Sim_Mean + (3*cef_importance$Sim_StD), y=c(1:14), cex=1.6, pch='|')
segments(x0=cef_importance$Sim_Mean - (3*cef_importance$Sim_StD), y0=c(1:14), 
         x1=cef_importance$Sim_Mean + (3*cef_importance$Sim_StD), y1=c(1:14))
mtext('A', side=2, line=2, las=2, adj=0.5, padj=-10, cex=1.3)

#------------------#

# B - Clindamycin
par(mar=c(4,3,1,1), yaxs='i', xaxs='i')
dotchart(clinda_importance$Metabolite_score, labels=clinda_importance$Compound_name, 
         xlab='Metabolite Importance Score', xlim=c(-14,14), pch=19,
         col=c('black','black','black','black','red'))
segments(x0=rep(-14, 5), y0=c(1:5), x1=rep(14, 5), y1=c(1:5), lty=2)
abline(v=0)

# Add simulated means and standard deviation
#points(x=clinda_importance$Sim_Mean, y=c(1:5), cex=2, pch=4) # mean
points(x=clinda_importance$Sim_Mean - (2*clinda_importance$Sim_StD), y=c(1:5), cex=1.6, pch='|')
points(x=clinda_importance$Sim_Mean + (2*clinda_importance$Sim_StD), y=c(1:5), cex=1.6, pch='|')
points(x=clinda_importance$Sim_Mean - (3*clinda_importance$Sim_StD), y=c(1:5), cex=1.6, pch='|')
points(x=clinda_importance$Sim_Mean + (3*clinda_importance$Sim_StD), y=c(1:5), cex=1.6, pch='|')
segments(x0=clinda_importance$Sim_Mean - (3*clinda_importance$Sim_StD), 
         y0=c(1:5), x1=clinda_importance$Sim_Mean + (3*clinda_importance$Sim_StD), y1=c(1:5))
mtext('B', side=2, line=2, las=2, adj=0.5, padj=-10, cex=1.3)

#------------------#

# C - Streptomycin
par(mar=c(4,3,1,1), yaxs='i', xaxs='i')
dotchart(strep_importance$Metabolite_score, labels=strep_importance$Compound_name, 
         xlab='Metabolite Importance Score', xlim=c(-14,14), pch=19)
segments(x0=rep(-14, 17), y0=c(1:17), x1=rep(14, 17), y1=c(1:17), lty=2)
abline(v=0)

# Add simulated means and standard deviation
#points(x=strep_importance$Sim_Mean, y=c(1:17), cex=2, pch=4) # mean
points(x=strep_importance$Sim_Mean - (2*strep_importance$Sim_StD), y=c(1:17), cex=1.6, pch='|')
points(x=strep_importance$Sim_Mean + (2*strep_importance$Sim_StD), y=c(1:17), cex=1.6, pch='|')
points(x=strep_importance$Sim_Mean - (3*strep_importance$Sim_StD), y=c(1:17), cex=1.6, pch='|')
points(x=strep_importance$Sim_Mean + (3*strep_importance$Sim_StD), y=c(1:17), cex=1.6, pch='|')
segments(x0=strep_importance$Sim_Mean - (3*strep_importance$Sim_StD), y0=c(1:17), 
         x1=strep_importance$Sim_Mean + (3*strep_importance$Sim_StD), y1=c(1:17))
mtext('C', side=2, line=2, las=2, adj=0.5, padj=-10, cex=1.3)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
#rm(cef_importance, clinda_importance, strep_importance, gf_importance, plot_file)


