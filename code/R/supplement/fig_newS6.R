
# Start with clean environment
rm(list=ls())
gc()

# Load dependencies
deps <- c('vegan', 'ggplot2', 'shape', 'wesanderson', 'matrixStats', 'flux');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
set.seed(42)

# Define function for mormatting growth curves for statistical analysis
format_curve <- function(raw_exp_data, exp_group, raw_control_data){
  formatted_data <- c()
  control_data <- c()
  for (time in 1:nrow(raw_exp_data)){
    temp_exp <- cbind(exp_group, time, raw_exp_data[time,])
    formatted_data <- rbind(formatted_data, temp_exp)
    temp_control <- cbind('control', time, raw_control_data[time,])
    control_data <- rbind(control_data, temp_control)
  }
  formatted_data <- as.data.frame(rbind(control_data, formatted_data))
  colnames(formatted_data) <- c('substrate','time','od')
  formatted_data$od <- as.numeric(as.character(formatted_data$od))
  
  return(formatted_data)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in substrate importance data
cef_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/importances.tsv'
clinda_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/clindamycin_630.bipartite.files/importances.tsv'
strep_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/streptomycin_630.bipartite.files/importances.tsv'

cef_importance <- read.delim(cef_importance_file, header=TRUE, sep='\t', row.names=1)
clinda_importance <- read.delim(clinda_importance_file, header=TRUE, sep='\t', row.names=1)
strep_importance <- read.delim(strep_importance_file, header=TRUE, sep='\t', row.names=1)
rm(cef_importance_file, clinda_importance_file, strep_importance_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format metabolite importance scores
cef_importance$Metabolite_score <- as.numeric(as.character(cef_importance$Metabolite_score))
clinda_importance$Metabolite_score <- as.numeric(as.character(clinda_importance$Metabolite_score))
strep_importance$Metabolite_score <- as.numeric(as.character(strep_importance$Metabolite_score))

# Remove non-significant scores
cef_importance <- as.data.frame(subset(cef_importance, (cef_importance[,3] != 'n.s.')))
clinda_importance <- as.data.frame(subset(clinda_importance, (clinda_importance[,3] != 'n.s.')))
strep_importance <- as.data.frame(subset(strep_importance, (strep_importance[,3] != 'n.s.')))

# Sort for least important, subset, and format names
cef_least_importance <- cef_importance[order(cef_importance$Metabolite_score),][c(1:10),]
clinda_least_importance <- clinda_importance[order(clinda_importance$Metabolite_score),][c(1:10),]
strep_least_importance <- strep_importance[order(strep_importance$Metabolite_score),][c(1:10),]
cef_least_importance$Compound_name <- gsub('_',' ',cef_least_importance$Compound_name)
clinda_least_importance$Compound_name <- gsub('_',' ',clinda_least_importance$Compound_name)
strep_least_importance$Compound_name <- gsub('_',' ',strep_least_importance$Compound_name)

# Sort for most important
cef_importance <- cef_importance[order(-cef_importance$Metabolite_score),]
clinda_importance <- clinda_importance[order(-clinda_importance$Metabolite_score),]
strep_importance <- strep_importance[order(-strep_importance$Metabolite_score),]

# Take top 50 scores
cef_importance <- cef_importance[c(1:50),]
clinda_importance <- clinda_importance[c(1:50),]
strep_importance <- strep_importance[c(1:50),]

# find shared important metabolites
shared_importance <- as.data.frame(subset(cef_importance, (cef_importance[,1] %in% clinda_importance[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% clinda_importance[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% strep_importance[,1])))
shared_importance <- shared_importance$Compound_name
shared_cef <- as.data.frame(subset(cef_importance, (cef_importance[,1] %in% shared_importance)))
shared_cef <- shared_cef[order(shared_cef$Compound_name),]
shared_cef$Metabolite_score <- as.numeric(as.character(shared_cef$Metabolite_score))
shared_clinda <- as.data.frame(subset(clinda_importance, (clinda_importance[,1] %in% shared_importance)))
shared_clinda <- shared_clinda[order(shared_clinda$Compound_name),]
shared_clinda$Metabolite_score <- as.numeric(as.character(shared_clinda$Metabolite_score))
shared_strep <- as.data.frame(subset(strep_importance, (strep_importance[,1] %in% shared_importance)))
shared_strep <- shared_strep[order(shared_strep$Compound_name),]
shared_strep$Metabolite_score <- as.numeric(as.character(shared_strep$Metabolite_score))

score_median <- as.data.frame(apply(cbind(shared_cef$Metabolite_score, shared_clinda$Metabolite_score, shared_strep$Metabolite_score), 1, median))
shared_importance <- cbind(shared_cef$Compound_name, score_median)
rownames(shared_importance) <- rownames(shared_cef)
colnames(shared_importance) <- c('Compound_name','Metabolite_score')
shared_importance <- shared_importance[order(shared_importance$Metabolite_score),]
rm(shared_cef, shared_clinda, shared_strep, score_median)

# Subset to most important, distinct metabolites
cef_importance <- cef_importance[c(1:25),]
clinda_importance <- clinda_importance[c(1:25),]
strep_importance <- strep_importance[c(1:25),]
cef_only_importance <- as.data.frame(subset(cef_importance, !(cef_importance[,1] %in% clinda_importance[,1])))
cef_only_importance <- as.data.frame(subset(cef_only_importance, !(cef_only_importance[,1] %in% strep_importance[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_importance, !(clinda_importance[,1] %in% cef_importance[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_only_importance, !(clinda_only_importance[,1] %in% strep_importance[,1])))
strep_only_importance <- as.data.frame(subset(strep_importance, !(strep_importance[,1] %in% clinda_importance[,1])))
strep_only_importance <- as.data.frame(subset(strep_only_importance, !(strep_only_importance[,1] %in% cef_importance[,1])))
rm(cef_importance, clinda_importance, strep_importance)
cef_only_importance <- cef_only_importance[order(cef_only_importance$Metabolite_score),]
clinda_only_importance <- clinda_only_importance[order(clinda_only_importance$Metabolite_score),]
strep_only_importance <- strep_only_importance[order(strep_only_importance$Metabolite_score),]
cef_only_importance$abx <- 'Cefoperazone (SPF)'
clinda_only_importance$abx <- 'Clindamycin (SPF)'
strep_only_importance$abx <- 'Streptomycin (SPF)'
top_importances <- rbind(cef_only_importance, clinda_only_importance, strep_only_importance)
top_importances$abx <- as.factor(top_importances$abx)
top_importances$abx <- ordered(top_importances$abx, levels=c('Streptomycin (SPF)', 'Cefoperazone (SPF)', 'Clindamycin (SPF)'))
rm(cef_only_importance, clinda_only_importance, strep_only_importance)

# Format names to look better for the plot
top_importances$Compound_name <- gsub('_',' ',top_importances$Compound_name)
top_importances$Compound_name <- gsub('monophosphate','p',top_importances$Compound_name)
top_importances$Compound_name <- gsub('phosphate','p',top_importances$Compound_name)
top_importances$Compound_name <- gsub('alpha-','',top_importances$Compound_name)
top_importances$Compound_name <- gsub('beta-','',top_importances$Compound_name)
top_importances$Compound_name <- gsub('\\(R\\)\\-','',top_importances$Compound_name)
top_importances$Compound_name <- gsub('\\(S\\)\\-','',top_importances$Compound_name)
top_importances$Compound_name[top_importances$Compound_name == '1-(5\'-Phosphoribosyl)-5-amino-4-(N-succinocarboxamide)-imidazole'] <- 'SAICAR' # shorten a long name
top_importances$Compound_name[top_importances$Compound_name == '5,6,7,8-Tetrahydromethanopterin'] <- 'THMPT' # shorten a long name
shared_importance$Compound_name <- gsub('_',' ',shared_importance$Compound_name)
shared_importance$Compound_name <- gsub('phosphate','p',shared_importance$Compound_name)
shared_importance$Compound_name <- gsub('alpha,alpha\'-','',shared_importance$Compound_name)
shared_importance$Compound_name <- gsub('\\(R\\)\\-','',shared_importance$Compound_name)
shared_importance$Compound_name[shared_importance$Compound_name == 'CO2'] <- expression(CO[2])

# Remove generic KEGG anotations
top_importances <- as.data.frame(subset(top_importances, !(top_importances$Compound_name %in% c('Malonyl-[acp] methyl ester','Glutaryl-[acp] methyl ester'))))

# Remove overlaps from distinct lists
top_importances <- as.data.frame(subset(top_importances, !(top_importances$Compound_name %in% shared_importance$Compound_name)))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Metabolome 

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
mock_metabolome <- subset(metabolome, infection == 'mock')
mock_metabolome$infection <- NULL
rm(metadata)

# Subset metabolomics - untreated vs mock infected comparison
# Saccharides
sialyllactose_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('3-sialyllactose','6\'-sialyllactose')))]
sialyllactose_mock$abx <- factor(sialyllactose_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
fructose_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('fructose')))]
fructose_mock$abx <- factor(fructose_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
sucrose_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('sucrose')))]
sucrose_mock$abx <- factor(sucrose_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
galactitol_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('galactitol_(dulcitol)')))]
galactitol_mock$abx <- factor(galactitol_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
mannose_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('mannose')))]
mannose_mock$abx <- factor(mannose_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
glucose_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('glucose','glucose_6-phosphate')))]
glucose_mock$abx <- factor(glucose_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))

# Amino acids
glycine_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('glycine')))]
glycine_mock$abx <- factor(glycine_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
hydroxyproline_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('trans-4-hydroxyproline')))]
hydroxyproline_mock$abx <- factor(hydroxyproline_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
prolylglycine_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('prolylglycine')))]
prolylglycine_mock$abx <- factor(prolylglycine_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))

# Aminoglycans
acetylmuramate_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('N-acetylmuramate')))]
acetylmuramate_mock$abx <- factor(acetylmuramate_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
acetylglucosamine_sulfate_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('N-acetylglucosamine_6-sulfate')))]
acetylglucosamine_sulfate_mock$abx <- factor(acetylglucosamine_sulfate_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
acetylglucosaminylasparagine_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('N-acetylglucosaminylasparagine')))]
acetylglucosaminylasparagine_mock$abx <- factor(acetylglucosaminylasparagine_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
acetylbetaglucosaminylamine_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('N-acetyl-beta-glucosaminylamine')))]
acetylbetaglucosaminylamine_mock$abx <- factor(acetylbetaglucosaminylamine_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
rm(mock_metabolome)

#--------------------#

# Subset metabolites - mock vs infected comparison
# Saccharides



# Amino acids



# Aminoglycans
acetylmuramate <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-acetylmuramate')))]
acetylmuramate_strep <- subset(acetylmuramate, abx == 'streptomycin')
acetylmuramate_strep$abx <- NULL
acetylmuramate_strep$infection <- factor(acetylmuramate_strep$infection, levels=c('mock','630'))
acetylmuramate_cef <- subset(acetylmuramate, abx == 'cefoperazone')
acetylmuramate_cef$abx <- NULL
acetylmuramate_cef$infection <- factor(acetylmuramate_cef$infection, levels=c('mock','630'))
acetylmuramate_clinda <- subset(acetylmuramate, abx == 'clindamycin')
acetylmuramate_clinda$abx <- NULL
acetylmuramate_clinda$infection <- factor(acetylmuramate_clinda$infection, levels=c('mock','630'))
acetylmuramate_gf <- subset(acetylmuramate, abx == 'germfree')
acetylmuramate_gf$abx <- NULL
acetylmuramate_gf$infection <- factor(acetylmuramate_gf$infection, levels=c('mock','630'))
rm(acetylmuramate)








#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S6.pdf'
pdf(file=plot_file, width=10, height=8)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))

#---------------------------------------#

# Shared metabolite importances
par(mar=c(4,5,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(shared_importance$Metabolite_score, labels=shared_importance$Compound_name, lcolor=NA, cex=1.2, color='black',
         xlab='Median Importance Score', xlim=c(0,10), pch=19, lwd=3)
segments(x0=rep(0,10), y0=c(1:27), x1=rep(12,10), y1=c(1:27), lty=2) # Dotted lines
mtext('A', side=2, line=2, las=2, adj=4, padj=-18, cex=1.5)

#---------------------------------------#

# Distinct metabolite importances
par(mar=c(4,4,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(top_importances$Metabolite_score, labels=top_importances$Compound_name,
         lcolor=NA, cex=1.2, groups=top_importances$abx, color='black',
         xlab='Importance Score', xlim=c(0,10), pch=19, lwd=3,
         gcolor=c(wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5]))
mtext('B', side=2, line=2, las=2, adj=2.3, padj=-18, cex=1.5)
segments(x0=rep(0, 15), y0=c(1:4, 7:8, 11:13), 
         x1=rep(12, 15), y1=c(1:4, 7:8, 11:13), lty=2) # Dotted lines

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()
