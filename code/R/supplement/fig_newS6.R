
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

#-------------------------------------------------------------------------------------------------------------------------------------#

# Subset metabolomics - untreated vs mock infected comparison
# Saccharides
sucrose_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('sucrose')))]
sucrose_mock$abx <- factor(sucrose_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
galactitol_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('galactitol_(dulcitol)')))]
galactitol_mock$abx <- factor(galactitol_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
mannose_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('mannose')))]
mannose_mock$abx <- factor(mannose_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))
sialyllactose_mock <- mock_metabolome[, c(1,which(colnames(mock_metabolome) %in% c('3-sialyllactose','6\'-sialyllactose')))]
sialyllactose_mock$abx <- factor(sialyllactose_mock$abx, levels=c('none','streptomycin','cefoperazone','clindamycin','germfree'))

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
fructose <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-fructose')))]
fructose_strep <- subset(fructose, abx == 'streptomycin')
fructose_strep$abx <- NULL
fructose_strep$infection <- factor(fructose_strep$infection, levels=c('mock','630'))
fructose_cef <- subset(fructose, abx == 'cefoperazone')
fructose_cef$abx <- NULL
fructose_cef$infection <- factor(fructose_cef$infection, levels=c('mock','630'))
fructose_clinda <- subset(fructose, abx == 'clindamycin')
fructose_clinda$abx <- NULL
fructose_clinda$infection <- factor(fructose_clinda$infection, levels=c('mock','630'))
fructose_gf <- subset(fructose, abx == 'germfree')
fructose_gf$abx <- NULL
fructose_gf$infection <- factor(fructose_gf$infection, levels=c('mock','630'))
rm(fructose)
sucrose <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-sucrose')))]
sucrose_strep <- subset(sucrose, abx == 'streptomycin')
sucrose_strep$abx <- NULL
sucrose_strep$infection <- factor(sucrose_strep$infection, levels=c('mock','630'))
sucrose_cef <- subset(sucrose, abx == 'cefoperazone')
sucrose_cef$abx <- NULL
sucrose_cef$infection <- factor(sucrose_cef$infection, levels=c('mock','630'))
sucrose_clinda <- subset(sucrose, abx == 'clindamycin')
sucrose_clinda$abx <- NULL
sucrose_clinda$infection <- factor(sucrose_clinda$infection, levels=c('mock','630'))
sucrose_gf <- subset(sucrose, abx == 'germfree')
sucrose_gf$abx <- NULL
sucrose_gf$infection <- factor(sucrose_gf$infection, levels=c('mock','630'))
rm(sucrose)
galactitol <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-galactitol')))]
galactitol_strep <- subset(galactitol, abx == 'streptomycin')
galactitol_strep$abx <- NULL
galactitol_strep$infection <- factor(galactitol_strep$infection, levels=c('mock','630'))
galactitol_cef <- subset(galactitol, abx == 'cefoperazone')
galactitol_cef$abx <- NULL
galactitol_cef$infection <- factor(galactitol_cef$infection, levels=c('mock','630'))
galactitol_clinda <- subset(galactitol, abx == 'clindamycin')
galactitol_clinda$abx <- NULL
galactitol_clinda$infection <- factor(galactitol_clinda$infection, levels=c('mock','630'))
galactitol_gf <- subset(galactitol, abx == 'germfree')
galactitol_gf$abx <- NULL
galactitol_gf$infection <- factor(galactitol_gf$infection, levels=c('mock','630'))
rm(galactitol)
mannose <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-mannose')))]
mannose_strep <- subset(mannose, abx == 'streptomycin')
mannose_strep$abx <- NULL
mannose_strep$infection <- factor(mannose_strep$infection, levels=c('mock','630'))
mannose_cef <- subset(mannose, abx == 'cefoperazone')
mannose_cef$abx <- NULL
mannose_cef$infection <- factor(mannose_cef$infection, levels=c('mock','630'))
mannose_clinda <- subset(mannose, abx == 'clindamycin')
mannose_clinda$abx <- NULL
mannose_clinda$infection <- factor(mannose_clinda$infection, levels=c('mock','630'))
mannose_gf <- subset(mannose, abx == 'germfree')
mannose_gf$abx <- NULL
mannose_gf$infection <- factor(mannose_gf$infection, levels=c('mock','630'))
rm(mannose)
sialyllactose <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('3-sialyllactose','6\'-sialyllactose')))]
sialyllactose_strep <- subset(sialyllactose, abx == 'streptomycin')
sialyllactose_strep$abx <- NULL
sialyllactose_strep$infection <- factor(sialyllactose_strep$infection, levels=c('mock','630'))
sialyllactose_cef <- subset(sialyllactose, abx == 'cefoperazone')
sialyllactose_cef$abx <- NULL
sialyllactose_cef$infection <- factor(sialyllactose_cef$infection, levels=c('mock','630'))
sialyllactose_clinda <- subset(sialyllactose, abx == 'clindamycin')
sialyllactose_clinda$abx <- NULL
sialyllactose_clinda$infection <- factor(sialyllactose_clinda$infection, levels=c('mock','630'))
sialyllactose_gf <- subset(sialyllactose, abx == 'germfree')
sialyllactose_gf$abx <- NULL
sialyllactose_gf$infection <- factor(sialyllactose_gf$infection, levels=c('mock','630'))
rm(sialyllactose)


# Amino acids
glycine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-glycine')))]
glycine_strep <- subset(glycine, abx == 'streptomycin')
glycine_strep$abx <- NULL
glycine_strep$infection <- factor(glycine_strep$infection, levels=c('mock','630'))
glycine_cef <- subset(glycine, abx == 'cefoperazone')
glycine_cef$abx <- NULL
glycine_cef$infection <- factor(glycine_cef$infection, levels=c('mock','630'))
glycine_clinda <- subset(glycine, abx == 'clindamycin')
glycine_clinda$abx <- NULL
glycine_clinda$infection <- factor(glycine_clinda$infection, levels=c('mock','630'))
glycine_gf <- subset(glycine, abx == 'germfree')
glycine_gf$abx <- NULL
glycine_gf$infection <- factor(glycine_gf$infection, levels=c('mock','630'))
rm(glycine)
hydroxyproline <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-hydroxyproline')))]
hydroxyproline_strep <- subset(hydroxyproline, abx == 'streptomycin')
hydroxyproline_strep$abx <- NULL
hydroxyproline_strep$infection <- factor(hydroxyproline_strep$infection, levels=c('mock','630'))
hydroxyproline_cef <- subset(hydroxyproline, abx == 'cefoperazone')
hydroxyproline_cef$abx <- NULL
hydroxyproline_cef$infection <- factor(hydroxyproline_cef$infection, levels=c('mock','630'))
hydroxyproline_clinda <- subset(hydroxyproline, abx == 'clindamycin')
hydroxyproline_clinda$abx <- NULL
hydroxyproline_clinda$infection <- factor(hydroxyproline_clinda$infection, levels=c('mock','630'))
hydroxyproline_gf <- subset(hydroxyproline, abx == 'germfree')
hydroxyproline_gf$abx <- NULL
hydroxyproline_gf$infection <- factor(hydroxyproline_gf$infection, levels=c('mock','630'))
rm(hydroxyproline)
prolylglycine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-prolylglycine')))]
prolylglycine_strep <- subset(prolylglycine, abx == 'streptomycin')
prolylglycine_strep$abx <- NULL
prolylglycine_strep$infection <- factor(prolylglycine_strep$infection, levels=c('mock','630'))
prolylglycine_cef <- subset(prolylglycine, abx == 'cefoperazone')
prolylglycine_cef$abx <- NULL
prolylglycine_cef$infection <- factor(prolylglycine_cef$infection, levels=c('mock','630'))
prolylglycine_clinda <- subset(prolylglycine, abx == 'clindamycin')
prolylglycine_clinda$abx <- NULL
prolylglycine_clinda$infection <- factor(prolylglycine_clinda$infection, levels=c('mock','630'))
prolylglycine_gf <- subset(prolylglycine, abx == 'germfree')
prolylglycine_gf$abx <- NULL
prolylglycine_gf$infection <- factor(prolylglycine_gf$infection, levels=c('mock','630'))
rm(prolylglycine)

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
acetylglucosamine_sulfate <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-acetylglucosamine_6-sulfate')))]
acetylglucosamine_sulfate_strep <- subset(acetylglucosamine_sulfate, abx == 'streptomycin')
acetylglucosamine_sulfate_strep$abx <- NULL
acetylglucosamine_sulfate_strep$infection <- factor(acetylglucosamine_sulfate_strep$infection, levels=c('mock','630'))
acetylglucosamine_sulfate_cef <- subset(acetylglucosamine_sulfate, abx == 'cefoperazone')
acetylglucosamine_sulfate_cef$abx <- NULL
acetylglucosamine_sulfate_cef$infection <- factor(acetylglucosamine_sulfate_cef$infection, levels=c('mock','630'))
acetylglucosamine_sulfate_clinda <- subset(acetylglucosamine_sulfate, abx == 'clindamycin')
acetylglucosamine_sulfate_clinda$abx <- NULL
acetylglucosamine_sulfate_clinda$infection <- factor(acetylglucosamine_sulfate_clinda$infection, levels=c('mock','630'))
acetylglucosamine_sulfate_gf <- subset(acetylglucosamine_sulfate, abx == 'germfree')
acetylglucosamine_sulfate_gf$abx <- NULL
acetylglucosamine_sulfate_gf$infection <- factor(acetylglucosamine_sulfate_gf$infection, levels=c('mock','630'))
rm(acetylglucosamine_sulfate)
acetylglucosaminylasparagine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-acetylglucosaminylasparagine')))]
acetylglucosaminylasparagine_strep <- subset(acetylglucosaminylasparagine, abx == 'streptomycin')
acetylglucosaminylasparagine_strep$abx <- NULL
acetylglucosaminylasparagine_strep$infection <- factor(acetylglucosaminylasparagine_strep$infection, levels=c('mock','630'))
acetylglucosaminylasparagine_cef <- subset(acetylglucosaminylasparagine, abx == 'cefoperazone')
acetylglucosaminylasparagine_cef$abx <- NULL
acetylglucosaminylasparagine_cef$infection <- factor(acetylglucosaminylasparagine_cef$infection, levels=c('mock','630'))
acetylglucosaminylasparagine_clinda <- subset(acetylglucosaminylasparagine, abx == 'clindamycin')
acetylglucosaminylasparagine_clinda$abx <- NULL
acetylglucosaminylasparagine_clinda$infection <- factor(acetylglucosaminylasparagine_clinda$infection, levels=c('mock','630'))
acetylglucosaminylasparagine_gf <- subset(acetylglucosaminylasparagine, abx == 'germfree')
acetylglucosaminylasparagine_gf$abx <- NULL
acetylglucosaminylasparagine_gf$infection <- factor(acetylglucosaminylasparagine_gf$infection, levels=c('mock','630'))
rm(acetylglucosaminylasparagine)
acetylbetaglucosaminylamine <- metabolome[, c(1,2,which(colnames(metabolome) %in% c('N-acetyl-beta-glucosaminylamine')))]
acetylbetaglucosaminylamine_strep <- subset(acetylbetaglucosaminylamine, abx == 'streptomycin')
acetylbetaglucosaminylamine_strep$abx <- NULL
acetylbetaglucosaminylamine_strep$infection <- factor(acetylbetaglucosaminylamine_strep$infection, levels=c('mock','630'))
acetylbetaglucosaminylamine_cef <- subset(acetylbetaglucosaminylamine, abx == 'cefoperazone')
acetylbetaglucosaminylamine_cef$abx <- NULL
acetylbetaglucosaminylamine_cef$infection <- factor(acetylbetaglucosaminylamine_cef$infection, levels=c('mock','630'))
acetylbetaglucosaminylamine_clinda <- subset(acetylbetaglucosaminylamine, abx == 'clindamycin')
acetylbetaglucosaminylamine_clinda$abx <- NULL
acetylbetaglucosaminylamine_clinda$infection <- factor(acetylbetaglucosaminylamine_clinda$infection, levels=c('mock','630'))
acetylbetaglucosaminylamine_gf <- subset(acetylbetaglucosaminylamine, abx == 'germfree')
acetylbetaglucosaminylamine_gf$abx <- NULL
acetylbetaglucosaminylamine_gf$infection <- factor(acetylbetaglucosaminylamine_gf$infection, levels=c('mock','630'))
rm(acetylbetaglucosaminylamine)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Calculate significant differences



#--------------------#

# Prep untreated groups for plotting
fructose_mock <- subset(fructose_mock, abx == 'none')
sucrose_mock <- subset(sucrose_mock, abx == 'none')
galactitol_mock <- subset(galactitol_mock, abx == 'none')
mannose_mock <- subset(mannose_mock, abx == 'none')
glucose_mock <- subset(glucose_mock, abx == 'none')
glycine_mock <- subset(glycine_mock, abx == 'none')
hydroxyproline_mock <- subset(hydroxyproline_mock, abx == 'none')
prolylglycine_mock <- subset(prolylglycine_mock, abx == 'none')
acetylmuramate_mock <- subset(acetylmuramate_mock, abx == 'none')
acetylglucosamine_sulfate_mock <- subset(acetylglucosamine_sulfate_mock, abx == 'none')
acetylglucosaminylasparagine_mock <- subset(acetylglucosaminylasparagine_mock, abx == 'none')
acetylbetaglucosaminylamine_mock <- subset(acetylbetaglucosaminylamine_mock, abx == 'none')

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S6.pdf'
pdf(file=plot_file, width=15, height=8)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE))

#-------------------------#

# Shared metabolite importances
par(mar=c(4,5,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(shared_importance$Metabolite_score, labels=shared_importance$Compound_name, lcolor=NA, cex=1.2, color='black',
         xlab='Median Importance Score', xlim=c(0,10), pch=19, lwd=3)
segments(x0=rep(0,10), y0=c(1:27), x1=rep(12,10), y1=c(1:27), lty=2) # Dotted lines
mtext('A', side=2, line=2, las=2, adj=4, padj=-18, cex=1.5)

#-------------------------#

# Distinct metabolite importances
par(mar=c(4,4,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(top_importances$Metabolite_score, labels=top_importances$Compound_name,
         lcolor=NA, cex=1.2, groups=top_importances$abx, color='black',
         xlab='Importance Score', xlim=c(0,10), pch=19, lwd=3,
         gcolor=c(wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5]))
mtext('B', side=2, line=2, las=2, adj=2.3, padj=-18, cex=1.5)
segments(x0=rep(0, 15), y0=c(1:7, 10:23, 26:37), 
         x1=rep(12, 15), y1=c(1:7, 10:23, 26:37), lty=2) # Dotted lines

#------------------------------------------------------#




#-------------------------#






#-------------------------#




#-------------------------#



#-------------------------#





#-------------------------#





#-------------------------#



#-------------------------#





#-------------------------#



dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()
