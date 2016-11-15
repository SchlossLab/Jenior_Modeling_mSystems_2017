
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
cef_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/cefoperazone_630.importance_score.tsv'
clinda_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/clindamycin_630.bipartite.files/clindamycin_630.importance_score.tsv'
strep_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/streptomycin_630.bipartite.files/streptomycin_630.importance_score.tsv'
gf_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/germfree_630.bipartite.files/germfree_630.importance_score.tsv'

cef_importance <- read.delim(cef_importance_file, header=TRUE, sep='\t', row.names=1)
clinda_importance <- read.delim(clinda_importance_file, header=TRUE, sep='\t', row.names=1)
strep_importance <- read.delim(strep_importance_file, header=TRUE, sep='\t', row.names=1)
gf_importance <- read.delim(gf_importance_file, header=TRUE, sep='\t', row.names=1)
rm(cef_importance_file, clinda_importance_file, strep_importance_file, gf_importance_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format metabolite importance scores
cef_importance$Metabolite_score <- as.numeric(as.character(cef_importance$Metabolite_score))
cef_importance <- cef_importance[order(-cef_importance$Metabolite_score),]
clinda_importance$Metabolite_score <- as.numeric(as.character(clinda_importance$Metabolite_score))
clinda_importance <- clinda_importance[order(-clinda_importance$Metabolite_score),]
strep_importance$Metabolite_score <- as.numeric(as.character(strep_importance$Metabolite_score))
strep_importance <- strep_importance[order(-strep_importance$Metabolite_score),]
gf_importance$Metabolite_score <- as.numeric(as.character(gf_importance$Metabolite_score))
gf_importance <- gf_importance[order(-gf_importance$Metabolite_score),]

# Remove non-significant scores
cef_importance <- as.data.frame(subset(cef_importance, (cef_importance[,3] != 'n.s.')))
clinda_importance <- as.data.frame(subset(clinda_importance, (clinda_importance[,3] != 'n.s.')))
strep_importance <- as.data.frame(subset(strep_importance, (strep_importance[,3] != 'n.s.')))
gf_importance <- as.data.frame(subset(gf_importance, (gf_importance[,3] != 'n.s.')))

# Write significant ranked data to supplementary table
write.table(cef_importance, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S3cef.tsv', quote=FALSE, sep='\t', row.names=TRUE)
write.table(clinda_importance, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S3clinda.tsv', quote=FALSE, sep='\t', row.names=TRUE)
write.table(strep_importance, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S3strep.tsv', quote=FALSE, sep='\t', row.names=TRUE)
write.table(gf_importance, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S3gf.tsv', quote=FALSE, sep='\t', row.names=TRUE)
# Assemble into multi-paneled Excel table downstream

# Take top 50 scores
cef_importance <- cef_importance[c(1:50),]
clinda_importance <- clinda_importance[c(1:50),]
strep_importance <- strep_importance[c(1:50),]
gf_importance <- gf_importance[c(1:50),]

# find shared important metabolites
shared_importance <- as.data.frame(subset(cef_importance, (cef_importance[,1] %in% clinda_importance[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% clinda_importance[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% strep_importance[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% gf_importance[,1])))
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
shared_gf <- as.data.frame(subset(gf_importance, (gf_importance[,1] %in% shared_importance)))
shared_gf <- shared_gf[order(shared_gf$Compound_name),]
shared_gf$Metabolite_score <- as.numeric(as.character(shared_gf$Metabolite_score))

score_median <- as.data.frame(apply(cbind(shared_cef$Metabolite_score, shared_clinda$Metabolite_score, shared_strep$Metabolite_score, shared_gf$Metabolite_score), 1, median))
shared_importance <- cbind(shared_cef$Compound_name, score_median)
rownames(shared_importance) <- rownames(shared_cef)
colnames(shared_importance) <- c('Compound_name','Metabolite_score')
shared_importance <- shared_importance[order(shared_importance$Metabolite_score),]
rm(shared_cef, shared_clinda, shared_strep, shared_gf, score_median)

# Subset to most important, distinct metabolites
cef_importance <- cef_importance[c(1:25),]
clinda_importance <- clinda_importance[c(1:25),]
strep_importance <- strep_importance[c(1:25),]
gf_importance <- gf_importance[c(1:25),]
cef_only_importance <- as.data.frame(subset(cef_importance, !(cef_importance[,1] %in% clinda_importance[,1])))
cef_only_importance <- as.data.frame(subset(cef_only_importance, !(cef_only_importance[,1] %in% strep_importance[,1])))
cef_only_importance <- as.data.frame(subset(cef_only_importance, !(cef_only_importance[,1] %in% gf_importance[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_importance, !(clinda_importance[,1] %in% cef_importance[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_only_importance, !(clinda_only_importance[,1] %in% strep_importance[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_only_importance, !(clinda_only_importance[,1] %in% gf_importance[,1])))
strep_only_importance <- as.data.frame(subset(strep_importance, !(strep_importance[,1] %in% clinda_importance[,1])))
strep_only_importance <- as.data.frame(subset(strep_only_importance, !(strep_only_importance[,1] %in% cef_importance[,1])))
strep_only_importance <- as.data.frame(subset(strep_only_importance, !(strep_only_importance[,1] %in% gf_importance[,1])))
gf_only_importance <- as.data.frame(subset(gf_importance, !(gf_importance[,1] %in% clinda_importance[,1])))
gf_only_importance <- as.data.frame(subset(gf_only_importance, !(gf_only_importance[,1] %in% strep_importance[,1])))
gf_only_importance <- as.data.frame(subset(gf_only_importance, !(gf_only_importance[,1] %in% cef_importance[,1])))
rm(cef_importance, clinda_importance, strep_importance, gf_importance)
cef_only_importance <- cef_only_importance[order(cef_only_importance$Metabolite_score),]
clinda_only_importance <- clinda_only_importance[order(clinda_only_importance$Metabolite_score),]
strep_only_importance <- strep_only_importance[order(strep_only_importance$Metabolite_score),]
gf_only_importance <- gf_only_importance[order(gf_only_importance$Metabolite_score),]
cef_only_importance$abx <- 'Cefoperazone (SPF)'
clinda_only_importance$abx <- 'Clindamycin (SPF)'
strep_only_importance$abx <- 'Streptomycin (SPF)'
gf_only_importance$abx <- 'No Antibiotics (GF)'
top_importances <- rbind(cef_only_importance, clinda_only_importance, strep_only_importance, gf_only_importance)
top_importances$abx <- as.factor(top_importances$abx)
top_importances$abx <- ordered(top_importances$abx, levels=c('Streptomycin (SPF)', 'Cefoperazone (SPF)', 'Clindamycin (SPF)', 'No Antibiotics (GF)'))
rm(cef_only_importance, clinda_only_importance, strep_only_importance, gf_only_importance)

# Format names to look better for the plot
top_importances$Compound_name <- gsub('_',' ',top_importances$Compound_name)
top_importances$Compound_name <- gsub('monophosphate','p',top_importances$Compound_name)
top_importances$Compound_name <- gsub('phosphate','p',top_importances$Compound_name)
top_importances$Compound_name <- gsub('alpha-','',top_importances$Compound_name)
top_importances$Compound_name <- gsub('beta-','',top_importances$Compound_name)
top_importances$Compound_name <- gsub('\\(R\\)\\-','',top_importances$Compound_name)
top_importances$Compound_name <- gsub('\\(S\\)\\-','',top_importances$Compound_name)
top_importances$Compound_name[top_importances$Compound_name == '1-(5\'-Phosphoribosyl)-5-amino-4-(N-succinocarboxamide)-imidazole'] <- 'SAICAR' # shorten a long name
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

# Read in growth rate data
# Define variables
growth_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/cd630_growth.tsv'

# Read in data
growth <- read.delim(growth_file, sep='\t', header=TRUE, row.names=1)
growth <- as.data.frame(t(growth))
rm(growth_file)

# Seperate to groups of each growth substrate and format
trehalose <- cbind(growth$trehalose_1, growth$trehalose_2, growth$trehalose_3) - growth$trehalose_blank
trehalose[trehalose < 0] <- 0
starch <- cbind(growth$starch_1, growth$starch_2, growth$starch_2) - growth$starch_blank
starch[starch < 0] <- 0
fructose <- cbind(growth$fructose_1, growth$fructose_2, growth$fructose_3) - growth$fructose_blank
fructose[fructose < 0] <- 0
mannitol <- cbind(growth$mannitol_1, growth$mannitol_2, growth$mannitol_3) - growth$mannitol_blank
mannitol[mannitol < 0] <- 0
salicin <- cbind(growth$salicin_1, growth$salicin_2, growth$salicin_3) - growth$salicin_blank
salicin[salicin < 0] <- 0
acetate <- cbind(growth$acetate_1, growth$acetate_2, growth$acetate_3) - growth$acetate_blank
acetate[acetate < 0] <- 0
acetylglucosamine <- cbind(growth$acetylglucosamine_1, growth$acetylglucosamine_2, growth$acetylglucosamine_3) - growth$acetylglucosamine_blank
acetylglucosamine[acetylglucosamine < 0] <- 0
acetylneuraminate <- cbind(growth$acetylneuraminate_1, growth$acetylneuraminate_2, growth$acetylneuraminate_3) - growth$acetylneuraminate_blank
acetylneuraminate[acetylneuraminate < 0] <- 0
no_carb <- cbind(growth$noCarb_1, growth$noCarb_2, growth$noCarb_3) - growth$noCarb_blank
no_carb[no_carb < 0] <- 0
no_aa <- cbind(growth$noAA_1, growth$noAA_2, growth$noAA_3) - growth$noAA_blank
no_aa[no_aa < 0] <- 0
bhi <- cbind(growth$bhi_1, growth$bhi_2, growth$bhi_3) - growth$bhi_blank
bhi[bhi < 0] <- 0
rm(growth)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format growth curves

# Find medians 
starch_median <- apply(starch, 1, median)
fructose_median <-  apply(fructose, 1, median)
trehalose_median <-  apply(trehalose, 1, median)
mannitol_median <- apply(mannitol, 1, median)
salicin_median <- apply(salicin, 1, median)
acetate_median <- apply(acetate, 1, median)
acetylglucosamine_median <- apply(acetylglucosamine, 1, median)
acetylneuraminate_median <- apply(acetylneuraminate, 1, median)
no_carb_median <- apply(no_carb, 1, median)
bhi_median <- apply(bhi, 1, median)
no_aa_median <- apply(no_aa, 1, median)

# Standard deviations
starch_sd <- rowSds(starch)
acetylneuraminate_sd <- rowSds(acetylneuraminate)
fructose_sd <-  rowSds(fructose)
trehalose_sd <-  rowSds(trehalose)
mannitol_sd <- rowSds(mannitol)
salicin_sd <- rowSds(salicin)
acetate_sd <- rowSds(acetate)
acetylglucosamine_sd <- rowSds(acetylglucosamine)
no_carb_sd <- rowSds(no_carb)
bhi_sd <- rowSds(bhi)
no_aa_sd <- rowSds(no_aa)

# Compile results
growth_medians <- as.data.frame(cbind(acetate_median, acetylglucosamine_median, trehalose_median,
                                      acetylneuraminate_median, starch_median, fructose_median, mannitol_median, salicin_median, no_carb_median, no_aa_median))
growth_sds <- as.data.frame(cbind(acetate_sd, acetylglucosamine_sd, acetylneuraminate_sd, trehalose_sd,
                                  starch_sd, fructose_sd, mannitol_sd, salicin_sd, no_carb_sd, no_aa_sd))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Analyze growth curves
substrates <- c('acetate','acetylglucosamine','acetylneuraminate','trehalose','fructose', 'mannitol','salicin','bhi','no_amino_acids','no_carbohydrates')

# Maximum growth rate
max_rate <- round(c(diff(acetate_median)[which.max(diff(acetate_median))],
                    diff(acetylglucosamine_median)[which.max(diff(acetylglucosamine_median))], 
                    diff(acetylneuraminate_median)[which.max(diff(acetylneuraminate_median))], diff(trehalose_median)[which.max(diff(trehalose_median))],
                    diff(fructose_median)[which.max(diff(fructose_median))], 
                    diff(mannitol_median)[which.max(diff(mannitol_median))], diff(salicin_median)[which.max(diff(salicin_median))], 
                    diff(bhi_median)[which.max(diff(bhi_median))], diff(no_aa_median)[which.max(diff(no_aa_median))], diff(no_carb_median)[which.max(diff(no_carb_median))]), digits=3)

# Time of maximum growth rate
time_max_rate <- round(c((which.max(diff(acetate_median)) * 0.5), (which.max(diff(acetylglucosamine_median)) * 0.5),
                         (which.max(diff(acetylneuraminate_median)) * 0.5), (which.max(diff(trehalose_median)) * 0.5),
                         (which.max(diff(fructose_median)) * 0.5), (which.max(diff(mannitol_median)) * 0.5), (which.max(diff(salicin_median)) * 0.5), 
                         (which.max(diff(bhi_median)) * 0.5), (which.max(diff(no_aa_median)) * 0.5), (which.max(diff(no_carb_median)) * 0.5)), digits=3) - 0.5
# Maximum OD
max_od <- round(c(max(acetate_median), max(acetylglucosamine_median), max(acetylneuraminate_median), max(trehalose_median), 
                  max(fructose_median), max(mannitol_median), max(salicin_median), max(bhi_median), max(no_aa_median), max(no_carb_median)), digits=3)

# Time of max OD
time_max_od <- round(c((which.max(acetate_median) * 0.5), (which.max(acetylglucosamine_median) * 0.5),
                       (which.max(acetylneuraminate_median) * 0.5), (which.max(trehalose_median) * 0.5), 
                       (which.max(fructose_median) * 0.5), (which.max(mannitol_median) * 0.5), (which.max(salicin_median) * 0.5), 
                       (which.max(bhi_median) * 0.5), (which.max(no_aa_median) * 0.5), (which.max(no_carb_median) * 0.5)), digits=3) - 0.5

# Growth rate at 24 hours
rate_24_hrs <- round(c(diff(acetate_median)[length(diff(acetate_median))],
                       diff(acetylglucosamine_median)[length(diff(acetylglucosamine_median))], diff(acetylneuraminate_median)[length(diff(acetylneuraminate_median))], 
                       diff(trehalose_median)[length(diff(trehalose_median))], diff(fructose_median)[length(diff(fructose_median))], 
                       diff(mannitol_median)[length(diff(mannitol_median))], diff(salicin_median)[length(diff(salicin_median))], 
                       diff(bhi_median)[length(diff(bhi_median))], diff(no_aa_median)[length(diff(no_aa_median))], diff(no_carb_median)[length(diff(no_carb_median))]), digits=3)

# Mean growth rate
mean_rate <- round(c(mean(diff(acetate_median)), mean(diff(acetylglucosamine_median)), mean(diff(acetylneuraminate_median)), mean(diff(trehalose_median)),
                     mean(diff(fructose_median)), mean(diff(mannitol_median)), mean(diff(salicin_median)), mean(diff(bhi_median)), mean(diff(no_aa_median)), mean(diff(no_carb_median))), digits=3)

# Area under curve
area_under <- round(c(auc(acetate_median, seq(1,49,1)), auc(acetylglucosamine_median, seq(1,49,1)),
                      auc(acetylneuraminate_median, seq(1,49,1)), auc(trehalose_median, seq(1,49,1)), auc(fructose_median, seq(1,49,1)), 
                      auc(mannitol_median, seq(1,49,1)), auc(salicin_median, seq(1,49,1)), auc(bhi_median, seq(1,49,1)), auc(no_aa_median, seq(1,49,1)), auc(no_carb_median, seq(1,49,1))), digits=3)

# Assemble the table
growth_summary <- cbind(substrates, max_rate, time_max_rate, max_od, time_max_od, rate_24_hrs, mean_rate, area_under)
colnames(growth_summary) <- c('Substrate', 'Max_Growth_Rate', 'Time_of_Max_Rate_in_Hours', 'Max_OD', 'Time_of_Max_OD_in_Hours', 'Rate_at_24_hours', 'Mean_Rate', 'AUC')

rm(substrates, max_rate, time_max_rate, max_od, time_max_od, rate_24_hrs, mean_rate, area_under)
rm(acetate_median, acetylglucosamine_median, acetylneuraminate_median, trehalose_median, starch_median, fructose_median, mannitol_median, salicin_median, no_carb_median, no_aa_median)
rm(acetate, acetylglucosamine, acetylneuraminate, trehalose, starch, fructose, mannitol, salicin, no_carb, no_aa)
rm(acetate_sd, acetylglucosamine_sd, acetylneuraminate_sd, trehalose_sd, starch_sd, fructose_sd, mannitol_sd, salicin_sd, no_carb_sd, no_aa_sd)

# Write growth summary data to supplementary table
table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/Table_S4.tsv'
write.table(growth_summary, file=table_file, quote=FALSE, sep='\t', row.names=FALSE)
rm(table_file, growth_summary)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_6.pdf'
pdf(file=plot_file, width=14, height=19)
layout(matrix(c(1,1,1,1,2,2,2,2,
                1,1,1,1,2,2,2,2,
                1,1,1,1,2,2,2,2,
                1,1,1,1,2,2,2,2,
                1,1,1,1,2,2,2,2,
                3,3,3,3,3,3,3,3,
                3,3,3,3,3,3,3,3,
                3,3,3,3,3,3,3,3,
                3,3,3,3,3,3,3,3), nrow=9, ncol=8, byrow=TRUE))

#---------------------------------------#

# Shared metabolite importances
par(mar=c(4,4,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(shared_importance$Metabolite_score, labels=shared_importance$Compound_name, lcolor=NA, cex=1.4, color='black',
         xlab=expression(paste('Importance Score (',Log[2],')')), xlim=c(0,10), pch=19, lwd=3)
segments(x0=rep(0,10), y0=c(1:10), x1=rep(12,10), y1=c(1:10), lty=2) # Dotted lines
mtext('a', side=2, line=2, las=2, adj=2.9, padj=-20, cex=1.8, font=2)

#---------------------------------------#

# Distinct metabolite importances
par(mar=c(4,4,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(top_importances$Metabolite_score, labels=top_importances$Compound_name,
         lcolor=NA, cex=1.4, groups=top_importances$abx, color='black',
         xlab=expression(paste('Importance Score (',Log[2],')')), xlim=c(0,10), pch=19, lwd=3,
         gcolor=c(wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5],'forestgreen'))
mtext('b', side=2, line=2, las=2, adj=2.4, padj=-20, cex=1.8, font=2)
segments(x0=rep(0, 15), y0=c(1:17, 20:23, 26:30, 33:35), x1=rep(12, 15), y1=c(1:17, 20:23, 26:30, 33:35), lty=2) # Dotted lines

#---------------------------------------#

# Growth on important carbohydrates
par(mar=c(7,7,1.5,2), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE, mgp=c(4,2,0))
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,50), ylim=c(-0.03,1.0), lwd=3, pch=15, xlab='Time (hours)', ylab=expression(OD[600]), cex=2.3)
abline(h=seq(0,1.0,0.1), lty=3, col='gray68') # adding gridlines
abline(v=seq(1,50,2), lty=3, col='gray68') # adding gridlines
axis(1, at=seq(1,49,4), labels=seq(0,24,2), tck=-0.018)
axis(2, at=seq(0.0,1.0,0.2), labels=c('0.0','0.2','0.4','0.6','0.8','1.0'), tck=-0.018)
mtext('c', side=2, line=2, las=2, adj=3.5, padj=-17, cex=1.8, font=2)

# Control
lines(growth_medians$no_carb_median, type='o', lwd=3, pch=16, cex=0, col='gray45')
segments(x0=seq(1,49,1), y0=growth_medians$no_carb_median+growth_sds$no_carb_sd, x1=seq(1,49,1), y1=growth_medians$no_carb_median-growth_sds$no_carb_sd, lwd=2.5, cex=2, col='gray45')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$no_carb_median+growth_sds$no_carb_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$no_carb_median+growth_sds$no_carb_sd, lwd=2.5, col='gray45')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$no_carb_median-growth_sds$no_carb_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$no_carb_median-growth_sds$no_carb_sd, lwd=2.5, col='gray45')

# Shared
lines(growth_medians$acetylglucosamine_median, type='o', col='black', lwd=2.5, pch=6, cex=2.3)
segments(x0=seq(1,49,1), y0=growth_medians$acetylglucosamine_median+growth_sds$acetylglucosamine_sd, x1=seq(1,49,1), y1=growth_medians$acetylglucosamine_median-growth_sds$acetylglucosamine_sd, lwd=2.5, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetylglucosamine_median+growth_sds$acetylglucosamine_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetylglucosamine_median+growth_sds$acetylglucosamine_sd, lwd=2.5, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetylglucosamine_median-growth_sds$acetylglucosamine_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetylglucosamine_median-growth_sds$acetylglucosamine_sd, lwd=2.5, col='black')
lines(growth_medians$trehalose_median, type='o', col='black', bg='black', lwd=2.5, pch=25, cex=2.3)
segments(x0=seq(1,49,1), y0=growth_medians$trehalose_median+growth_sds$trehalose_sd, x1=seq(1,49,1), y1=growth_medians$trehalose_median-growth_sds$trehalose_sd, lwd=2.5, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$trehalose_median+growth_sds$trehalose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$trehalose_median+growth_sds$trehalose_sd, lwd=2.5, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$trehalose_median-growth_sds$trehalose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$trehalose_median-growth_sds$trehalose_sd, lwd=2.5, col='black')

# Streptomycin
lines(growth_medians$fructose_median, type='o', col=wes_palette('FantasticFox')[1], lwd=2.5, pch=0, cex=2.6)
segments(x0=seq(1,49,1), y0=growth_medians$fructose_median+growth_sds$fructose_sd, x1=seq(1,49,1), y1=growth_medians$fructose_median-growth_sds$fructose_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$fructose_median+growth_sds$fructose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$fructose_median+growth_sds$fructose_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$fructose_median-growth_sds$fructose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$fructose_median-growth_sds$fructose_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])

# Cefoperazone
lines(growth_medians$mannitol_median, type='o', col=wes_palette('FantasticFox')[3], lwd=2.5, pch=1, cex=2.7)
segments(x0=seq(1,49,1), y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, x1=seq(1,49,1), y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$mannitol_median+growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$mannitol_median-growth_sds$mannitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])

# Clindamycin
lines(growth_medians$salicin_median, type='o', col=wes_palette('FantasticFox')[5], lwd=2.5, pch=2, cex=2.5)
segments(x0=seq(1,49,1), y0=growth_medians$salicin_median+growth_sds$salicin_sd, x1=seq(1,49,1), y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$salicin_median+growth_sds$salicin_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$salicin_median+growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$salicin_median-growth_sds$salicin_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])

# Gnotobiotic
lines(growth_medians$acetylneuraminate_median, type='o', col='forestgreen', lwd=2.5, pch=5, cex=2.5)
segments(x0=seq(1,49,1), y0=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd, x1=seq(1,49,1), y1=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd, lwd=2.5, col='forestgreen')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd, lwd=2.5, col='forestgreen')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd, lwd=2.5, col='forestgreen')

legend('topleft', legend=c('No Carbohydrates','N-Acetyl-D-glucosamine','Trehalose','D-Fructose','Mannitol','Salicin','N-Acetylneuriminate'), 
       col=c('gray45','black','black',wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5],'forestgreen'), 
       pch=c(16,6,25,0,1,2,5), cex=2.4, pt.cex=c(0,3.4,3.2,3.4,3.4,3.4,3.4), pt.bg=c('white','white','black','white','white','white','white'), lwd=3, bg='white')

dev.off()

#---------------------------------------#

plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/Figure_S3.pdf'
pdf(file=plot_file, width=14, height=10)
par(mar=c(7,7,1.5,2), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE, mgp=c(4,2,0))
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,50), ylim=c(-0.03,0.83), lwd=2, pch=15, xlab='Time (hours)', ylab=expression(OD[600]), cex=2.3)
abline(h=seq(0,0.9,0.1), lty=3, col='gray68') # adding gridlines
abline(v=seq(1,50,2), lty=3, col='gray68') # adding gridlines
axis(1, at=seq(1,49,4), labels=seq(0,24,2), tck=-0.018)
axis(2, at=seq(0.0,0.8,0.2), labels=c('0.0','0.2','0.4','0.6','0.8'), tck=-0.018)

# Control
lines(growth_medians$no_carb_median, type='o', lwd=3, pch=16, cex=0, col='gray45')
segments(x0=seq(1,49,1), y0=growth_medians$no_carb_median+growth_sds$no_carb_sd, x1=seq(1,49,1), y1=growth_medians$no_carb_median-growth_sds$no_carb_sd, lwd=2.5, cex=2, col='gray45')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$no_carb_median+growth_sds$no_carb_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$no_carb_median+growth_sds$no_carb_sd, lwd=2.5, col='gray45')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$no_carb_median-growth_sds$no_carb_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$no_carb_median-growth_sds$no_carb_sd, lwd=2.5, col='gray45')

lines(growth_medians$starch_median, type='o', col=wes_palette('FantasticFox')[1], lwd=2.5, pch=0, cex=2)
segments(x0=seq(1,49,1), y0=growth_medians$starch_median+growth_sds$starch_sd, x1=seq(1,49,1), y1=growth_medians$starch_median-growth_sds$starch_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$starch_median+growth_sds$starch_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$starch_median+growth_sds$starch_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$starch_median-growth_sds$starch_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$starch_median-growth_sds$starch_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])

lines(growth_medians$no_aa_median, type='o', lwd=2.5, pch=15, cex=1, col='gray45')
segments(x0=seq(1,49,1), y0=growth_medians$no_aa_median+growth_sds$no_aa_sd, x1=seq(1,49,1), y1=growth_medians$no_aa_median-growth_sds$no_aa_sd, lwd=2.5, cex=2, col='gray45')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$no_aa_median+growth_sds$no_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$no_aa_median+growth_sds$no_aa_sd, lwd=2.5, col='gray45')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$no_aa_median-growth_sds$no_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$no_aa_median-growth_sds$no_aa_sd, lwd=2.5, col='gray45')

lines(bhi_median, type='o', col='darkorchid3', lwd=2.5, pch=15, cex=0.8)
segments(x0=seq(1,49,1), y0=bhi_median+bhi_sd, x1=seq(1,49,1), y1=bhi_median-bhi_sd, lwd=2.5, col='darkorchid3')
segments(x0=seq(1,49,1)-0.2, y0=bhi_median+bhi_sd, x1=seq(1,49,1)+0.2, y1=bhi_median+bhi_sd, lwd=2.5, col='darkorchid3')
segments(x0=seq(1,49,1)-0.2, y0=bhi_median-bhi_sd, x1=seq(1,49,1)+0.2, y1=bhi_median-bhi_sd, lwd=2.5, col='darkorchid3')

lines(growth_medians$acetate_median, type='o', col='forestgreen', lwd=2.5, pch=18, cex=1.5)
segments(x0=seq(1,49,1), y0=growth_medians$acetate_median+growth_sds$acetate_sd, x1=seq(1,49,1), y1=growth_medians$acetate_median-growth_sds$acetate_sd, lwd=2.5, col='forestgreen')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetate_median+growth_sds$acetate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetate_median+growth_sds$acetate_sd, lwd=2.5, col='forestgreen')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetate_median-growth_sds$acetate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetate_median-growth_sds$acetate_sd, lwd=2.5, col='forestgreen')

legend('topleft', legend=c('No Carbohydrates','No Amino acids','BHI','Starch','Acetate'), 
       col=c('gray45','gray45','darkorchid3',wes_palette('FantasticFox')[1],'forestgreen'), 
       pch=c(16,16,15,15,18), cex=1.3, pt.cex=c(0,1,1,2.5,2.5), bg='white', lwd=2.5)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
rm(plot_file, top_importances, shared_importance, growth_medians, growth_sds, format_curve, bhi, bhi_median, bhi_sd)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()



