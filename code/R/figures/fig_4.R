
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
gf_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/germfree_630.bipartite.files/importances.tsv'

cef_importance <- read.delim(cef_importance_file, header=TRUE, sep='\t', row.names=1)
clinda_importance <- read.delim(clinda_importance_file, header=TRUE, sep='\t', row.names=1)
strep_importance <- read.delim(strep_importance_file, header=TRUE, sep='\t', row.names=1)
gf_importance <- read.delim(gf_importance_file, header=TRUE, sep='\t', row.names=1)
rm(cef_importance_file, clinda_importance_file, strep_importance_file, gf_importance_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format metabolite importance scores
cef_importance$Metabolite_score <- as.numeric(as.character(cef_importance$Metabolite_score))
clinda_importance$Metabolite_score <- as.numeric(as.character(clinda_importance$Metabolite_score))
strep_importance$Metabolite_score <- as.numeric(as.character(strep_importance$Metabolite_score))
gf_importance$Metabolite_score <- as.numeric(as.character(gf_importance$Metabolite_score))

# Remove non-significant scores
cef_importance <- as.data.frame(subset(cef_importance, (cef_importance[,3] != 'n.s.')))
cef_importance$p_value <- NULL
clinda_importance <- as.data.frame(subset(clinda_importance, (clinda_importance[,3] != 'n.s.')))
clinda_importance$p_value <- NULL
strep_importance <- as.data.frame(subset(strep_importance, (strep_importance[,3] != 'n.s.')))
strep_importance$p_value <- NULL
gf_importance <- as.data.frame(subset(gf_importance, (gf_importance[,3] != 'n.s.')))
gf_importance$p_value <- NULL

# Sort for most important
cef_importance <- cef_importance[order(-cef_importance$Metabolite_score),]
clinda_importance <- clinda_importance[order(-clinda_importance$Metabolite_score),]
strep_importance <- strep_importance[order(-strep_importance$Metabolite_score),]
gf_importance <- gf_importance[order(-gf_importance$Metabolite_score),]

# Write significant ranked data to supplementary table
write.table(cef_importance, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S3cef.tsv', quote=FALSE, sep='\t', row.names=TRUE)
write.table(clinda_importance, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S3clinda.tsv', quote=FALSE, sep='\t', row.names=TRUE)
write.table(strep_importance, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S3strep.tsv', quote=FALSE, sep='\t', row.names=TRUE)
write.table(gf_importance, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S3gf.tsv', quote=FALSE, sep='\t', row.names=TRUE)
# Assemble into multi-paneled Excel table downstream

# Take top 50 scores
cef_importance_top <- cef_importance[c(1:40),]
clinda_importance_top <- clinda_importance[c(1:40),]
strep_importance_top <- strep_importance[c(1:40),]
gf_importance_top <- gf_importance[c(1:40),]

# find shared important metabolites
shared_importance <- as.data.frame(subset(cef_importance_top, (cef_importance_top[,1] %in% clinda_importance_top[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% clinda_importance_top[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% strep_importance_top[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% gf_importance_top[,1])))
shared_importance <- shared_importance$Compound_name
shared_cef <- as.data.frame(subset(cef_importance_top, (cef_importance_top[,1] %in% shared_importance)))
shared_cef <- shared_cef[order(shared_cef$Compound_name),]
shared_cef$Metabolite_score <- as.numeric(as.character(shared_cef$Metabolite_score))
shared_clinda <- as.data.frame(subset(clinda_importance_top, (clinda_importance_top[,1] %in% shared_importance)))
shared_clinda <- shared_clinda[order(shared_clinda$Compound_name),]
shared_clinda$Metabolite_score <- as.numeric(as.character(shared_clinda$Metabolite_score))
shared_strep <- as.data.frame(subset(strep_importance_top, (strep_importance_top[,1] %in% shared_importance)))
shared_strep <- shared_strep[order(shared_strep$Compound_name),]
shared_strep$Metabolite_score <- as.numeric(as.character(shared_strep$Metabolite_score))
shared_gf <- as.data.frame(subset(gf_importance_top, (gf_importance_top[,1] %in% shared_importance)))
shared_gf <- shared_gf[order(shared_gf$Compound_name),]
shared_gf$Metabolite_score <- as.numeric(as.character(shared_gf$Metabolite_score))
antilog_scores <- 2 ^ cbind(shared_cef$Metabolite_score, shared_clinda$Metabolite_score, shared_strep$Metabolite_score, shared_gf$Metabolite_score)
score_sd <- log2(as.data.frame(apply(antilog_scores, 1, sd)))
score_median <- log2(as.data.frame(apply(antilog_scores, 1, median)))
shared_importance <- cbind(shared_cef$Compound_name, score_median)
rownames(shared_importance) <- rownames(shared_cef)
colnames(shared_importance) <- c('Compound_name','Metabolite_score')
shared_importance$abx <- 'Median Shared Importance'
shared_importance <- shared_importance[order(-shared_importance$Metabolite_score),]
#if (nrow(shared_importance) > 10) {shared_importance <- shared_importance[1:10,]}

shared_importance <- shared_importance[order(shared_importance$Metabolite_score),]
shared_importance$Compound_name <- gsub('-D-','', shared_importance$Compound_name)
rm(shared_cef, shared_clinda, shared_strep, shared_gf, score_median)

# Subset to most important, distinct metabolites
cef_only_importance <- as.data.frame(subset(cef_importance_top, !(cef_importance_top[,1] %in% clinda_importance_top[,1])))
cef_only_importance <- as.data.frame(subset(cef_only_importance, !(cef_only_importance[,1] %in% strep_importance_top[,1])))
cef_only_importance <- as.data.frame(subset(cef_only_importance, !(cef_only_importance[,1] %in% gf_importance_top[,1])))
cef_only_importance <- as.data.frame(subset(cef_only_importance, !(cef_only_importance$Compound_name %in% c('[Enzyme]-cysteine')))) 
#if (nrow(cef_only_importance) > 10) {cef_only_importance <- cef_only_importance[1:10,]}

clinda_only_importance <- as.data.frame(subset(clinda_importance_top, !(clinda_importance_top[,1] %in% cef_importance_top[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_only_importance, !(clinda_only_importance[,1] %in% strep_importance_top[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_only_importance, !(clinda_only_importance[,1] %in% gf_importance_top[,1])))
#if (nrow(clinda_only_importance) > 10) {clinda_only_importance <- clinda_only_importance[1:10,]}

strep_only_importance <- as.data.frame(subset(strep_importance_top, !(strep_importance_top[,1] %in% clinda_importance_top[,1])))
strep_only_importance <- as.data.frame(subset(strep_only_importance, !(strep_only_importance[,1] %in% cef_importance_top[,1])))
strep_only_importance <- as.data.frame(subset(strep_only_importance, !(strep_only_importance[,1] %in% gf_importance_top[,1])))
#if (nrow(strep_only_importance) > 10) {strep_only_importance <- strep_only_importance[1:10,]}

gf_only_importance <- as.data.frame(subset(gf_importance_top, !(gf_importance_top[,1] %in% clinda_importance_top[,1])))
gf_only_importance <- as.data.frame(subset(gf_only_importance, !(gf_only_importance[,1] %in% strep_importance_top[,1])))
gf_only_importance <- as.data.frame(subset(gf_only_importance, !(gf_only_importance[,1] %in% cef_importance_top[,1])))
gf_only_importance <- as.data.frame(subset(gf_only_importance, !(gf_only_importance$Compound_name %in% c('Peptide')))) 
#if (nrow(gf_only_importance) > 10) {gf_only_importance <- gf_only_importance[1:10,]}

rm(cef_importance, clinda_importance, strep_importance, gf_importance,
   cef_importance_top, clinda_importance_top, strep_importance_top, gf_importance_top)
cef_only_importance <- cef_only_importance[order(cef_only_importance$Metabolite_score),]
clinda_only_importance <- clinda_only_importance[order(clinda_only_importance$Metabolite_score),]
strep_only_importance <- strep_only_importance[order(strep_only_importance$Metabolite_score),]
gf_only_importance <- gf_only_importance[order(gf_only_importance$Metabolite_score),]
cef_only_importance$abx <- 'Cefoperazone-pretreated'
clinda_only_importance$abx <- 'Clindamycin-pretreated'
strep_only_importance$abx <- 'Streptomycin-pretreated'
gf_only_importance$abx <- 'ex-Germfree'
top_importances <- rbind(cef_only_importance, clinda_only_importance, strep_only_importance, gf_only_importance)
top_importances$abx <- as.factor(top_importances$abx)
top_importances$abx <- ordered(top_importances$abx, levels=c('Streptomycin-pretreated', 'Cefoperazone-pretreated', 'Clindamycin-pretreated', 'ex-Germfree'))
rm(cef_only_importance, clinda_only_importance, strep_only_importance, gf_only_importance)

# Format names to look better for the plot
top_importances$Compound_name <- gsub('_',' ',top_importances$Compound_name)
top_importances$Compound_name <- gsub('-phosphate','P',top_importances$Compound_name)
top_importances$Compound_name <- gsub('Salicin ','Salicin-',top_importances$Compound_name)
top_importances$Compound_name <- gsub('D-Ribose ','D-Ribose-',top_importances$Compound_name)
top_importances$Compound_name <- gsub('-monophosphate','P',top_importances$Compound_name)
top_importances$Compound_name <- gsub('2-Hydroxy-3-m','3-M',top_importances$Compound_name)
top_importances$Compound_name <- gsub('2-Hydroxy-4-h','4-H',top_importances$Compound_name)
top_importances$Compound_name <- gsub('2-Hydroxy-3-c','3-C',top_importances$Compound_name)
top_importances$Compound_name[top_importances$Compound_name == '1-(5\'-Phosphoribosyl)-5-amino-4-(N-succinocarboxamide)-imidazole'] <- 'SAICAR' # shorten a long name
top_importances$Compound_name[top_importances$Compound_name == 'Nicotinamide-beta-riboside'] <- expression(Nicotinamide- ~ beta ~ -riboside)
shared_importance$Compound_name <- gsub('_',' ',shared_importance$Compound_name)
shared_importance$Compound_name[shared_importance$Compound_name == 'CO2'] <- expression(CO[2])

# Combine Shared and top hits
importances <- rbind(shared_importance, top_importances)
importances$Metabolite_score <- as.numeric(as.character(importances$Metabolite_score))
importances$abx <- ordered(importances$abx, levels=c('Median Shared Importance', 'Streptomycin-pretreated', 'Cefoperazone-pretreated', 'Clindamycin-pretreated', 'ex-Germfree'))
rm(shared_importance, top_importances)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in growth rate data
# Define variables
growth_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/cd630_growth.tsv'

# Read in data
growth <- read.delim(growth_file, sep='\t', header=TRUE, row.names=1)
growth <- as.data.frame(t(growth))
rm(growth_file)

# Seperate to groups of each growth substrate and format
fructose <- cbind(growth$fructose_1, growth$fructose_2, growth$fructose_3) - growth$fructose_blank
fructose[fructose < 0] <- 0
mannitol <- cbind(growth$mannitol_1, growth$mannitol_2, growth$mannitol_3) - growth$mannitol_blank
mannitol[mannitol < 0] <- 0
salicin <- cbind(growth$salicin_1, growth$salicin_2, growth$salicin_3) - growth$salicin_blank
salicin[salicin < 0] <- 0
galactitol <- cbind(growth$galactitol_1, growth$galactitol_2, growth$galactitol_3) - growth$galactitol_blank
galactitol[galactitol < 0] <- 0
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
sorbitol <- cbind(growth$sorbitol_1, growth$sorbitol_2, growth$sorbitol_3) - growth$sorbitol_blank
sorbitol[sorbitol < 0] <- 0
rm(growth)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format growth curves

# Find medians 
fructose_median <-  apply(fructose, 1, median)
mannitol_median <- apply(mannitol, 1, median)
salicin_median <- apply(salicin, 1, median)
acetylglucosamine_median <- apply(acetylglucosamine, 1, median)
acetylneuraminate_median <- apply(acetylneuraminate, 1, median)
galactitol_median <- apply(galactitol, 1, median)
no_carb_median <- apply(no_carb, 1, median)
bhi_median <- apply(bhi, 1, median)
no_aa_median <- apply(no_aa, 1, median)
sorbitol_median <- apply(sorbitol, 1, median)

# Standard deviations
acetylneuraminate_sd <- rowSds(acetylneuraminate)
fructose_sd <-  rowSds(fructose)
mannitol_sd <- rowSds(mannitol)
salicin_sd <- rowSds(salicin)
acetylglucosamine_sd <- rowSds(acetylglucosamine)
no_carb_sd <- rowSds(no_carb)
bhi_sd <- rowSds(bhi)
no_aa_sd <- rowSds(no_aa)
sorbitol_sd <- rowSds(sorbitol)
galactitol_sd <- rowSds(galactitol)

# Compile results
growth_medians <- as.data.frame(cbind(acetylglucosamine_median, sorbitol_median, galactitol_median, 
                                      acetylneuraminate_median, fructose_median, mannitol_median, salicin_median, no_carb_median, no_aa_median))
growth_sds <- as.data.frame(cbind(acetylglucosamine_sd, acetylneuraminate_sd, sorbitol_sd, galactitol_sd, 
                                  fructose_sd, mannitol_sd, salicin_sd, no_carb_sd, no_aa_sd))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Analyze growth curves
substrates <- c('acetylglucosamine','acetylneuraminate','sorbitol', 'galactitol','mannitol','salicin','bhi','no_amino_acids','no_carbohydrates')

# Maximum growth rate
max_rate <- round(c(diff(acetylglucosamine_median)[which.max(diff(acetylglucosamine_median))], 
                    diff(acetylneuraminate_median)[which.max(diff(acetylneuraminate_median))],
                    diff(sorbitol_median)[which.max(diff(sorbitol_median))], diff(galactitol_median)[which.max(diff(galactitol_median))],
                    diff(mannitol_median)[which.max(diff(mannitol_median))], diff(salicin_median)[which.max(diff(salicin_median))], 
                    diff(bhi_median)[which.max(diff(bhi_median))], diff(no_aa_median)[which.max(diff(no_aa_median))], 
                    diff(no_carb_median)[which.max(diff(no_carb_median))]), digits=3)

# Time of maximum growth rate
time_max_rate <- round(c((which.max(diff(acetylglucosamine_median)) * 0.5),
                         (which.max(diff(acetylneuraminate_median)) * 0.5), (which.max(diff(sorbitol_median)) * 0.5),
                         (which.max(diff(galactitol_median)) * 0.5), (which.max(diff(mannitol_median)) * 0.5), (which.max(diff(salicin_median)) * 0.5), 
                         (which.max(diff(bhi_median)) * 0.5), (which.max(diff(no_aa_median)) * 0.5), (which.max(diff(no_carb_median)) * 0.5)), digits=3) - 0.5
# Maximum OD
max_od <- round(c(max(acetylglucosamine_median), max(acetylneuraminate_median), max(sorbitol_median),
                  max(galactitol_median), max(mannitol_median), max(salicin_median), max(bhi_median), max(no_aa_median), max(no_carb_median)), digits=3)

# Time of max OD
time_max_od <- round(c((which.max(acetylglucosamine_median) * 0.5),
                       (which.max(acetylneuraminate_median) * 0.5), (which.max(sorbitol_median) * 0.5),
                       (which.max(galactitol_median) * 0.5), (which.max(mannitol_median) * 0.5), (which.max(salicin_median) * 0.5), 
                       (which.max(bhi_median) * 0.5), (which.max(no_aa_median) * 0.5), (which.max(no_carb_median) * 0.5)), digits=3) - 0.5

# Growth rate at 24 hours
rate_24_hrs <- round(c(diff(acetylglucosamine_median)[length(diff(acetylglucosamine_median))], diff(acetylneuraminate_median)[length(diff(acetylneuraminate_median))], 
                       diff(sorbitol_median)[length(diff(sorbitol_median))], diff(galactitol_median)[length(diff(galactitol_median))], 
                       diff(mannitol_median)[length(diff(mannitol_median))], diff(salicin_median)[length(diff(salicin_median))], 
                       diff(bhi_median)[length(diff(bhi_median))], diff(no_aa_median)[length(diff(no_aa_median))], diff(no_carb_median)[length(diff(no_carb_median))]), digits=3)

# Mean growth rate
mean_rate <- round(c(mean(diff(acetylglucosamine_median)), mean(diff(acetylneuraminate_median)), mean(diff(sorbitol_median)),
                     mean(diff(galactitol_median)), mean(diff(mannitol_median)), mean(diff(salicin_median)), mean(diff(bhi_median)), mean(diff(no_aa_median)), mean(diff(no_carb_median))), digits=3)

# Area under curve
area_under <- round(c(auc(acetylglucosamine_median, seq(1,37,1)),
                      auc(acetylneuraminate_median, seq(1,37,1)), auc(sorbitol_median, seq(1,37,1)),  auc(galactitol_median, seq(1,37,1)),
                      auc(mannitol_median, seq(1,37,1)), auc(salicin_median, seq(1,37,1)), auc(bhi_median, seq(1,37,1)), auc(no_aa_median, seq(1,37,1)), auc(no_carb_median, seq(1,37,1))), digits=3)

# Assemble the table
growth_summary <- cbind(substrates, max_rate, time_max_rate, max_od, time_max_od, rate_24_hrs, mean_rate, area_under)
colnames(growth_summary) <- c('Substrate', 'Max_Growth_Rate', 'Time_of_Max_Rate_in_Hours', 'Max_OD', 'Time_of_Max_OD_in_Hours', 'Rate_at_18_hours', 'Mean_Rate', 'AUC')

# Write growth summary data to supplementary table
table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/Table_S4.tsv'
write.table(growth_summary, file=table_file, quote=FALSE, sep='\t', row.names=FALSE)
rm(table_file, growth_summary)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Cut data down to 18 hours for plot
growth_medians <- growth_medians[1:37,]
growth_sds <- growth_sds[1:37,]

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_4.pdf'
pdf(file=plot_file, width=6, height=12)
layout(matrix(c(1,
                1,
                2), nrow=3, ncol=1, byrow=TRUE))

#---------------------------------------#

# Metabolite importances
par(mar=c(3,3,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(importances$Metabolite_score, labels=importances$Compound_name,
         lcolor=NA, cex=0.85, groups=importances$abx, color='black',
         xlab='Importance Score', xlim=c(0,10), pch=19, lwd=3, 
         gcolor=c('darkmagenta',wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5],'forestgreen'))
mtext('A', side=2, line=2, las=2, adj=0.6, padj=-18, cex=1.7)
segments(x0=rep(0, 40), y0=c(1:21, 24:27, 30:35, 38:44, 47:52), 
         x1=rep(12, 40), y1=c(1:21, 24:27, 30:35, 38:44, 47:52), lty=2) # Dotted lines

#---------------------------------------#

# Growth on important carbohydrates
par(mar=c(3.5,4.5,1,1), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE, mgp=c(2.5,1,0), xaxs='r')
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(1,37), ylim=c(-0.03,1.0), pch=15, xlab='', ylab='')
abline(h=seq(0,1.0,0.1), lty=3, col='gray68') # adding gridlines
abline(v=seq(1,37,2), lty=3, col='gray68') # adding gridlines
axis(1, at=seq(1,37,4), labels=seq(0,18,2), tck=-0.018, cex.axis=1.2)
mtext('Time (hours)', side=1, at=19, padj=2.5, cex=)
axis(2, at=seq(0.0,1.0,0.2), labels=c('0.0','0.2','0.4','0.6','0.8','1.0'), tck=-0.018, cex.axis=1.2)
text(x=-4, y=0.5, labels=expression(OD[600]), cex=1.5, xpd=TRUE, srt=90)
mtext('B', side=2, line=2, las=2, adj=1.5, padj=-8, cex=1.7)
legend('topleft', legend=c('No Carbohydrates','No Amino Acids','N-Acetylglucosamine','D-Sorbitol', 'Mannitol','Salicin','N-Acetylneuraminate'), 
       col=c('black','black','darkmagenta',wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5],'forestgreen'), 
       pch=c(15,19,6,0,1,2,5), cex=1.2, pt.cex=c(1.2,1.2,2,2.4,2.4,2,2), lwd=3, bg='white')

# Controls
lines(growth_medians$no_carb_median, type='o', lwd=3, pch=15, cex=1.1, col='black') # Carbohydrates
segments(x0=seq(1,37,1), y0=growth_medians$no_carb_median+growth_sds$no_carb_sd, 
         x1=seq(1,37,1), y1=growth_medians$no_carb_median-growth_sds$no_carb_sd, lwd=2.5, cex=2, col='black')
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$no_carb_median+growth_sds$no_carb_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$no_carb_median+growth_sds$no_carb_sd, lwd=2.5, col='black')
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$no_carb_median-growth_sds$no_carb_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$no_carb_median-growth_sds$no_carb_sd, lwd=2.5, col='black')
lines(growth_medians$no_aa_median, type='o', lwd=3, pch=19, cex=1.1, col='black') # Amino Acids
segments(x0=seq(1,37,1), y0=growth_medians$no_aa_median+growth_sds$no_aa_sd, 
         x1=seq(1,37,1), y1=growth_medians$no_aa_median-growth_sds$no_aa_sd, lwd=2.5, cex=2, col='black')
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$no_aa_median+growth_sds$no_aa_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$no_aa_median+growth_sds$no_aa_sd, lwd=2.5, col='black')
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$no_aa_median-growth_sds$no_aa_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$no_aa_median-growth_sds$no_aa_sd, lwd=2.5, col='black')

# Shared
lines(growth_medians$acetylglucosamine_median, type='o', col='darkmagenta', lwd=2.5, pch=6, cex=1.3)
segments(x0=seq(1,37,1), y0=growth_medians$acetylglucosamine_median+growth_sds$acetylglucosamine_sd, 
         x1=seq(1,37,1), y1=growth_medians$acetylglucosamine_median-growth_sds$acetylglucosamine_sd, lwd=2.5, col='darkmagenta')
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$acetylglucosamine_median+growth_sds$acetylglucosamine_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$acetylglucosamine_median+growth_sds$acetylglucosamine_sd, lwd=2.5, col='darkmagenta')
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$acetylglucosamine_median-growth_sds$acetylglucosamine_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$acetylglucosamine_median-growth_sds$acetylglucosamine_sd, lwd=2.5, col='darkmagenta')

# Streptomycin
lines(growth_medians$sorbitol_median, type='o', col=wes_palette('FantasticFox')[1], lwd=2.5, pch=0, cex=1.6)
segments(x0=seq(1,37,1), y0=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, 
         x1=seq(1,37,1), y1=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])

# Cefoperazone
lines(growth_medians$mannitol_median, type='o', col=wes_palette('FantasticFox')[3], lwd=2.5, pch=1, cex=1.7)
segments(x0=seq(1,37,1), y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, 
         x1=seq(1,37,1), y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$mannitol_median+growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$mannitol_median-growth_sds$mannitol_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])

# Clindamycin
lines(growth_medians$salicin_median, type='o', col=wes_palette('FantasticFox')[5], lwd=2.5, pch=2, cex=1.5)
segments(x0=seq(1,37,1), y0=growth_medians$salicin_median+growth_sds$salicin_sd, 
         x1=seq(1,37,1), y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$salicin_median+growth_sds$salicin_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$salicin_median+growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$salicin_median-growth_sds$salicin_sd, 
         x1=seq(1,37,1)+0.2, y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])

# Gnotobiotic 
lines(growth_medians$acetylneuraminate_median, type='o', col='forestgreen', lwd=2.5, pch=5, cex=1.5) 
segments(x0=seq(1,37,1), y0=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd,  
         x1=seq(1,37,1), y1=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd, lwd=2.5, col='forestgreen') 
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd,  
         x1=seq(1,37,1)+0.2, y1=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd, lwd=2.5, col='forestgreen') 
segments(x0=seq(1,37,1)-0.2, y0=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd,  
         x1=seq(1,37,1)+0.2, y1=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd, lwd=2.5, col='forestgreen') 

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()
