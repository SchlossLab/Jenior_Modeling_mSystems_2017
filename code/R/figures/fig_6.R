
# Load dependencies
deps <- c('vegan', 'ggplot2', 'shape', 'wesanderson', 'flux', 'matrixStats');
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

# Read in growth rate data
# Define variables
growth_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/formatted_growth.tsv'

# Read in data
growth <- read.delim(growth_file, sep='\t', header=TRUE, row.names=1)
growth <- as.data.frame(t(growth))
rm(growth_file)

# Seperate to groups of each growth substrate and format

# add acetylglucosamine

sorbitol <- cbind(growth$B9, growth$B10, growth$B11) - growth$B8
sorbitol[sorbitol < 0] <- 0
fructose <- cbind(growth$E9, growth$E10, growth$E11) - growth$E8
fructose[fructose < 0] <- 0
mannitol <- cbind(growth$F9, growth$F10, growth$F11) - growth$F8
mannitol[mannitol < 0] <- 0
salicin <- cbind(growth$G9, growth$G10, growth$G11) - growth$G8
salicin[salicin < 0] <- 0
acetylglucosamine <- cbind(growth$D9, growth$D10, growth$D11) - growth$D8
acetylglucosamine[acetylglucosamine < 0] <- 0
acetylneuraminate <- cbind(growth$C9, growth$C10, growth$C11) - growth$C8
acetylneuraminate[acetylneuraminate < 0] <- 0
y_glucose_y_aa <- cbind(growth$B3, growth$B4, growth$B5) - growth$B2
y_glucose_y_aa[y_glucose_y_aa < 0] <- 0
n_glucose_y_aa <- cbind(growth$D3, growth$D4, growth$D5) - growth$D2
n_glucose_y_aa[n_glucose_y_aa < 0] <- 0
y_glucose_n_aa <- cbind(growth$C3, growth$C4, growth$C5) - growth$C2
y_glucose_n_aa[y_glucose_n_aa < 0] <- 0
n_glucose_n_aa <- cbind(growth$E3, growth$E4, growth$E5) - growth$E2
n_glucose_n_aa[n_glucose_n_aa < 0] <- 0
bhi <- cbind(growth$F3, growth$F4, growth$F5) - growth$F2
bhi[bhi < 0] <- 0

#-------------------------------------------------------------------------------------------------------------------------------------#

# Prepare data for statistical tests
fructose_test <- format_curve(fructose, 'fructose', n_glucose_y_aa)
sorbitol_test <- format_curve(sorbitol, 'sorbitol', n_glucose_y_aa)
mannitol_test <- format_curve(mannitol, 'mannitol', n_glucose_y_aa)
salicin_test <- format_curve(salicin, 'salicin', n_glucose_y_aa)
acetylneuraminate_test <- format_curve(acetylneuraminate, 'acetylneuraminate', n_glucose_y_aa)
y_glucose_y_aa_test <- format_curve(y_glucose_y_aa, 'y_glucose_y_aa', n_glucose_y_aa)
y_glucose_n_aa_test <- format_curve(y_glucose_n_aa, 'y_glucose_n_aa', n_glucose_y_aa)
n_glucose_n_aa_test <- format_curve(n_glucose_n_aa, 'n_glucose_n_aa', n_glucose_y_aa)
bhi_test <- format_curve(bhi, 'bhi', n_glucose_y_aa)

# Calculate differences
sorbitol_sig <- aov(formula=od ~ substrate * time, data=sorbitol_test)
summary(sorbitol_sig) # p = 0.0494 *, corrected = 4.940e-01 n.s.
fructose_sig <- aov(formula=od ~ substrate * time, data=fructose_test)
summary(fructose_sig) # p < 2e-16 ***, corrected = 2.000e-15 ***
mannitol_sig <- aov(formula=od ~ substrate * time, data=mannitol_test)
summary(mannitol_sig) # p < 2e-16 ***, corrected = 2.000e-15 ***
salicin_sig <- aov(formula=od ~ substrate * time, data=salicin_test)
summary(salicin_sig) # p < 2e-16 ***, corrected = 2.000e-15 ***
acetylneuraminate_sig <- aov(formula=od ~ substrate * time, data=acetylneuraminate_test)
summary(acetylneuraminate_sig) # p < 2e-16 ***, corrected = 2.000e-15 ***
y_glucose_y_aa_sig <- aov(formula=od ~ substrate * time, data=y_glucose_y_aa_test)
summary(y_glucose_y_aa_sig) # p < 2e-16 ***, corrected = 2.000e-15 ***
y_glucose_n_aa_sig <- aov(formula=od ~ substrate * time, data=y_glucose_y_aa_test)
summary(y_glucose_n_aa_sig) # p < 2e-16 ***, corrected = 2.000e-15 ***
n_glucose_n_aa_sig <- aov(formula=od ~ substrate * time, data=y_glucose_y_aa_test)
summary(n_glucose_n_aa_sig) # p < 2e-16 ***, corrected = 2.000e-15 ***
bhi_sig <- aov(formula=od ~ substrate * time, data=bhi_test)
summary(bhi_sig) # p < 2e-16 ***, corrected = 2.000e-15 ***

c('acetylneuraminate','sorbitol', 'fructose', 'mannitol','salicin','y_glucose_y_aa','n_glucose_y_aa','y_glucose_n_aa','n_glucose_n_aa','bhi')
p_values <- c(2e-16, 0.0494, 2e-16, 2e-16, 2e-16, 2e-16, 2e-16, 2e-16, 2e-16)
corrected_p_values <- as.character(p.adjust(p_values, method='bonferroni'))
corrected_p_values <- append(corrected_p_values, 'NA', after=5) 

# Clean up
rm(p_values)
rm(fructose_test, sorbitol_test, mannitol_test, salicin_test, y_glucose_y_aa_test, y_glucose_n_aa_test, n_glucose_n_aa_test, bhi_test, acetylneuraminate_test)
rm(fructose_sig, sorbitol_sig, mannitol_sig, salicin_sig, y_glucose_y_aa_sig, y_glucose_n_aa_sig, n_glucose_n_aa_sig, bhi_sig, acetylneuraminate_sig)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format growth curves (individual)

# Find medians of treatement groups
sorbitol_median <- apply(sorbitol, 1, median)
fructose_median <-  apply(fructose, 1, median)
mannitol_median <- apply(mannitol, 1, median)
salicin_median <- apply(salicin, 1, median)
acetylneuraminate_median <- apply(acetylneuraminate, 1, median)
y_glucose_y_aa_median <- apply(y_glucose_y_aa, 1, median)
n_glucose_y_aa_median <- apply(n_glucose_y_aa, 1, median)
y_glucose_n_aa_median <- apply(y_glucose_n_aa, 1, median)
n_glucose_n_aa_median <- apply(n_glucose_n_aa, 1, median)
bhi_median <- apply(bhi, 1, median)
growth_medians <- as.data.frame(rbind(acetylneuraminate_median, sorbitol_median, fructose_median, mannitol_median, salicin_median, 
                                      y_glucose_y_aa_median, n_glucose_y_aa_median, y_glucose_n_aa_median, n_glucose_n_aa_median, bhi_median))

# Determine some features of the 12 hour growth curves
substrates <- c('acetylneuraminate','sorbitol', 'fructose', 'mannitol','salicin','y_glucose_y_aa','n_glucose_y_aa','y_glucose_n_aa','n_glucose_n_aa','bhi')

# Maximum growth rate
max_rate <- round(c(diff(acetylneuraminate_median)[which.max(diff(acetylneuraminate_median))],diff(sorbitol_median)[which.max(diff(sorbitol_median))], 
                    diff(fructose_median)[which.max(diff(fructose_median))], diff(mannitol_median)[which.max(diff(mannitol_median))], diff(salicin_median)[which.max(diff(salicin_median))],
                    diff(y_glucose_y_aa_median)[which.max(diff(y_glucose_y_aa_median))], diff(n_glucose_y_aa_median)[which.max(diff(n_glucose_y_aa_median))], diff(y_glucose_n_aa_median)[which.max(diff(y_glucose_n_aa_median))],
                    diff(n_glucose_n_aa_median)[which.max(diff(n_glucose_n_aa_median))], diff(bhi_median)[which.max(diff(bhi_median))]), digits=3)

# Time of maximum growth rate
time_max_rate <- round(c((which.max(diff(acetylneuraminate_median)) * 0.5), (which.max(diff(sorbitol_median)) * 0.5), 
                         (which.max(diff(fructose_median)) * 0.5), (which.max(diff(mannitol_median)) * 0.5), (which.max(diff(salicin_median)) * 0.5),
                         (which.max(diff(y_glucose_y_aa_median)) * 0.5), (which.max(diff(n_glucose_y_aa_median)) * 0.5), (which.max(diff(y_glucose_n_aa_median)) * 0.5),
                         (which.max(diff(n_glucose_n_aa_median)) * 0.5), (which.max(diff(bhi_median)) * 0.5)), digits=3) - 0.5
# Maximum OD
max_od <- round(c(max(acetylneuraminate_median), max(sorbitol_median), max(fructose_median), max(mannitol_median), max(salicin_median), 
                  max(y_glucose_y_aa_median), max(n_glucose_y_aa_median), max(y_glucose_n_aa_median), max(n_glucose_n_aa_median), max(bhi_median)), digits=3)

# Time of max OD
time_max_od <- round(c((which.max(acetylneuraminate_median) * 0.5), (which.max(sorbitol_median) * 0.5), 
                       (which.max(fructose_median) * 0.5), (which.max(mannitol_median) * 0.5), (which.max(salicin_median) * 0.5),
                       (which.max(y_glucose_y_aa_median) * 0.5), (which.max(n_glucose_y_aa_median) * 0.5), (which.max(y_glucose_n_aa_median) * 0.5),
                       (which.max(n_glucose_n_aa_median) * 0.5), (which.max(bhi_median) * 0.5)), digits=3) - 0.5

# Growth rate at 24 hours
rate_24_hrs <- round(c(diff(acetylneuraminate_median)[length(diff(acetylneuraminate_median))], diff(sorbitol_median)[length(diff(sorbitol_median))], 
                       diff(fructose_median)[length(diff(fructose_median))], diff(mannitol_median)[length(diff(mannitol_median))], diff(salicin_median)[length(diff(salicin_median))],
                       diff(y_glucose_y_aa_median)[length(diff(y_glucose_y_aa_median))], diff(n_glucose_y_aa_median)[length(diff(n_glucose_y_aa_median))], diff(y_glucose_n_aa_median)[length(diff(y_glucose_n_aa_median))],
                       diff(n_glucose_n_aa_median)[length(diff(n_glucose_n_aa_median))], diff(bhi_median)[length(diff(bhi_median))]), digits=3)

# Mean growth rate
mean_rate <- round(c(mean(diff(acetylneuraminate_median)), mean(diff(sorbitol_median)), mean(diff(fructose_median)), mean(diff(mannitol_median)),
                     mean(diff(salicin_median)), mean(diff(n_glucose_y_aa_median)), mean(diff(y_glucose_y_aa_median)), mean(diff(y_glucose_n_aa_median)), 
                     mean(diff(n_glucose_n_aa_median)), mean(diff(bhi_median))), digits=3)

# Area under curve
area_under <- round(c(auc(acetylneuraminate_median, seq(1,49,1)), auc(sorbitol_median, seq(1,49,1)), 
                      auc(fructose_median, seq(1,49,1)), auc(mannitol_median, seq(1,49,1)), auc(salicin_median, seq(1,49,1)),
                      auc(y_glucose_y_aa_median, seq(1,49,1)), auc(n_glucose_y_aa_median, seq(1,49,1)), auc(y_glucose_n_aa_median, seq(1,49,1)), 
                      auc(n_glucose_n_aa_median, seq(1,49,1)), auc(bhi_median, seq(1,49,1)), digits=3))

# Assemble the table
growth_summary <- cbind(substrates, max_rate, time_max_rate, max_od, time_max_od, rate_24_hrs, mean_rate, area_under, corrected_p_values)
colnames(growth_summary) <- c('Substrate', 'Max_Growth_Rate', 'Time_of_Max_Rate_in_Hours', 'Max_OD', 'Time_of_Max_OD_in_Hours', 'Rate_at_24_hours', 'Mean_Rate', 'AUC', 'Corrected_AoV_pvalue')
rm(substrates, max_rate, time_max_rate, max_od, time_max_od, rate_24_hrs, mean_rate, area_under, corrected_p_values)

# Write growth summary data to supplementary table
table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S6.tsv'
write.table(growth_summary, file=table_file, quote=FALSE, sep='\t', row.names=FALSE)
rm(table_file, growth_summary)

# Standard deviations
sorbitol_sd <- rowSds(sorbitol, na.rm=TRUE)
acetylneuraminate_sd <- rowSds(acetylneuraminate, na.rm=TRUE)
fructose_sd <-  rowSds(fructose, na.rm=TRUE)
mannitol_sd <- rowSds(mannitol, na.rm=TRUE)
salicin_sd <- rowSds(salicin, na.rm=TRUE)
y_glucose_y_aa_sd <- rowSds(y_glucose_y_aa, na.rm=TRUE)
n_glucose_y_aa_sd <- rowSds(n_glucose_y_aa, na.rm=TRUE)
y_glucose_n_aa_sd <- rowSds(y_glucose_n_aa, na.rm=TRUE)
n_glucose_n_aa_sd <- rowSds(n_glucose_n_aa, na.rm=TRUE)
bhi_sd <- rowSds(bhi, na.rm=TRUE)
growth_sds <- as.data.frame(rbind(acetylneuraminate_sd, sorbitol_sd, fructose_sd, mannitol_sd, salicin_sd, y_glucose_y_aa_sd, n_glucose_y_aa_sd, y_glucose_n_aa_sd, n_glucose_n_aa_sd, bhi_sd))
rm(acetylneuraminate_median, sorbitol_median, fructose_median, mannitol_median, salicin_median, y_glucose_y_aa_median, n_glucose_y_aa_median, y_glucose_n_aa_median, n_glucose_n_aa_median, bhi_median)
rm(acetylneuraminate, sorbitol, fructose, mannitol, salicin, y_glucose_y_aa, y_glucose_n_aa, n_glucose_n_aa, bhi)
rm(acetylneuraminate_sd, sorbitol_sd, fructose_sd, mannitol_sd, salicin_sd, y_glucose_y_aa_sd, n_glucose_y_aa_sd, y_glucose_n_aa_sd, n_glucose_n_aa_sd, bhi_sd)
rm(growth)

# Subset for first 12 hours of assay
growth_medians <- as.data.frame(t(growth_medians[,1:49]))
growth_sds <- as.data.frame(t(growth_sds[,1:49]))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format growth curves (combined)

# Read in growth rate data
# Define variables
growth_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/pairwise-carb-combos.formatted.tsv'

# Read in data
growth <- read.delim(growth_file, sep='\t', header=TRUE, row.names=1)
growth <- as.data.frame(t(growth))
rm(growth_file)

# Seperate to groups of each growth substrate 
fructose_sorbitol <- cbind(growth$C3, growth$C4, growth$C5) - growth$C2
fructose_sorbitol[fructose_sorbitol < 0] <- 0
fructose_sorbitol_median <- rowMeans(fructose_sorbitol, na.rm=TRUE)
fructose_sorbitol_sd <- rowSds(fructose_sorbitol, na.rm=TRUE)
fructose_salicin <- cbind(growth$D3, growth$D4, growth$D5) - growth$D2
fructose_salicin[fructose_salicin < 0] <- 0
fructose_salicin_median <- rowMeans(fructose_salicin, na.rm=TRUE)
fructose_salicin_median[fructose_salicin_median < 0] <- 0
fructose_salicin_sd <- rowSds(fructose_salicin, na.rm=TRUE)
fructose_mannitol <- cbind(growth$E3, growth$E4, growth$E5) - growth$E2
fructose_mannitol[fructose_mannitol < 0] <- 0
fructose_mannitol_median <- rowMeans(fructose_mannitol, na.rm=TRUE)
fructose_mannitol_median[fructose_mannitol_median < 0] <- 0
fructose_mannitol_sd <- rowSds(fructose_mannitol, na.rm=TRUE)
fructose_acetylneuriminate <- cbind(growth$F3, growth$F4, growth$F5) - growth$F2
fructose_acetylneuriminate[fructose_acetylneuriminate < 0] <- 0
fructose_acetylneuriminate_median <- rowMeans(fructose_acetylneuriminate, na.rm=TRUE)
fructose_acetylneuriminate_median[fructose_acetylneuriminate_median < 0] <- 0
fructose_acetylneuriminate_sd <- rowSds(fructose_acetylneuriminate, na.rm=TRUE)
sorbitol_salicin <- cbind(growth$G3, growth$G4, growth$G5) - growth$G2
sorbitol_salicin[sorbitol_salicin < 0] <- 0
sorbitol_salicin_median <- rowMeans(sorbitol_salicin, na.rm=TRUE)
sorbitol_salicin_median[sorbitol_salicin_median < 0] <- 0
sorbitol_salicin_sd <- rowSds(sorbitol_salicin, na.rm=TRUE)
sorbitol_mannitol <- cbind(growth$B8, growth$B9, growth$B10) - growth$B7
sorbitol_mannitol[sorbitol_mannitol < 0] <- 0
sorbitol_mannitol_median <- rowMeans(sorbitol_mannitol, na.rm=TRUE)
sorbitol_mannitol_median[sorbitol_mannitol_median < 0] <- 0
sorbitol_mannitol_sd <- rowSds(sorbitol_mannitol, na.rm=TRUE)
sorbitol_mannitol <- cbind(growth$C8, growth$C9, growth$C10) - growth$C7
sorbitol_mannitol[sorbitol_mannitol < 0] <- 0
sorbitol_acetylneuriminate_median <- rowMeans(sorbitol_acetylneuriminate, na.rm=TRUE)
sorbitol_acetylneuriminate_median[sorbitol_acetylneuriminate_median < 0] <- 0
sorbitol_acetylneuriminate_sd <- rowSds(sorbitol_acetylneuriminate, na.rm=TRUE)
salicin_mannitol <- cbind(growth$D8, growth$D9, growth$D10) - growth$D7
salicin_mannitol[salicin_mannitol < 0] <- 0
salicin_mannitol_median <- rowMeans(salicin_mannitol, na.rm=TRUE)
salicin_mannitol_median[salicin_mannitol_median < 0] <- 0
salicin_mannitol_sd <- rowSds(salicin_mannitol, na.rm=TRUE)
salicin_acetylneuriminate <- cbind(growth$E8, growth$E9, growth$E10) - growth$E7
salicin_acetylneuriminate[salicin_acetylneuriminate < 0] <- 0
salicin_acetylneuriminate_median <- rowMeans(salicin_acetylneuriminate, na.rm=TRUE)
salicin_acetylneuriminate_median[salicin_acetylneuriminate_median < 0] <- 0
salicin_acetylneuriminate_sd <- rowSds(salicin_acetylneuriminate, na.rm=TRUE)
mannitol_acetylneuriminate <- cbind(growth$F8, growth$F9, growth$F10) - growth$F7
mannitol_acetylneuriminate[mannitol_acetylneuriminate < 0] <- 0

mannitol_acetylneuriminate_median <- rowMeans(mannitol_acetylneuriminate, na.rm=TRUE)
mannitol_acetylneuriminate_median[mannitol_acetylneuriminate_median < 0] <- 0
mannitol_acetylneuriminate_sd <- rowSds(mannitol_acetylneuriminate, na.rm=TRUE)

combo_medians <- t(as.data.frame(rbind(bhi_median, fructose_sorbitol_median, fructose_salicin_median, fructose_mannitol_median, 
                                       fructose_acetylneuriminate_median, sorbitol_salicin_median, sorbitol_mannitol_median, 
                                       sorbitol_acetylneuriminate_median, salicin_mannitol_median, salicin_acetylneuriminate_median, 
                                       mannitol_acetylneuriminate_median)))
combo_sds <- t(as.data.frame(rbind(bhi_sd, fructose_sorbitol_sd, fructose_salicin_sd, fructose_mannitol_sd, 
                                   fructose_acetylneuriminate_sd, sorbitol_salicin_sd, sorbitol_mannitol_sd, 
                                   sorbitol_acetylneuriminate_sd, salicin_mannitol_sd, salicin_acetylneuriminate_sd, 
                                   mannitol_acetylneuriminate_sd)))

fructose_sorbitol_test <- format_curve(fructose_sorbitol, 'fructose_sorbitol', n_glucose_y_aa)
fructose_salicin_test <- format_curve(fructose_salicin, 'fructose_salicin', n_glucose_y_aa)
fructose_mannitol_test <- format_curve(fructose_mannitol, 'fructose_mannitol', n_glucose_y_aa)
fructose_acetylneuriminate_test <- format_curve(fructose_acetylneuriminate, 'fructose_acetylneuriminate', n_glucose_y_aa)
sorbitol_salicin_test <- format_curve(sorbitol_salicin, 'sorbitol_salicin', n_glucose_y_aa)
sorbitol_mannitol_test <- format_curve(sorbitol_mannitol, 'sorbitol_mannitol', n_glucose_y_aa)
sorbitol_acetylneuriminate_test <- format_curve(sorbitol_acetylneuriminate, 'sorbitol_acetylneuriminate', n_glucose_y_aa)
salicin_mannitol_test <- format_curve(salicin_mannitol, 'salicin_mannitol', n_glucose_y_aa)
salicin_acetylneuriminate_test <- format_curve(salicin_acetylneuriminate, 'salicin_acetylneuriminate', n_glucose_y_aa)
mannitol_acetylneuriminate_test <- format_curve(mannitol_acetylneuriminate, 'mannitol_acetylneuriminate', n_glucose_y_aa)

rm(growth, bhi, fructose_sorbitol, fructose_salicin, fructose_mannitol, 
   fructose_acetylneuriminate, sorbitol_salicin, sorbitol_mannitol, 
   sorbitol_acetylneuriminate, salicin_mannitol, salicin_acetylneuriminate, 
   mannitol_acetylneuriminate, n_glucose_y_aa) 

rm(bhi_median, fructose_sorbitol_median, fructose_salicin_median, fructose_mannitol_median, 
   fructose_acetylneuriminate_median, sorbitol_salicin_median, sorbitol_mannitol_median, 
   sorbitol_acetylneuriminate_median, salicin_mannitol_median, salicin_acetylneuriminate_median, 
   mannitol_acetylneuriminate_median) 

rm(bhi_sd, fructose_sorbitol_sd, fructose_salicin_sd, fructose_mannitol_sd, 
   fructose_acetylneuriminate_sd, sorbitol_salicin_sd, sorbitol_mannitol_sd, 
   sorbitol_acetylneuriminate_sd, salicin_mannitol_sd, salicin_acetylneuriminate_sd, 
   mannitol_acetylneuriminate_sd)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_6.pdf'
pdf(file=plot_file, width=16, height=9)
layout(matrix(c(1,1,2,2,3,
                1,1,2,2,3), nrow=2, ncol=5, byrow=TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Growth on important compounds (separate)
par(mar=c(5,6,1,1), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE, mgp=c(3,1,0))
plot(0, type='n', xaxt='n', xlim=c(0,50), ylim=c(-0.03,1.0), lwd=2, pch=15, xlab='Time (hours)', ylab=expression(OD[600]), cex=2.3)
abline(h=seq(0,1,0.1), lty=3, col='gray68') # adding gridlines
abline(v=seq(1,50,2), lty=3, col='gray68') # adding gridlines
axis(1, at=seq(1,49,4), labels=seq(0,24,2))

# controls
lines(growth_medians$n_glucose_y_aa_median, type='o', lwd=2, pch=17, cex=2.7, col='gray50')
segments(x0=seq(1,49,1), y0=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, x1=seq(1,49,1), y1=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, lwd=3, cex=2, col='gray50')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, lwd=3, col='gray50')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, lwd=3, col='gray50')
lines(growth_medians$y_glucose_n_aa_median, type='o', lwd=2, pch=16, cex=2.7, col='gray50')
segments(x0=seq(1,49,1), y0=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, x1=seq(1,49,1), y1=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, lwd=3, cex=2, col='gray50')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, lwd=3, col='gray50')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, lwd=3, col='gray50')

# needs to be replaced by n-acetylglucosamine
lines(growth_medians$y_glucose_y_aa_median, type='o', lwd=2, pch=1, cex=2, col='black')
segments(x0=seq(1,49,1), y0=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, x1=seq(1,49,1), y1=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, lwd=3, cex=2, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, lwd=3, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, lwd=3, col='black')

lines(growth_medians$fructose_median, type='o', col=wes_palette('FantasticFox')[1], lwd=3, pch=0, cex=2)
segments(x0=seq(1,49,1), y0=growth_medians$fructose_median+growth_sds$fructose_sd, x1=seq(1,49,1), y1=growth_medians$fructose_median-growth_sds$fructose_sd, lwd=3, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$fructose_median+growth_sds$fructose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$fructose_median+growth_sds$fructose_sd, lwd=3, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$fructose_median-growth_sds$fructose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$fructose_median-growth_sds$fructose_sd, lwd=3, col=wes_palette('FantasticFox')[1])
#lines(growth_medians$sorbitol_median, type='o', col=wes_palette('FantasticFox')[1], lwd=3, pch=1, cex=2.5)
#segments(x0=seq(1,49,1), y0=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, x1=seq(1,49,1), y1=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, lwd=3, col=wes_palette('FantasticFox')[1])
#segments(x0=seq(1,49,1)-0.2, y0=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, lwd=3, col=wes_palette('FantasticFox')[1])
#segments(x0=seq(1,49,1)-0.2, y0=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, lwd=3, col=wes_palette('FantasticFox')[1])

lines(growth_medians$mannitol_median, type='o', col=wes_palette('FantasticFox')[3], lwd=3, pch=2, cex=2)
segments(x0=seq(1,49,1), y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, x1=seq(1,49,1), y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=3, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$mannitol_median+growth_sds$mannitol_sd, lwd=3, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$mannitol_median-growth_sds$mannitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=3, col=wes_palette('FantasticFox')[3])

lines(growth_medians$salicin_median, type='o', col=wes_palette('FantasticFox')[5], lwd=3, pch=5, cex=2)
segments(x0=seq(1,49,1), y0=growth_medians$salicin_median+growth_sds$salicin_sd, x1=seq(1,49,1), y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=3, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$salicin_median+growth_sds$salicin_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$salicin_median+growth_sds$salicin_sd, lwd=3, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$salicin_median-growth_sds$salicin_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=3, col=wes_palette('FantasticFox')[5])

lines(growth_medians$acetylneuraminate_median, type='o', col='forestgreen', lwd=3, pch=6, cex=2)
segments(x0=seq(1,49,1), y0=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd, x1=seq(1,49,1), y1=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd, lwd=3, col='forestgreen')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd, lwd=3, col='forestgreen')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd, lwd=3, col='forestgreen')

#legend('topleft', legend=c('No Carbohydrates','No Amino acids','N-Acetylglucosamine (GlcNAc)','D-Fructose','D-Sorbitol','Mannitol','Salicin','N-Acetylneuriminate (Neu5Ac)'), 
#       col=c('black','black','black',wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5],'forestgreen'), 
#       pch=c(17,16,6,0,1,2,5,6), cex=1.3, pt.cex=2.5, bg='white', lwd=3)
legend('topleft', legend=c('No Carbohydrates','No Amino acids','N-Acetylglucosamine (GlcNAc)','D-Fructose','Mannitol','Salicin','N-Acetylneuriminate (Neu5Ac)'), 
       col=c('gray50','gray50','black',wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5],'forestgreen'), 
       pch=c(17,16,1,0,2,5,6), cex=1.4, pt.cex=2.5, bg='white', lwd=3)

mtext('a', side=2, line=2, las=2, adj=2, padj=-21, cex=1.6, font=2)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Growth on important compounds (combinations)
par(mar=c(5,5,1,0), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE, mgp=c(3,1,0))
plot(0, type='n', xaxt='n', xlim=c(0,50), ylim=c(-0.03,1.0), lwd=2, pch=15, xlab='Time (hours)', ylab=expression(OD[600]), cex=2.3)
abline(h=seq(0,1,0.1), lty=3, col='gray68') # adding gridlines
abline(v=seq(1,50,2), lty=3, col='gray68') # adding gridlines
axis(1, at=seq(1,49,4), labels=seq(0,24,2))
mtext('b', side=2, line=2, las=2, adj=2, padj=-21, cex=1.6, font=2)


# Use wes_palette("Cavalcanti") for colors
par(mar=c(0,0,0,0))
plot(0, type='n', axes=F, xlab='', ylab='') # Empty plot


legend('left', legend=c('Fructose + Sorbitol','Fructose + Salicin','Fructose + Mannitol','Fructose + GlcNAc','Fructose + Neu5Ac',
                        'Sorbitol + Salicin','Sorbitol + Mannitol','Sorbitol + GlcNAc','Sorbitol + Neu5Ac',
                        'Salicin + Mannitol','Salicin + GlcNAc','Salicin + Neu5Ac',
                        'Mannitol + GlcNAc','Mannitol + Neu5Ac',
                        'GlcNAc + Neu5Ac'), 
       pch=c(0,1,2,5,6), cex=2, pt.cex=3, bg='white', lwd=3,
       col=c('black'))

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

# Control supplementary plot (Glucose and BHI controls)

plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S6.pdf'
pdf(file=plot_file, width=12, height=9)

par(mar=c(5,6,1,1), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE, mgp=c(3,1,0))
plot(0, type='n', xaxt='n', xlim=c(0,50), ylim=c(-0.03,1.0), lwd=2, pch=15, xlab='Time (hours)', ylab=expression(OD[600]), cex=2.3)
abline(h=seq(0,1,0.1), lty=3, col='gray68') # adding gridlines
abline(v=seq(1,50,2), lty=3, col='gray68') # adding gridlines
axis(1, at=seq(1,49,4), labels=seq(0,24,2))

lines(growth_medians$y_glucose_y_aa_median, type='o', col='black', lwd=2, pch=19, cex=2)
segments(x0=seq(1,49,1), y0=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, x1=seq(1,49,1), y1=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, lwd=3, cex=2)
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, lwd=3)
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, lwd=3)
lines(growth_medians$n_glucose_y_aa_median, type='o', lwd=2, pch=17, cex=2.5)
segments(x0=seq(1,49,1), y0=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, x1=seq(1,49,1), y1=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, lwd=3, cex=2)
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, lwd=3)
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, lwd=3)
lines(growth_medians$y_glucose_n_aa_median, type='o', lwd=2, pch=15, cex=2)
segments(x0=seq(1,49,1), y0=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, x1=seq(1,49,1), y1=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, lwd=3, cex=2)
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, lwd=3)
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, lwd=3)
lines(growth_medians$n_glucose_n_aa_median, type='o', lwd=2, pch=18, cex=2.5)
segments(x0=seq(1,49,1), y0=growth_medians$n_glucose_n_aa_median+growth_sds$n_glucose_n_aa_sd, x1=seq(1,49,1), y1=growth_medians$n_glucose_n_aa_median-growth_sds$n_glucose_n_aa_sd, lwd=3, cex=2)
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$n_glucose_n_aa_median+growth_sds$n_glucose_n_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$n_glucose_n_aa_median+growth_sds$n_glucose_n_aa_sd, lwd=3)
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$n_glucose_n_aa_median-growth_sds$n_glucose_n_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$n_glucose_n_aa_median-growth_sds$n_glucose_n_aa_sd, lwd=3)

lines(growth_medians$bhi_median, type='o', lwd=2, pch=19, cex=2.5, col='azure4')
segments(x0=seq(1,49,1), y0=growth_medians$bhi_median+growth_sds$bhi_sd, x1=seq(1,49,1), y1=growth_medians$bhi_median-growth_sds$bhi_sd, lwd=3, cex=2, col='azure4')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$bhi_median+growth_sds$bhi_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$bhi_median+growth_sds$bhi_sd, lwd=3, col='azure4')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$bhi_median-growth_sds$bhi_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$bhi_median-growth_sds$bhi_sd, lwd=3, col='azure4')

legend('topleft', legend=c('CDMM +Glucose +Amino acids','CDMM -Glucose +Amino acids','CDMM +Glucose -Amino acids','CDMM -Glucose -Amino acids','BHI'), 
       col=c('black','black','black','black','azure4'), pch=c(19,17,15,18,19), bg='white', cex=1.2, pt.cex=2, lwd=3)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

rm(plot_file, growth_sds, growth_medians)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

