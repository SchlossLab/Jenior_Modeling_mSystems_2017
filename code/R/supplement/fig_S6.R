
# Load dependencies
deps <- c('vegan', 'wesanderson', 'matrixStats', 'flux');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in growth rate data
# Define variables
growth_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/formatted_growth.tsv'

# Read in data
growth <- read.delim(growth_file, sep='\t', header=TRUE, row.names=1)
growth <- as.data.frame(t(growth))
rm(growth_file)

# Seperate to groups of each growth substrate
sorbitol <- cbind(growth$B9, growth$B10, growth$B11)
galactitol <- cbind(growth$C9, growth$C10, growth$C11)
starch <- cbind(growth$D9, growth$D10, growth$D11)
fructose <- cbind(growth$E9, growth$E10, growth$E11)
combination <- cbind(growth$G3, growth$G4, growth$G5)
mannitol <- cbind(growth$F9, growth$F10, growth$F11)
salicin <- cbind(growth$G9, growth$G10, growth$G11)
y_glucose_y_aa <- cbind(growth$B3, growth$B4, growth$B5)
n_glucose_y_aa <- cbind(growth$D3, growth$D4, growth$D5)
y_glucose_n_aa <- cbind(growth$C3, growth$C4, growth$C5)
n_glucose_n_aa <- cbind(growth$E3, growth$E4, growth$E5)
bhi <- cbind(growth$F3, growth$F4, growth$F5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Medians of treatement groups and subtract blanks
sorbitol_median <- rowMedians(sorbitol, na.rm=TRUE) - growth$B8
sorbitol_median[sorbitol_median < 0] <- 0
galactitol_median <- rowMedians(galactitol, na.rm=TRUE) - growth$C8
galactitol_median[galactitol_median < 0] <- 0
starch_median <- rowMedians(starch, na.rm=TRUE) - growth$D8
starch_median[starch_median < 0] <- 0
fructose_median <-  rowMedians(fructose, na.rm=TRUE) - growth$E8
fructose_median[fructose_median < 0] <- 0
combination_median <- rowMedians(combination, na.rm=TRUE) - growth$G2
combination_median[combination_median < 0] <- 0
mannitol_median <- rowMedians(mannitol, na.rm=TRUE) - growth$F8
mannitol_median[mannitol_median < 0] <- 0
salicin_median <- rowMedians(salicin, na.rm=TRUE) - growth$G7
salicin_median[salicin_median < 0] <- 0
y_glucose_y_aa_median <- rowMedians(y_glucose_y_aa, na.rm=TRUE) - growth$B2
y_glucose_y_aa_median[y_glucose_y_aa_median < 0] <- 0
n_glucose_y_aa_median <- rowMedians(n_glucose_y_aa, na.rm=TRUE) - growth$D2
n_glucose_y_aa_median[n_glucose_y_aa_median < 0] <- 0
y_glucose_n_aa_median <- rowMedians(y_glucose_n_aa, na.rm=TRUE) - growth$C2
y_glucose_n_aa_median[y_glucose_n_aa_median < 0] <- 0
n_glucose_n_aa_median <- rowMedians(n_glucose_n_aa, na.rm=TRUE) - growth$E2
n_glucose_n_aa_median[n_glucose_n_aa_median < 0] <- 0
bhi_median <- rowMedians(bhi, na.rm=TRUE) - growth$F2
bhi_median[bhi_median < 0] <- 0
growth_medians <- as.data.frame(cbind(sorbitol_median, galactitol_median, starch_median, fructose_median, combination_median, mannitol_median, salicin_median, y_glucose_y_aa_median, n_glucose_y_aa_median, y_glucose_n_aa_median, n_glucose_n_aa_median, bhi_median))

# Standard deviations
sorbitol_sd <- rowSds(sorbitol, na.rm=TRUE)
galactitol_sd <- rowSds(galactitol, na.rm=TRUE)
starch_sd <-  rowSds(starch, na.rm=TRUE)
fructose_sd <-  rowSds(fructose, na.rm=TRUE)
combination_sd <-  rowSds(combination, na.rm=TRUE)
mannitol_sd <- rowSds(mannitol, na.rm=TRUE)
salicin_sd <- rowSds(salicin, na.rm=TRUE)
y_glucose_y_aa_sd <- rowSds(y_glucose_y_aa, na.rm=TRUE)
n_glucose_y_aa_sd <- rowSds(n_glucose_y_aa, na.rm=TRUE)
y_glucose_n_aa_sd <- rowSds(y_glucose_n_aa, na.rm=TRUE)
n_glucose_n_aa_sd <- rowSds(n_glucose_n_aa, na.rm=TRUE)
bhi_sd <- rowSds(bhi, na.rm=TRUE)
growth_sds <- as.data.frame(cbind(sorbitol_sd, galactitol_sd, starch_sd, fructose_sd, combination_sd, mannitol_sd, salicin_sd, y_glucose_y_aa_sd, n_glucose_y_aa_sd, y_glucose_n_aa_sd, n_glucose_n_aa_sd, bhi_sd))

# Clean up
rm(sorbitol_median, galactitol_median, starch_median, fructose_median, combination_median, mannitol_median, salicin_median, y_glucose_y_aa_median, n_glucose_y_aa_median, y_glucose_n_aa_median, n_glucose_n_aa_median, bhi_median)
rm(sorbitol, galactitol, starch, fructose, combination, mannitol, salicin, y_glucose_y_aa, n_glucose_y_aa, y_glucose_n_aa, n_glucose_n_aa, bhi)
rm(sorbitol_sd, galactitol_sd, starch_sd, fructose_sd, combination_sd, mannitol_sd, salicin_sd, y_glucose_y_aa_sd, n_glucose_y_aa_sd, y_glucose_n_aa_sd, n_glucose_n_aa_sd, bhi_sd)
rm(growth)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plot file
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S6.pdf'
pdf(file=plot_file, width=35, height=14)
layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Plot the intersting 24 curves
par(mar=c(5,6,1,1), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE)
plot(0, type='n', xaxt='n', xlim=c(0,50), ylim=c(-0.03,1), lwd=2, pch=19, xlab='Hours Postinoculation', ylab=expression(OD[600]), cex=2.3)
axis(1, at=seq(1,50,4), labels=seq(0,24,2))
abline(h=seq(0,1,0.1), lty=3, col='gray68')
abline(v=seq(1,50,2), lty=3, col='gray68')

lines(growth_medians$y_glucose_y_aa_median, type='o', lwd=2, pch=19, cex=2, col='black')
segments(x0=seq(1,49,1), y0=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, x1=seq(1,49,1), y1=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, lwd=2.5, cex=2, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, lwd=2, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, lwd=2, col='black')

lines(growth_medians$bhi_median, type='o', lwd=2, pch=19, cex=2, col=wes_palette('Darjeeling2')[2])
segments(x0=seq(1,49,1), y0=growth_medians$bhi_median+growth_sds$bhi_sd, x1=seq(1,49,1), y1=growth_medians$bhi_median-growth_sds$bhi_sd, lwd=2.5, cex=2, col=wes_palette('Darjeeling2')[2])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$bhi_median+growth_sds$bhi_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$bhi_median+growth_sds$bhi_sd, lwd=2, col=wes_palette('Darjeeling2')[2])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$bhi_median-growth_sds$bhi_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$bhi_median-growth_sds$bhi_sd, lwd=2, col=wes_palette('Darjeeling2')[2])

lines(growth_medians$salicin_median, type='o', col=wes_palette('FantasticFox')[5], lwd=2.5, pch=19, cex=2)
segments(x0=seq(1,49,1), y0=growth_medians$salicin_median+growth_sds$salicin_sd, x1=seq(1,49,1), y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$salicin_median+growth_sds$salicin_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$salicin_median+growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$salicin_median-growth_sds$salicin_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])

legend('topleft', legend=c('+Glucose +AA','Salicin','BHI'), pch=19, cex=2.5, pt.cex=3.5, bg='white',
       col=c('black',wes_palette('FantasticFox')[5],wes_palette('Darjeeling2')[2]))

mtext('A', side=2, line=2, las=2, adj=2, padj=-26.5, cex=2)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Plot the stationary curves
par(mar=c(5,6,1,1), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE)
plot(0, type='n', xaxt='n', xlim=c(0,50), ylim=c(-0.03,1), lwd=2, pch=19, xlab='Hours Postinoculation', ylab=expression(OD[600]), cex=2.3)
axis(1, at=seq(1,50,4), labels=seq(0,24,2))
abline(h=seq(0,1,0.1), lty=3, col='gray68')
abline(v=seq(1,50,2), lty=3, col='gray68')

lines(growth_medians$n_glucose_y_aa_median, type='o', lwd=2, pch=19, cex=2, col='black')
segments(x0=seq(1,49,1), y0=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, x1=seq(1,49,1), y1=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, lwd=2.5, cex=2, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, lwd=2, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, lwd=2, col='black')

lines(growth_medians$y_glucose_n_aa_median, type='o', lwd=2, pch=19, cex=2, col='black')
segments(x0=seq(1,49,1), y0=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, x1=seq(1,49,1), y1=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, lwd=2.5, cex=2, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, lwd=2, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, lwd=2, col='black')

lines(growth_medians$n_glucose_n_aa_median, type='o', lwd=2, pch=19, cex=2, col='black')
segments(x0=seq(1,49,1), y0=growth_medians$n_glucose_n_aa_median+growth_sds$n_glucose_n_aa_sd, x1=seq(1,49,1), y1=growth_medians$n_glucose_n_aa_median-growth_sds$n_glucose_n_aa_sd, lwd=2.5, cex=2, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$n_glucose_n_aa_median+growth_sds$n_glucose_n_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$n_glucose_n_aa_median+growth_sds$n_glucose_n_aa_sd, lwd=2, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$n_glucose_n_aa_median-growth_sds$n_glucose_n_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$n_glucose_n_aa_median-growth_sds$n_glucose_n_aa_sd, lwd=2, col='black')

lines(growth_medians$sorbitol_median, type='o', lwd=2, pch=19, cex=2.7, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1), y0=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, x1=seq(1,49,1), y1=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, lwd=2.5, cex=2, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, lwd=2, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, lwd=2, col=wes_palette('FantasticFox')[1])

lines(growth_medians$galactitol_median, type='o', lwd=2, pch=19, cex=2.7, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1), y0=growth_medians$galactitol_median+growth_sds$galactitol_sd, x1=seq(1,49,1), y1=growth_medians$galactitol_median-growth_sds$galactitol_sd, lwd=2.5, cex=2, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$galactitol_median+growth_sds$galactitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$galactitol_median+growth_sds$galactitol_sd, lwd=2, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$galactitol_median-growth_sds$galactitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$galactitol_median-growth_sds$galactitol_sd, lwd=2, col=wes_palette('FantasticFox')[1])

lines(growth_medians$starch_median, type='o', lwd=2, pch=19, cex=2.7, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1), y0=growth_medians$starch_median+growth_sds$starch_sd, x1=seq(1,49,1), y1=growth_medians$starch_median-growth_sds$starch_sd, lwd=2.5, cex=2, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$starch_median+growth_sds$starch_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$starch_median+growth_sds$starch_sd, lwd=2, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$starch_median-growth_sds$starch_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$starch_median-growth_sds$starch_sd, lwd=2, col=wes_palette('FantasticFox')[1])

lines(growth_medians$fructose_median, type='o', lwd=2, pch=19, cex=2.7, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1), y0=growth_medians$fructose_median+growth_sds$fructose_sd, x1=seq(1,49,1), y1=growth_medians$fructose_median-growth_sds$fructose_sd, lwd=2.5, cex=2, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$fructose_median+growth_sds$fructose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$fructose_median+growth_sds$fructose_sd, lwd=2, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$fructose_median-growth_sds$fructose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$fructose_median-growth_sds$fructose_sd, lwd=2, col=wes_palette('FantasticFox')[1])

lines(growth_medians$mannitol_median, type='o', lwd=2, pch=19, cex=2.7, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,49,1), y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, x1=seq(1,49,1), y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=2.5, cex=2, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$mannitol_median+growth_sds$mannitol_sd, lwd=2, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$mannitol_median-growth_sds$mannitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=2, col=wes_palette('FantasticFox')[3])

legend('topleft', legend=c('-Glucose +AA','+Glucose -AA','-Glucose -AA','D-Sorbitol','Galactitol','Starch','D-Fructose','Mannitol'), 
       col=c('black','black','black',wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3]), 
       pch=c(16,17,18,15,16,18,17,19), cex=2.5, pt.cex=3.5, bg='white')

mtext('B', side=2, line=2, las=2, adj=2, padj=-26.5, cex=2)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Plot the stationary curves
par(mar=c(5,6,1,1), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE)
plot(0, type='n', xaxt='n', xlim=c(0,50), ylim=c(-0.03,1), lwd=2, pch=19, xlab='Hours Postinoculation', ylab=expression(OD[600]), cex=2.3)
axis(1, at=seq(1,50,4), labels=seq(0,24,2))
abline(h=seq(0,1,0.1), lty=3, col='gray68')
abline(v=seq(1,50,2), lty=3, col='gray68')

lines(growth_medians$combination_median, type='o', lwd=2, pch=19, cex=2.7, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1), y0=growth_medians$combination_median+growth_sds$combination_sd, x1=seq(1,49,1), y1=growth_medians$combination_median-growth_sds$combination_sd, lwd=2.5, cex=2, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$combination_median+growth_sds$combination_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$combination_median+growth_sds$combination_sd, lwd=2, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$combination_median-growth_sds$combination_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$combination_median-growth_sds$combination_sd, lwd=2, col=wes_palette('FantasticFox')[1])

legend('topleft', legend='Equal concentrations of D-Sorbitol, Galactitol, Starch, and D-Fructose', pch=19, cex=2.5, pt.cex=3.5, bg='white', col=wes_palette('FantasticFox')[1])

mtext('C', side=2, line=2, las=2, adj=2, padj=-26.5, cex=2)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
rm(plot_file, growth_sds, growth_medians)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

