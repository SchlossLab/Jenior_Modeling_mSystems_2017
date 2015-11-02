

# Read in data
cfu_sporeVtoxin <- read.delim('', sep='\t', header=T)
cfu_sporeVtoxin$Spores_log <- log10(cfu_sporeVtoxin$Spores)


# Denote which points correspond to each treament
clindamycin <- cfu_sporeVtoxin[cfu_sporeVtoxin$Treatment=='clindamycin',c(2,4)]
streptomycin <- cfu_sporeVtoxin[cfu_sporeVtoxin$Treatment=='streptomycin',c(2,4)]
cefoperazone <- cfu_sporeVtoxin[cfu_sporeVtoxin$Treatment=='cefoperazone',c(2,4)]
germfree <- cfu_sporeVtoxin[cfu_sporeVtoxin$Treatment=='germfree',c(2,4)]


# Calculate strength of correlation
cor.test(cfu_sporeVtoxin$Toxin, cfu_sporeVtoxin$Spores, method='spearman')  # S = 3579.569, p-value = 0.0006885, rho = 0.539309 
# 0.539309 ** 2 = 0.2908542


# Plot s
par(mar = rep(5, 4) , las = 1)
plot(cfu_sporeVtoxin$Toxin, cfu_sporeVtoxin$Spores_log, pch=16, xlab=substitute(paste(italic('C. diffile'), ' Toxin Titer')), ylim=c(1,7), yaxt='n', ylab=substitute(paste(italic('C. diffile'), ' Spores per gram cecal content')))
labelsY=parse(text=paste(rep(10,7), "^", seq(1,7,1), sep=""))
axis(side = 2, at = 1.301030:7.301030, labelsY, tick = TRUE)
abline(lm(log10(cfu_sporeVtoxin$Spores) ~ cfu_sporeVtoxin$Toxin), col="black", lwd=2) # fit line
#abline(h=2.301030, col='black', lty=2) # Limit of detection 1
#abline(v=2, col='black', lty=2) # Limit of detection 2


# Color points on plot
points(clindamycin, col='black', pch=16)
points(streptomycin, col='green', pch=16)
points(cefoperazone, col='red', pch=16)
points(germfree, col='blue', pch=16)


# Add legend and other correlation metrics
legend(x='bottomright',legend=c('Clindamyin', 'Streptomycin', 'Cefoperzone', 'Germfree'), col=c('black', 'green', 'red', 'blue'), pch=16)
text(2.915, 2.55, substitute(paste(italic('p value'), ' < 0.001')), cex=0.8)
text(2.857, 2.35, parse(text=paste('R', '^', '2', sep="")), cex=0.8)
text(2.905, 2.35, '= 0.291', cex=0.8)
