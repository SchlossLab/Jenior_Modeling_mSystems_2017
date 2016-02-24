
# Creates plots for CFU comparisons

# Define variable
cfu_file <- '/Users/pschloss/Desktop/Jenior_812/data/cfu.dat'
vegetative_file <- '/Users/pschloss/Desktop/vegetative.cfu.pdf'
spore_file <- '/Users/pschloss/Desktop/spore.cfu.pdf'
percent_file <- '/Users/pschloss/Desktop/percent.cfu.pdf'

# Read, format, and separate the data
cfu_raw <- read.delim(cfu_file, sep='\t', header=T)

cfu <- cfu_raw
cfu[cfu == 0] <- 100
cfu$cfu <- log10(cfu$cfu)
cfu$mouse <- NULL
cfu <- subset(cfu, cage < 4)
cfu$cage <- NULL
cfu <- subset(cfu, treatment != 'conventional')

cfu_raw[cfu_raw == 0] <- 100
cfu_raw$mouse <- NULL
cfu_raw <- subset(cfu_raw, cage < 4)
cfu_raw$cage <- NULL
cfu_raw <- subset(cfu_raw, treatment != 'conventional')

vegetative_cfu <- subset(cfu, type == 'vegetative')
vegetative_cfu$type <- NULL
vegetative_cfu$treatment <- factor(vegetative_cfu$treatment, levels = c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree'))
vegetative_cfu <- droplevels(vegetative_cfu)
spore_cfu <- subset(cfu, type == 'spore')
spore_cfu$type <- NULL
spore_cfu$treatment <- factor(spore_cfu$treatment, levels = c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree'))
spore_cfu <- droplevels(spore_cfu)

vegetative_cfu_raw <- subset(cfu_raw, type == 'vegetative')
vegetative_cfu_raw$type <- NULL
vegetative_cfu_raw$treatment <- factor(vegetative_cfu_raw$treatment, levels = c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree'))
vegetative_cfu_raw <- droplevels(vegetative_cfu_raw)
spore_cfu_raw <- subset(cfu_raw, type == 'spore')
spore_cfu_raw$type <- NULL
spore_cfu_raw$treatment <- factor(spore_cfu_raw$treatment, levels = c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree'))
spore_cfu_raw <- droplevels(spore_cfu_raw)


cef_iqr_vege <- quantile(vegetative_cfu[vegetative_cfu$treatment == 'Cefoperazone', 2])
clinda_iqr_vege <- quantile(vegetative_cfu[vegetative_cfu$treatment == 'Clindamycin', 2])
strep_iqr_vege <- quantile(vegetative_cfu[vegetative_cfu$treatment == 'Streptomycin', 2])
gf_iqr_vege <- quantile(vegetative_cfu[vegetative_cfu$treatment == 'Germfree', 2])

cef_iqr_spore <- quantile(spore_cfu[spore_cfu$treatment == 'Cefoperazone', 2])
clinda_iqr_spore <- quantile(spore_cfu[spore_cfu$treatment == 'Clindamycin', 2])
strep_iqr_spore <- quantile(spore_cfu[spore_cfu$treatment == 'Streptomycin', 2])
gf_iqr_spore <- quantile(spore_cfu[spore_cfu$treatment == 'Germfree', 2])

#-------------------------------------------------------------------------------------------------------------------------------------------------------#

# Plot formatted data - vegetative
pdf(file=vegetative_file, width=10, height=7)
par(las=1, mar=c(3.5,4,1,1), mgp=c(2.5,0.7,0))
stripchart(cfu~treatment, data=vegetative_cfu, bg='firebrick', 
        ylim=c(5,9), xaxt='n', yaxt='n', pch=21, vertical=T, method='jitter',
        jitter=0.2, cex=3, ylab='Vegetative CFU/g Cecal Content', cex.lab=1.5)
axis(side=1, at=c(1:4), c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree'), tick = FALSE, font=2, cex.axis=1.4)
mtext(c('0.5 mg/ml DW', '0.5 mg/ml DW', '10 mg/kg IP', ''), side=1, at=c(1:4), cex=0.9, padj=3.8)
labelsY <- parse(text=paste(rep(10,5), '^', seq(5,9,1), sep=''))
axis(side=2, at=c(5:9), labelsY, tick=TRUE, cex.axis=1.2, las=1)

# Plot medians
segments(0.6, cef_iqr_vege[3], 1.4, cef_iqr_vege[3], lwd=8) # cefoperazone
segments(1.6, strep_iqr_vege[3], 2.4, strep_iqr_vege[3], lwd=8) # streptomycin
segments(2.6, clinda_iqr_vege[3], 3.4, clinda_iqr_vege[3], lwd=8) # clindamycin
segments(3.6, gf_iqr_vege[3], 4.4, gf_iqr_vege[3], lwd=8) # germfree
dev.off()


# Plot formatted data - spores
pdf(file=spore_file, width=10, height=7)
par(las=1, mar=c(3.5,4,1,1), mgp=c(2.5,0.7,0))
stripchart(cfu~treatment, data=spore_cfu, bg='blue2', 
           ylim=c(1,7), xaxt='n', yaxt='n', pch=21, vertical=T, method='jitter',
           jitter=0.2, cex=3, ylab='Spore CFU/g Cecal Content', cex.lab=1.5)
axis(side=1, at=c(1:4), c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree'), tick = FALSE, font=2, cex.axis=1.4)
mtext(c('0.5 mg/ml DW', '0.5 mg/ml DW', '10 mg/kg IP', ''), side=1, at=c(1:4), cex=0.9, padj=3.8)
labelsY <- parse(text=paste(rep(10,7), '^', seq(1,7,1), sep=''))
axis(side=2, at=c(1:7), labelsY, tick=TRUE, cex.axis=1.2)

# Plot medians
segments(0.6, cef_iqr_spore[3], 1.4, cef_iqr_spore[3], lwd=8) # cefoperazone
segments(1.6, strep_iqr_spore[3], 2.4, strep_iqr_spore[3], lwd=8) # streptomycin
segments(2.6, clinda_iqr_spore[3], 3.4, clinda_iqr_spore[3], lwd=8) # clindamycin
segments(3.6, gf_iqr_spore[3], 4.4, gf_iqr_spore[3], lwd=8) # germfree

# Plot limit of detection and significance
abline(h=2, col="black", lty=2, lwd=4)
text(4, 7, '***', cex=3, font=2)
dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------------------------#

# Calculate differences - vegetative
cef_vegetative <- subset(vegetative_cfu_raw, treatment=='Cefoperazone')[,2]
strep_vegetative <- subset(vegetative_cfu_raw, treatment=='Streptomycin')[,2]
clinda_vegetative <- subset(vegetative_cfu_raw, treatment=='Clindamycin')[,2]
gf_vegetative <- subset(vegetative_cfu_raw, treatment=='Germfree')[,2]

shapiro.test(cef_vegetative) # p-value = 0.4805
shapiro.test(strep_vegetative) # p-value = 0.07191
shapiro.test(clinda_vegetative) # p-value = 0.3428
shapiro.test(gf_vegetative) # p-value = 0.0768

wilcox.test(cef_vegetative, strep_vegetative, exact=F) # p-value = 0.563
wilcox.test(cef_vegetative, clinda_vegetative, exact=F) # p-value = 0.2868
wilcox.test(cef_vegetative, gf_vegetative, exact=F) # p-value = 0.01197  *
wilcox.test(strep_vegetative, clinda_vegetative, exact=F) # p-value = 0.6908
wilcox.test(strep_vegetative, gf_vegetative, exact=F) # p-value = 0.5644
wilcox.test(clinda_vegetative, gf_vegetative, exact=F) # p-value = 0.0035  **


# Calculate differences - spore
cef_spore <- subset(spore_cfu_raw, treatment=='Cefoperazone')[,2]
strep_spore <- subset(spore_cfu_raw, treatment=='Streptomycin')[,2]
clinda_spore <- subset(spore_cfu_raw, treatment=='Clindamycin')[,2]
gf_spore <- subset(spore_cfu_raw, treatment=='Germfree')[,2]

shapiro.test(cef_spore) # p-value = 0.001587
shapiro.test(strep_spore) # p-value = 0.002147
shapiro.test(clinda_spore) # p-value = 0.008843
shapiro.test(gf_spore) # p-value = 0.5463

wilcox.test(cef_spore, strep_spore, exact=F) # p-value = 0.5617
wilcox.test(cef_spore, clinda_spore, exact=F) # p-value = 0.4241
wilcox.test(cef_spore, gf_spore, exact=F) # p-value = 0.001892  ***
wilcox.test(strep_spore, clinda_spore, exact=F) # p-value = 1
wilcox.test(strep_spore, gf_spore, exact=F) # p-value = 0.001086  ***
wilcox.test(clinda_spore, gf_spore, exact=F) # p-value = 0.0005699  ***

#-------------------------------------------------------------------------------------------------------------------------------------------------------#


# Plotting spores as a percentage of total CFU
percentage <- (spore_cfu_raw$cfu / vegetative_cfu_raw$cfu) * 100
treatment <- as.vector(vegetative_cfu_raw$treatment)
cfu_percent <- as.data.frame(cbind(treatment, percentage))
cfu_percent$percentage <- as.numeric(as.character(cfu_percent$percentage))
cfu_percent$treatment <- factor(cfu_percent$treatment, levels = c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree'))


# Plot the percentages
pdf(file=percent_file, width=10, height=7)
par(las=1, mar=c(4,3.5,1,1), mgp=c(2.5,0.7,0))
stripchart(percentage~treatment, data=cfu_percent, vertical=T, method='jitter', 
           jitter=.1, pch=16, ylim=c(0,8), yaxt='n', xaxt='n', cex=1.4, 
           ylab=substitute(paste('Percent Spores of Cultureable ', italic('C. diffile'))))
axis(side=1, at=c(1:4), c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Germfree'), 
     tick = FALSE, font=2, cex.axis=1.4)
mtext(c('0.5 mg/ml DW', '0.5 mg/ml DW', '10 mg/kg IP', ''), side=1, at=c(1:4), 
      cex=0.9, padj=3.5)
labelsY <- c('0%', '', '2%', '', '4%', '', '6%', '', '8%')
axis(side=2, at=0:8, labelsY, tick=TRUE, cex.axis=1.2)
abline(h=0, col='black', lty=2)
abline(h=2, col='black', lty=2)
abline(h=4, col='black', lty=2)
abline(h=6, col='black', lty=2)
abline(h=8, col='black', lty=2)
text(4, 4.2, '***', cex=2, font=2) 

# Plot medians
cef_median <- median(subset(cfu_percent, treatment=='Cefoperazone')[,2])
strep_median <- median(subset(cfu_percent, treatment=='Clindamycin')[,2])
clinda_median <- median(subset(cfu_percent, treatment=='Streptomycin')[,2])
gf_median <- median(subset(cfu_percent, treatment=='Germfree')[,2])

segments(x0=0.7, x1=1.3, y0=cef_median, y1=cef_median, col="black", lwd=3)
segments(x0=1.7, x1=2.3, y0=strep_median, y1=strep_median, col="black", lwd=3)
segments(x0=2.7, x1=3.3, y0=clinda_median, y1=clinda_median, col="black", lwd=3)
segments(x0=3.7, x1=4.3, y0=gf_median, y1=gf_median, col="black", lwd=3)


# Significance
text(4, 4.2, '***', cex=2, font=2)
dev.off()

# Calculate differences
wilcox.test(subset(cfu_percent, treatment=='Cefoperazone')[,2], 
            subset(cfu_percent, treatment=='Clindamycin')[,2], exact=F) # p-value = 0.2891
wilcox.test(subset(cfu_percent, treatment=='Cefoperazone')[,2], 
            subset(cfu_percent, treatment=='Streptomycin')[,2], exact=F) # p-value = 0.4011
wilcox.test(subset(cfu_percent, treatment=='Cefoperazone')[,2], 
            subset(cfu_percent, treatment=='Germfree')[,2], exact=F) # p-value = 0.01337  *
wilcox.test(subset(cfu_percent, treatment=='Clindamycin')[,2], 
            subset(cfu_percent, treatment=='Streptomycin')[,2], exact=F) # p-value = 0.4268
wilcox.test(subset(cfu_percent, treatment=='Clindamycin')[,2], 
            subset(cfu_percent, treatment=='Germfree')[,2], exact=F) # p-value = 0.004718  ***
wilcox.test(subset(cfu_percent, treatment=='Streptomycin')[,2], 
            subset(cfu_percent, treatment=='Germfree')[,2], exact=F) # p-value = 0.008071  ***


