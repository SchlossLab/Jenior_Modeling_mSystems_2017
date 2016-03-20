
# Define variables
toxin_file <- '/Users/pschloss/Desktop/Repositories/Jenior_Transcriptomics_2015/data/raw/toxin_titer.dat'
cfu_file <- '/Users/pschloss/Desktop/Repositories/Jenior_Transcriptomics_2015/data/raw/cfu.dat'
figure_file <- '/Users/pschloss/Desktop/toxin_spore_correlation.pdf'

# Read in data
toxin <- read.delim(toxin_file, sep='\t', header=T)
cfu <- read.delim(cfu_file, sep='\t', header=T)

# Format data seperately
toxin$mouse <- NULL
toxin$cage <- NULL
toxin <- droplevels(toxin)
toxin <- toxin[with(toxin, order(treatment)), ]

cfu <- subset(cfu, treatment != 'conventional')
cfu <- subset(cfu, cage < 4)
cfu$mouse <- NULL
cfu$cage <- NULL
cfu[cfu == 0] <- 100
cfu <- subset(cfu, type == 'spore')
cfu$type <- NULL
cfu <- droplevels(cfu)
cfu <- cfu[with(cfu, order(treatment)), ]

# Combine the datasets
combined_toxin_spores <- cbind(toxin, cfu[,2], deparse.level = 0)
colnames(combined_toxin_spores) <- c('treatment', 'toxin', 'spores')

# Calculate strength of correlation
glm_out <- glm(combined_toxin_spores$toxin ~ combined_toxin_spores$spores)
summary(glm_out)
# p = 3.91e-05 ***
# AIC: 26.423

# Plot the combined data
pdf(file=figure_file, width=8, height=7)
par(las=1, mar=c(4,4,1,1))
plot(jitter(combined_toxin_spores$toxin, factor=0.75), log10(combined_toxin_spores$spores), 
     xlim=c(1.75,3.25), ylim=c(1, 7), yaxt='n', pch=21, cex=1.8, cex.lab=1.2,
     xlab='Toxin Titer', ylab='CFU Spores per gram ceal content', 
     bg=c('firebrick2', 'blue2', 'chartreuse4', 'black')[combined_toxin_spores$treatment])
labelsY <- parse(text=paste(rep(10,7), '^', seq(1,7,1), sep=''))
axis(side=2, at=c(1:7), labelsY, tick=TRUE, cex.axis=1.2, las=1)
legend('bottomright', legend=c('Cefoperzone', 'Clindamycin', 'Streptomycin', 'Germfree'), 
       col=c('firebrick2', 'black', 'blue2', 'chartreuse4'), pch=15, cex=1.5, pt.cex=2.5, bty = "n")

#abline(glm_out, col="red")

dev.off()

