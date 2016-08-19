
# Load colors
if ('wesanderson' %in% installed.packages()[,"Package"] == FALSE){
  install.packages(as.character('wesanderson'), quiet=TRUE);
} 
library('wesanderson', verbose=FALSE, character.only=TRUE)

# Define files
figure_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S1.pdf'
toxin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/toxin_titer.dat'
cfu_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/cfu.dat'

#---------------------------------------------------------------------------------------------------------------#


# Read in the data
toxin <- read.delim(toxin_file, sep='\t', header=T)
cfu <- read.delim(cfu_file, sep='\t', header=T)
cfu$treatment <- gsub('cefoperazone','Cefoperazone',cfu$treatment)
cfu$treatment <- gsub('streptomycin','Streptomycin',cfu$treatment)
cfu$treatment <- gsub('clindamycin','Clindamycin',cfu$treatment)
cfu$treatment <- gsub('conventional','Conventional',cfu$treatment)
cfu$treatment <- gsub('germfree','Germfree',cfu$treatment)

# Format datasets separately
toxin$mouse <- NULL
toxin$cage <- NULL
toxin <- subset(toxin, treatment != 'Conventional')
toxin <- droplevels(toxin)
toxin <- toxin[order(toxin$treatment), ]
cfu <- subset(cfu, treatment != 'Conventional')
cfu <- subset(cfu, cage < 4)
cfu$mouse <- NULL
cfu$cage <- NULL
cfu[cfu == 0] <- 100
cfu$cfu_vegetative <- NULL
cfu <- cfu[order(cfu$treatment), ]
color_palette <-c(rep(wes_palette("FantasticFox")[3],9),rep(wes_palette("FantasticFox")[5],9),rep('black',9),rep(wes_palette("FantasticFox")[1],9))

all_data <- as.data.frame(cbind(toxin$titer, cfu$cfu_spore, color_palette))
all_data <- cbind(toxin$treatment, all_data)
colnames(all_data) <- c('treatment', 'toxin', 'cfu', 'color')
all_data$treatment <- as.character(all_data$treatment)
all_data$toxin <- as.numeric(as.character(all_data$toxin))
all_data$cfu <- as.numeric(as.character(all_data$cfu))
all_data$color <- as.character(all_data$color)

#---------------------------------------------------------------------------------------------------------------#

# Set up the plot area
pdf(file=figure_file, width=21, height=7)
layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow=TRUE))

#---------------------------------------------------------------------------------------------------------------#

# A











mtext('A', side=2, line=2, las=2, adj=2, padj=-12, cex=2)

#---------------------------------------------------------------------------------------------------------------#

# B
# Remove gnotobiotic
part_B <- subset(all_data, treatment != 'Germfree')

# Calculate strength of correlation
cor.test(part_A$toxin, part_A$cfu, method='spearman')

# Plot it
par(las=1, mar=c(4,5,1,1))
plot(jitter(part_B$toxin, factor=0.75), log10(part_B$cfu), 
     xlim=c(1.75,3.25), ylim=c(1, 7), yaxt='n', pch=21, cex=1.8, cex.lab=1.2,
     xlab='Toxin Titer (Log10)', ylab='CFU Spores per gram ceal content', bg=part_B$color)
labelsY <- parse(text=paste(rep(10,7), '^', seq(1,7,1), sep=''))
axis(side=2, at=c(1:7), labelsY, tick=TRUE, cex.axis=1.2, las=1)
legend('bottomright', legend=c('Streptomycin', 'Cefoperzone', 'Clindamycin'), 
       pt.bg=c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5]), 
       pch=21, cex=1.5, pt.cex=2.5, bty='n')
legend('topleft', legend=c('rho = 0.07755254', 'p-value = 0.7006'), bty='n', cex=1.5)
abline(lm(log10(part_B$cfu)~part_A$toxin))
mtext('B', side=2, line=2, las=2, adj=2, padj=-12, cex=2)

#---------------------------------------------------------------------------------------------------------------#

# B
# Calculate strength of correlation
cor.test(all_data$toxin, all_data$cfu, method='spearman')

# Plot it
par(las=1, mar=c(4,5,1,1))
plot(jitter(all_data$toxin, factor=0.75), log10(all_data$cfu), 
     xlim=c(1.75,3.25), ylim=c(1, 7), yaxt='n', pch=21, cex=1.8, cex.lab=1.2,
     xlab='Toxin Titer (Log10)', ylab='CFU Spores per gram ceal content', bg=all_data$color)
labelsY <- parse(text=paste(rep(10,7), '^', seq(1,7,1), sep=''))
axis(side=2, at=c(1:7), labelsY, tick=TRUE, cex.axis=1.2, las=1)
legend('bottomright', legend=c('Streptomycin', 'Cefoperzone', 'Clindamycin', 'Gnotobiotic'), 
       pt.bg=c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'black'), 
       pch=21, cex=1.5, pt.cex=2.5, bty='n')
legend('topleft', legend=c('rho = 0.539309 ', 'p-value = 0.0006885'), bty='n', cex=1.5)
abline(lm(log10(all_data$cfu)~all_data$toxin))
mtext('C', side=2, line=2, las=2, adj=2, padj=-12, cex=2)

#---------------------------------------------------------------------------------------------------------------#

dev.off()

