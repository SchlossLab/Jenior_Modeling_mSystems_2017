



# Define variables
toxin_file <- '/Users/pschloss/Desktop/Repositories/Jenior_Transcriptomics_2015/data/raw/toxin_titer.dat'
cfu_file <- '/Users/pschloss/Desktop/Jenior_812/data/cfu.dat'
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
colnames(combined_toxin_spores) <- c('treatement', 'toxin', 'spores')

# Calculate strength of correlation
#cor.test(combined_toxin_spores$toxin, combined_toxin_spores$spores, method='spearman')
# S = 3579.6, p-value = 0.0006885, rho = 0.539309 

# Plot the combined data
plot(combined_toxin_spores$toxin, combined_toxin_spores$spores)
abline(lm(combined_toxin_spores$spores ~ combined_toxin_spores$toxin), lwd=4)

