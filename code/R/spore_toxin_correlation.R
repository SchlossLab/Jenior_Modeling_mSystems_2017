

# Read in data
cfu <- read.delim('~/Desktop/Jenior_Transcriptomics_2015/data/raw/cfu.dat', sep='\t', header=T)
cfu <- subset(cfu, (type=='spore') & (treatment!='conventional') & (mouse<10))
toxin <- read.delim('~/Desktop/Jenior_Transcriptomics_2015/data/raw/toxin_titer.dat', sep='\t', header=T)
cfu$toxin_titer <- c(toxin$Cefoperazone, toxin$Clindamycin, toxin$Streptomycin, toxin$Germfree)

# Calculate strength of correlation
cor.test(cfu$cfu, cfu$toxin_titer, method='spearman')  
# data:  cfu$cfu and cfu$toxin_titer
# S = 3579.6, p-value = 0.0006885
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
# rho 
# 0.539309 

