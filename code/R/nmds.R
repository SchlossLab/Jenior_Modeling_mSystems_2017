
# Define variables
nmds_file <- '/Users/schloss/Desktop/nmds.pdf'

# Load in data
abx.nmds.axes <- read.delim('/Users/schloss/Desktop/Jenior_Transcriptomics_2015/data/processed/otu/allabx.0.03.thetayc.0.03.lt.ave.nmds.axes', 
                            sep='\t', header=T, row.names=1)
abx.nmds.axes <- abx.nmds.axes[!rownames(abx.nmds.axes) %in% c('CefC5M2'), ] # Remove contaminated sample

metadata <- read.delim('/Users/schloss/Desktop/Jenior_Transcriptomics_2015/data/metadata.txt', 
                       sep='\t', header=T, row.names=1)

# Combine the metadata and axes
abx.nmds.axes <- merge(metadata, abx.nmds.axes, by='row.names')

# Select points corresponding to each treatment grops
conventional <- abx.nmds.axes[abx.nmds.axes$abx=='none', c(7,8)]
cefoperzone <- abx.nmds.axes[abx.nmds.axes$abx=='cefoperazone', c(7,8)]
clindamycin <- abx.nmds.axes[abx.nmds.axes$abx=='clindamycin', c(7,8)]
streptomycin <- abx.nmds.axes[abx.nmds.axes$abx=='streptomycin', c(7,8)]

# Plot the unlabeled data
pdf(file=nmds_file, width=7, height=6)
plot(abx.nmds.axes$axis1, abx.nmds.axes$axis2, pch='.', 
     xlim=c(-0.8,0.8), ylim=c(-0.8,0.8), 
     xlab='NMDS Axis 1', ylab='NMDS Axis 2')

# Color points by grouping
points(conventional, bg='black', pch=21, cex=1.4)
points(cefoperzone, bg='firebrick2', pch=21, cex=1.4)
points(clindamycin, bg='blue2', pch=21, cex=1.4)
points(streptomycin, bg='chartreuse4', pch=21, cex=1.4)

# Add legend
legend('bottomright', legend=c('Conventional', 'Cefoperzone', 'Clindamycin', 'Streptomycin'), 
       col=c('black', 'firebrick2', 'blue2', 'chartreuse4'), pch=15, pt.cex=1.8, bty = "n")
dev.off()
