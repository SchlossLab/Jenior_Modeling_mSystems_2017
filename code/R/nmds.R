
# Define variables
nmds_file <- '/Users/pschloss/Desktop/Repositories/Jenior_Transcriptomics_2015/data/processed/otu/allabx.0.03.thetayc.0.03.lt.ave.nmds.axes'
metadata_file <- '/Users/pschloss/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.txt'
figure_file <- '/Users/pschloss/Desktop/nmds.pdf'

# Load in data
abx.nmds.axes <- read.delim(nmds_file, sep='\t', header=T, row.names=1)
abx.nmds.axes <- abx.nmds.axes[!rownames(abx.nmds.axes) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)

# Combine the metadata and axes
abx.nmds.axes <- merge(metadata, abx.nmds.axes, by='row.names')

# Remove uninfected cages
pick.abx.nmds.axes <- subset(abx.nmds.axes, treatment != 630)

pdf(file=nmds_file, width=7, height=6)

# Plot the unlabeled data
par(las=1, mar=c(3,3,1,1), mgp=c(2,0.8,0))
plot(pick.abx.nmds.axes$axis1, pick.abx.nmds.axes$axis2, pch=21, cex=2,
     bg=c('black', 'firebrick2', 'blue2', 'chartreuse4')[pick.abx.nmds.axes$abx],
     xlim=c(-0.8,0.8), ylim=c(-0.8,0.8), 
     xlab='NMDS Axis 1', ylab='NMDS Axis 2')

# Add legend
legend('bottomright', legend=c('Conventional', 'Cefoperzone', 'Clindamycin', 'Streptomycin'), 
       col=c('blue2', 'black', 'firebrick2', 'chartreuse4'), pch=15, pt.cex=1.8, bty = "n")

# Plot centroids
points(mean(subset(pick.abx.nmds.axes, abx == 'none')$axis1), mean(subset(pick.abx.nmds.axes, abx == 'none')$axis2), bg='blue2', pch=21, cex=4)
points(mean(subset(pick.abx.nmds.axes, abx == 'cefoperazone')$axis1), mean(subset(pick.abx.nmds.axes, abx == 'cefoperazone')$axis2), bg='black', pch=21, cex=4)
points(mean(subset(pick.abx.nmds.axes, abx == 'clindamycin')$axis1), mean(subset(pick.abx.nmds.axes, abx == 'clindamycin')$axis2), bg='firebrick2', pch=21, cex=4)
points(mean(subset(pick.abx.nmds.axes, abx == 'streptomycin')$axis1), mean(subset(pick.abx.nmds.axes, abx == 'streptomycin')$axis2), bg='chartreuse4', pch=21, cex=4)







#dev.off()