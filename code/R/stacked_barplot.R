
# Define variables
figure_file <- '/Users/schloss/Desktop/family_barplot.pdf'
legend_file <- '/Users/schloss/Desktop/family_barplot.legend.pdf'
  
# Load in and format data
metadata <- read.delim('/Users/schloss/Desktop/Jenior_Transcriptomics_2015/data/metadata.txt', 
                       sep='\t', header=T, row.names=1)
taxonomy <- read.delim('/Users/schloss/Desktop/Jenior_Transcriptomics_2015/data/processed/phylotype/allabx.2.cons.family.format.taxonomy', 
                       sep='\t', header=T, row.names=1)
shared <- read.delim('/Users/schloss/Desktop/Jenior_Transcriptomics_2015/data/processed/phylotype/allabx.2.subsample.shared', 
                     sep='\t', header=T, row.names=2)
shared$label <- NULL
shared$numOtus <- NULL
shared <- shared[!rownames(shared) %in% c('CefC5M2'), ] # Remove contaminated sample
rel_abund <- (shared / rowSums(shared)) * 100

# Rename each OTU by it's corresponding taxa
name_shared <- merge(taxonomy, t(rel_abund), by='row.names')
rownames(name_shared) <- name_shared$Taxonomy
name_shared$Taxonomy <- NULL
name_shared$Size <- NULL
name_shared$Size <- NULL
name_shared$Row.names <- NULL
name_shared <- t(name_shared)

# Label with metadata and filter shared file
labeled_shared <- merge(metadata, name_shared, by='row.names')
filtered_shared <- subset(labeled_shared, labeled_shared$treatment != '630')
rownames(filtered_shared) <- filtered_shared$Row.names
filtered_shared$Row.names <- NULL

# Remove unused features
filtered_shared$cage <- NULL
filtered_shared$mouse <- NULL
filtered_shared$gender <- NULL
filtered_shared$treatment <- NULL

# Find the mean of OTUs with each treatment group
final_shared <- aggregate(.~abx, data=filtered_shared, mean)
rownames(final_shared) <- c('Cefoperazone', 'Clindamycin', 'Conventional', 'Streptomycin')

# Reorder rows to make more sense for plotting
final_shared <- final_shared[match(c('none', 'cefoperazone', 'clindamycin', 'streptomycin'), final_shared$abx),]
final_shared$abx <- NULL

# Bin lowly abundant OTUs into an 'Other' category
final_shared[final_shared < 1] <- 0
final_shared <- final_shared[, colSums(final_shared != 0) > 0]
final_shared$Other <- 100 - rowSums(final_shared)

# Plot the final formatted table
pdf(file=figure_file, width=12, height=9)
par(las=1, mar=c(2.2,3,1,1), mgp=c(2,0.8,0), cex=2)
colors <- rainbow(ncol(final_shared))
barplot(t(final_shared), col=colors, yaxt='n',
        ylim=c(-0.5,100.5), ylab='% Relative Abundance')
axis(side=2, at=seq(0,100,20), tick=TRUE)
abline(h=c(0, 100))
abline(h=seq(20,80,20), lty=2)
segments(x0=4.98, y0=0, x1=4.98, y1=100)
dev.off()

# Create a figure legend on a blank plot
pdf(file=legend_file, width=7, height=6)
plot(0, type='n', axes=F, xlab='', ylab='')
taxa <- gsub('_', ', ', rownames(t(final_shared)))
legend('center', legend=taxa, pt.bg=colors, pch=22, pt.cex=1.8)
dev.off()
