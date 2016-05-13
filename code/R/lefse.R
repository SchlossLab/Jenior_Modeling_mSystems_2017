

cef_lefse_file <- '/Users/mattjenior/Desktop/lefse/cefoperazone.final.0.03.Cdifficile630.0.03.lefse_summary'
clinda_lefse_file <- '/Users/mattjenior/Desktop/lefse/clindamycin.final.0.03.Cdifficile630.0.03.lefse_summary'
strep_lefse_file <- '/Users/mattjenior/Desktop/lefse/streptomycin.final.0.03.Cdifficile630.0.03.subsample.0.03.lefse_summary'
taxonomy_file <- '/Users/mattjenior/Desktop/lefse/allabx.final.0.03.cons.last.format.taxonomy'
plot_file <- '/Users/mattjenior/Desktop/lefse.pdf'

# Load in data
cef_lefse <- read.delim(cef_lefse_file, sep='\t', header=T, row.names=1)
clinda_lefse <- read.delim(clinda_lefse_file, sep='\t', header=T, row.names=1)
strep_lefse <- read.delim(strep_lefse_file, sep='\t', header=T, row.names=1)
taxonomy <- read.delim(taxonomy_file, sep='\t', header=T, row.names=1)

# Remove NAs
cef_lefse <- subset(cef_lefse, !is.na(LDA))
clinda_lefse <- subset(clinda_lefse, !is.na(LDA))
strep_lefse <- subset(strep_lefse, !is.na(LDA))

# Format data
cef_lefse <- cef_lefse[order(cef_lefse$Class, cef_lefse$LDA),]
clinda_lefse <- clinda_lefse[order(clinda_lefse$Class, clinda_lefse$LDA),]
strep_lefse <- strep_lefse[order(strep_lefse$Class, strep_lefse$LDA),]
taxonomy$Size <- NULL

# Merge with taxa names
cef_lefse <- merge(cef_lefse, taxonomy, by='row.names')
clinda_lefse <- merge(clinda_lefse, taxonomy, by='row.names')
strep_lefse <- merge(strep_lefse, taxonomy, by='row.names')

# Plot it
pdf(file=plot_file, width=24, height=10)
layout(matrix(1:3, nrow=1, byrow=T))
par(mar=c(5,17,3,1))
barplot(cef_lefse$LDA, horiz=TRUE, xlim=c(0,4), xlab='LDA', main='Cefoperazone', col='firebrick', names.arg=cef_lefse$Taxonomy, las=1, cex.main=3, cex.lab=2, cex.axis=1.5, cex.names=1.2)
box()
par(mar=c(5,17,3,1))
barplot(clinda_lefse$LDA, horiz=TRUE, xlim=c(0,5), xlab='LDA', main='Clindamycin', col=c('firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','darkblue','darkblue','darkblue'), names.arg=clinda_lefse$Taxonomy, las=1, cex.main=3, cex.lab=2, cex.axis=1.5, cex.names=1.2)
box()
par(mar=c(5,17,3,1), xpd=TRUE)
barplot(strep_lefse$LDA, horiz=TRUE, xlim=c(0,6), xlab='LDA', main='Streptomycin', col=c('firebrick','firebrick','firebrick','darkblue','darkblue','darkblue'), names.arg=strep_lefse$Taxonomy, las=1, cex.main=3, cex.lab=2, cex.axis=1.5, cex.names=1.2)
box()
dev.off()
