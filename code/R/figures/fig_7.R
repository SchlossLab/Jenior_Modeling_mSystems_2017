
# Load dependencies
deps <- c('wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Select files
concentrations <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/ms_substrates.tsv'
metadata <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'

# Read in data
concentrations <- read.delim(concentrations, sep='\t', header=T, row.names=1)
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$type <- NULL

# Subset data and dd metadata column
shared <- t(concentrations[1:2,])
shared <- merge(shared, metadata, by='row.names')
rownames(shared) <- shared$Row.names
shared$Row.names <- NULL
shared <- subset(shared, infection != '630')
shared$infection <- NULL
strep <- t(concentrations[3:4,c(1:9,46:63)])
strep <- merge(strep, metadata, by='row.names')
rownames(strep) <- strep$Row.names
strep$Row.names <- NULL
strep <- subset(strep, infection != '630')
strep$infection <- NULL
cef <- t(concentrations[3,c(1:9,10:27)])
cef <- merge(cef, metadata, by='row.names')
rownames(cef) <- cef$Row.names
cef$Row.names <- NULL
cef <- subset(cef, infection != '630')
cef$infection <- NULL
clinda <- t(concentrations[5,c(1:9,28:45)])
clinda <- merge(clinda, metadata, by='row.names')
rownames(clinda) <- clinda$Row.names
clinda$Row.names <- NULL
clinda <- subset(clinda, infection != '630')
clinda$infection <- NULL
gf <- t(concentrations[6:12,c(1:9,64:81)])
gf <- merge(gf, metadata, by='row.names')
rownames(gf) <- gf$Row.names
gf$Row.names <- NULL
gf <- subset(gf, infection != '630')
gf$infection <- NULL
rm(concentrations, metadata)

# Calculate medians








#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_7.pdf'
select_palette <- c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black')
pdf(file=plot_file, width=7.5, height=9)
layout(matrix(c(1,2,3,
                4,5,5),
              nrow=2, ncol=3, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# A. Shared
acetylglucosamine <- shared[,c(1,3)]
proline <- shared[,2:3]
par(las=1, mar=c(0.7,4,1,1), mgp=c(2.5,0.7,0), yaxs='i')
stripchart(Nacetylglucosamine_Nacetylgalactosamine~abx, data=acetylglucosamine, vertical=T, pch=1, lwd=2.5, col=select_palette,
           xaxt='n', cex=2, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(proline~abx, data=proline, vertical=T, pch=1, lwd=2.5, col=select_palette,
           xaxt='n', cex=2, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2)



axis()


segments(0.6, vege_medians[1], 1.4, vege_medians[1], lwd=3) # cefoperazone
text(c(2.5,3,3.5), c(3.16,3.31,3.46), labels=c('*','*','*'), cex=2.2)
mtext('A', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.1)









#--------------------------------#

# B. Streptomycin-only
par(las=1, mar=c(0.7,4,1,1), mgp=c(2.5,0.7,0), yaxs='i')
stripchart(mannitol_sorbitol~abx, data=strep, vertical=T, pch=1, lwd=2.5,
           cex=2, col=wes_palette("FantasticFox")[1], ylim=c(0,5),
           ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2)



axis()


segments(0.6, vege_medians[1], 1.4, vege_medians[1], lwd=3)
text(c(2.5,3,3.5), c(3.16,3.31,3.46), labels=c('*','*','*'), cex=2.2)

mtext('B', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.1)

#--------------------------------#

# C. Cefoperazone-only
par(las=1, mar=c(0.7,4,1,1), mgp=c(2.5,0.7,0), yaxs='i')
stripchart(mannitol_sorbitol~abx, data=cef, vertical=T, pch=1, lwd=2.5,
           ylim=c(-5,50), xaxt='n', cex=2, col=wes_palette("FantasticFox")[3],
           ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2)




axis()



segments(0.6, vege_medians[1], 1.4, vege_medians[1], lwd=3)
text(c(2.5,3,3.5), c(3.16,3.31,3.46), labels=c('*','*','*'), cex=2.2)

mtext('C', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.1)

#--------------------------------#

# D. Clindamycin-only
par(las=1, mar=c(0.7,4,1,1), mgp=c(2.5,0.7,0), yaxs='i')
stripchart(~abx, data=clinda, vertical=T, pch=1, lwd=2.5,
           ylim=c(1,9), xaxt='n', cex=2, col=wes_palette("FantasticFox")[5],
           ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2)


axis()


segments(0.6, vege_medians[1], 1.4, vege_medians[1], lwd=3)
text(c(2.5,3,3.5), c(3.16,3.31,3.46), labels=c('*','*','*'), cex=2.2)

mtext('D', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.1)

#--------------------------------#

# E. Germ free-only
par(las=1, mar=c(0.7,4,1,1), mgp=c(2.5,0.7,0), yaxs='i')
stripchart(~abx, data=gf, vertical=T, pch=1, lwd=2.5,
           ylim=c(1,9), xaxt='n', cex=2, col='forestgreen',
           ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2)




axis()



segments(0.6, vege_medians[1], 1.4, vege_medians[1], lwd=3)
text(c(2.5,3,3.5), c(3.16,3.31,3.46), labels=c('*','*','*'), cex=2.2)

mtext('E', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.1)

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
dev.off()
rm(plot_file, select_palette)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

