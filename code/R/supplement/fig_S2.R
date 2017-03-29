
# Load dependencies
deps <- c('wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Plot logarithmic tick marks on axes
minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){
  
  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]
  
  major.ticks <- pretty(lims,n=5)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)
  
  major.ticks <- c(mn:mx)
  labelsY <- c(0,parse(text=paste(rep(10,mx), '^', seq(mn+1,mx,1), sep='')))
  axis(ax,at=major.ticks, labels=labelsY, las=1)
  
  n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]
  
  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]
  
  
  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}

# Agreggate c diff read abundances
collect_cdiff <- function(shared, antibiotic) {
  
  inf_shared <- subset(shared, infection == '630')
  inf_shared$infection <- NULL
  inf_shared$abx <- NULL
  mock_shared <- subset(shared, infection == 'mock')
  mock_shared$infection <- NULL
  mock_shared$abx <- NULL
  
  cdiff <- rowSums(cbind(inf_shared$Otu0004, inf_shared$Otu0309))
  inf_shared$Otu0004 <- NULL
  inf_shared$Otu0309 <- NULL
  other <- rowSums(inf_shared)
  inf_shared <- as.data.frame(cbind(cdiff, other))
  colnames(inf_shared) <- c('cdifficile', 'other')
  inf_shared$infection <- rep('630', nrow(inf_shared))
  inf_shared$abx <- rep(antibiotic, nrow(inf_shared))
  
  cdiff <- rowSums(cbind(mock_shared$Otu0004, mock_shared$Otu0309))
  mock_shared$Otu0004 <- NULL
  mock_shared$Otu0309 <- NULL
  other <- rowSums(mock_shared)
  mock_shared <- as.data.frame(cbind(cdiff, other))
  colnames(mock_shared) <- c('cdifficile', 'other')
  mock_shared$infection <- rep('mock', nrow(mock_shared))
  mock_shared$abx <- rep(antibiotic, nrow(mock_shared))
  
  final_data <- rbind(inf_shared, mock_shared)
  
  return(final_data)
}

# Define variables
nmds_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.unique_list.thetayc.0.03.lt.ave.nmds.axes'
summary_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.unique_list.groups.summary'
metadata_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'
shared_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
taxonomy_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.cons.taxonomy'

# Load in data
nmds <- read.delim(nmds_file, sep='\t', header=T, row.names=1)
nmds <- nmds[!rownames(nmds) %in% c('CefC5M2','StrepC4M1'), ] # Remove contaminated sample
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2','StrepC4M1'), ] # Remove contaminated sample
summary <- read.delim(summary_file, sep='\t', header=T, row.names=2)
summary <- summary[!rownames(summary) %in% c('CefC5M2','StrepC4M1'), ] # Remove contaminated sample
taxonomy <- read.delim(taxonomy_file, sep='\t', header=T, row.names=1)
shared <- read.delim(shared_file, sep='\t', header=T, row.names=2)
shared <- shared[!rownames(shared) %in% c('CefC5M2','StrepC4M1'), ]  # Remove contaminated sample
shared$numOtus <- NULL
shared$label <- NULL
rm(nmds_file, summary_file, metadata_file, shared_file, taxonomy_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data

# Combine the metadata with axes and format
nmds <- merge(metadata, nmds, by='row.names')
mock_axes <- subset(nmds, infection == 'mock')
mock_axes <- subset(mock_axes, type != 'germfree')
infected_axes <- subset(nmds, infection == '630')
infected_axes <- subset(infected_axes, type != 'germfree')

# Combine the metadata with summary, format, and get median invsimpson diversity
summary <- merge(metadata, summary, by='row.names')
summary <- subset(summary, type != 'germfree')
summary$abx <- factor(summary$abx, levels=c('streptomycin', 'cefoperazone', 'clindamycin', 'none'))
strep_div <- as.numeric(median(summary[summary$abx == 'streptomycin', 13]))
cef_div <- as.numeric(median(summary[summary$abx == 'cefoperazone', 13]))
clinda_div <- as.numeric(median(summary[summary$abx == 'clindamycin', 13]))
conv_div <- as.numeric(median(summary[summary$abx == 'none', 13]))

# Calculate C. diffile reads vs the rest of the community in each condition
shared <- merge(metadata, shared, by='row.names')
rownames(shared) <- shared$Row.names
shared$Row.names <- NULL
shared$cage <- NULL
shared$mouse <- NULL
shared$gender <- NULL
shared <- subset(shared, type != 'germfree')
shared$type <- NULL
strep_shared <- subset(shared, abx == 'streptomycin')
strep_shared <- collect_cdiff(strep_shared, 'streptomycin') 
strep_630 <- subset(strep_shared, infection == '630')
strep_mock <- subset(strep_shared, infection == 'mock')
cef_shared <- subset(shared, abx == 'cefoperazone')
cef_shared <- collect_cdiff(cef_shared, 'cefoperazone') 
cef_630 <- subset(cef_shared, infection == '630')
cef_mock <- subset(cef_shared, infection == 'mock')
clinda_shared <- subset(shared, abx == 'clindamycin')
clinda_shared <- collect_cdiff(clinda_shared, 'clindamycin') 
clinda_630 <- subset(clinda_shared, infection == '630')
clinda_mock <- subset(clinda_shared, infection == 'mock')
conv_shared <- subset(shared, abx == 'none')
conv_mock <- collect_cdiff(conv_shared, 'none') 
rm(metadata, strep_shared, cef_shared, clinda_shared, conv_shared)

shared_630 <- rbind(strep_630, cef_630, clinda_630)
shared_630$cdifficile <- as.numeric(as.character(shared_630$cdifficile))
shared_630$cdifficile <- log10(shared_630$cdifficile + 1)
shared_630$other <- as.numeric(as.character(shared_630$other))
shared_630$other <- log10(shared_630$other + 1)
shared_630$abx <- factor(shared_630$abx, levels=c('streptomycin', 'cefoperazone', 'clindamycin'))
shared_mock <- rbind(strep_mock, cef_mock, clinda_mock, conv_mock)
shared_mock$cdifficile <- as.numeric(as.character(shared_mock$cdifficile))
shared_mock$cdifficile <- log10(shared_mock$cdifficile + 1)
shared_mock$other <- as.numeric(as.character(shared_mock$other))
shared_mock$other <- log10(shared_mock$other + 1)
shared_mock$abx <- factor(shared_mock$abx, levels=c('streptomycin', 'cefoperazone', 'clindamycin', 'none'))

# Percentages of C diff reads to all others
percent_diff <- c(paste('~', as.character(round((mean(strep_630[,1]) / (mean(quantile(strep_630[,2]))+mean(strep_630[,1]))) * 100, digits=3)), '%', sep=''), 
                  paste('~', as.character(round((mean(strep_mock[,1]) / (mean(quantile(strep_mock[,2]))+mean(strep_mock[,1]))) * 100, digits=3)), '%', sep=''), 
                  paste('~', as.character(round((mean(cef_630[,1]) / (mean(quantile(cef_630[,2]))+mean(strep_630[,1]))) * 100, digits=3)), '%', sep=''), 
                  paste('~', as.character(round((mean(cef_mock[,1]) / (mean(quantile(cef_mock[,2]))+mean(cef_mock[,1]))) * 100, digits=3)), '%', sep=''),
                  paste('~', as.character(round((mean(clinda_630[,1]) / (mean(quantile(strep_630[,2]))+mean(clinda_630[,1]))) * 100, digits=3)), '%', sep=''), 
                  paste('~', as.character(round((mean(clinda_mock[,1]) / (mean(quantile(strep_mock[,2]))+mean(clinda_mock[,1]))) * 100, digits=3)), '%', sep=''),
                  paste('~', as.character(round((mean(conv_mock[,1]) / (mean(quantile(conv_mock[,2]))+mean(conv_mock[,1]))) * 100, digits=3)), '%', sep=''))
rm(strep_630, cef_630, clinda_630, strep_mock, cef_mock, clinda_mock, conv_mock)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S2.pdf'
pdf(file=plot_file, width=10, height=7)
layout(matrix(c(1,1,2,2,
                3,3,3,3), nrow=2, ncol=4, byrow=TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# NMDS of treatment groups
par(las=1, mar=c(4,4,1,1), mgp=c(2.5,0.7,0))
plot(nmds$axis1, nmds$axis2, pch=21, cex=0,
     xlim=c(-0.8,0.8), ylim=c(-0.8,0.8), cex.lab=1.2,
     xlab='NMDS Axis 1', ylab='NMDS Axis 2')

# add mock points
points(x=mock_axes$axis1, y=mock_axes$axis2, 
       col=c(wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black', wes_palette("FantasticFox")[1])[mock_axes$abx], 
       pch=1, lwd=2, cex=2)
# add mock infected points
points(x=infected_axes$axis1, y=infected_axes$axis2, 
       col=c(wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black', wes_palette("FantasticFox")[1])[infected_axes$abx], 
       pch=1, lwd=2, cex=2)

# Add legends
legend('bottomleft', legend=c('Streptomycin-treated', 'Cefoperzone-treated', 'Clindamycin-treated', 'No Antibiotics'), 
       col=c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'black'), 
       pch=15, pt.cex=1.7, bty='n')

mtext('A', side=2, line=2, las=2, adj=1.3, padj=-7.5, cex=1.5)

#-----------------------#

# Stripchart of Inverse simpson diversity
par(las=1, mar=c(2,4,1,1), mgp=c(2.5,0.7,0))
stripchart(invsimpson~abx, data=summary, vertical=T, pch=19, 
           ylim=c(0,20), xaxt='n', cex=1.5,
           col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5],'black'),
           ylab='Inv. Simpson Diversity', method='jitter', jitter=0.15, cex.lab=1.2)
axis(side=1, at=c(1:4), labels=c('Streptomycin','Cefoperazone','Clindamycin','No Antibiotics'), cex.axis=1.2, tick=FALSE)

# add medians
segments(x0=0.7, y0=strep_div, x1=1.3, y1=strep_div, lwd=2.5)
segments(x0=1.7, y0=cef_div, x1=2.3, y1=cef_div, lwd=2.5)
segments(x0=2.7, y0=clinda_div, x1=3.3, y1=clinda_div, lwd=2.5)
segments(x0=3.7, y0=conv_div, x1=4.3, y1=conv_div, lwd=2.5)
text(x=4, y=19, '*', font=2, cex=2) # add significance

mtext('B', side=2, line=2, las=2, adj=1.3, padj=-8.2, cex=1.5)

#-----------------------#

# C. difficile read abundance relative to community
par(las=1, mar=c(4,4,1,1), mgp=c(2.5,0.7,0))
boxplot(cdifficile~abx, data=shared_630, cex=0, lwd=2.5, xlim=c(0.5,23.5), ylim=c(0,4), at=c(1,8,15),
           col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5]),
           ylab='16S V4 Region Amplicon Read Abundance', staplewex=0.8, boxwex=1, lty=1, medlwd=2.5, xaxt='n', yaxt='n')
boxplot(other~abx, data=shared_630, cex=0, lwd=2.5, xlim=c(0.5,23.5), ylim=c(0,4), at=c(2,9,16),
           col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5]),
           ylab='16S V4 Region Amplicon Read Abundance', staplewex=0.8, boxwex=1, lty=1, medlwd=2.5, xaxt='n', yaxt='n', add=TRUE)
boxplot(cdifficile~abx, data=shared_mock, cex=0, lwd=2.5, xlim=c(0.5,23.5), ylim=c(0,4), at=c(4,11,18,22),
        col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5],'gray50'),
        ylab='16S V4 Region Amplicon Read Abundance', staplewex=0.8, boxwex=1, lty=1, medlwd=2.5, xaxt='n', yaxt='n', add=TRUE)
boxplot(other~abx, data=shared_mock, cex=0, lwd=2.5, xlim=c(0.5,23.5), ylim=c(0,4), at=c(5,12,19,23),
        col=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5],'gray50'),
        ylab='16S V4 Region Amplicon Read Abundance', staplewex=0.8, boxwex=1, lty=1, medlwd=2.5, xaxt='n', yaxt='n', add=TRUE)
mtext(c('630','Other','630','Other','630','Other','630','Other','630','Other','630','Other','630','Other'), # define symbols
      side=1, at=c(1,2,4,5,8,9,11,12,15,16,18,19,22,23), padj=0.2, cex=0.7)
mtext(c('Infected','Mock','Infected','Mock','Infected','Mock','Mock'), side=1, at=c(1.5,4.5,8.5,11.5,15.5,18.5,22.5), padj=1.8, cex=0.8)
mtext(c('Streptomycin','Cefoperazone','Clindamycin','No Antibiotics'), side=1, at=c(3,10,17,22.5), padj=3)
minor.ticks.axis(2, 10, mn=0, mx=4)
abline(v=c(6.5,13.5,20.5), lty=2)
text(x=c(1.5,4.5,8.5,11.5,15.5,18.5,22.5), y=4, labels=percent_diff, font=2)
text(x=c(2,9,16), y=c(3.6,3.6,3.6), labels=c('*','*','*'), font=2, cex=2)
mtext('C', side=2, line=2, las=2, adj=1.3, padj=-7.5, cex=1.5)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps) {
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only=TRUE)
}
rm(list=ls())
gc()
