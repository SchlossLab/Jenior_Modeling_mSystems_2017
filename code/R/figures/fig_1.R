
# Load dependencies
deps <- c('shape', 'plotrix', 'grid', 'gridExtra');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in and format antibiotic table
abx_table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/abx.tsv'
abx_table <- read.delim(abx_table_file, sep='\t', header=TRUE)

abx_table[,2] <- gsub('_', ' ', abx_table[,2])
abx_table[,3] <- gsub('_', ' ', abx_table[,3])
abx_table[,4] <- gsub('_', ' ', abx_table[,4])
abx_table[,5] <- gsub('_', ' ', abx_table[,5])
abx_table[,6] <- gsub('_', ' ', abx_table[,6])

abx_table$Target <- sapply(lapply(abx_table$Target, strwrap, width=40), paste, collapse="\n")
abx_table$Activity <- sapply(lapply(abx_table$Activity, strwrap, width=40), paste, collapse="\n")

#-------------------------------------------------------------------------------------------------------------------------------------#

plot_file1 <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_1A.pdf'
pdf(file=plot_file1, width=7, height=4)

# Timeline of mouse experiments

# Create an empty plot
par(mar=c(0,3,0,0))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-4.8,4), ylim=c(-2,5))

# Abx in drinking water timeline
rect(xleft=-4, ybottom=2.8, xright=0, ytop=3.2, col='darkorchid2', border='black')
Arrows(x0=-4, y0=3, x1=3.5, y1=3, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,0,2,2.75), y0=c(3.5,3.5,3.5,3.5), x1=c(-4,0,2,2.75), y1=c(2.5,2.5,2.5,2.5), lwd=4)
segments(x0=c(-4,-3,-2,-1,1), y0=c(3.25,3.25,3.25,3.25,3.25), x1=c(-4,-3,-2,-1,1), y1=c(2.75,2.75,2.75,2.75,2.75), lwd=2)
points(x=c(2,2.75), y=c(4,4), pch=25, bg=c('white','black'), col='black', cex=2.5)
text(x=c(-4,0,2,2.75), y=c(2.2,2.2,2.2,2.2), c('Day -7', 'Day -2', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.6, y=3.2, 'Cefoperazone', cex=0.7)
text(x=-4.6, y=2.95, 'or', font=2, cex=0.8)
text(x=-4.6, y=2.7, 'Streptomycin', cex=0.7)

# IP injection abx timeline
Arrows(x0=-4, y0=0, x1=-1.5, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,-3,-2.25), y0=c(-0.5,-0.5,-0.5), x1=c(-4,-3,-2.25), y1=c(0.5,0.5,0.5), lwd=4)
points(x=c(-4,-3,-2.25), y=c(1,1,1), pch=c(25,25,25), bg=c('gray60','white','black'), col='black', cex=2.5)
text(x=c(-4,-3,-2.25), y=c(-0.8,-0.8,-0.8), c('Day -1', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.6, y=0, 'Clindamycin', cex=0.7)

# Legend
legend(x=0, y=1.3, legend=expression('Antibiotic in Drinking Water', 'IP Injection of Antibiotic', paste(italic('C. difficile'), ' Spore Gavage'), 'Euthanize & Necropsy'), 
       pt.bg=c('darkorchid2','gray60','white','black'), pch=c(22,25,25,25), pt.cex=c(2.5,2,2,2), bty='n')

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

plot_file2 <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_1B.pdf'
pdf(file=plot_file2, width=15, height=3)

grid.table(abx_table, rows=NULL)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

# COMBINE THE TIMELINE AND ABX TABLE OUTIDE OF R

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
rm(abx_table_file, abx_table, plot_file1, plot_file2)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()
