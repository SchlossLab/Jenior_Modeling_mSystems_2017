
# Load dependencies
deps <- c('shape', 'plotrix');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_1.pdf'
pdf(file=plot_file, width=8, height=4)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Timeline of mouse experiments
par(mar=c(0,3,0,0))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-4.8,4), ylim=c(-2,5))

# Abx in drinking water timeline
rect(xleft=-4, ybottom=2.8, xright=0, ytop=3.2, col='darkorchid2', border='black')
Arrows(x0=-4, y0=3, x1=3.5, y1=3, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,0,2,2.75), y0=c(3.5,3.5,3.5,3.5), x1=c(-4,0,2,2.75), y1=c(2.5,2.5,2.5,2.5), lwd=4)
segments(x0=c(-4,-3,-2,-1,1), y0=c(3.25,3.25,3.25,3.25,3.25), x1=c(-4,-3,-2,-1,1), y1=c(2.75,2.75,2.75,2.75,2.75), lwd=2)
points(x=c(2,2.75), y=c(4,4), pch=25, bg=c('white','black'), col='black', cex=2.5)
text(x=c(-4,0,2,2.75), y=c(2.2,2.2,2.2,2.2), c('Day -7', 'Day -2', 'Day 0', '18 hrs'), cex=0.9)
#text(x=-4.6, y=3.2, 'Cefoperazone', cex=0.7)
#text(x=-4.6, y=2.95, 'or', font=2, cex=0.8)
#text(x=-4.6, y=2.7, 'Streptomycin', cex=0.7)
text(x=-4.6, y=2.95, 'a', cex=2)

# IP injection abx timeline
Arrows(x0=-4, y0=0, x1=-1.5, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,-3,-2.25), y0=c(-0.5,-0.5,-0.5), x1=c(-4,-3,-2.25), y1=c(0.5,0.5,0.5), lwd=4)
points(x=c(-4,-3,-2.25), y=c(1,1,1), pch=c(25,25,25), bg=c('chocolate1','white','black'), col='black', cex=2.5)
text(x=c(-4,-3,-2.25), y=c(-0.8,-0.8,-0.8), c('Day -1', 'Day 0', '18 hrs'), cex=0.9)
#text(x=-4.6, y=0, 'Clindamycin', cex=0.7)
text(x=-4.6, y=0, 'b', cex=2)

# Legend
legend(x=0, y=1.1, legend=expression('Antibiotic in Drinking Water', 'IP Injection of Antibiotic', paste(italic('C. difficile'), ' Spore Gavage'), 'Euthanize & Necropsy'), 
       pt.bg=c('darkorchid2','chocolate1','white','black'), pch=c(22,25,25,25), pt.cex=c(2.5,2,2,2), bty='n')

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
rm(plot_file)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()
