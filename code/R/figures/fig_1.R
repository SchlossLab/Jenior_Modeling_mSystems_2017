
# Load dependencies
deps <- c('shape', 'plotrix', 'wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_1.pdf'
pdf(file=plot_file, width=7, height=10)

# Create layout for multi-plot
layout(mat=matrix(c(1,
                    2,
                    3), nrow=3, ncol=1, byrow=TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

par(mar=c(0,0,0,0))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-5,5), ylim=c(-2,2)) # Empty plot

# Strep in drinking water timeline
rect(xleft=-4, ybottom=-0.4, xright=1, ytop=0.4, col=wes_palette("FantasticFox")[1], border='black')
Arrows(x0=-4, y0=0, x1=4.8, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,1,3,3.75), y0=c(0.5,0.5,0.5,0.5), x1=c(-4,1,3,3.75), y1=c(-0.5,-0.5,-0.5,-0.5), lwd=4)
segments(x0=c(-4,-3,-2,-1,0,2), y0=c(0.25,0.25,0.25,0.25,0.25), x1=c(-4,-3,-2,-1,0,2), y1=c(-0.25,-0.25,-0.25,-0.25,-0.25), lwd=2)
points(x=c(0,3,3.75), y=c(0.9,0.9,0.9), pch=25, bg=c('white','darkorchid2','black'), col='black', cex=3.4)
text(x=c(-4,1,3,3.75), y=c(-0.8,-0.8,-0.8,-0.8), c('Day -7', 'Day -2', 'Day 0', '18 Hrs'), cex=1.3)
text(x=-4.6, y=0, 'a', cex=2, font=2)

#----------------------------#

par(mar=c(0,0,0,0))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-5,5), ylim=c(-2,2)) # Empty plot

# Cef in drinking water timeline
rect(xleft=-4, ybottom=-0.4, xright=1, ytop=0.4, col=wes_palette("FantasticFox")[3], border='black')
Arrows(x0=-4, y0=0, x1=4.8, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,1,3,3.75), y0=c(0.5,0.5,0.5,0.5), x1=c(-4,1,3,3.75), y1=c(-0.5,-0.5,-0.5,-0.5), lwd=4)
segments(x0=c(-4,-3,-2,-1,0,2), y0=c(0.25,0.25,0.25,0.25,0.25), x1=c(-4,-3,-2,-1,0,2), y1=c(-0.25,-0.25,-0.25,-0.25,-0.25), lwd=2)
points(x=c(0,3,3.75), y=c(0.9,0.9,0.9), pch=25, bg=c('white','darkorchid2','black'), col='black', cex=3.4)
text(x=c(-4,1,3,3.75), y=c(-0.8,-0.8,-0.8,-0.8), c('Day -7', 'Day -2', 'Day 0', '18 Hrs'), cex=1.3)
text(x=-4.6, y=0, 'b', cex=2, font=2)

#----------------------------#

par(mar=c(0,0,0,0))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-5,5), ylim=c(-2,2)) # Empty plot

# Clinda IP injection abx timeline
Arrows(x0=-3, y0=0, x1=-0.2, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-3,-2,-1.25), y0=c(-0.5,-0.5,-0.5), x1=c(-3,-2,-1.25), y1=c(0.5,0.5,0.5), lwd=4)
points(x=c(-3,-2,-1.25), y=c(1,1,1), pch=c(25,25,25), bg=c(wes_palette("FantasticFox")[5],'darkorchid3','black'), col='black', cex=3.4)
text(x=c(-3,-2,-1.25), y=c(-0.8,-0.8,-0.8), c('Day -1', 'Day 0', '18 Hrs'), cex=1.3)
text(x=-3.6, y=0, 'c', cex=2, font=2)

# Legend
legend(x=0.9, y=1.4, legend=expression('Streptomycin in drinking water', 'Cefoperazone in drinking water', 'Clindamycin IP injection', 'Switch to untreated drinking water', paste(italic('C. difficile'), ' str. 630 spore gavage'), 'Euthanize & necropsy'), 
       pt.bg=c(wes_palette("FantasticFox")[1],wes_palette("FantasticFox")[3],wes_palette("FantasticFox")[5],'white','darkorchid2','black'), cex=1.5,  pch=c(22,22,25,25,25,25), pt.cex=c(3.5,3.5,3,3,3,3), bty='n')

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
dev.off()
rm(plot_file, abx_table_file, abx_table)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()
