
# Load dependencies
deps <- c('shape', 'plotrix', 'wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S1.pdf'
pdf(file=plot_file, width=6, height=9)

# Create layout for multi-plot
layout(mat=matrix(c(1,
                    2,
                    3,
                    4,
                    5), nrow=5, ncol=1, byrow=TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

par(mar=c(0,0,0,0))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-5,5), ylim=c(-2,2)) # Empty plot

# Legend
legend('center', legend=expression('Untreated water', paste(italic('C. difficile'), ' gavage'), 'Euthanize & Necropsy'), 
       pt.bg=c(wes_palette("Chevalier")[3],'darkorchid2','black'), cex=1.5,  pch=c(22,25,25), pt.cex=c(4,3,3), bty='n', ncol=3)

#----------------------------#

plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-5,5), ylim=c(-2,2)) # Empty plot

# Strep in drinking water timeline
rect(xleft=-4, ybottom=-0.4, xright=1, ytop=0.4, col=wes_palette("FantasticFox")[1], border='black')
rect(xleft=1, ybottom=-0.4, xright=3.75, ytop=0.4, col=wes_palette("Chevalier")[3], border='black')
Arrows(x0=-4, y0=0, x1=4.8, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,1,3,3.75), y0=c(0.5,0.5,0.5,0.5), x1=c(-4,1,3,3.75), y1=c(-0.5,-0.5,-0.5,-0.5), lwd=4)
segments(x0=c(-3,-2,-1,0,2), y0=c(0.25,0.25,0.25,0.25), x1=c(-3,-2,-1,0,2), y1=c(-0.25,-0.25,-0.25,-0.25), lwd=2)
text(x=c(-4,3), y=c(-0.8,-0.8), c('Day -7', 'Day 0'), cex=1.5)
legend('topleft', legend='Streptomycin (SPF)', pt.cex=0, cex=1.5, bty='n')
text(x=-4.5, y=0, 'A', cex=1.9)

points(x=c(3,3.75), y=c(1,1), pch=25, bg=c('darkorchid2','black'), col='black', cex=3.4)

#----------------------------#

plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-5,5), ylim=c(-2,2)) # Empty plot

# Cef in drinking water timeline
rect(xleft=-4, ybottom=-0.4, xright=1, ytop=0.4, col=wes_palette("FantasticFox")[3], border='black')
rect(xleft=1, ybottom=-0.4, xright=3.75, ytop=0.4, col=wes_palette("Chevalier")[3], border='black')
Arrows(x0=-4, y0=0, x1=4.8, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,1,3,3.75), y0=c(0.5,0.5,0.5,0.5), x1=c(-4,1,3,3.75), y1=c(-0.5,-0.5,-0.5,-0.5), lwd=4)
segments(x0=c(-4,-3,-2,-1,0,2), y0=c(0.25,0.25,0.25,0.25,0.25), x1=c(-4,-3,-2,-1,0,2), y1=c(-0.25,-0.25,-0.25,-0.25,-0.25), lwd=2)
text(x=c(-4,3), y=c(-0.8,-0.8), c('Day -7', 'Day 0'), cex=1.5)
legend('topleft', legend='Cefoperazone (SPF)', pt.cex=0, cex=1.5, bty='n')

text(x=-4.5, y=0, 'B', cex=1.9)

points(x=c(3,3.75), y=c(1,1), pch=25, bg=c('darkorchid2','black'), col='black', cex=3.4)

#----------------------------#

plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-5,5), ylim=c(-2,2)) # Empty plot

# Clinda IP injection abx timeline
rect(xleft=-4, ybottom=-0.4, xright=3.75, ytop=0.4, col=wes_palette("Chevalier")[3], border='black')
Arrows(x0=-4, y0=0, x1=4.8, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,2,3,3.75), y0=c(0.5,0.5,0.5,0.5), x1=c(-4,2,3,3.75), y1=c(-0.5,-0.5,-0.5,-0.5), lwd=4)
segments(x0=c(-4,-3,-2,-1,0,1), y0=c(0.25,0.25,0.25,0.25,0.25), x1=c(-4,-3,-2,-1,0,1), y1=c(-0.25,-0.25,-0.25,-0.25,-0.25), lwd=2)
points(x=2, y=1, pch=25, bg=c(wes_palette("FantasticFox")[5]), col='black', cex=3.4)
text(x=c(-4,3), y=c(-0.8,-0.8), c('Day -7', 'Day 0'), cex=1.5)
legend('topleft', legend='Clindamycin IP Injection (SPF)', pt.cex=0, cex=1.5, bty='n')
text(x=-4.5, y=0, 'C', cex=1.9)

points(x=c(3,3.75), y=c(1,1), pch=25, bg=c('darkorchid2','black'), col='black', cex=3.4)

#----------------------------#

plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-5,5), ylim=c(-2,2)) # Empty plot

# No antibiotics and germ free timeline
rect(xleft=-4, ybottom=-0.4, xright=3.75, ytop=0.4, col=wes_palette("Chevalier")[3], border='black')
Arrows(x0=-4, y0=0, x1=4.8, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,3,3.75), y0=c(0.5,0.5,0.5), x1=c(-4,3,3.75), y1=c(-0.5,-0.5,-0.5), lwd=4)
segments(x0=c(-4,-3,-2,-1,0,1,2), y0=c(0.25,0.25,0.25,0.25,0.25,0.25,0.25), x1=c(-4,-3,-2,-1,0,1,2), y1=c(-0.25,-0.25,-0.25,-0.25,-0.25,-0.25,-0.25), lwd=2)
text(x=c(-4,3), y=c(-0.8,-0.8), c('Day -7', 'Day 0'), cex=1.5)
legend('topleft', legend='No Antibiotics (SPF & GF)', pt.cex=0, cex=1.5, bty='n')
text(x=-4.5, y=0, 'D', cex=1.9)

points(x=c(3,3.75), y=c(1,1), pch=25, bg=c('darkorchid2','black'), col='black', cex=3.4)

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
