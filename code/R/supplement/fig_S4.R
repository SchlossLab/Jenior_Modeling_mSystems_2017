

pdf(file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S3.pdf', width=10, height=10)

par(mar=c(0,0,0,0))
triplot(x=combined_mapping[,1], y=combined_mapping[,2], z=combined_mapping[,3], 
        frame=TRUE, label=c('Cefoperazone','Clindamycin','Streptomycin'), grid=FALSE, cex=0)

# 50% lines
lines(x=c(-0.577,0.288), y=c(-0.333,0.1665))
lines(x=c(0,0), y=c(-0.333,0.665))
lines(x=c(-0.288,0.577), y=c(0.1665,-0.333))

dev.off()

# The rest of the figure is done in an image editor


