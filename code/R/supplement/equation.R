
pdf(file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/equation.pdf', width=10, height=6)
plot(1, type='n', axes=F, xlab='', ylab='', xlim=c(-5,5), ylim=c(-3,3))
text(0, 0, expression(Importance == paste(log[2],'( ',frac(Sigma * t[i], e[o]),' ','-',' ',frac(Sigma * t[o], e[i]),' )')), cex = 1.9) # Importance algorithm
dev.off()
