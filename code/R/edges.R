# Read in data
degree_file <- '/Users/schloss/Desktop/select_degree.tsv'
degrees <- read.delim(degree_file, sep='\t', row.names=1, header=TRUE)

# Subset each group
both <- subset(degrees, group == 'both')
both$group <- NULL
rownames(both) <- c('Pyruvate', 'Acetyl-CoA', 'Fructose-6-P', 'Oxoglutarate', 'Glyceraldehyde-3-P', 'Glutamate')
input <- subset(degrees, group == 'input')
input$group <- NULL
rownames(input) <- c('Tyrosine', 'P-Glycerate', 'Homoserine', 'Mannose', 'Glycerol-3-P', 'Glucose')
output <- subset(degrees, group == 'output')
output$group <- NULL
rownames(output) <- c('DNA', 'Malonyl-CoA', 'Glycerone-P', 'Acetyl-P', 'Glyceroyl-P', 'Isopropylmaleate', 'Uracil')

# Figure variables
both_plot <- '/Users/pschloss/Desktop/both_select.pdf'
input_plot <- '/Users/pschloss/Desktop/input_select.pdf'
output_plot <- '/Users/pschloss/Desktop/output_select.pdf'
legend_plot <- '/Users/pschloss/Desktop/legend_degree.pdf'

# Plot each seperately
pdf(file=both_plot, width=10, height=6)
par(las=1, mar=c(6,4,1,1), mgp=c(3.5,1,0))
barplot(t(both), col=c('red2', 'blue2'), beside=TRUE, xlim=c(0.8,18.3), ylim=c(0,100), xaxt='n', yaxt='n')
text(c(1,4,7,10,13,16.5), par("usr")[3] - 10, labels = rownames(both), srt = 45, pos = 1, xpd = TRUE, cex=0.8, font=2)
axis(side=2, at=seq(0,100,20), seq(0,100,20), tick=TRUE, cex=1.2)
mtext('Number of Edges', side=2, at=50, cex=1.5, padj=-3, las=3)
abline(h=c(20,40,60,80), lty=2)
box()
legend('topright', legend=c('Reaction Output', 'Reaction Input'), col=c('red2', 'blue2'), pch=15, cex=2, pt.cex=4)
dev.off()

pdf(file=input_plot, width=10, height=6)
par(las=1, mar=c(6,4,1,1), mgp=c(3.5,1,0))
barplot(t(input), col=c('red2', 'blue2'), beside=TRUE, xlim=c(0.8,18.3), ylim=c(0,25), xaxt='n', yaxt='n')
text(c(1.4,4.3,7.2,10.3,13.3,16.5), par("usr")[3] - 2, labels = rownames(input), srt = 45, pos = 1, xpd = TRUE, cex=0.8, font=2)
axis(side=2, at=seq(0,25,5), seq(0,25,5), tick=TRUE, cex=1.2)
mtext('Number of Edges', side=2, at=12.5, cex=1.5, padj=-3, las=3)
abline(h=c(10,20), lty=2)
box()
legend('topright', legend=c('Reaction Output', 'Reaction Input'), col=c('red2', 'blue2'), pch=15, cex=2, pt.cex=4)
dev.off()

pdf(file=output_plot, width=10, height=6)
par(las=1, mar=c(6,4,1,1), mgp=c(3.5,1,0))
barplot(t(output), col=c('red2', 'blue2'), beside=TRUE, xlim=c(0.8,21), ylim=c(0,50), xaxt='n', yaxt='n')
text(c(1.7,4.5,7.5,10.5,13,16,20), par("usr")[3] - 5, labels = rownames(output), srt = 45, pos = 1, xpd = TRUE, cex=0.8, font=2)
axis(side=2, at=seq(0,50,10), seq(0,50,10), tick=TRUE, cex=1.2)
mtext('Number of Edges', side=2, at=25, cex=1.5, padj=-3, las=3)
abline(h=c(10,20,30,40), lty=2)
box()
legend('topright', legend=c('Reaction Output', 'Reaction Input'), col=c('red2', 'blue2'), pch=15, cex=2, pt.cex=4)
dev.off()






