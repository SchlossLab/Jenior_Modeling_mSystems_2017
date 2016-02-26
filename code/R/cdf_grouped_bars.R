

expression_file <- '/Users/pschloss/Desktop/Jenior_812/cdf_genes.txt'

expression <- read.delim(expression_file, header = T, sep='\t', row.names=1)
pick_expression <- pick_expression[,3:6]
pick_expression [pick_expression  == 0] <- 1
log_pick_expression <- log10(pick_expression)
pick_log_pick_expression <- log_pick_expression[which(rowSums(log_pick_expression) > 1),]

par(las=1, mar=c(9,3,1,3), mgp=c(1,0,0))
barplot(t(pick_log_pick_expression), col=c('red2', 'blue2', 'forestgreen', 'gold2'), 
        beside=TRUE, xlim=c(0,750), ylim=c(0.5,4), xaxt='n', yaxt='n', 
        ylab='KEGG Pathway', cex.lab=1.5)

