
rep1 <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/scaria_2013/cd630_sucrose_rep1.final.pool.norm.tsv', sep='\t', header=FALSE)
colnames(rep1) <- c('gene', 'description', 'reads')
rep1$description <- NULL
rep2 <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/scaria_2013/cd630_sucrose_rep2.final.pool.norm.tsv', sep='\t', header=FALSE)
colnames(rep2) <- c('gene', 'description', 'reads')
rep2$description <- NULL
ko <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/kegg/cdf_ko.list', sep='\t', header=FALSE)
colnames(ko) <- c('gene', 'ko')

rep1 <- merge(ko, rep1, by='gene', all.y=TRUE)
rep1 <- rep1[complete.cases(rep1),]
rep1$gene <- NULL
rep2 <- merge(ko, rep2, by='gene', all.y=TRUE)
rep2 <- rep2[complete.cases(rep2),]
rep2$gene <- NULL
rm(ko)

write.table(rep1, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/scaria_2013/sucrose_rep1_ko.tsv', sep='\t', row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(rep2, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/scaria_2013/sucrose_rep2_ko.tsv', sep='\t', row.names=FALSE, quote=FALSE, col.names=FALSE)
rm(rep1, rep2)

#----------------------------------------------------------------------------------------------#

scores1 <- read.delim('/home/matt/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/scaria_2013/sucrose_rep1.bipartite.files/importances.tsv', sep='\t', header=TRUE, row.names=1)
scores1$p_value <- NULL
scores2 <- read.delim('/home/matt/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/scaria_2013/sucrose_rep2.bipartite.files/importances.tsv', sep='\t', header=TRUE, row.names=1)
scores2$Metabolite_name <- NULL
scores2$p_value <- NULL
scores <- merge(scores1, scores2, by='row.names')
scores$Row.names <- NULL
rm(scores1, scores2)
colnames(scores) <- c('metabolite','rep1','rep2')
scores$difference <- abs(scores$rep1 - scores$rep2)
quantile(scores$difference)
var(scores$difference)
scores$rep1 <- 2 ^ abs(scores$rep1)
scores$rep2 <- 2 ^ abs(scores$rep2)
scores$difference <- abs(scores$rep1 - scores$rep2)
quantile(scores$difference)
var(scores$difference)




