
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
