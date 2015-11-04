
coverage <- read.delim('~/Desktop/cefoperazone_630.mapped2genes.txt', sep='\t', header=F)
colnames(coverage) <- c('ref_name','ref_length','mapped','unmapped')
coverage$unmapped <- NULL
coverage <- coverage[order(coverage$mapped, decreasing=T),]
coverage$read_len <- 50
coverage$coverage <- (coverage$mapped / coverage$ref_length) * coverage$read_len
rownames(coverage) <- coverage$ref_name
coverage$ref_name <- NULL
coverage$mapped <- NULL
coverage$read_len <- NULL

write.table(coverage, file='~/Desktop/cef_630.coverage.txt', sep='\t', quote=F, row.names=coverage$ref_name)
