
# Import read-to-gene mapping data
coverage <- read.delim('~/Desktop/cefoperazone_630.mapped2genes.txt', sep='\t', header=F)
colnames(coverage) <- c('ref_name','ref_length','mapped','unmapped')
coverage$unmapped <- NULL
coverage <- coverage[order(coverage$mapped, decreasing=T),]
coverage$read_len <- 50
coverage$norm_coverage <- round((coverage$mapped / coverage$ref_length) * coverage$read_len)
rownames(coverage) <- coverage$ref_name
coverage$ref_name <- NULL
coverage$mapped <- NULL
coverage$read_len <- NULL
coverage$ref_length <- NULL

# Import BLAST hits
blast <- read.delim('~/Desktop/Cefoperazone.protVprot.out', sep='\t', header=F, 
                         row.names=1, colClasses = c('character', 'character', 'numeric', rep('NULL', 7), 'numeric', 'NULL'))
colnames(blast) <- c('gene', 'percent_id', 'evalue')

# Import reference tables
genes <- read.delim('~/Desktop/ko_genes.list', sep='\t', header=F) # Gene to KO translation
colnames(genes) <- c('ko', 'gene')
definition <- read.delim('~/Desktop/ko_definition.list', sep='\t', header=F) # KO annotations
colnames(definition) <- c('ko', 'definition')

# Merge the tables to annotate genes
annotated <- merge(blast, coverage, by='row.names')
annotated <- merge(annotated, genes, by='gene')
annotated <- merge(annotated, definition, by='ko')
annotated <- annotated[order(annotated$norm_coverage, decreasing=T),]


# Write the new table to a file
write.table(annotated, file='~/Desktop/cef_630.coverage.txt', sep='\t', quote=F, col.names=NA)
