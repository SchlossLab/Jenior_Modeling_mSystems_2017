# Data variables
cefoperazone_file <- '/Users/schloss/Desktop/Repos/presentations/annotated_cdf/cefoperazone_630.mapped2cdf630.annotated.txt'
streptomycin_file <- '/Users/schloss/Desktop/Repos/presentations/annotated_cdf/streptomycin_630.mapped2cdf630.annotated.txt'
clindamycin_file <- '/Users/schloss/Desktop/Repos/presentations/annotated_cdf/clindamycin_630.mapped2cdf630.annotated.txt'
germfree_file <- '/Users/schloss/Desktop/Repos/presentations/annotated_cdf/germfree.mapped2cdf630.annotated.txt'
output_file <- '/Users/schloss/Desktop/Repos/presentations/annotated_cdf/expression_figure.pdf'

# Read in data
cefoperazone <- read.delim(cefoperazone_file, sep='\t', row.names=2, header=F)
colnames(cefoperazone) <- c('transcripts', 'gene', 'ko', 'pathway_code', 'pathway_name', 'metadata1', 'metadata3')
streptomycin <- read.delim(streptomycin_file, sep='\t', row.names=2, header=F)
colnames(streptomycin) <- c('transcripts', 'gene', 'ko', 'pathway_code', 'pathway_name', 'metadata1', 'metadata3')
clindamycin <- read.delim(clindamycin_file, sep='\t', row.names=2, header=F)
colnames(clindamycin) <- c('transcripts', 'gene', 'ko', 'pathway_code', 'pathway_name', 'metadata1', 'metadata3')
germfree <- read.delim(germfree_file, sep='\t', row.names=2, header=F)
colnames(germfree) <- c('transcripts', 'gene', 'ko', 'pathway_code', 'pathway_name', 'metadata1', 'metadata3')

# Subset the data
cefoperazone_pathways <- subset(cefoperazone, pathway_name != 'metadata_unannotated')
streptomycin_pathways <- subset(streptomycin, pathway_name != 'metadata_unannotated')
clindamycin_pathways <- subset(clindamycin, pathway_name != 'metadata_unannotated')
germfree_pathways <- subset(germfree, pathway_name != 'metadata_unannotated')

# Pool data seperately
cefoperazone_pooled_pathways <- aggregate(transcripts ~ pathway_name, FUN = sum, data = cefoperazone_pathways)  
rownames(cefoperazone_pooled_pathways) <- cefoperazone_pooled_pathways[,1]
cefoperazone_pooled_pathways[,1] <- NULL
colnames(cefoperazone_pooled_pathways) <- 'Cefoperazone'
streptomycin_pooled_pathways <- aggregate(transcripts ~ pathway_name, FUN = sum, data = streptomycin_pathways)  
rownames(streptomycin_pooled_pathways) <- streptomycin_pooled_pathways[,1]
streptomycin_pooled_pathways[,1] <- NULL
colnames(streptomycin_pooled_pathways) <- 'Streptomycin'
clindamycin_pooled_pathways <- aggregate(transcripts ~ pathway_name, FUN = sum, data = clindamycin_pathways)  
rownames(clindamycin_pooled_pathways) <- clindamycin_pooled_pathways[,1]
clindamycin_pooled_pathways[,1] <- NULL
colnames(clindamycin_pooled_pathways) <- 'Clindamycin'
germfree_pooled_pathways <- aggregate(transcripts ~ pathway_name, FUN = sum, data = germfree_pathways)  
rownames(germfree_pooled_pathways) <- germfree_pooled_pathways[,1]
germfree_pooled_pathways[,1] <- NULL
colnames(germfree_pooled_pathways) <- 'Germfree'

# Merge datasets
abx_pathways <- merge(cefoperazone_pooled_pathways, streptomycin_pooled_pathways, by='row.names')
rownames(abx_pathways) <- abx_pathways$Row.names
abx_pathways$Row.names <- NULL
abx_pathways <- merge(abx_pathways, clindamycin_pooled_pathways, by='row.names')
rownames(abx_pathways) <- abx_pathways$Row.names
abx_pathways$Row.names <- NULL
abx_pathways <- merge(abx_pathways, germfree_pooled_pathways, by='row.names')
rownames(abx_pathways) <- abx_pathways$Row.names
abx_pathways$Row.names <- NULL

# Screen low abundance pathways and transform data
pick_abx_pathways <- abx_pathways[ which(rowSums(abx_pathways) > 150),]
pick_abx_pathways <- abx_pathways[ which(rowSums(pick_abx_pathways) < 50000),]
pick_abx_pathways[pick_abx_pathways == 0] <- 1
transformed_pick_abx_pathways <- log10(pick_abx_pathways)

# Filter into groups for manageable plots




# Plot the data
#pdf(file=spore_file, width=30, height=7)
par(las=1, mar=c(9,3,1,3), mgp=c(1,0,0))
barplot(t(transformed_pick_abx_pathways), col=c('red2', 'blue2', 'forestgreen', 'gold2'), 
        beside=TRUE, xlim=c(0,170), ylim=c(1,5), xaxt='n', yaxt='n', 
        ylab='KEGG Pathway', cex.lab=1.5)




#rotate 60 degrees, srt=60
text(seq(1.5,end_point,by=2), par("usr")[3]-0.25, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = paste(rownames(mtcars)), cex=0.65)

labelsY <- parse(text=paste(rep(10,5), '^', seq(5,9,1), sep=''))
axis(side=2, at=c(5:9), labelsY, tick=TRUE, cex.axis=1.2, las=1)




legend()

#dev.off()