
library(vegan)
library(stringr)
RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

cef_reads <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/cefoperazone_630.RNA_reads2cdf630.norm.ko.pick.txt', 
                    sep='\t', header=F)
colnames(cef_reads) <- c('ko', 'cefoperazone')
cef_reads <- aggregate(cefoperazone~ko,data=cef_reads,FUN=sum)
clinda_reads <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/clindamycin_630.RNA_reads2cdf630.norm.ko.pick.txt', 
                        sep='\t', header=F)
colnames(clinda_reads) <- c('ko', 'clindamycin')
clinda_reads <- aggregate(clindamycin~ko,data=clinda_reads,FUN=sum)
strep_reads <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/streptomycin_630.RNA_reads2cdf630.norm.ko.pick.txt', 
                        sep='\t', header=F)
colnames(strep_reads) <- c('ko', 'streptomycin')
strep_reads <- aggregate(streptomycin~ko,data=strep_reads,FUN=sum)
gf_reads <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/expression/germfree.RNA_reads2cdf630.norm.ko.pick.txt', 
                        sep='\t', header=F)
colnames(gf_reads) <- c('ko', 'germfree')
gf_reads <- aggregate(germfree~ko,data=gf_reads,FUN=sum)
reads <- merge(cef_reads, clinda_reads, by='ko')
reads <- merge(reads, strep_reads, by='ko')
reads <- merge(reads, gf_reads, by='ko')
rm(cef_reads, clinda_reads, strep_reads, gf_reads)
sub_size <- ceiling(min(colSums(reads[,2:5])) * 0.95)
reads[,2] <- as.vector(rrarefy(floor(reads[,2]), sample=sub_size))
reads[,3] <- as.vector(rrarefy(floor(reads[,3]), sample=sub_size))
reads[,4] <- as.vector(rrarefy(floor(reads[,4]), sample=sub_size))
reads[,5] <- as.vector(rrarefy(floor(reads[,5]), sample=sub_size))
rm(sub_size)

# K02469 = gyrA - DNA gyrase subunit A
# K03040 = rpoA - DNA-directed RNA polymerase subunit alpha
# K10670 = glycine_reductase_[EC:1.21.4.2]
# K10793 = D-proline_reductase_(dithiol)_PrdA_[EC:1.21.4.1]
# K02406 = flagellin
# K06334 = spore_coat_protein_JC
test <- reads
test[,1] <- NULL
reads$variance <- RowVar(test)
reads <- reads[order(-reads$variance),]
defs <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/kegg/ko_definition.list', 
                       sep='\t', header=FALSE, quote='')
colnames(defs) <- c('ko', 'definition')
defs$definition <- gsub(' ', '_', defs$definition)
reads <- merge(reads, defs, by='ko')
test <- reads[reads$ko %in% c('K02469','K03040','K10670','K10793','K02406','K06334'),]
write.table(test, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/variance_ko.tsv', row.names=FALSE, quote=FALSE, sep='\t')
rm(test)
reads$variance <- NULL
rm(defs)
cdf <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/kegg/cdf_ko.list', 
                   sep='\t', header=FALSE)
colnames(cdf) <- c('CD630_gene', 'ko')
cdf$ko <- cdf$ko %>% str_replace('ko:', '')
cdf <- aggregate(CD630_gene~ko, cdf, paste, collapse=",")
reads <- merge(reads, cdf, by='ko', all.x=TRUE, all.y=FALSE)
rm(cdf)


ko_path <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/kegg/ko.list', 
                      sep='\t', header=FALSE)
ko_path$V3 <- NULL
colnames(ko_path) <- c('pathway','ko')
ko_path$pathway <- ko_path$pathway %>% str_replace('path:ko', '')
ko_path$ko <- ko_path$ko %>% str_replace('ko:', '')
pathway <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/kegg/pathway.list', 
                  sep='\t', header=FALSE)
colnames(pathway) <- c('pathway','pathway_name')
pathway$pathway <- pathway$pathway %>% str_replace('path:', '')
pathway <- merge(pathway, ko_path, by='pathway', all.x=TRUE, all.y=FALSE)
rm(ko_path)
pathway$pathway <- NULL
pathway <- aggregate(pathway_name~ko, pathway, paste, collapse=",")
pathway <- pathway[grep('K', pathway$ko),]
reads <- merge(reads, pathway, by='ko', all.x=TRUE, all.y=FALSE)
rm(pathway)
reaction <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/kegg/ko_reaction.list', 
                  sep='\t', header=FALSE)
colnames(reaction) <- c('ko', 'reaction')
reaction$ko <- reaction$ko %>% str_replace('ko:', '')
reaction$reaction <- reaction$reaction %>% str_replace('rn:', '')
reaction_map <- read.delim('~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/kegg/reaction_mapformula.lst', 
                           sep=':', header=FALSE)
reaction_map$V2 <- NULL
colnames(reaction_map) <- c('reaction', 'formulas')
reaction_map$formulas <- gsub(' ', '', reaction_map$formulas)
reaction_map <- unique(reaction_map)
reaction_map <- aggregate(formulas~reaction, reaction_map, paste, collapse=",")
reaction <- merge(reaction, reaction_map, by='reaction', all.x=TRUE, all.y=FALSE)
rm(reaction_map)
reaction <- reaction[complete.cases(reaction),]
reaction <- aggregate(cbind(reaction,formulas)~ko, reaction, paste, collapse=",")
reads <- merge(reads, reaction, by='ko', all.x=TRUE, all.y=FALSE)
rm(reaction)
colnames(reads) <- c('KO','cefoperazone-pretreated','clindamycin-pretreated','streptomycin-pretreated',
                     'ex-germfree','gene_name','CD630_genes','pathways','enzymatic_reactions','chemical_formulas')
write.table(reads, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/kegg/ko_mapping.tsv', row.names=FALSE, quote=FALSE, sep='\t')
rm(reads)

