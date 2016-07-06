
deps <- c('vegan', 'wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

# Import and format transcript mapping data
cefoperazone_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/cefoperazone_630.RNA_reads2cdf630.norm.annotated.txt'
clindamycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/clindamycin_630.RNA_reads2cdf630.norm.annotated.txt'
streptomycin_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/streptomycin_630.RNA_reads2cdf630.norm.annotated.txt'
germfree_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/cdifficile630/germfree.RNA_reads2cdf630.norm.annotated.txt'

# Load in data
cefoperazone <- read.delim(cefoperazone_file, sep='\t', header=FALSE, row.names=1)
colnames(cefoperazone) <- c('Cefoperazone', 'ko', 'gene', 'pathway')
clindamycin <- read.delim(clindamycin_file, sep='\t', header=FALSE, row.names=1)
colnames(clindamycin) <- c('Clindamycin', 'ko', 'gene', 'pathway')
streptomycin <- read.delim(streptomycin_file, sep='\t', header=FALSE, row.names=1)
colnames(streptomycin) <- c('Streptomycin', 'ko', 'gene', 'pathway')
germfree <- read.delim(germfree_file, sep='\t', header=FALSE, row.names=1)
colnames(germfree) <- c('Germfree', 'ko', 'gene', 'pathway')
rm(cefoperazone_file, clindamycin_file, streptomycin_file, germfree_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data for merging
cefoperazone$ko <- NULL
cefoperazone$gene <- NULL
cefoperazone$pathway <- NULL
clindamycin$ko <- NULL
clindamycin$gene <- NULL
clindamycin$pathway <- NULL
streptomycin$ko <- NULL
streptomycin$gene <- NULL
streptomycin$pathway <- NULL

# Merge tables
combined_mapping <- merge(cefoperazone, clindamycin, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
combined_mapping <- merge(combined_mapping, streptomycin, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
combined_mapping <- merge(combined_mapping, germfree, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
rm(cefoperazone, clindamycin, streptomycin, germfree)

# Rarefy mappings
sub_size <- round(min(colSums(combined_mapping[,1:4])) * 0.9) # 30645
combined_mapping$Cefoperazone <- t(rrarefy(combined_mapping$Cefoperazone, sample=sub_size))
combined_mapping$Clindamycin <- t(rrarefy(combined_mapping$Clindamycin, sample=sub_size))
combined_mapping$Streptomycin <- t(rrarefy(combined_mapping$Streptomycin, sample=sub_size))
combined_mapping$Germfree <- t(rrarefy(combined_mapping$Germfree, sample=sub_size))
rm(sub_size)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Select genes associated with lifecycle and virulence

# Sporulation associated genes
spo <- rbind(subset(combined_mapping, grepl('spo.....;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('spo....;', combined_mapping$gene)),
             subset(combined_mapping, grepl('spo...;', combined_mapping$gene)),
             subset(combined_mapping, grepl('spo..;', combined_mapping$gene)))
sig <- subset(combined_mapping, grepl('sig.;', combined_mapping$gene))
sip <- subset(combined_mapping, grepl('sip.;', combined_mapping$gene))
spm <- subset(combined_mapping, grepl('spm.;', combined_mapping$gene))
oxa <- subset(combined_mapping, grepl('oxa..;', combined_mapping$gene))
cot <- subset(combined_mapping, grepl('cot...;', combined_mapping$gene))
bcl <- subset(combined_mapping, grepl('bcl..;', combined_mapping$gene))
ssp <- subset(combined_mapping, grepl('ssp.;', combined_mapping$gene))
CD630_01250 <- subset(combined_mapping, grepl('CD630_01250', rownames(combined_mapping)))
CD630_01290 <- subset(combined_mapping, grepl('CD630_01290', rownames(combined_mapping)))
CD630_05960 <- subset(combined_mapping, grepl('CD630_05960', rownames(combined_mapping)))
CD630_10450 <- subset(combined_mapping, grepl('CD630_10450', rownames(combined_mapping)))
CD630_11680 <- subset(combined_mapping, grepl('CD630_11680', rownames(combined_mapping)))
CD630_12210 <- subset(combined_mapping, grepl('CD630_12210', rownames(combined_mapping)))
CD630_12900 <- subset(combined_mapping, grepl('CD630_12900', rownames(combined_mapping)))
CD630_12980 <- subset(combined_mapping, grepl('CD630_12980', rownames(combined_mapping)))
CD630_13210 <- subset(combined_mapping, grepl('CD630_13210', rownames(combined_mapping)))
CD630_14330 <- subset(combined_mapping, grepl('CD630_14330', rownames(combined_mapping)))
CD630_15110 <- subset(combined_mapping, grepl('CD630_15110', rownames(combined_mapping)))
CD630_16130 <- subset(combined_mapping, grepl('CD630_16130', rownames(combined_mapping)))
CD630_21440 <- subset(combined_mapping, grepl('CD630_21440', rownames(combined_mapping)))
CD630_23760 <- subset(combined_mapping, grepl('CD630_23760', rownames(combined_mapping)))
CD630_26410 <- subset(combined_mapping, grepl('CD630_26410', rownames(combined_mapping)))
CD630_26850 <- subset(combined_mapping, grepl('CD630_26850', rownames(combined_mapping)))
CD630_34640 <- subset(combined_mapping, grepl('CD630_34640', rownames(combined_mapping)))
CD630_35630 <- subset(combined_mapping, grepl('CD630_35630', rownames(combined_mapping)))
CD630_35690 <- subset(combined_mapping, grepl('CD630_35690', rownames(combined_mapping)))
sporulation_genes <- rbind(spo, sig, sip, spm, oxa, cot, bcl, ssp, CD630_01250, CD630_01290, 
                           CD630_05960, CD630_10450, CD630_11680, CD630_12210, CD630_12900, 
                           CD630_12980, CD630_13210, CD630_14330, CD630_15110, CD630_16130, 
                           CD630_21440, CD630_23760, CD630_26410, CD630_26850, CD630_34640, 
                           CD630_35630, CD630_35690)
sporulation_sums <- colSums(sporulation_genes[,c(1:4)])
sporulation_genes$pathway <- NULL
sporulation_genes$ko <- NULL
rm(spo, sig, sip, spm, oxa, cot, bcl, ssp, CD630_01250, CD630_01290, 
      CD630_05960, CD630_10450, CD630_11680, CD630_12210, CD630_12900, 
      CD630_12980, CD630_13210, CD630_14330, CD630_15110, CD630_16130, 
      CD630_21440, CD630_23760, CD630_26410, CD630_26850, CD630_34640, 
      CD630_35630, CD630_35690)


# Toxin associated genes
tcd <- subset(combined_mapping, grepl('tcd.;', combined_mapping$gene))
cdt <- subset(combined_mapping, grepl('cdt.;', combined_mapping$gene))
toxin_genes <- rbind(tcd, cdt)
toxin_sums <- colSums(toxin_genes[,c(1:4)])
toxin_genes$pathway <- NULL
toxin_genes$ko <- NULL
rm(tcd, cdt, combined_mapping)

# Assemble full gene table
select_mapping <- rbind(sporulation_genes, toxin_genes)
select_mapping[select_mapping == 0] <- 1
transformed_mapping <- log10(select_mapping[,c(1:4)])
transformed_mapping$gene <- select_mapping$gene
rm(select_mapping, sporulation_genes, toxin_genes)

# Assemble summary table
sum_mapping <- log10(rbind(sporulation_sums, toxin_sums))
rm(sporulation_sums, toxin_sums)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plot
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_3.pdf'
darjeeling <- wes_palette("Darjeeling")
pdf(file=plot_file, width=14, height=12)
layout(matrix(c(1,2), nrow=2, ncol=1, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Plot all genes
par(las=1, mar=c(8,4,1,1))
barplot(t(transformed_mapping[,c(1:4)]), col=c(darjeeling[1],darjeeling[2],darjeeling[4],darjeeling[5]), 
        beside=TRUE, xaxt='n', yaxt='n', ylab='Transcript Abundance (Log10)', ylim=c(0,3))
box()
axis(side=2, at=c(1:4), parse(text=paste(rep(10,3), '^', seq(1,3,1), sep='')), tick=TRUE, las=1)
abline(h=c(1:4), lty=2)
legend('topleft', legend=c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Gnotobiotic'), pt.cex=2, bty='n',
       pch=22, col='black', pt.bg=c(darjeeling[1],darjeeling[2],darjeeling[4],darjeeling[5]), ncol=2)
mtext('A', side=2, line=2, las=2, adj=2, padj=-6.2, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Plot catagory sums
par(las=1, mar=c(3,4,1,1))
barplot(t(sum_mapping), col=c(darjeeling[1],darjeeling[2],darjeeling[4],darjeeling[5]), 
        beside=TRUE, xaxt='n', yaxt='n', ylab='Transcript Abundance (Log10)', ylim=c(0,4))
box()
axis(side=2, at=c(1:4), parse(text=paste(rep(10,4), '^', seq(1,4,1), sep='')), tick=TRUE, las=1)
abline(h=c(1:4), lty=2)
legend('topleft', legend=c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Gnotobiotic'), pt.cex=2, bty='n',
       pch=22, col='black', pt.bg=c(darjeeling[1],darjeeling[2],darjeeling[4],darjeeling[5]), ncol=2)
text(x=c(3.5,8.5), y=c(-0.2,-0.2), labels=c('Sporulation','Toxin Production'), adj=1, xpd=TRUE, cex=1.4)
mtext('B', side=2, line=2, las=2, adj=2, padj=-6.2, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
rm(transformed_mapping, sum_mapping, darjeeling, plot_file)

