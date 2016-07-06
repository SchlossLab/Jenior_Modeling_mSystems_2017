
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

#-------------------------------------------------------------------------------------------------------------------------#

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

#-------------------------------------------------------------------------------------------------------------------------#

# Select genes associated with lifecycle and virulence

# Sporulation associated genes
spo <- rbind(subset(combined_mapping, grepl('spo.....;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('spo....;', combined_mapping$gene)),
             subset(combined_mapping, grepl('spo...;', combined_mapping$gene)),
             subset(combined_mapping, grepl('spo..;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('spo.;', combined_mapping$gene)), 
             subset(combined_mapping, grepl('spo;', combined_mapping$gene)))
sig <- rbind(subset(combined_mapping, grepl('sig.;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('sig;', combined_mapping$gene)))
sod <- rbind(subset(combined_mapping, grepl('sod.;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('sod;', combined_mapping$gene)))
sip <- rbind(subset(combined_mapping, grepl('sip.;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('sip;', combined_mapping$gene)))
spm <- rbind(subset(combined_mapping, grepl('spm.;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('spm;', combined_mapping$gene)))
yqf <- rbind(subset(combined_mapping, grepl('yqf.;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('yqf;', combined_mapping$gene)))
yab <- rbind(subset(combined_mapping, grepl('yab.;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('yab;', combined_mapping$gene)))
oxa <- rbind(subset(combined_mapping, grepl('oxa.;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('oxa;', combined_mapping$gene)))
cot <- rbind(subset(combined_mapping, grepl('cot.;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('cot;', combined_mapping$gene)))
bcl <- rbind(subset(combined_mapping, grepl('bcl.;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('bcl;', combined_mapping$gene)))
ssp <- rbind(subset(combined_mapping, grepl('ssp.;', combined_mapping$gene)), # 
             subset(combined_mapping, grepl('ssp;', combined_mapping$gene)))

# Toxin associated genes
tcd <- rbind(subset(combined_mapping, grepl('tcd.;', combined_mapping$gene)), # Toxin operon
             subset(combined_mapping, grepl('tcd;', combined_mapping$gene)))
cdt <- rbind(subset(combined_mapping, grepl('cdt.;', combined_mapping$gene)), # Binary toxin
             subset(combined_mapping, grepl('cdt;', combined_mapping$gene)))
rm(combined_mapping)

# Assemble final table
select_mapping <- rbind(spo, sig, sod, sip, spm, yqf, yab, oxa, cot, bcl, ssp, tcd, cdt)
rm(spo, sig, sod, sip, spm, yqf, yab, oxa, cot, bcl, ssp, tcd, cdt)
rownames(select_mapping) <- c('spo', 'sig', 'sod', 'sip', 'spm', 'yqf', 'yab', 'oxa', 'cot', 'bcl', 'ssp', 'tcd', 'cdt')
select_mapping[select_mapping == 0] <- 1
transformed_mapping <- log10(select_mapping)
rm(select_mapping)

#-------------------------------------------------------------------------------------------------------------------------#

# Set up plot
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_3.pdf'
darjeeling <- wes_palette("Darjeeling")
pdf(file=plot_file, width=7, height=14)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Plot the data
par(las=1, mar=c(8,4,1,1))
barplot(t(transformed_mapping), col=c(darjeeling[1],darjeeling[2],darjeeling[4],darjeeling[5]), 
        beside=TRUE, xaxt='n', yaxt='n', ylab='Transcript Abundance (Log10)', ylim=c(0,5))
box()
axis(side=2, at=c(1:4), parse(text=paste(rep(10,4), '^', seq(1,4,1), sep='')), tick=TRUE, las=1)
abline(h=c(1:4), lty=2)
legend('topleft', legend=c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Gnotobiotic'), pt.cex=2, bty='n',
       pch=22, col='black', pt.bg=c(darjeeling[1],darjeeling[2],darjeeling[4],darjeeling[5]), ncol=2)
text(x=seq(4,59,5), y=par()$usr[3]-0.03*(par()$usr[4]-par()$usr[3]),
     labels=rownames(transformed_mapping), srt=45, adj=1, xpd=TRUE, cex=0.8)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
rm(transformed_mapping, darjeeling, plot_file)

