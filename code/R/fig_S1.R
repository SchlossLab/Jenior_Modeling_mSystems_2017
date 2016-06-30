
# Set up environment
deps <- c('qrqc');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

read_file <- '~/Desktop/48670_AGGCAGAA-TATCCTCT_L001_R1_001.fastq.gz'
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figure_S1.pdf'

# Read in data
reads <- readSeqFile(read_file, type='fastq', quality='sanger')

# Generate figure
pdf(file=plot_file, width=12, height=10)
qualPlot(reads)

#Clean up
dev.off()
rm(read_file, plot_file, reads)
