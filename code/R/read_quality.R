
install.packages('qrqc')
library('qrqc')
mock1 <- readSeqFile('~/Desktop/48670_AGGCAGAA-TATCCTCT_L001_R1_001.fastq.gz', type='fastq', quality='sanger')
qualPlot(mock1)


