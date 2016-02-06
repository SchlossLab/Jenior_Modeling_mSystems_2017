
# Load igraph package
source("https://bioconductor.org/biocLite.R")
biocLite("pathview")





test_dist <- read.delim('~/Desktop/cefoperazone_630.mapped2cdf630.ko_expression.lst', sep='\t', header=FALSE, row.names=1)
