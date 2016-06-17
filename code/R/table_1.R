deps <- c('pander');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

file_name <- '~/Desktop/Jenior_812/bipartite_graph.txt'

abx_table <- read.table(file_name, header = T, sep = '\t')
