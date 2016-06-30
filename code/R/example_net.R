

deps <- c('vegan', 'igraph', 'ggplot2');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Define file variables for network plot
network <- matrix(c('K01', 'C01',
                    'C01', 'K02',
                    'C01', 'K03'), nrow=3, ncol=2, byrow=TRUE)
node_size <- matrix(c('K01', '1',
                      'K02', '2',
                      'K03', '2.4',
                      'C01', '1',
                      'C02', '1',
                      'C03', '1'), nrow=3, ncol=2, byrow=TRUE)

# Format directed graph
network <- graph.data.frame(network, directed=TRUE)

# Scale points by number of transcripts mapped
node_size <- node_size[match(V(network)$name, node_size[,1]),]
node_size <- setNames(as.numeric(node_size[,2]), as.character(node_size[,1]))
V(network)$size <- as.matrix(node_size)
rm(node_size)

# Color graph
V(network)$color <- ifelse(grepl('K', V(network)$name), adjustcolor('firebrick3', alpha.f=0.6), 'blue3') # Color nodes
E(network)$color <- 'gray15' # Color edges

# Calculate optimal layout
optimal <- layout.graphopt(graph=network, niter=1000, charge=0.001, mass=50, spring.length=0, spring.constant=1)

#------------------------------------------------------------------------------------------------------------------#

plot(network, vertex.label=NA, layout=optimal,
     edge.arrow.size=0.5, edge.arrow.width=0.8, vertex.frame.color='black')



