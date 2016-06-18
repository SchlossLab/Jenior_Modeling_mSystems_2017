
deps <- c('vegan', 'igraph');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

# Define and check color palette
palette_plot <- function(col, border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}

# Define variables
network_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/meabolic_models/'
ko_score_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/meabolic_models/'
substrate_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/meabolic_models/'

# Read in data
network <- read.table(network_file, header=FALSE, sep='\t')
ko_score <- read.table(ko_score_file, header=FALSE, sep='\t')
substrate_file <- as.vector(read.table(substrate_file, header=FALSE, sep='\t')$V1)
rm(network_file, ko_score_file, substrate_file)

# Format directed graph
raw_graph <- graph.data.frame(network, directed=TRUE)

# Remove loops and multiple edges to make visualzation easier
simple_graph <- simplify(raw_graph)

# Get largest component
largest_component <- which.max(sapply(simple_graph, vcount))
largest_simple_graph <- simple.graph[[largest_component]]

# Format data for plotting
V(largest.simple.graph)$size <- # need these to scale by expression....where is that other script????
V(largest.simple.graph)$color <- ifelse(V(largest.simple.graph)$name %in% node_group_2, colors[2], colors[1]) # Color nodes
E(largest.simple.graph)$color <- 'gray15' # Color edges





#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up ployying environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_4.pdf'
pdf(file=plot_file, width=12, height=6)
layout(matrix(c(1,2,
                3,4), 
              nrow=2, ncol=2, byrow = TRUE))

# Plot the large component of th graph
par(mar=c(0,0,0,0))
plot(largest_simple_graph, vertex.label=NA, layout=layout.graphopt,
     edge.arrow.size=0.5, edge.arrow.width=0.8, vertex.frame.color='black')
legend('bottomleft', legend=c(nodes_1_label, nodes_2_label), 
       pt.bg=c('red', 'blue'), col='black', pch=21, pt.cex=3, cex=1.5, bty = "n")




















dev.off()







