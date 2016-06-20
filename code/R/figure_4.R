
deps <- c('vegan', 'igraph', 'ggplot2');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

# Define variables
network_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/bipartite_graph.txt'
ko_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/cefoperazone_630.RNA_reads2cdf630.norm.ko.txt'
substrate_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/'

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in metabolic network data
network <- read.table(network_file, header=FALSE, sep='\t')
ko <- read.table(ko_file, header=FALSE, sep='\t')
rm(network_file, ko_file)

# Format directed graph
raw_graph <- graph.data.frame(network, directed=TRUE)
rm(network)

# Remove loops and multiple edges to make visualzation easier
simple_graph <- simplify(raw_graph)
rm(raw_graph)

# Decompose graph
decomp_simple_graph <- decompose.graph(simple_graph)
rm(simple_graph)

# Get largest component
largest_component <- which.max(sapply(decomp_simple_graph, vcount))
largest_simple_graph <- decomp_simple_graph[[largest_component]]
rm(decomp_simple_graph, largest_component)

# Remove zeroes so transformation doesn't return negative infinity
ko[,2][ko[,2] == 0] <- 1
ko[,2] <- log10(ko[,2])
ko[,2][ko[,2] < 0.4] <- 0.4

# Scale points by number of transcripts mapped
ko_simp <- as.vector(grep('K', V(largest_simple_graph)$name, value=TRUE))
substrate_simp <- as.vector(grep('C', V(largest_simple_graph)$name, value=TRUE))
nodes <- c(ko_simp, substrate_simp)
ko <- as.data.frame(subset(ko, ko[,1] %in% ko_simp))
ko <- ko[match(ko_simp, ko$V1),]
mappings <- c(as.vector(ko[,2] * 5), rep(2, length(substrate_simp)))
node_size <- cbind.data.frame(nodes, mappings)
node_size <- node_size[match(V(largest_simple_graph)$name, node_size$nodes),]
node_size <- setNames(as.numeric(node_size[,2]), as.character(node_size[,1]))
V(largest_simple_graph)$size <- as.matrix(node_size)
rm(ko, ko_simp, substrate_simp, node_size, mappings, nodes)
#V(largest_simple_graph)$size <- degree(largest_simple_graph) * 5 # Scale by degree!

# Color graph
V(largest_simple_graph)$color <- ifelse(grepl('K', V(largest_simple_graph)$name), adjustcolor('firebrick3', alpha.f=0.6), 'blue3') # Color nodes
E(largest_simple_graph)$color <- 'gray15' # Color edges

# Calculate optimal layout
#optimal <- layout.graphopt(graph=largest_simple_graph, niter=1000, charge=0.001, mass=50, spring.length=0, spring.constant=1)
#optimal <- layout.kamada.kawai(graph=largest_simple_graph)
#optimal <- layout.fruchterman.reingold(graph=largest_simple_graph)
optimal <- layout_nicely(graph=largest_simple_graph, dim=2)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_4.pdf'
pdf(file=plot_file, width=6, height=12)
layout(matrix(c(1,
                2), 
              nrow=2, ncol=1, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 4A

# Plot the large component of the graph
par(mar=c(0,0,0,0))
plot(largest_simple_graph, vertex.label=NA, layout=optimal,
     edge.arrow.size=0.5, edge.arrow.width=0.8, vertex.frame.color='black')
legend('bottomright', legend=c('KEGG Ortholog', 'Enzyme Substrate'), 
       pt.bg=c('firebrick3', 'blue3'), col='black', pch=21, pt.cex=2.3)
text(x=-1.5, y=1, labels='A', cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 4B

par(mar=c(0,0,0,0))

x <- mtcars[order(mtcars$mpg),] # sort by mpg
x$cyl <- factor(x$cyl) # it must be a factor
x$color[x$cyl==4] <- "red"
x$color[x$cyl==6] <- "blue"
x$color[x$cyl==8] <- "darkgreen"	
dotchart(x$mpg,labels=row.names(x),cex=.7,groups= x$cyl,
         main="Gas Milage for Car Models\ngrouped by cylinder",
         xlab="Miles Per Gallon", gcolor="black", color=x$color)


# Compound importance dotplot  , all 3 on same plot - 3 different colors of dots






dev.off()







