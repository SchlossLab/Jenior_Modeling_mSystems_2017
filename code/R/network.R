
# Load igraph package
#install.packages("igraph")
library(igraph)

# Define variables
file_name <- '~/Desktop/cefoperazone_630.bipartite.graph'
nodes_1 <- '~/Desktop/cefoperazone_630.compound.lst'
nodes_1_label <- 'Substrate'
nodes_2 <- '~/Desktop/cefoperazone_630.enzyme.lst'
nodes_2_label <- 'KEGG Ortholog'
figure_file <- '~/Desktop/Cdifficile630.bipartite.scc.pdf'

# Read in data
graph.file <- read.table(file_name, header = F, sep = '\t')
node_group_1 <- as.vector(read.table(nodes_1, header = F, sep = '\t')$V1)
node_group_2 <- as.vector(read.table(nodes_2, header = F, sep = '\t')$V1)

# Format directed graph
raw.graph <- graph.data.frame(graph.file, directed = T)

# Remove loops and multiple edges to make visualzation easier
simple.graph <- simplify(raw.graph)

# Decompose graph
all.simple.graph <- decompose.graph(simple.graph)

# Get largest component
largest <- which.max(sapply(all.simple.graph, vcount))
largest.simple.graph <- all.simple.graph[[largest]]

# Format data for plotting
V(largest.simple.graph)$size <- 3 # Node size
V(largest.simple.graph)$color <- ifelse(V(largest.simple.graph)$name %in% node_group_2, "blue", "red") # Color nodes
E(largest.simple.graph)$color <- 'gray15' # Color edges

# Plot the network
pdf(file=figure_file, width=10, height=7)
par(mar=c(0,0,0,0), font=2)
plot(largest.simple.graph, vertex.label = NA, layout = layout.graphopt,
     edge.arrow.size = 0.5, edge.arrow.width = 0.8, vertex.frame.color = 'black')
legend('bottomleft', legend=c(nodes_2_label, nodes_1_label), 
       pt.bg=c('blue', 'red'), col='black', pch=21, pt.cex=3, cex=1.5, bty = "n")
dev.off()

