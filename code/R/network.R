
# Load igraph package
#install.packages("igraph")
library(igraph)

# Define variables
file_name <- '~/Desktop/Matt/data/cdf_bipartite/cefoperazone_630.bipartite.files/cefoperazone_630.bipartite.graph'
nodes_1 <- '~/Desktop/Matt/data/cdf_bipartite/cefoperazone_630.bipartite.files/cefoperazone_630.compound.lst'
nodes_1_label <- 'Substrate'
nodes_2 <- '~/Desktop/Matt/data/cdf_bipartite/cefoperazone_630.bipartite.files/cefoperazone_630.enzyme.lst'
nodes_2_label <- 'KEGG Ortholog'
figure_file <- '~/Desktop/Matt/data/cdf_bipartite/cefoperazone_630.bipartite.files/Cdifficile630.bipartite.scc.pdf'

/Users/pschloss/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.txt

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

# When needed, use this pallete I made
final_colors <- c("gold1", "orangered1", "aquamarine3", "firebrick", "forestgreen", "blue3", 
                  "mediumorchid2", "violetred4", "mediumpurple4", "dodgerblue3", "goldenrod3", "chartreuse3")

# Pick to different colors at random
colors <- final_colors[sample(1:length(color_palette), 2, replace=F)]

# Format data for plotting
V(largest.simple.graph)$size <- 3 # Node size
V(largest.simple.graph)$color <- ifelse(V(largest.simple.graph)$name %in% node_group_2, color_palette[2], color_palette[1]) # Color nodes
E(largest.simple.graph)$color <- 'gray15' # Color edges

# Plot the network
pdf(file=figure_file, width=10, height=7)
par(mar=c(0,0,0,0), font=2)
plot(largest.simple.graph, vertex.label = NA, layout = layout.graphopt,
     edge.arrow.size = 0.5, edge.arrow.width = 0.8, vertex.frame.color = 'black')
legend('bottomleft', legend=c(nodes_1_label, nodes_2_label), 
       pt.bg=c('red', 'blue'), col='black', pch=21, pt.cex=3, cex=1.5, bty = "n")
dev.off()







