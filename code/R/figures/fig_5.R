
# Load dependencies
deps <- c('vegan', 'igraph', 'ggplot2', 'shape');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
set.seed(42)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Define file variables for network plot
network_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/bipartite_graph.tsv'
ko_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/cefoperazone_630.mapping.tsv'
layout_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/optimal_layout.tsv'

# Read in metabolic network data
network <- read.delim(network_file, header=FALSE, sep='\t')
ko <- read.delim(ko_file, header=TRUE, sep='\t')
optimal_layout1 <- as.matrix(read.delim(layout_file, header=TRUE, sep='\t', row.names=1))
rm(network_file, ko_file, layout_file)

# Format directed graph
raw_graph <- graph.data.frame(network, directed=TRUE)
rm(network)

# Decompose graph
decomp_whole_graph <- decompose.graph(raw_graph)

# Get largest component and get node information
largest_component <- which.max(sapply(decomp_whole_graph, vcount))
largest_whole_graph <- decomp_whole_graph[[largest_component]]
print(length(as.vector(grep('K', V(largest_whole_graph)$name, value=TRUE)))) # 404
print(length(as.vector(grep('C', V(largest_whole_graph)$name, value=TRUE)))) # 666

# Find strongly-connected components
largest_scc <- rownames(as.data.frame(clusters(largest_whole_graph, mode='strong')[1]))
rm(decomp_whole_graph, largest_component)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Determine some statistics about graph

# Print a summary of nodes and edges for entire graph
summary(raw_graph)
print(length(as.vector(grep('K', V(raw_graph)$name, value=TRUE))))
print(length(as.vector(grep('C', V(raw_graph)$name, value=TRUE))))

# Find degrees of nodes
graph_indegree <- as.data.frame(degree(raw_graph, v=V(raw_graph), mode='in'))
graph_outdegree <- as.data.frame(degree(raw_graph, v=V(raw_graph), mode='out'))
graph_undirected <- as.data.frame(degree(raw_graph, v=V(raw_graph), mode='all'))

# Calculate betweensness of entrire graph
graph_betweenness <- as.data.frame(betweenness(raw_graph))

# Calculate closeness of entrire graph
graph_closeness <- as.data.frame(closeness(raw_graph, vids=V(raw_graph), mode='all'))

# Calculate closeness of largest, strongly-connected component
occi <- as.data.frame(closeness(largest_whole_graph, vids=largest_scc, mode='all')) # Ma (2003)

# Merge characteristic tables
graph_topology <- merge(graph_outdegree, graph_indegree, by='row.names')
rownames(graph_topology) <- graph_topology$Row.names
graph_topology$Row.names <- NULL
graph_topology <- merge(graph_topology, graph_undirected, by='row.names')
rownames(graph_topology) <- graph_topology$Row.names
graph_topology$Row.names <- NULL
graph_topology <- merge(graph_topology, graph_betweenness, by='row.names')
rownames(graph_topology) <- graph_topology$Row.names
graph_topology$Row.names <- NULL
graph_topology$HBLC <- graph_topology[,4] / graph_topology[,3] # Calculate betweensness topolgy ratio, Joy et. al. (2005)
graph_topology[,3] <- NULL
graph_topology <- merge(graph_topology, graph_closeness, by='row.names')
rownames(graph_topology) <- graph_topology$Row.names
graph_topology$Row.names <- NULL
graph_topology <- merge(graph_topology, occi, by='row.names', all=TRUE)
graph_topology[is.na(graph_topology)] <- 0
colnames(graph_topology) <- c('KEGG_ID','Outdegree','Indegree','Betweenness', 'HBLC', 'Closeness', 'OCCI')
rm(graph_indegree, graph_outdegree, graph_undirected, graph_betweenness, graph_closeness, occi, largest_scc, largest_whole_graph)

# Read in KEGG code translation files
kegg_substrate_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/kegg/compound_names.tsv'
kegg_enzyme_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/kegg/ko_names.tsv'
kegg_substrate <- read.delim(kegg_substrate_file, header=FALSE, sep='\t', row.names=1)
kegg_enzyme <- read.delim(kegg_enzyme_file, header=FALSE, sep='\t', quote='', row.names=1)
rm(kegg_substrate_file, kegg_enzyme_file)

# Subset for Enzymes and Substrates
substrate_topology <- subset(graph_topology, grepl('C', graph_topology$KEGG_ID))
rownames(substrate_topology) <- substrate_topology$KEGG_ID
substrate_topology <- merge(substrate_topology, kegg_substrate, by='row.names')
substrate_topology$Row.names <- NULL
colnames(substrate_topology)[8] <- 'Common_name'
substrate_topology <- substrate_topology[order(-substrate_topology$Betweenness),]
enzyme_topology <- subset(graph_topology, grepl('K', graph_topology$KEGG_ID))
rownames(enzyme_topology) <- enzyme_topology$KEGG_ID
enzyme_topology <- merge(enzyme_topology, kegg_enzyme, by='row.names')
enzyme_topology$Row.names <- NULL
colnames(enzyme_topology)[8] <- 'Common_name'
enzyme_topology <- enzyme_topology[order(-enzyme_topology$Betweenness),]
rm(graph_topology, kegg_substrate, kegg_enzyme)

# Write tables to files, ranked by betweenness
table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/substrate_topology.tsv'
write.table(substrate_topology, file=table_file, quote=FALSE, sep='\t', row.names=FALSE)
rm(table_file, substrate_topology)
table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/enzyme_topology.tsv'
write.table(enzyme_topology, file=table_file, quote=FALSE, sep='\t', row.names=FALSE)
rm(table_file, enzyme_topology)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Transform largest component graph for plotting

# Remove loops and multiple edges to make visualzation easier
simple_graph <- simplify(raw_graph)
decomp_simple_graph <- decompose.graph(simple_graph)
largest_component <- which.max(sapply(decomp_simple_graph, vcount))
largest_simple_graph <- decomp_simple_graph[[largest_component]]
ko_simp <- as.vector(grep('K', V(largest_simple_graph)$name, value=TRUE))
substrate_simp <- as.vector(grep('C', V(largest_simple_graph)$name, value=TRUE))
nodes <- c(ko_simp, substrate_simp)
summary(largest_simple_graph)
print(length(ko_simp))
print(length(substrate_simp))

# Remove zeroes so transformation doesn't return negative infinity
ko[,2][ko[,2] == 0] <- 1
ko[,2] <- log10(ko[,2])
ko[,2][ko[,2] < 0.4] <- 0.4

# Scale points by number of transcripts mapped
ko <- as.data.frame(subset(ko, ko[,1] %in% ko_simp))
ko <- ko[match(ko_simp, ko$KO_code),]
mappings <- c(as.vector(ko[,2] * 5), rep(2, length(substrate_simp)))
node_size <- cbind.data.frame(nodes, mappings)
node_size <- node_size[match(V(largest_simple_graph)$name, node_size$nodes),]
node_size <- setNames(as.numeric(node_size[,2]), as.character(node_size[,1]))
V(largest_simple_graph)$size <- as.matrix(node_size)
rm(raw_graph, ko, ko_simp, substrate_simp, node_size, mappings, nodes, simple_graph, decomp_simple_graph, largest_component)

# Color graph
V(largest_simple_graph)$color <- ifelse(grepl('K', V(largest_simple_graph)$name), adjustcolor('firebrick3', alpha.f=0.6), 'blue3') # Color nodes
E(largest_simple_graph)$color <- 'gray15' # Color edges

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up example network

# Define file variables for network plot and layout
network <- matrix(c('K01', 'C01',
                    'C01', 'K02',
                    'C01', 'K03'), nrow=3, ncol=2, byrow=TRUE)
node_size <- matrix(c('K01', '12.5',
                      'K02', '45',
                      'K03', '50',
                      'C01', '20'), nrow=4, ncol=2, byrow=TRUE)
optimal_layout2 <- matrix(c(-21.09826017, 22.1407060,
                            0.09077637, 0.2154631,
                            -8.32243732, -29.0949351,
                            29.67130628, 7.6231375), nrow=4, ncol=2, byrow=TRUE)

# Format directed graph
network <- graph.data.frame(network, directed=TRUE)

# Assign nodes sizes
node_size <- node_size[match(V(network)$name, node_size[,1]),]
node_size <- setNames(as.numeric(node_size[,2]), as.character(node_size[,1]))
V(network)$size <- as.matrix(node_size)
rm(node_size)

# Color graph
V(network)$color <- ifelse(grepl('K', V(network)$name), 'firebrick3', 'blue3') # Color nodes
E(network)$color <- 'gray15' # Color edges

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_5.pdf'
pdf(file=plot_file, width=10, height=7)
layout(matrix(c(1,2,2,
                3,2,2), nrow=2, ncol=3, byrow=TRUE))

#---------------------------------------#

# Large component of C. difficile 630 graph
par(mar=c(0,1,1,1))
plot(largest_simple_graph, vertex.label=NA, layout=optimal_layout1,
     edge.arrow.size=0.25, edge.arrow.width=0.4, vertex.frame.color='black')
# Draw box
rect(xleft=0.8, ybottom=0, xright=1.15, ytop=0.3, lwd=2)
text(x=-1, y=1, 'a', cex=2, font=2)

#---------------------------------------#

# Example network and importance calculation
par(mar=c(0,0,0,1))
plot(network, vertex.label=NA, layout=optimal_layout2, vertex.frame.color='black', xlim=c(-1.2,1.2), ylim=c(-1.2,1.2))
text(x=-1.5, y=1.5, 'b', cex=2, font=2)

text(x=-0.95, y=1.11, labels='dAdo Aminohydrolase', font=2, cex=1.3) # Enzyme 1 name
text(x=-1, y=1, labels='7', col='white', cex=1.3) # Enzyme 1 transcription
text(x=-0.5, y=-1.3, labels='ATP:dAdo 5\'-Phosphotransferase', font=2, cex=1.3) # Enzyme 2 name
text(x=-0.5, y=-1, labels='94', col='white', cex=2.4) # Enzyme 2 transcription
text(x=1, y=0.75, labels=expression(bold(paste('dAdo:', PO[4]^-3, ' Ribosyltransferase'))), cex=1.3) # Enzyme 3
text(x=0.99, y=0.44, labels='115', col='white', cex=2.6) # Enzyme 3 transcription
text(x=-0.165, y=0.145, 'm', col='white', cex=2.1) # Substrate node label

text(x=c(-0.8,-0.8), y=c(0.15,0.05), labels=c('Deoxyadenosine (dAdo)','Importance = 6.554'), cex=1.5, font=c(2,1)) # Compound & calculated importance
segments(x0=-1.15, y0=0, x1=-0.45, y1=0, lwd=2)

legend(x=0.81, y=1.3, legend=c('Enzyme node', 'Metabolite node'),
       pt.bg=c('firebrick3', 'blue3'), col='black', pch=21, pt.cex=3, cex=1.5)
text(x=-0.5, y=-1.12, expression(t[i]), col='white', cex=2) # labeled transcription for input reactions
text(x=0.99, y=0.32, expression(t[i]), col='white', cex=2)
text(x=-1, y=0.87, expression(t[o]), col='black', cex=2) # labeled transcription for output reactions
text(x=-0.6, y=0.7, expression(e[i]), col='black', cex=2) # labeled indegree
text(x=-0.4, y=-0.3, expression(e[o]), col='black', cex=2) # labeled outdegree
text(x=0.3, y=0.33, expression(e[o]), col='black', cex=2)

Arrows(x0=0.63, y0=-0.5, x1=0.12, y1=-0.5, lwd=3, arr.type='triangle', arr.length=0.5, arr.width=0.2) # Score explanation line
Arrows(x0=0.63, y0=-0.5, x1=1.14, y1=-0.5, lwd=3, arr.type='triangle', arr.length=0.5, arr.width=0.2)
segments(x0=0.63, y0=-0.45, x1=0.63, y1=-0.55, lwd=2)
text(x=0.63, y=-0.63, '0', cex=1.8)
text(x=0.12, y=-0.61, '-', cex=2.2)
text(x=1.14, y=-0.61, '+', cex=2.2)
text(x=1.15, y=-0.4, 'More likely consumed', cex=1.35)
text(x=0.14, y=-0.4, 'More likely released', cex=1.35)
text(x=0.63, y=-0.75, 'Importance Score', cex=1.5, font=2)

# Continuation lines
Arrows(x0=-1.38, y0=1, x1=-1.12, y1=1, arr.type='curved', arr.length=0.45, arr.width=0.25, col='gray15')
segments(x0=1.25, y0=0.44, x1=1.56, y1=0.44, col='gray15', lwd=1.1)
segments(x0=-1.38, y0=-1, x1=-0.725, y1=-1, col='gray15', lwd=1.1)

# Draw box
rect(xleft=-1.38, ybottom=-1.4, xright=1.56, ytop=1.3, lwd=2)

#---------------------------------------#

# Importance algorithm

plot.new()
text(x=0.28, y=0.9, expression(mu[i] == paste(bgroup('(',frac(Sigma * t[i], italic(n) * (e[o])),')'))), cex = 1.8)
text(x=0.8, y=0.9, expression(mu[o] == paste(bgroup('(',frac(Sigma * t[o], italic(n) * (e[i])),')'))), cex = 1.8)
text(x=0.58, y=0.6, expression(Importance(m) == paste(Log[2], '(', mu[i], ' - ', mu[o], ')')), cex = 1.8)
text(x=c(0.05,0.56,0.07), y=c(0.9,0.9,0.6), labels=c('(i)','(ii)','(iii)'), cex = 1.8, font=2)                                                                  

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
rm(optimal_layout1, optimal_layout2, largest_simple_graph, network, plot_file)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

