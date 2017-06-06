
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
network_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/graph.tsv'
ko_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/KO_mapping.tsv'
layout_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/optimal_layout.tsv'
example_sim_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/test_distribution.txt'

# Read in metabolic network data
network <- read.delim(network_file, header=FALSE, sep='\t')
ko <- read.delim(ko_file, header=TRUE, sep='\t')
optimal_layout1 <- as.matrix(read.delim(layout_file, header=TRUE, sep='\t', row.names=1))
rm(network_file, ko_file, layout_file)

# Format directed graph
raw_graph <- graph.data.frame(network, directed=TRUE)
rm(network)

# Simplify graph
simple_graph <- simplify(raw_graph)

# Decompose graph
decomp_whole_graph <- decompose.graph(simple_graph)

# Get largest component and get node information
largest_component <- which.max(sapply(decomp_whole_graph, vcount))
largest_whole_graph <- decomp_whole_graph[[largest_component]]
print(length(as.vector(grep('K', V(largest_whole_graph)$name, value=TRUE)))) # 404
print(length(as.vector(grep('C', V(largest_whole_graph)$name, value=TRUE)))) # 666

# Find strongly-connected components
largest_scc <- rownames(as.data.frame(clusters(largest_whole_graph, mode='strong')[1]))
rm(decomp_whole_graph, largest_component)

# Read in simulated score distribution
example_sim <- read.delim(example_sim_file, header=FALSE)[,1]
rm(example_sim_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format simulated distribution
score_density <- density(example_sim)
score_median <- median(example_sim)
example_sim <- sort(unique(example_sim))
score_n <- length(example_sim)
score_q <- 0.5
score_nq <- score_n * score_q
score_range <- 1.96 * sqrt(score_n * score_q * (1 - score_q))
score_j <- ceiling(score_nq - score_range)
score_k <- ceiling(score_nq + score_range)
lower_95 <- example_sim[score_j]
upper_95 <- example_sim[score_k]

#-------------------------------------------------------------------------------------------------------------------------------------#

# Determine some statistics about graph

# Print a summary of nodes and edges for entire graph
summary(simple_graph)
print(length(as.vector(grep('K', V(simple_graph)$name, value=TRUE))))
print(length(as.vector(grep('C', V(simple_graph)$name, value=TRUE))))

# Find degree centrality
graph_degree <- as.data.frame(degree(simple_graph, v=V(simple_graph), mode='all'))

# Calculate betweenness centrality
graph_betweenness <- as.data.frame(betweenness(simple_graph))

# Calculate closeness centrality
graph_closeness_in <- as.data.frame(closeness(simple_graph, vids=V(simple_graph), mode='in'))
graph_closeness_out <- as.data.frame(closeness(simple_graph, vids=V(simple_graph), mode='out'))
graph_closeness_total <- as.data.frame(closeness(simple_graph, vids=V(simple_graph), mode='total'))

# Merge characteristic tables
graph_topology <- merge(graph_degree, graph_betweenness, by='row.names')
rownames(graph_topology) <- graph_topology$Row.names
graph_topology$Row.names <- NULL
graph_topology <- merge(graph_topology, graph_closeness_total, by='row.names')
rownames(graph_topology) <- graph_topology$Row.names
colnames(graph_topology) <- c('KEGG_ID','Degree_centrality','Betweenness_centrality','Closeness_centrality')

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
colnames(substrate_topology)[5] <- 'Compound_name'
enzyme_topology <- subset(graph_topology, grepl('K', graph_topology$KEGG_ID))
rownames(enzyme_topology) <- enzyme_topology$KEGG_ID
enzyme_topology <- merge(enzyme_topology, kegg_enzyme, by='row.names')
enzyme_topology$Row.names <- NULL
colnames(enzyme_topology)[5] <- 'Compound_name'
rm(kegg_substrate, kegg_enzyme)

# Get top ranking substrates from each
top_enzyme_DC <- as.vector(enzyme_topology[order(-enzyme_topology$Degree_centrality),][1:10,5])
top_enzyme_BC <- as.vector(enzyme_topology[order(-enzyme_topology$Betweenness_centrality),][1:10,5])
top_enzyme_CC <- as.vector(enzyme_topology[order(-enzyme_topology$Closeness_centrality),][1:10,5])
top_substrate_DC <- as.vector(substrate_topology[order(-substrate_topology$Degree_centrality),][1:10,5])
top_substrate_BC <- as.vector(substrate_topology[order(-substrate_topology$Betweenness_centrality),][1:10,5])
top_substrate_CC <- as.vector(substrate_topology[order(-substrate_topology$Closeness_centrality),][1:10,5])
shared_enzyme <- intersect(top_enzyme_BC, top_enzyme_CC)
shared_substrate <- intersect(top_substrate_BC, top_substrate_CC)

# Write tables to files, ranked by betweenness
table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/Table_S2_enzyme.tsv'
write.table(enzyme_topology, file=table_file, quote=FALSE, sep='\t', row.names=FALSE)
table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/Table_S2_substrate.tsv'
write.table(substrate_topology, file=table_file, quote=FALSE, sep='\t', row.names=FALSE)
rm(table_file, graph_topology, enzyme_topology, substrate_topology)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Transform largest component graph for plotting

# Get nodes for visualization
ko_simp <- as.vector(grep('K', V(largest_whole_graph)$name, value=TRUE))
substrate_simp <- as.vector(grep('C', V(largest_whole_graph)$name, value=TRUE))
nodes <- c(ko_simp, substrate_simp)
summary(largest_whole_graph)
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
node_size <- node_size[match(V(largest_whole_graph)$name, node_size$nodes),]
node_size <- setNames(as.numeric(node_size[,2]), as.character(node_size[,1]))
V(largest_whole_graph)$size <- as.matrix(node_size)
rm(raw_graph, ko, ko_simp, substrate_simp, node_size, mappings, nodes, simple_graph)

# Color graph
V(largest_whole_graph)$color <- ifelse(grepl('K', V(largest_whole_graph)$name), adjustcolor('firebrick3', alpha.f=0.6), 'blue3') # Color nodes
E(largest_whole_graph)$color <- 'gray15' # Color edges

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
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_3.pdf'
pdf(file=plot_file, width=7, height=12)
layout(matrix(c(1,2,
                3,3,
                3,3,
                4,4), nrow=4, ncol=2, byrow=TRUE))

#---------------------------------------#

# Large component of C. difficile 630 graph
par(mar=c(0,1,1,1), xpd=TRUE)
plot(largest_whole_graph, vertex.label=NA, layout=optimal_layout1,
     edge.arrow.size=0.25, edge.arrow.width=0.4, vertex.frame.color='black')
rect(xleft=0.4, ybottom=-0.7, xright=0.75, ytop=-0.4, lwd=2) # box

text(x=-1, y=1, 'A', cex=2)

#---------------------------------------#

# Metabolite score calculation
plot.new()
text(x=0.2, y=0.65, expression(paste(bold('( I ) '), mu[i]) == paste(bgroup('(',frac(Sigma * t[i], italic(n) * (e[o])),')'))), cex = 1.8)
text(x=0.7, y=0.65, expression(paste(bold('( II ) '), mu[o]) == paste(bgroup('(',frac(Sigma * t[o], italic(n) * (e[i])),')'))), cex = 1.8)
text(x=0.46, y=0.35, expression(paste(bold('( III ) '), Score(m)) == paste(log[2], '(', mu[i], ' - ', mu[o], ')')), cex = 1.8)

#---------------------------------------#

# Example network and score calculation
par(mar=c(0,0,0,1), xpd=FALSE)
plot(network, vertex.label=NA, layout=optimal_layout2, vertex.frame.color='black', xlim=c(-1.2,1.2), ylim=c(-1.2,1.2))
text(x=-1.5, y=1.35, 'B', cex=2)

text(x=-0.93, y=1.11, labels='Sucrose Glucohydrolase', font=2, cex=1.3) # Enzyme 1 name
text(x=-1, y=1, labels='7', col='white', cex=1.3) # Enzyme 1 transcription
text(x=-0.5, y=-1.3, labels='Fructose Phosphotransferase', font=2, cex=1.3) # Enzyme 2 name
text(x=-0.5, y=-1, labels='98', col='white', cex=2.4) # Enzyme 2 transcription
text(x=1, y=0.75, labels='Fructose Dehydrogenase', font=2, cex=1.3) # Enzyme 3
text(x=0.99, y=0.44, labels='165', col='white', cex=2.6) # Enzyme 3 transcription
text(x=-0.165, y=0.145, 'm', col='white', cex=2.1) # Substrate node label

text(x=c(-0.8,-0.8), y=c(0.15,0.05), labels=c('D-Fructose','Score = 7.0'), cex=1.5, font=c(2,1))
segments(x0=-1.15, y0=0, x1=-0.45, y1=0, lwd=2)

legend(x=0.6, y=1.28, legend=c('Enzyme node', 'Metabolite node'),
       pt.bg=c('firebrick3', 'blue3'), col='black', pch=21, pt.cex=3.3, cex=1.7, bty='n')
text(x=-0.5, y=-1.14, expression(t[i]), col='white', cex=2) # labeled transcription for input reactions
text(x=0.99, y=0.28, expression(t[i]), col='white', cex=2)
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
text(x=1.13, y=-0.4, 'More likely consumed', cex=1.35)
text(x=0.15, y=-0.4, 'More likely produced', cex=1.35)
text(x=0.63, y=-0.75, 'Metabolite Score', cex=1.5, font=2)

# Continuation arrows
Arrows(x0=-1.38, y0=1, x1=-1.12, y1=1, arr.type='curved', arr.length=0.45, arr.width=0.25, col='gray15')
Arrows(x0=1.25, y0=0.44, x1=1.5, y1=0.44, arr.type='curved', arr.length=0.45, arr.width=0.25, col='gray15')
Arrows(x0=-0.725, y0=-1, x1=-1.32, y1=-1, arr.type='curved', arr.length=0.45, arr.width=0.25, col='gray15')

# Draw box
rect(xleft=-1.38, ybottom=-1.4, xright=1.56, ytop=1.3, lwd=2)

#---------------------------------------#

# Example simulated distribution

par(mar=c(5,6,1,3), las=1, mgp=c(3.5,1,0), xaxs='i', yaxs='i', xpd=FALSE)
plot(score_density, xlim=c(-12,12), ylim=c(0,0.14), main='', xaxt='n', cex.lab=1.7, cex.axis=1.1,
     xlab='Simulated Metabolite Score', ylab='Score Density') 
axis(side=1, at=seq(-12,12,4), labels=seq(-12,12,4), cex.axis=1.2, tck=-0.03)
axis(side=1, at=c(-12:12), tck=-0.015, labels=FALSE)
polygon(score_density, col='gray80', border='black', lwd=1.5) 
abline(v=score_median, lty=2, lwd=2, col='black') # Median
abline(v=c(lower_95,upper_95), lty=2, lwd=2, col='red') # 0.95 Confidence Interval
abline(v=7, lwd=2, col='blue') # Measured Score
legend('topleft', legend=c('Simulated Median','Simulated 95% Confidence Interval','Measured Metabolite Score'), 
       col=c('black','red','blue'), lty=c(2,2,1), lwd=2, bg='white', cex=1.1)

mtext('C', side=2, line=2, las=2, adj=2.5, padj=-6.5, cex=1.3)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

