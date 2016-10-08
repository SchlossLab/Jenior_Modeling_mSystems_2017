
# Load dependencies
deps <- c('vegan', 'igraph', 'ggplot2', 'shape', 'wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
set.seed(42)

# Define function for mormatting growth curves for statistical analysis
format_curve <- function(raw_exp_data, exp_group, raw_control_data){
  formatted_data <- c()
  control_data <- c()
  for (time in 1:49){
    temp_exp <- cbind(exp_group, time, raw_exp_data[time,])
    formatted_data <- rbind(formatted_data, temp_exp)
    temp_control <- cbind('control', time, raw_control_data[time,])
    control_data <- rbind(control_data, temp_control)
  }
  formatted_data <- as.data.frame(rbind(control_data, formatted_data))
  colnames(formatted_data) <- c('substrate','time','od')
  formatted_data$od <- as.numeric(as.character(formatted_data$od))
  
  return(formatted_data)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

# Define file variables for network plot
network_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/bipartite_graph.txt'
ko_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/cefoperazone_630.mapping.txt'
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
#V(largest_simple_graph)$size <- degree(largest_simple_graph) * 5 # Scale by degree!

# Color graph
V(largest_simple_graph)$color <- ifelse(grepl('K', V(largest_simple_graph)$name), adjustcolor('firebrick3', alpha.f=0.6), 'blue3') # Color nodes
E(largest_simple_graph)$color <- 'gray15' # Color edges

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up example network

# Define file variables for network plot and layout
network <- matrix(c('K01', 'C01',
                    'C01', 'K02',
                    'C01', 'K03'), nrow=3, ncol=2, byrow=TRUE)
node_size <- matrix(c('K01', '10',
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

# Read in substrate importance data
cef_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/cefoperazone_630.monte_carlo.score.txt'
clinda_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/clindamycin_630.bipartite.files/clindamycin_630.monte_carlo.score.txt'
strep_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/streptomycin_630.bipartite.files/streptomycin_630.monte_carlo.score.txt'
gf_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/germfree_630.bipartite.files/germfree_630.monte_carlo.score.txt'

cef_importance <- read.delim(cef_importance_file, header=TRUE, sep='\t', row.names=1)
clinda_importance <- read.delim(clinda_importance_file, header=TRUE, sep='\t', row.names=1)
strep_importance <- read.delim(strep_importance_file, header=TRUE, sep='\t', row.names=1)
gf_importance <- read.delim(gf_importance_file, header=TRUE, sep='\t', row.names=1)
rm(cef_importance_file, clinda_importance_file, strep_importance_file, gf_importance_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format metabolite importance scores
cef_importance <- cef_importance[order(-cef_importance$Metabolite_score),]
clinda_importance <- clinda_importance[order(-clinda_importance$Metabolite_score),]
strep_importance <- strep_importance[order(-strep_importance$Metabolite_score),]
gf_importance <- gf_importance[order(-gf_importance$Metabolite_score),]

cef_importance <- cef_importance[c(1:50),]
clinda_importance <- clinda_importance[c(1:50),]
strep_importance <- strep_importance[c(1:50),]
gf_importance <- gf_importance[c(1:50),]

shared_importance <- as.data.frame(subset(cef_importance, (cef_importance[,1] %in% clinda_importance[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% clinda_importance[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% strep_importance[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% gf_importance[,1])))

shared_importance <- shared_importance$Compound_name
shared_cef <- as.data.frame(subset(cef_importance, (cef_importance[,1] %in% shared_importance)))
shared_cef <- shared_cef[order(shared_cef$Compound_name),] 
shared_cef$Metabolite_score <- as.numeric(as.character(shared_cef$Metabolite_score))
shared_cef$Sim_Median <- as.numeric(as.character(shared_cef$Sim_Median))
shared_clinda <- as.data.frame(subset(clinda_importance, (clinda_importance[,1] %in% shared_importance)))
shared_clinda <- shared_clinda[order(shared_clinda$Compound_name),] 
shared_clinda$Metabolite_score <- as.numeric(as.character(shared_clinda$Metabolite_score))
shared_clinda$Sim_Median <- as.numeric(as.character(shared_clinda$Sim_Median))
shared_strep <- as.data.frame(subset(strep_importance, (strep_importance[,1] %in% shared_importance)))
shared_strep <- shared_strep[order(shared_strep$Compound_name),] 
shared_strep$Metabolite_score <- as.numeric(as.character(shared_strep$Metabolite_score))
shared_strep$Sim_Median <- as.numeric(as.character(shared_strep$Sim_Median))
shared_gf <- as.data.frame(subset(gf_importance, (gf_importance[,1] %in% shared_importance)))
shared_gf <- shared_gf[order(shared_gf$Compound_name),] 
shared_gf$Metabolite_score <- as.numeric(as.character(shared_gf$Metabolite_score))
shared_gf$Sim_Median <- as.numeric(as.character(shared_gf$Sim_Median))

shared_score <- cbind(shared_cef$Metabolite_score, shared_clinda$Metabolite_score, shared_strep$Metabolite_score, shared_gf$Metabolite_score)
rownames(shared_score) <- shared_cef$Compound_name
shared_sim <- cbind(shared_cef$Sim_Median, shared_clinda$Sim_Median, shared_strep$Sim_Median, shared_gf$Sim_Median)
rownames(shared_sim) <- shared_cef$Compound_name
shared_importance <- as.data.frame(cbind(apply(shared_score, 1, median) , apply(shared_sim, 1, median)))
colnames(shared_importance) <- c('Metabolite_score', 'Sim_Median')
shared_importance$Compound_name <- rownames(shared_importance)
shared_importance <- shared_importance[order(shared_importance$Metabolite_score),]
shared_importance$Compound_name <- gsub('_',' ',shared_importance$Compound_name)
shared_importance$Compound_name <- gsub('phosphate','p',shared_importance$Compound_name)
rm(shared_cef, shared_clinda, shared_strep, shared_gf, shared_score, shared_sim)

cef_importance <- cef_importance[c(1:25),]
clinda_importance <- clinda_importance[c(1:25),]
strep_importance <- strep_importance[c(1:25),]
gf_importance <- gf_importance[c(1:25),]

cef_only_importance <- as.data.frame(subset(cef_importance, !(cef_importance[,1] %in% clinda_importance[,1])))
cef_only_importance <- as.data.frame(subset(cef_only_importance, !(cef_only_importance[,1] %in% strep_importance[,1])))
cef_only_importance <- as.data.frame(subset(cef_only_importance, !(cef_only_importance[,1] %in% gf_importance[,1])))
#cef_only_importance <- as.data.frame(subset(cef_only_importance, !(cef_only_importance[,1] %in% shared_importance$Compound_name)))
clinda_only_importance <- as.data.frame(subset(clinda_importance, !(clinda_importance[,1] %in% cef_importance[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_only_importance, !(clinda_only_importance[,1] %in% strep_importance[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_only_importance, !(clinda_only_importance[,1] %in% gf_importance[,1])))
#clinda_only_importance <- as.data.frame(subset(clinda_only_importance, !(clinda_only_importance[,1] %in% shared_importance$Compound_name)))
strep_only_importance <- as.data.frame(subset(strep_importance, !(strep_importance[,1] %in% clinda_importance[,1])))
strep_only_importance <- as.data.frame(subset(strep_only_importance, !(strep_only_importance[,1] %in% cef_importance[,1])))
strep_only_importance <- as.data.frame(subset(strep_only_importance, !(strep_only_importance[,1] %in% gf_importance[,1])))
#strep_only_importance <- as.data.frame(subset(strep_only_importance, !(strep_only_importance[,1] %in% shared_importance$Compound_name)))
gf_only_importance <- as.data.frame(subset(gf_importance, !(gf_importance[,1] %in% clinda_importance[,1])))
gf_only_importance <- as.data.frame(subset(gf_only_importance, !(gf_only_importance[,1] %in% strep_importance[,1])))
gf_only_importance <- as.data.frame(subset(gf_only_importance, !(gf_only_importance[,1] %in% cef_importance[,1])))
#gf_only_importance <- as.data.frame(subset(gf_only_importance, !(gf_only_importance[,1] %in% shared_importance$Compound_name)))
rm(cef_importance, clinda_importance, strep_importance, gf_importance)

cef_only_importance <- cef_only_importance[order(cef_only_importance$Metabolite_score),]
clinda_only_importance <- clinda_only_importance[order(clinda_only_importance$Metabolite_score),]
strep_only_importance <- strep_only_importance[order(strep_only_importance$Metabolite_score),]
gf_only_importance <- gf_only_importance[order(gf_only_importance$Metabolite_score),]

cef_only_importance$abx <- 'Cefoperazone'
cef_only_importance$color <- wes_palette("FantasticFox")[3]
clinda_only_importance$abx <- 'Clindamycin'
clinda_only_importance$color <- wes_palette("FantasticFox")[5]
strep_only_importance$abx <- 'Streptomycin'
strep_only_importance$color <- wes_palette("FantasticFox")[1]
gf_only_importance$abx <- 'Gnotobiotic'
gf_only_importance$color <- 'forestgreen'

top_importances <- rbind(cef_only_importance[,c(1,2,3,4,6,7)], clinda_only_importance[,c(1,2,3,4,6,7)], 
                         strep_only_importance[,c(1,2,3,4,6,7)], gf_only_importance[,c(1,2,3,4,6,7)])
top_importances$abx <- as.factor(top_importances$abx)
top_importances$abx <- ordered(top_importances$abx, levels=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Gnotobiotic'))
#top_importances$Sim_StD <- top_importances$Sim_StD * 1.645 # 90% confidence interval
top_importances$Compound_name <- gsub('_',' ',top_importances$Compound_name)
top_importances$Compound_name <- gsub('mono', '', top_importances$Compound_name)
top_importances$Compound_name <- gsub('phosphate','p',top_importances$Compound_name)
top_importances$Compound_name[top_importances$Compound_name == '5,6,7,8-Tetrahydromethanopterin'] <- 'Tetrahydromethanopterin'
top_importances$Compound_name[top_importances$Compound_name == '1-(5\'-Phosphoribosyl)-5-amino-4-(N-succinocarboxamide)-imidazole'] <- 'SAICAR'
top_importances <- subset(top_importances, rownames(top_importances) != 'C00012') # Remove generic Peptide 
rm(cef_only_importance, clinda_only_importance, strep_only_importance, gf_only_importance)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_5.pdf'
pdf(file=plot_file, width=14, height=14)
layout(matrix(c(1,2,3,4,4,4,4,5,
                6,7,7,4,4,4,4,8,
                9,7,7,4,4,4,4,10,
                11,12,13,4,4,4,4,14,
                15,15,15,15,16,16,16,16,
                15,15,15,15,16,16,16,16,
                15,15,15,15,16,16,16,16,
                15,15,15,15,16,16,16,16), nrow=8, ncol=8, byrow=TRUE))

plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot

#-------------------------------------------------------------------------------------------------------------------------------------#

# A - Example network and importance calculation
par(mar=c(0,1,0,0))
plot(network, vertex.label=NA, layout=optimal_layout2, vertex.frame.color='black', xlim=c(-1.2,1.2), ylim=c(-1.2,1.2))

text(-0.96, 0.08, expression(Importance(m) == paste(log[2], 
                                                   bgroup('(',frac(Sigma * t[i], Sigma * e[o]),''),' - ', 
                                                   bgroup('',frac(Sigma * t[o], Sigma * e[i]),')'))), cex = 1.6) # Importance algorithm

text(x=-0.95, y=1.14, labels='Tetrathionate reductase', font=2, cex=1.6) # Enzyme 1 name
text(x=-1, y=1, labels='7', col='white', cex=1.3) # Enzyme 1 transcription
text(x=-0.5, y=-1.3, labels='Sulfate reductase', font=2, cex=1.6) # Enzyme 2 name
text(x=-0.5, y=-1, labels='94', col='white', cex=2.4) # Enzyme 2 transcription
text(x=1, y=0.77, labels='Thiosulfate Oxidase', font=2, cex=1.6) # Enzyme 3
text(x=0.99, y=0.44, labels='115', col='white', cex=2.6) # Enzyme 3 transcription
text(x=-0.165, y=0.145, 'm', col='white', cex=2.1) # Substrate node label
text(x=c(0.1,0.1), y=c(-0.02,-0.12), labels=c('Thiosulfate','= 6.554'), cex=1.7, font=c(2,1)) # Compound & calculated importance
segments(x0=-0.05, y0=-0.18, x1=0.25, y1=-0.18, lwd=2)
legend(x=0.7, y=1.3, legend=c('Enzyme node', 'Metabolite node'), 
       pt.bg=c('firebrick3', 'blue3'), col='black', pch=21, pt.cex=3, cex=1.7)
text(x=-0.5, y=-1.12, expression(t[i]), col='white', cex=2) # labeled transcription for input reactions
text(x=0.99, y=0.32, expression(t[i]), col='white', cex=2)
text(x=-1, y=0.87, expression(t[o]), col='black', cex=2) # labeled transcription for output reactions
text(x=-0.6, y=0.7, expression(e[i]), col='black', cex=2) # labeled indegree
text(x=-0.4, y=-0.3, expression(e[o]), col='black', cex=2) # labeled outdegree
text(x=0.3, y=0.33, expression(e[o]), col='black', cex=2)
Arrows(x0=0.63, y0=-0.7, x1=0.12, y1=-0.7, lwd=3, arr.type='triangle', arr.length=0.6, arr.width=0.3) # Score explanation line
Arrows(x0=0.63, y0=-0.7, x1=1.14, y1=-0.7, lwd=3, arr.type='triangle', arr.length=0.6, arr.width=0.3)
segments(x0=0.63, y0=-0.65, x1=0.63, y1=-0.75, lwd=2)
text(x=0.63, y=-0.83, '0', cex=1.8) 
text(x=0.12, y=-0.81, '-', cex=2.2)
text(x=1.14, y=-0.81, '+', cex=2.2)
text(x=1.15, y=-0.6, 'More likely consumed', cex=1.3)
text(x=0.14, y=-0.6, 'More likely released', cex=1.3)
text(x=0.63, y=-0.95, 'Importance Score', cex=1.6, font=2) 

plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot

# Large component of C. difficile 630 graph
par(mar=c(1,3,1,1))
plot(largest_simple_graph, vertex.label=NA, layout=optimal_layout1,
     edge.arrow.size=0.2, edge.arrow.width=0.4, vertex.frame.color='black')
mtext('a', side=2, line=2, las=2, adj=-2, padj=-7, cex=1.9, font=2)

plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot

# Boxes in A are drawn in seperate software (Gimp)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Shared metabolite importances
par(mar=c(4,4,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(shared_importance$Metabolite_score, labels=shared_importance$Compound_name, lcolor=NA, cex=1.4, color='black', 
         xlab='Median Metabolite Importance Score', xlim=c(-4,10), 
         pch=c(19,19,1,19,1,1,19,1,1,19,19,19,19,19,19,19,19))
segments(x0=rep(-4, 16), y0=c(1:17), x1=rep(13, 16), y1=c(1:17), lty=2)
abline(v=0, col='gray68', lwd=1.7)
points(x=shared_importance$Sim_Median, y=c(1:17), cex=2.5, col='black', pch='|') # Add simulated medians
mtext('b', side=2, line=2, las=2, adj=2, padj=-12.5, cex=1.9, font=2)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Unique metabolite importances
par(mar=c(4,4,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(top_importances$Metabolite_score, labels=top_importances$Compound_name, 
         lcolor=NA, cex=1.4, groups=top_importances$abx, color='black', 
         xlab='Metabolite Importance Score', xlim=c(-4,10), 
         gcolor=c(wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5],'forestgreen'), 
         pch=c(19,1,1,19,1,1,19,19,19,1,19,19,19,19,19,19,19,19,19,19,19), lwd=3)
segments(x0=rep(-4, 14), y0=c(1:14, 17:18, 21, 24:26), x1=rep(10, 14), y1=c(1:14, 17:18, 21, 24:26), lty=2)
abline(v=0, col='gray68', lwd=1.7)

# Add simulated medians
points(x=top_importances[c(20:7),3], y=c(1:14), cex=2.5, col='black', pch='|') # Gnotobiotic
#points(x=top_importances[c(20:7),3]+top_importances[c(20:7),4], y=c(1:14), cex=1.8, col='black', pch='|')
#points(x=top_importances[c(20:7),3]-top_importances[c(20:7),4], y=c(1:14), cex=1.8, col='black', pch='|')
points(x=top_importances[c(2:3),3], y=c(17:18), cex=2.5, col='black', pch='|') # Clindamycin
#points(x=top_importances[c(2:3),3]+top_importances[c(2:3),4], y=c(17:18), cex=1.8, col='black', pch='|')
#points(x=top_importances[c(2:3),3]-top_importances[c(2:3),4], y=c(17:18), cex=1.8, col='black', pch='|')
points(x=top_importances[1,3], y=21, cex=2.5, col='black', pch='|') # Cefoperazone
#points(x=top_importances[1,3]+top_importances[1,4], y=21, cex=1.8, col='black', pch='|')
#points(x=top_importances[1,3]-top_importances[1,4], y=21, cex=1.8, col='black', pch='|')
points(x=top_importances[c(4:6),3], y=c(24:26), cex=2.5, col='black', pch='|') # Streptomycin
#points(x=top_importances[c(4:6),3]+top_importances[c(4:6),4], y=c(24:26), cex=1.8, col='black', pch='|')
#points(x=top_importances[c(4:6),3]-top_importances[c(4:6),4], y=c(24:26), cex=1.8, col='black', pch='|')

mtext('c', side=2, line=2, las=2, adj=2, padj=-12.5, cex=1.9, font=2)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()











rm(optimal_layout1, optimal_layout2, largest_simple_graph, network, plot_file, top_importances, shared_importance)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

