
# Load dependencies
deps <- c('vegan', 'igraph', 'ggplot2', 'shape', 'matrixStats');
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

# Read in substrate importance data
invitro_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/invitro_glucose_630.bipartite.files/invitro_glucose_630.importance_score.tsv'
invitro_importance <- read.delim(invitro_importance_file, header=TRUE, sep='\t', row.names=1)
rm(invitro_importance_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format metabolite importance scores
invitro_importance$Metabolite_score <- as.numeric(as.character(invitro_importance$Metabolite_score))
invitro_importance$Sim_Median <- as.numeric(as.character(invitro_importance$Sim_Median))
invitro_importance$Sim_Lower_95_Confidence <- as.numeric(as.character(invitro_importance$Sim_Lower_95_Confidence))
invitro_importance$Sim_Upper_95_Confidence <- as.numeric(as.character(invitro_importance$Sim_Upper_95_Confidence))
invitro_importance <- invitro_importance[order(-invitro_importance$Metabolite_score),]
invitro_importance <- invitro_importance[c(1:25),]

# Format names to look better for the plot
invitro_importance$Compound_name <- gsub('_',' ',invitro_importance$Compound_name)
invitro_importance$Compound_name <- gsub('phosphate','p', invitro_importance$Compound_name)
invitro_importance[invitro_importance == '2-Methylpropanoyl-CoA'] <- 'Isobutyryl-CoA'

# Change point color based on significance
invitro_importances$entry_color <- ifelse(invitro_importances$Compound_name == 'D-Glucose', 'red', 'black')
invitro_importances$point_color <- ifelse(invitro_importances$Metabolite_score > invitro_importances$Sim_Upper_95_Confidence, 'black', 'white')

# set same minimums as plotting area
invitro_importance[invitro_importance < -2] <- -2

# Make really low confidence intervals visible
invitro_importance$Sim_Upper_95_Confidence[invitro_importance$Sim_Upper_95_Confidence == -2] <- -1.7

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in growth rate data
# Define variables
growth_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/cd630_growth.tsv'

# Read in data
growth <- read.delim(growth_file, sep='\t', header=TRUE, row.names=1)
growth <- as.data.frame(t(growth))
rm(growth_file)

# Seperate to groups of each growth substrate and format
glucose <- cbind(growth$glucose_1, growth$glucose_2, growth$glucose_3) - growth$glucose_blank
glucose[glucose < 0] <- 0
no_carb <- cbind(growth$noCarb_1, growth$noCarb_2, growth$noCarb_3) - growth$noCarb_blank
no_carb[no_carb < 0] <- 0
no_aa <- cbind(growth$noAA_1, growth$noAA_2, growth$noAA_3) - growth$noAA_blank
no_aa[no_aa < 0] <- 0
bhi <- cbind(growth$bhi_1, growth$bhi_2, growth$bhi_3) - growth$bhi_blank
bhi[bhi < 0] <- 0
rm(growth)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format growth curves

# Find medians 
glucose_median <- apply(glucose, 1, median)[1:25]
no_carb_median <- apply(no_carb, 1, median)[1:25]
no_aa_median <- apply(no_aa, 1, median)[1:25]
bhi_median <- apply(bhi, 1, median)[1:25]

# Standard deviations
glucose_sd <- rowSds(glucose)[1:25]
no_carb_sd <- rowSds(no_carb)[1:25]
no_aa_sd <- rowSds(no_aa)[1:25]
bhi_sd <- rowSds(bhi)[1:25]
rm(glucose, no_carb, no_aa, bhi)

# Time of maximum glucose growth rate
max_gluc_time <- which.max(diff(glucose_median)) # 7 hours
max_gluc_rate <- diff(glucose_median)[which.max(diff(glucose_median))]

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_5.pdf'
pdf(file=plot_file, width=28, height=10)
layout(matrix(c(1,2,2,2,3,3,3,3,4,4,4,
                5,2,2,2,3,3,3,3,4,4,4,
                6,2,2,2,3,3,3,3,4,4,4,
                7,2,2,2,3,3,3,3,4,4,4), nrow=4, ncol=11, byrow=TRUE))

plot.new()
text(x=0.9, y=0.1, 'a', cex=2.1, font=2)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Example network and importance calculation
par(mar=c(0,0,0,0))
plot(network, vertex.label=NA, layout=optimal_layout2, vertex.frame.color='black', xlim=c(-1.2,1.2), ylim=c(-1.2,1.2))

text(x=-0.95, y=1.11, labels='dAdo Aminohydrolase', font=2, cex=1.3) # Enzyme 1 name
text(x=-1, y=1, labels='7', col='white', cex=1.3) # Enzyme 1 transcription
text(x=-0.5, y=-1.3, labels='ATP:dAdo 5\'-Phosphotransferase', font=2, cex=1.3) # Enzyme 2 name
text(x=-0.5, y=-1, labels='94', col='white', cex=2.4) # Enzyme 2 transcription
text(x=1, y=0.75, labels=expression(bold(paste('dAdo:', PO[4]^-3, ' Ribosyltransferase'))), cex=1.3) # Enzyme 3
text(x=0.99, y=0.44, labels='115', col='white', cex=2.6) # Enzyme 3 transcription
text(x=-0.165, y=0.145, 'm', col='white', cex=2.1) # Substrate node label

text(x=c(-0.8,-0.8), y=c(0.15,0.05), labels=c('Deoxyadenosine (dAdo)','Importance = 6.554'), cex=1.5, font=c(2,1)) # Compound & calculated importance
segments(x0=-1.17, y0=0, x1=-0.43, y1=0, lwd=2)

legend(x=0.92, y=1.3, legend=c('Enzyme node', 'Metabolite node'),
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

# Growth on glucose
par(mar=c(7,7,1.5,2), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE, mgp=c(4,2,0))
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,26), ylim=c(-0.03,0.63), lwd=2, pch=15, xlab='Time (hours)', ylab=expression(OD[600]))
abline(h=seq(0,1,0.1), lty=3, col='gray68') # adding gridlines
abline(v=seq(1,26,2), lty=3, col='gray68') # adding gridlines
axis(1, at=seq(1,25,4), labels=seq(0,12,2), tck=-0.018)
axis(2, at=seq(0.0,1.0,0.2), labels=c('0.0','0.2','0.4','0.6','0.8','1.0'), tck=-0.018)
mtext('b', side=2, line=2, las=2, adj=3, padj=-24, cex=1.5, font=2)

lines(glucose_median, type='o', lwd=3, pch=15, cex=3, col='dodgerblue4')
segments(x0=seq(1,25,1), y0=glucose_median+glucose_sd, x1=seq(1,25,1), y1=glucose_median-glucose_sd, lwd=3, cex=2, col='dodgerblue4')
segments(x0=seq(1,25,1)-0.3, y0=glucose_median+glucose_sd, x1=seq(1,25,1)+0.3, y1=glucose_median+glucose_sd, lwd=3, col='dodgerblue4')
segments(x0=seq(1,25,1)-0.3, y0=glucose_median-glucose_sd, x1=seq(1,25,1)+0.3, y1=glucose_median-glucose_sd, lwd=3, col='dodgerblue4')

lines(no_carb_median, type='o', lwd=3, pch=16, cex=2, col='black')
segments(x0=seq(1,25,1), y0=no_carb_median+no_carb_sd, x1=seq(1,25,1), y1=no_carb_median-no_carb_sd, lwd=3, cex=2, col='black')
segments(x0=seq(1,25,1)-0.3, y0=no_carb_median+no_carb_sd, x1=seq(1,25,1)+0.3, y1=no_carb_median+no_carb_sd, lwd=3, col='black')
segments(x0=seq(1,25,1)-0.3, y0=no_carb_median-no_carb_sd, x1=seq(1,25,1)+0.3, y1=no_carb_median-no_carb_sd, lwd=3, col='black')

lines(no_aa_median, type='o', lwd=3, pch=17, cex=2, col='gray60')
segments(x0=seq(1,25,1), y0=no_aa_median+no_aa_sd, x1=seq(1,25,1), y1=no_aa_median-no_aa_sd, lwd=3, cex=2, col='gray60')
segments(x0=seq(1,25,1)-0.3, y0=no_aa_median+no_aa_sd, x1=seq(1,25,1)+0.3, y1=no_aa_median+no_aa_sd, lwd=3, col='gray60')
segments(x0=seq(1,25,1)-0.3, y0=no_aa_median-no_aa_sd, x1=seq(1,25,1)+0.3, y1=no_aa_median-no_aa_sd, lwd=3, col='gray60')

legend('topleft', legend=c('CDMM + D-Glucose','CDMM - Carbohydrates','CDMM - Amino acids'), 
       col=c('dodgerblue4','black','gray60'), pch=c(15,16,17), cex=2.7, pt.cex=4, bg='white', lwd=2.5)

# Time point for transcriptomic sequencing
Arrows(x0=max_gluc_time, y0=0.415, x1=max_gluc_time, y1=0.35, lwd=3, arr.type='triangle', arr.length=0.6, arr.width=0.3, col='firebrick')
text(x=max_gluc_time, y=0.43, 'Max rate', col='firebrick', cex=1.5, font=2)

#---------------------------------------#

# Shared metabolite importances
par(mar=c(4,4,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
plot(0, type='n', axes=F, xlab='', ylab='') # Empty plot
#dotchart(invitro_importance$Metabolite_score, labels=invitro_importance$Compound_name, lcolor=NA, cex=1.2, color=entry_color,
#         xlab='Median Metabolite Importance Score', xlim=c(-2,12), pch=19)
#segments(x0=rep(-2, 14), y0=c(1:14), x1=rep(12, 14), y1=c(1:14), lty=2) # Dotted lines
mtext('c', side=2, line=2, las=2, adj=2, padj=-25, cex=1.5, font=2)

#segments(x0=invitro_importance$Sim_Upper_95_Confidence, y0=seq(0.7,13.7,1), x1=invitro_importance$Sim_Upper_95_Confidence, y1=seq(1.3,14.3,1), lwd=1.5, col='gray75') # lower 95% confidence
#segments(x0=invitro_importance$Sim_Median, y0=seq(0.9,13.9,1), x1=invitro_importance$Sim_Median, y1=seq(1.1,14.1,1), lwd=1.5, col='gray75') # median
#segments(x0=invitro_importance$Sim_Upper_95_Confidence, y0=seq(0.7,13.7,1), x1=invitro_importance$Sim_Upper_95_Confidence, y1=seq(1.3,14.3,1), lwd=1.5, col='gray75') # upper 95% confidence
#segments(x0=invitro_importance$Sim_Lower_95_Confidence, y0=c(1:14), x1=invitro_importance$Sim_Upper_95_Confidence, y1=c(1:14), lwd=1.5, col='gray75') # connection

# Labeled significance
#points(x=shared_importance$Metabolite_score, y=c(1:14), pch=21, cex=2.75, lwd=2, col='black', bg=shared_importance$point_color)

#segments(x0=-2, y0=0, x1=-2, y1=18) # Left side of plot

#---------------------------------------#

# Add some network stats (Complete network, not just largest component)
plot.new()
text(x=0.19, y=0.45, 'C. difficile', cex=1.4, font=c(4,2))
text(x=0.64, y=0.45, '630 network', cex=1.4, font=2)
segments(x0=0.01, y0=0.38, x1=0.88, y1=0.38, lwd=2)
text(x=0.26, y=0.28, '- 447 Enzymes', cex=1.3)
text(x=0.295, y=0.15, '- 758 Metabolites', cex=1.3)
text(x=0.23, y=0, '- 2135 Edges', cex=1.3)

#---------------------------------------#

# Large component of C. difficile 630 graph
par(mar=c(0,1,0,0))
plot(largest_simple_graph, vertex.label=NA, layout=optimal_layout1,
     edge.arrow.size=0.25, edge.arrow.width=0.4, vertex.frame.color='black')

# Draw box
rect(xleft=0.4, ybottom=-0.4, xright=0.8, ytop=-0.75, lwd=2)

#---------------------------------------#

plot(0, type='n', axes=F, xlab='', ylab='') # Empty plot
dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
rm(max_gluc_time, max_gluc_rate, glucose_median, no_carb_median, no_aa_median, bhi_median, glucose_sd, no_carb_sd, no_aa_sd, bhi_sd)
rm(optimal_layout1, optimal_layout2, largest_simple_graph, network, plot_file, invitro_importances)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

