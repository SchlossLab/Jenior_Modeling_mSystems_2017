
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
network_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/bipartite_graph.txt'
ko_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/cefoperazone_630.RNA_reads2cdf630.norm.ko.txt'

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
optimal <- layout.graphopt(graph=largest_simple_graph, niter=1000, charge=0.001, mass=50, spring.length=0, spring.constant=1)
#optimal <- layout.kamada.kawai(graph=largest_simple_graph)
#optimal <- layout.fruchterman.reingold(graph=largest_simple_graph)
#optimal <- layout_nicely(graph=largest_simple_graph, dim=2)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in substrate importance data

cef_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/cefoperazone_630.monte_carlo.score.txt'
clinda_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/clindamycin_630.bipartite.files/clindamycin_630.monte_carlo.score.txt'
strep_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/streptomycin_630.bipartite.files/streptomycin_630.monte_carlo.score.txt'
gf_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/germfree.bipartite.files/germfree.monte_carlo.score.txt'

cef_importance <- read.table(cef_importance_file, header=TRUE, sep='\t', row.names=1)
clinda_importance <- read.table(clinda_importance_file, header=TRUE, sep='\t', row.names=1)
strep_importance <- read.table(strep_importance_file, header=TRUE, sep='\t', row.names=1)
gf_importance <- read.table(gf_importance_file, header=TRUE, sep='\t', row.names=1)
rm(cef_importance_file, clinda_importance_file, strep_importance_file, gf_importance_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format inout and output importances
cef_input <- cef_importance[,c(1,2,4)]
cef_input <- subset(cef_input, cef_input$StD_from_Sim_Input_Mean > 1)
cef_input <- cef_input[order(-cef_input$Input_metabolite_score),] 
cef_input <- cef_input[c(1:10),]
cef_input <- cef_input[order(cef_input$Input_metabolite_score),] 
clinda_input <- clinda_importance[,c(1,2,4)]
clinda_input <- subset(clinda_input, clinda_input$StD_from_Sim_Input_Mean > 1)
clinda_input <- clinda_input[order(-clinda_input$Input_metabolite_score),]
clinda_input <- clinda_input[c(1:10),]
clinda_input <- clinda_input[order(clinda_input$Input_metabolite_score),]
strep_input <- strep_importance[,c(1,2,4)]
strep_input <- subset(strep_input, strep_input$StD_from_Sim_Input_Mean > 1)
strep_input <- strep_input[order(-strep_input$Input_metabolite_score),]
strep_input <- strep_input[c(1:10),]
strep_input <- strep_input[order(strep_input$Input_metabolite_score),]
gf_input <- gf_importance[,c(1,2,4)]
gf_input <- subset(gf_input, gf_input$StD_from_Sim_Input_Mean > 1)
gf_input <- gf_input[order(-gf_input$Input_metabolite_score),]
gf_input <- gf_input[c(1:10),]
gf_input <- gf_input[order(gf_input$Input_metabolite_score),]

rm(cef_importance, clinda_importance, strep_importance, gf_importance)

# Merge tables
cef_input$abx <- 'Cefoperazone'
cef_input$abx <- factor(cef_input$abx)
cef_input$color <- 'forestgreen'
clinda_input$abx <- 'Clindamycin'
clinda_input$abx <- factor(clinda_input$abx)
clinda_input$color <- 'darkorange3'
strep_input$abx <- 'Streptomycin'
strep_input$abx <- factor(strep_input$abx)
strep_input$color <- 'dodgerblue4'
gf_input$abx <- 'Gnotobiotic'
gf_input$abx <- factor(gf_input$abx)
gf_input$color <- 'darkorchid4'
input_importance <- rbind(cef_input, clinda_input, strep_input, gf_input)

rm(cef_input, clinda_input, strep_input, gf_input)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_5.pdf'
pdf(file=plot_file, width=6, height=15)
layout(matrix(c(1, 2,2), nrow=3, ncol=1, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 4A

# Plot the large component of the graph
par(mar=c(1,3,1,1))
plot(largest_simple_graph, vertex.label=NA, layout=optimal,
     edge.arrow.size=0.5, edge.arrow.width=0.8, vertex.frame.color='black')
legend('topright', legend=c('KEGG Ortholog', 'Enzyme Substrate'), 
       pt.bg=c('firebrick3', 'blue3'), col='black', pch=21, pt.cex=2.3)
legend(x=0.7, y=-0.8, legend=c('Total nodes: 1070', 'KO nodes: 404', 'Substrate nodes: 666'), pt.cex=0, text.font=c(2,1,1), bty='n')
mtext('A', side=2, line=2, las=2, adj=0.5, padj=-12, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 4B

par(mar=c(4,3,1,1))
dotchart(input_importance$Input_metabolite_score, labels=gsub('_',' ',input_importance$Compound_name), 
         groups= input_importance$abx, color=input_importance$color,
         xlab="Input Metabolite Score", gcolor="black", pch=19)
mtext('B', side=2, line=2, las=2, adj=0.5, padj=-26, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

dev.off()




