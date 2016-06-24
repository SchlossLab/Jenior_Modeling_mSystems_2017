
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
ko_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/cefoperazone_630.original_mapping.txt'

# Read in metabolic network data
network <- read.table(network_file, header=FALSE, sep='\t')
ko <- read.table(ko_file, header=TRUE, sep='\t')
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
ko <- ko[match(ko_simp, ko$KO_code),]
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

# Format metabolite importance scores
cef_importance <- cef_importance[order(-cef_importance$Metabolite_score),]
clinda_importance <- clinda_importance[order(-clinda_importance$Metabolite_score),]
strep_importance <- strep_importance[order(-strep_importance$Metabolite_score),]
gf_importance <- gf_importance[order(-gf_importance$Metabolite_score),]

cef_importance <- cef_importance[c(1:25),]
clinda_importance <- clinda_importance[c(1:25),]
strep_importance <- strep_importance[c(1:25),]
gf_importance <- gf_importance[c(1:25),]

cef_only_importance <- as.data.frame(subset(cef_importance, !(cef_importance[,1] %in% clinda_importance[,1])))
cef_only_importance <- as.data.frame(subset(cef_only_importance, !(cef_only_importance[,1] %in% strep_importance[,1])))
cef_only_importance <- as.data.frame(subset(cef_only_importance, !(cef_only_importance[,1] %in% gf_importance[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_importance, !(clinda_importance[,1] %in% cef_importance[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_only_importance, !(clinda_only_importance[,1] %in% strep_importance[,1])))
clinda_only_importance <- as.data.frame(subset(clinda_only_importance, !(clinda_only_importance[,1] %in% gf_importance[,1])))
strep_only_importance <- as.data.frame(subset(strep_importance, !(strep_importance[,1] %in% clinda_importance[,1])))
strep_only_importance <- as.data.frame(subset(strep_only_importance, !(strep_only_importance[,1] %in% cef_importance[,1])))
strep_only_importance <- as.data.frame(subset(strep_only_importance, !(strep_only_importance[,1] %in% gf_importance[,1])))
gf_only_importance <- as.data.frame(subset(gf_importance, !(gf_importance[,1] %in% clinda_importance[,1])))
gf_only_importance <- as.data.frame(subset(gf_only_importance, !(gf_only_importance[,1] %in% strep_importance[,1])))
gf_only_importance <- as.data.frame(subset(gf_only_importance, !(gf_only_importance[,1] %in% cef_importance[,1])))
rm(cef_importance, clinda_importance, strep_importance, gf_importance)

cef_only_importance <- cef_only_importance[order(cef_only_importance$Metabolite_score),]
clinda_only_importance <- clinda_only_importance[order(clinda_only_importance$Metabolite_score),]
strep_only_importance <- strep_only_importance[order(strep_only_importance$Metabolite_score),]
gf_only_importance <- gf_only_importance[order(gf_only_importance$Metabolite_score),]

cef_only_importance$abx <- 'Cefoperazone'
cef_only_importance$color <- 'blue2'
clinda_only_importance$abx <- 'Clindamycin'
clinda_only_importance$color <- 'chartreuse4'
strep_only_importance$abx <- 'Streptomycin'
strep_only_importance$color <- 'firebrick3'
gf_only_importance$abx <- 'Gnotobiotic'
gf_only_importance$color <- 'darkorchid4'

top_importances <- rbind(cef_only_importance[,c(1,2,5,6)], clinda_only_importance[,c(1,2,5,6)], strep_only_importance[,c(1,2,5,6)], gf_only_importance[,c(1,2,5,6)])
top_importances$abx <- as.factor(top_importances$abx)
top_importances$abx <- ordered(top_importances$abx, levels=c('Cefoperazone', 'Streptomycin', 'Clindamycin', 'Gnotobiotic'))
top_importances$Compound_name <- gsub('_',' ',top_importances$Compound_name)
rm(cef_only_importance, clinda_only_importance, strep_only_importance, gf_only_importance)

top_importances <- top_importances[ !(rownames(top_importances) %in% c('C03688', 'C11436')), ]

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_5.pdf'
pdf(file=plot_file, width=16, height=8)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 4A - Large component of graph
par(mar=c(1,3,1,1))
plot(largest_simple_graph, vertex.label=NA, layout=optimal,
     edge.arrow.size=0.5, edge.arrow.width=0.8, vertex.frame.color='black')
legend(x=0.4, y=1.2, legend=c('KEGG Ortholog', 'Enzyme Substrate'), 
       pt.bg=c('firebrick3', 'blue3'), col='black', pch=21, pt.cex=2.3)
legend(x=0.4, y=-0.8, legend=c('Total nodes: 1070', 'Enzyme nodes: 404', 'Substrate nodes: 666'), pt.cex=0, text.font=c(2,1,1), bty='n')
mtext('A', side=2, line=2, las=2, adj=0.5, padj=-20.5, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 4B - Top compound importances
par(mar=c(4,3,1,1), xaxs='i')
dotchart(top_importances$Metabolite_score, labels=top_importances$Compound_name, 
         lcolor=NA, cex=1.5, groups=top_importances$abx, color=top_importances$color, 
         xlab='Metabolite Importance Score', xlim=c(0,14), gcolor="black", pch=19)
segments(x0=rep(0, 15), y0=c(1:6, 9:12, 15:18, 21), x1=rep(14, 15), y1=c(1:6, 9:12, 15:18, 21), lty=2)
mtext('B', side=2, line=2, las=2, adj=0.5, padj=-18, cex=1.5)

#-------------------------------------------------------------------------------------------------------------------------------------#

dev.off()




