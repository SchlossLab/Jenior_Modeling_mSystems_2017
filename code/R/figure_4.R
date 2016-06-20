
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
substrate_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/compound.lst'

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in metabolic network data
network <- read.table(network_file, header=FALSE, sep='\t')
ko <- read.table(ko_file, header=FALSE, sep='\t')
substrate <- as.vector(read.table(substrate_file, header=FALSE, sep='\t')$V1)
rm(network_file, ko_file, substrate_file)

# Format directed graph
raw_graph <- graph.data.frame(network, directed=TRUE)
rm(network)

# Scale points by number of transcripts mapped
nodes <- c(as.vector(ko[,1]), substrate)
ko[,2][ko[,2] == 0] <- 0.25
mappings <- c(as.vector(ko[,2] * 5), rep(2, length(substrate)))
node_size <- as.matrix(setNames(mappings, nodes))
V(raw_graph)$size <- node_size

# Color nodes
V(raw_graph)$color <- ifelse(grepl('K', V(raw_graph)$name), 'firebrick', 'blue3') # Color nodes

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
#V(net)$size=degree(net)*5 Scale by degree!

# Calculate optimal layout
optimal <- layout.graphopt(graph=largest_simple_graph, niter=1000, charge=0.01, mass=50, spring.length=0, spring.constant=1)




E(largest_simple_graph)$color <- 'gray15' # Color edges

node transparency!!!!!!! will let you see them all!

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_4.pdf'
pdf(file=plot_file, width=8, height=12)
layout(matrix(c(1,
                2), 
              nrow=2, ncol=1, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 4A

# Plot the large component of the graph
par(mar=c(0,0,0,0))
plot(largest_simple_graph, vertex.label=NA, layout=optimal,
     edge.arrow.size=0.5, edge.arrow.width=0.8, vertex.frame.color='black', vertex.size=node_size)



plot(g, 
     vertex.color = adjustcolor("SkyBlue2", alpha.f = .5), 
     vertex.label.color = adjustcolor("black", .5))



legend('bottomleft', legend=c('KEGG Ortholog', 'Enzyme Substrate'), 
       pt.bg=c('firebrick', 'blue3'), col='black', pch=21, pt.cex=2.3)

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







