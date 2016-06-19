
deps <- c('vegan', 'igraph');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

# Define and check color palette
palette_plot <- function(col, border='light gray', ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}

# Define variables
network_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/'
ko_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/'
substrate_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/'

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in metabolic network data
network <- read.table(network_file, header=FALSE, sep='\t')
ko <- read.table(ko_file, header=FALSE, sep='\t')
rm(network_file, ko_score_file, substrate_file)

# Format directed graph
raw_graph <- graph.data.frame(network, directed=TRUE)

# Remove loops and multiple edges to make visualzation easier
simple_graph <- simplify(raw_graph)

# Get largest component
largest_component <- which.max(sapply(simple_graph, vcount))
largest_simple_graph <- simple.graph[[largest_component]]

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in compound importace data
substrate <- as.vector(read.table(substrate_file, header=FALSE, sep='\t')$V1)







#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_4.pdf'
pdf(file=plot_file, width=12, height=6)
layout(matrix(c(1,2,
                3,4), 
              nrow=2, ncol=2, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 4A

# Scale points by number of transcripts mapped
mappings <- (log10(ko_score[,2])) / 2
V(largest_simple_graph)$size <- mappings
  
# Color nodes
V(largest_simple_graph)$color <- ifelse(V(largest_simple_graph)$name %in% substrate, 'firebrick', 'blue3') # Color nodes
E(largest_simple_graph)$color <- 'gray15' # Color edges

# Plot the large component of th graph
par(mar=c(0,0,0,0))
plot(largest_simple_graph, vertex.label=NA, layout=layout.graphopt,
     edge.arrow.size=0.5, edge.arrow.width=0.8, vertex.frame.color='black')
legend('bottomleft', legend=c('KEGG Ortholog', 'Enzyme Substrate'), 
       pt.bg=c('firebrick', 'blue3'), col='black', pch=21, pt.cex=2.5, bty = "n")

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







