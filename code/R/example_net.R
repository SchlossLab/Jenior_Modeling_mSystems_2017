

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
network <- matrix(c('K01', 'C01',
                    'C01', 'K02',
                    'C01', 'K03'), nrow=3, ncol=2, byrow=TRUE)
node_size <- matrix(c('K01', '5',
                      'K02', '45',
                      'K03', '50',
                      'C01', '20'), nrow=4, ncol=2, byrow=TRUE)

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

# Establish layout
optimal <- matrix(c(-21.09826017, 22.1407060,
                      0.09077637, 0.2154631,
                      -8.32243732, -29.0949351,
                      29.67130628, 7.6231375), nrow=4, ncol=2, byrow=TRUE)

#-------------------------------------------------------------------------------------------------------------------------------------#

par(mar=c(0,0,0,0))
plot(network, vertex.label=NA, layout=optimal_layout, vertex.frame.color='black', xlim=c(-1.2,1.2), ylim=c(-1.2,1.2))
text(0.6, -0.8, expression(I == paste(log[2],'( ',frac(Sigma * t[i], e[o]),' ','-',' ',frac(Sigma * t[o], e[i]),' )')), cex = 1.7) # Importance algorithm
text(x=-1, y=1.1, labels='Tetrathionate reductase', font=2) # Enzyme 1 name
text(x=-1, y=1, labels='5', col='white') # Enzyme 1 transcription
text(x=-0.5, y=-1.3, labels='Sulfane reductase', font=2) # Enzyme 2 name
text(x=-0.5, y=-1, labels='97', col='white', cex=2) # Enzyme 2 transcription
text(x=1, y=0.75, labels='Thiosulfate Oxidase', font=2) # Enzyme 3
text(x=0.99, y=0.44, labels='115', col='white', cex=2.1) # Enzyme 3 transcription
text(x=c(0.1,0.1), y=c(-0.02,-0.12), labels=c('Thiosulfate','= 6.658'), cex=1.3, font=c(2,1)) # Compound & calculated importance
segments(x0=-0.05, y0=-0.17, x1=0.25, y1=-0.17, lwd=2)










