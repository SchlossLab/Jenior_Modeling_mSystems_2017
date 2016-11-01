
# Load dependencies
deps <- c('vegan', 'igraph', 'ggplot2', 'shape', 'wesanderson', 'matrixStats');
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
  for (time in 1:nrow(raw_exp_data)){
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
cef_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/cefoperazone_630.importance_score.tsv'
clinda_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/clindamycin_630.bipartite.files/clindamycin_630.importance_score.tsv'
strep_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/streptomycin_630.bipartite.files/streptomycin_630.importance_score.tsv'
gf_importance_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/germfree_630.bipartite.files/germfree_630.importance_score.tsv'

cef_importance <- read.delim(cef_importance_file, header=TRUE, sep='\t', row.names=1)
clinda_importance <- read.delim(clinda_importance_file, header=TRUE, sep='\t', row.names=1)
strep_importance <- read.delim(strep_importance_file, header=TRUE, sep='\t', row.names=1)
gf_importance <- read.delim(gf_importance_file, header=TRUE, sep='\t', row.names=1)
rm(cef_importance_file, clinda_importance_file, strep_importance_file, gf_importance_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format metabolite importance scores
cef_importance$Metabolite_score <- as.numeric(as.character(cef_importance$Metabolite_score))
cef_importance$Sim_Median <- as.numeric(as.character(cef_importance$Sim_Median))
cef_importance$Sim_Lower_95_Confidence <- as.numeric(as.character(cef_importance$Sim_Lower_95_Confidence))
cef_importance$Sim_Upper_95_Confidence <- as.numeric(as.character(cef_importance$Sim_Upper_95_Confidence))
cef_importance <- cef_importance[order(-cef_importance$Metabolite_score),]
clinda_importance$Metabolite_score <- as.numeric(as.character(clinda_importance$Metabolite_score))
clinda_importance$Sim_Median <- as.numeric(as.character(clinda_importance$Sim_Median))
clinda_importance$Sim_Lower_95_Confidence <- as.numeric(as.character(clinda_importance$Sim_Lower_95_Confidence))
clinda_importance$Sim_Upper_95_Confidence <- as.numeric(as.character(clinda_importance$Sim_Upper_95_Confidence))
clinda_importance <- clinda_importance[order(-clinda_importance$Metabolite_score),]
strep_importance$Metabolite_score <- as.numeric(as.character(strep_importance$Metabolite_score))
strep_importance$Sim_Median <- as.numeric(as.character(strep_importance$Sim_Median))
strep_importance$Sim_Lower_95_Confidence <- as.numeric(as.character(strep_importance$Sim_Lower_95_Confidence))
strep_importance$Sim_Upper_95_Confidence <- as.numeric(as.character(strep_importance$Sim_Upper_95_Confidence))
strep_importance <- strep_importance[order(-strep_importance$Metabolite_score),]
gf_importance$Metabolite_score <- as.numeric(as.character(gf_importance$Metabolite_score))
gf_importance$Sim_Median <- as.numeric(as.character(gf_importance$Sim_Median))
gf_importance$Sim_Lower_95_Confidence <- as.numeric(as.character(gf_importance$Sim_Lower_95_Confidence))
gf_importance$Sim_Upper_95_Confidence <- as.numeric(as.character(gf_importance$Sim_Upper_95_Confidence))
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
shared_clinda <- as.data.frame(subset(clinda_importance, (clinda_importance[,1] %in% shared_importance)))
shared_clinda <- shared_clinda[order(shared_clinda$Compound_name),]
shared_clinda$Metabolite_score <- as.numeric(as.character(shared_clinda$Metabolite_score))
shared_strep <- as.data.frame(subset(strep_importance, (strep_importance[,1] %in% shared_importance)))
shared_strep <- shared_strep[order(shared_strep$Compound_name),]
shared_strep$Metabolite_score <- as.numeric(as.character(shared_strep$Metabolite_score))
shared_gf <- as.data.frame(subset(gf_importance, (gf_importance[,1] %in% shared_importance)))
shared_gf <- shared_gf[order(shared_gf$Compound_name),]
shared_gf$Metabolite_score <- as.numeric(as.character(shared_gf$Metabolite_score))

score_median <- as.data.frame(apply(cbind(shared_cef$Metabolite_score, shared_clinda$Metabolite_score, shared_strep$Metabolite_score, shared_gf$Metabolite_score), 1, median))
sim_median <- as.data.frame(apply(cbind(shared_cef$Sim_Median, shared_clinda$Sim_Median, shared_strep$Sim_Median, shared_gf$Sim_Median), 1, median))
sim_lower95 <- as.data.frame(apply(cbind(shared_cef$Sim_Lower_95_Confidence, shared_clinda$Sim_Lower_95_Confidence, shared_strep$Sim_Lower_95_Confidence, shared_gf$Sim_Lower_95_Confidence), 1, median))
sim_upper95 <- as.data.frame(apply(cbind(shared_cef$Sim_Upper_95_Confidence, shared_clinda$Sim_Upper_95_Confidence, shared_strep$Sim_Upper_95_Confidence, shared_gf$Sim_Upper_95_Confidence), 1, median))

shared_importance <- cbind(shared_cef$Compound_name, score_median, sim_median, sim_lower95, sim_upper95)
rownames(shared_importance) <- rownames(shared_cef)
colnames(shared_importance) <- c('Compound_name','Metabolite_score','Sim_Median','Sim_Lower95','Sim_Upper95')
shared_importance <- shared_importance[order(shared_importance$Metabolite_score),]
rm(shared_cef, shared_clinda, shared_strep, shared_gf, score_median, sim_median, sim_lower95, sim_upper95)

# Format names to look better for the plot
shared_importance$Compound_name <- gsub('_',' ',shared_importance$Compound_name)
shared_importance$Compound_name <- gsub('phosphate','p',shared_importance$Compound_name)
shared_importance[shared_importance == 'N-Acetyl-D-glucosamine'] <- 'N-Acetyl-D-glucosamine    '
shared_importance[shared_importance == '2-Methylpropanoyl-CoA'] <- 'Isobutyryl-CoA'

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

cef_only_importance$abx <- 'Cefoperazone (SPF)'
clinda_only_importance$abx <- 'Clindamycin (SPF)'
strep_only_importance$abx <- 'Streptomycin (SPF)'
gf_only_importance$abx <- 'No Antibiotics (GF)'
top_importances <- rbind(cef_only_importance, clinda_only_importance, strep_only_importance, gf_only_importance)
top_importances$abx <- as.factor(top_importances$abx)
top_importances$abx <- ordered(top_importances$abx, levels=c('Streptomycin (SPF)', 'Cefoperazone (SPF)', 'Clindamycin (SPF)', 'No Antibiotics (GF)'))

# Format names to look better for the plot
top_importances$Compound_name <- gsub('_',' ',top_importances$Compound_name)
top_importances$Compound_name <- gsub('mono', '', top_importances$Compound_name)
top_importances$Compound_name <- gsub('phosphate','p',top_importances$Compound_name)
top_importances$Compound_name[top_importances$Compound_name == '5,6,7,8-Tetrahydromethanopterin'] <- 'Tetrahydromethanopterin'
top_importances$Compound_name[top_importances$Compound_name == '1-(5\'-Phosphoribosyl)-5-amino-4-(N-succinocarboxamide)-imidazole'] <- 'SAICAR'
top_importances <- subset(top_importances, rownames(top_importances) != 'C00012') # Remove generic Peptide
rm(cef_only_importance, clinda_only_importance, strep_only_importance, gf_only_importance)

# Remove non-significant points
top_importances <- rbind(top_importances[1,], as.data.frame(subset(top_importances, (top_importances[,2] > top_importances[,5]))))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,2] > shared_importance[,5])))

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
sorbitol <- cbind(growth$sorbitol_1, growth$sorbitol_2, growth$sorbitol_2) - growth$sorbitol_blank
sorbitol[sorbitol < 0] <- 0
fructose <- cbind(growth$fructose_1, growth$fructose_2, growth$fructose_3) - growth$fructose_blank
fructose[fructose < 0] <- 0
mannitol <- cbind(growth$mannitol_1, growth$mannitol_2, growth$mannitol_3) - growth$mannitol_blank
mannitol[mannitol < 0] <- 0
salicin <- cbind(growth$salicin_1, growth$salicin_2, growth$salicin_3) - growth$salicin_blank
salicin[salicin < 0] <- 0
acetate <- cbind(growth$acetate_1, growth$acetate_2, growth$acetate_3) - growth$acetate_blank
acetate[acetate < 0] <- 0
hydroxybutanoate <- cbind(growth$hydroxybutanoate_1, growth$hydroxybutanoate_2, growth$hydroxybutanoate_3) - growth$hydroxybutanoate_blank
hydroxybutanoate[hydroxybutanoate < 0] <- 0
acetylglucosamine <- cbind(growth$acetylglucosamine_1, growth$acetylglucosamine_2, growth$acetylglucosamine_3) - growth$acetylglucosamine_blank
acetylglucosamine[acetylglucosamine < 0] <- 0
acetylneuraminate <- cbind(growth$acetylneuraminate_1, growth$acetylneuraminate_2, growth$acetylneuraminate_3) - growth$acetylneuraminate_blank
acetylneuraminate[acetylneuraminate < 0] <- 0
no_carb <- cbind(growth$noCarb_1, growth$noCarb_2, growth$noCarb_3) - growth$noCarb_blank
no_carb[no_carb < 0] <- 0
no_aa <- cbind(growth$noAA_1, growth$noAA_2, growth$noAA_3) - growth$noAA_blank
no_aa[no_aa < 0] <- 0
bhi <- cbind(growth$bhi_1, growth$bhi_2, growth$bhi_3) - growth$bhi_blank
bhi[bhi < 0] <- 0
rm(growth)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Prepare data for statistical tests
glucose_test <- format_curve(glucose, 'glucose', no_carb)
fructose_test <- format_curve(fructose, 'fructose', no_carb)
sorbitol_test <- format_curve(sorbitol, 'sorbitol', no_carb)
mannitol_test <- format_curve(mannitol, 'mannitol', no_carb)
salicin_test <- format_curve(salicin, 'salicin', no_carb)
acetate_test <- format_curve(acetate, 'salicin', no_carb)
hydroxybutanoate_test <- format_curve(hydroxybutanoate, 'hydroxybutanoate', no_carb)
acetylglucosamine_test <- format_curve(acetylglucosamine, 'acetylglucosamine', no_carb)
acetylneuraminate_test <- format_curve(acetylneuraminate, 'acetylneuraminate', no_carb)
bhi_test <- format_curve(bhi, 'bhi', no_carb)
no_aa_test <- format_curve(no_aa, 'no_amino_acids', no_carb)

# Calculate differences
summary(aov(formula=od ~ substrate * time, data=glucose_test)) # <2e-16
summary(aov(formula=od ~ substrate * time, data=sorbitol_test)) # 0.0414
summary(aov(formula=od ~ substrate * time, data=fructose_test)) # <2e-16
summary(aov(formula=od ~ substrate * time, data=mannitol_test)) # <2e-16
summary(aov(formula=od ~ substrate * time, data=salicin_test)) # <2e-16
summary(aov(formula=od ~ substrate * time, data=acetate_test)) # 0.0632
summary(aov(formula=od ~ substrate * time, data=hydroxybutanoate_test)) # 0.0721
summary(aov(formula=od ~ substrate * time, data=acetylglucosamine_test)) # <2e-16
summary(aov(formula=od ~ substrate * time, data=acetylneuraminate_test)) # <2e-16
summary(aov(formula=od ~ substrate * time, data=bhi_test)) # <2e-16
summary(aov(formula=od ~ substrate * time, data=no_aa_test)) # <2e-16

# Correct p-values
corrected_p_values <- as.character(p.adjust(c(2e-16, 0.0632, 0.0721, 2e-16, 2e-16, 0.0414, 2e-16, 2e-16, 2e-16, 2e-16, 2e-16), method='bonferroni'))
corrected_p_values <- append(corrected_p_values, 'NA') # need to fix once i add the new data

# Clean up
rm(glucose_test, fructose_test, sorbitol_test, mannitol_test, salicin_test, acetylneuraminate_test, acetate_test, hydroxybutanoate_test, acetylglucosamine_test, bhi_test, no_aa_test)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format growth curves

# Find medians 
glucose_median <- apply(glucose, 1, median)
sorbitol_median <- apply(sorbitol, 1, median)
fructose_median <-  apply(fructose, 1, median)
mannitol_median <- apply(mannitol, 1, median)
salicin_median <- apply(salicin, 1, median)
acetate_median <- apply(acetate, 1, median)
hydroxybutanoate_median <- apply(hydroxybutanoate, 1, median)
acetylglucosamine_median <- apply(acetylglucosamine, 1, median)
acetylneuraminate_median <- apply(acetylneuraminate, 1, median)
no_carb_median <- apply(no_carb, 1, median)
bhi_median <- apply(bhi, 1, median)
no_aa_median <- apply(no_aa, 1, median)

# Standard deviations
glucose_sd <- rowSds(glucose)
sorbitol_sd <- rowSds(sorbitol)
acetylneuraminate_sd <- rowSds(acetylneuraminate)
fructose_sd <-  rowSds(fructose)
mannitol_sd <- rowSds(mannitol)
salicin_sd <- rowSds(salicin)
acetate_sd <- rowSds(acetate)
hydroxybutanoate_sd <- rowSds(hydroxybutanoate)
acetylglucosamine_sd <- rowSds(acetylglucosamine)
no_carb_sd <- rowSds(no_carb)
bhi_sd <- rowSds(bhi)
no_aa_sd <- rowSds(no_aa)

# Compile results
growth_medians <- as.data.frame(cbind(glucose_median, acetate_median, hydroxybutanoate_median, acetylglucosamine_median, 
                                      acetylneuraminate_median, sorbitol_median, fructose_median, mannitol_median, salicin_median, no_carb_median, no_aa_median))
growth_sds <- as.data.frame(cbind(glucose_sd, acetate_sd, hydroxybutanoate_sd, acetylglucosamine_sd, acetylneuraminate_sd, 
                                  sorbitol_sd, fructose_sd, mannitol_sd, salicin_sd, no_carb_sd, no_aa_sd))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Analyze growth curves
substrates <- c('glucose','acetate','hydroxybutanoate','acetylglucosamine','acetylneuraminate','fructose', 'mannitol','salicin','bhi','no_amino_acids','no_carbohydrates')

# Maximum growth rate
max_rate <- round(c(diff(glucose_median)[which.max(diff(glucose_median))], diff(acetate_median)[which.max(diff(acetate_median))], diff(hydroxybutanoate_median)[which.max(diff(hydroxybutanoate_median))],
                    diff(acetylglucosamine_median)[which.max(diff(acetylglucosamine_median))], diff(acetylneuraminate_median)[which.max(diff(acetylneuraminate_median))],
                    diff(fructose_median)[which.max(diff(fructose_median))], 
                    diff(mannitol_median)[which.max(diff(mannitol_median))], diff(salicin_median)[which.max(diff(salicin_median))], 
                    diff(bhi_median)[which.max(diff(bhi_median))], diff(no_aa_median)[which.max(diff(no_aa_median))], diff(no_carb_median)[which.max(diff(no_carb_median))]), digits=3)

# Time of maximum growth rate
time_max_rate <- round(c((which.max(diff(glucose_median)) * 0.5), (which.max(diff(acetate_median)) * 0.5), (which.max(diff(hydroxybutanoate_median)) * 0.5), (which.max(diff(acetylglucosamine_median)) * 0.5),
                         (which.max(diff(acetylneuraminate_median)) * 0.5),
                         (which.max(diff(fructose_median)) * 0.5), (which.max(diff(mannitol_median)) * 0.5), (which.max(diff(salicin_median)) * 0.5), 
                         (which.max(diff(bhi_median)) * 0.5), (which.max(diff(no_aa_median)) * 0.5), (which.max(diff(no_carb_median)) * 0.5)), digits=3) - 0.5
# Maximum OD
max_od <- round(c(max(glucose_median), max(acetate_median), max(hydroxybutanoate_median), max(acetylglucosamine_median), max(acetylneuraminate_median), 
                  max(fructose_median), max(mannitol_median), max(salicin_median), max(bhi_median), max(no_aa_median), max(no_carb_median)), digits=3)

# Time of max OD
time_max_od <- round(c((which.max(glucose_median) * 0.5), (which.max(acetate_median) * 0.5), (which.max(hydroxybutanoate_median) * 0.5), (which.max(acetylglucosamine_median) * 0.5),
                       (which.max(acetylneuraminate_median) * 0.5), 
                       (which.max(fructose_median) * 0.5), (which.max(mannitol_median) * 0.5), (which.max(salicin_median) * 0.5), 
                       (which.max(bhi_median) * 0.5), (which.max(no_aa_median) * 0.5), (which.max(no_carb_median) * 0.5)), digits=3) - 0.5

# Growth rate at 24 hours
rate_24_hrs <- round(c(diff(glucose_median)[length(diff(glucose_median))], diff(acetate_median)[length(diff(acetate_median))], diff(hydroxybutanoate_median)[length(diff(hydroxybutanoate_median))],
                       diff(acetylglucosamine_median)[length(diff(acetylglucosamine_median))], diff(acetylneuraminate_median)[length(diff(acetylneuraminate_median))], 
                       diff(fructose_median)[length(diff(fructose_median))], 
                       diff(mannitol_median)[length(diff(mannitol_median))], diff(salicin_median)[length(diff(salicin_median))], 
                       diff(bhi_median)[length(diff(bhi_median))], diff(no_aa_median)[length(diff(no_aa_median))], diff(no_carb_median)[length(diff(no_carb_median))]), digits=3)

# Mean growth rate
mean_rate <- round(c(mean(diff(glucose_median)), mean(diff(acetate_median)), mean(diff(hydroxybutanoate_median)), mean(diff(acetylglucosamine_median)), mean(diff(acetylneuraminate_median)),
                     mean(diff(fructose_median)), mean(diff(mannitol_median)), mean(diff(salicin_median)), mean(diff(bhi_median)), mean(diff(no_aa_median)), mean(diff(no_carb_median))), digits=3)

# Area under curve
area_under <- round(c(auc(glucose_median, seq(1,49,1)), auc(acetate_median, seq(1,49,1)), auc(hydroxybutanoate_median, seq(1,49,1)), auc(acetylglucosamine_median, seq(1,49,1)),
                      auc(acetylneuraminate_median, seq(1,49,1)), auc(fructose_median, seq(1,49,1)), 
                      auc(mannitol_median, seq(1,49,1)), auc(salicin_median, seq(1,49,1)), auc(bhi_median, seq(1,49,1)), auc(no_aa_median, seq(1,49,1)), auc(no_carb_median, seq(1,49,1))), digits=3)

# Assemble the table
growth_summary <- cbind(substrates, max_rate, time_max_rate, max_od, time_max_od, rate_24_hrs, mean_rate, area_under, corrected_p_values)
colnames(growth_summary) <- c('Substrate', 'Max_Growth_Rate', 'Time_of_Max_Rate_in_Hours', 'Max_OD', 'Time_of_Max_OD_in_Hours', 'Rate_at_24_hours', 'Mean_Rate', 'AUC', 'Corrected_AoV_pvalue')
rm(substrates, max_rate, time_max_rate, max_od, time_max_od, rate_24_hrs, mean_rate, area_under, corrected_p_values)

rm(glucose_median, acetate_median, hydroxybutanoate_median, acetylglucosamine_median, acetylneuraminate_median, sorbitol_median, fructose_median, mannitol_median, salicin_median, no_carb_median, no_aa_median)
rm(glucose, acetate, hydroxybutanoate, acetylglucosamine, acetylneuraminate, sorbitol, fructose, mannitol, salicin, no_carb, no_aa)
rm(glucose_sd, acetate_sd, hydroxybutanoate_sd, acetylglucosamine_sd, acetylneuraminate_sd, sorbitol_sd, fructose_sd, mannitol_sd, salicin_sd, no_carb_sd, no_aa_sd)

# Write growth summary data to supplementary table
table_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/tables/table_S6.tsv'
write.table(growth_summary, file=table_file, quote=FALSE, sep='\t', row.names=FALSE)
rm(table_file, growth_summary)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_5.pdf'
pdf(file=plot_file, width=16, height=15)
layout(matrix(c(1,2,2,2,2,3,3,3,
                4,2,2,2,2,3,3,3,
                5,2,2,2,2,3,3,3,
                6,2,2,2,2,3,3,3,
                7,7,7,8,8,8,8,8,
                7,7,7,8,8,8,8,8,
                7,7,7,8,8,8,8,8,
                7,7,7,8,8,8,8,8), nrow=8, ncol=8, byrow=TRUE))

plot.new()
text(x=0.9, y=0.1, 'a', cex=2.1, font=2)

#---------------------------------------#

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
segments(x0=-1.1, y0=0, x1=-0.5, y1=0, lwd=2)

legend(x=0.91, y=1.3, legend=c('Enzyme node', 'Metabolite node'),
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

# Shared metabolite importances
par(mar=c(4,4,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(shared_importance$Metabolite_score, labels=shared_importance$Compound_name, lcolor=NA, cex=1.2, color='black',
         xlab='Median Metabolite Importance Score', xlim=c(0,12), pch=19)
segments(x0=rep(0, 13), y0=c(1:13), x1=rep(12, 13), y1=c(1:13), lty=2) # Dotted lines
mtext('b', side=2, line=2, las=2, adj=2.5, padj=-18, cex=1.4, font=2)

#---------------------------------------#

plot.new()
text(x=0.22, y=0.46, 'C. difficile', cex=1.1, font=c(4,2))
text(x=0.734, y=0.46, '630 network', cex=1.1, font=2)
segments(x0=0.02, y0=0.39, x1=0.99, y1=0.39, lwd=2)
text(x=0.29, y=0.29, '- 447 Enzymes')
text(x=0.325, y=0.16, '- 758 Metabolites')
text(x=0.26, y=0.01, '- 2135 Edges')

#---------------------------------------#

# Large component of C. difficile 630 graph
par(mar=c(0,1,0,0))
plot(largest_simple_graph, vertex.label=NA, layout=optimal_layout1,
     edge.arrow.size=0.25, edge.arrow.width=0.4, vertex.frame.color='black')

# Draw box
rect(xleft=0.4, ybottom=-0.4, xright=0.8, ytop=-0.75, lwd=2)

plot(0, type='n', axes=F, xlab='', ylab='') # Empty plot

#---------------------------------------#

# Unique metabolite importances
par(mar=c(4,4,1,1), xaxs='i', xpd=FALSE, mgp=c(2,1,0))
dotchart(top_importances$Metabolite_score, labels=top_importances$Compound_name,
         lcolor=NA, cex=1.2, groups=top_importances$abx, color='black',
         xlab='Metabolite Importance Score', xlim=c(0,12), pch=19, lwd=3,
         gcolor=c(wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5],'forestgreen'))
mtext('c', side=2, line=2, las=2, adj=2.5, padj=-18, cex=1.4, font=2)
segments(x0=rep(0, 15), y0=c(1:11, 14:15, 18, 21), x1=rep(12, 15), y1=c(1:11, 14:15, 18, 21), lty=2) # Dotted lines

#---------------------------------------#

# Growth on important compounds (separate)
par(mar=c(7,7,1.5,2), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE, mgp=c(4,2,0))
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,50), ylim=c(-0.03,1.0), lwd=2, pch=15, xlab='Time (hours)', ylab=expression(OD[600]), cex=2.3)
abline(h=seq(0,1.0,0.1), lty=3, col='gray68') # adding gridlines
abline(v=seq(1,50,2), lty=3, col='gray68') # adding gridlines
axis(1, at=seq(1,49,4), labels=seq(0,24,2), tck=-0.018)
axis(2, at=seq(0.0,1.0,0.2), labels=c('0.0','0.2','0.4','0.6','0.8','1.0'), tck=-0.018)
mtext('d', side=2, line=2, las=2, adj=3.7, padj=-18, cex=1.4, font=2)

# Control
lines(growth_medians$no_carb_median, type='o', lwd=2.5, pch=16, cex=0, col='gray60')
segments(x0=seq(1,49,1), y0=growth_medians$no_carb_median+growth_sds$no_carb_sd, x1=seq(1,49,1), y1=growth_medians$no_carb_median-growth_sds$no_carb_sd, lwd=2.5, cex=2, col='gray60')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$no_carb_median+growth_sds$no_carb_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$no_carb_median+growth_sds$no_carb_sd, lwd=2.5, col='gray60')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$no_carb_median-growth_sds$no_carb_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$no_carb_median-growth_sds$no_carb_sd, lwd=2.5, col='gray60')

lines(growth_medians$acetylglucosamine_median, type='o', col='black', lwd=2.5, pch=6, cex=2.5)
segments(x0=seq(1,49,1), y0=growth_medians$acetylglucosamine_median+growth_sds$acetylglucosamine_sd, x1=seq(1,49,1), y1=growth_medians$acetylglucosamine_median-growth_sds$acetylglucosamine_sd, lwd=2.5, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetylglucosamine_median+growth_sds$acetylglucosamine_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetylglucosamine_median+growth_sds$acetylglucosamine_sd, lwd=2.5, col='black')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetylglucosamine_median-growth_sds$acetylglucosamine_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetylglucosamine_median-growth_sds$acetylglucosamine_sd, lwd=2.5, col='black')

lines(growth_medians$fructose_median, type='o', col=wes_palette('FantasticFox')[1], lwd=2.5, pch=0, cex=2.5)
segments(x0=seq(1,49,1), y0=growth_medians$fructose_median+growth_sds$fructose_sd, x1=seq(1,49,1), y1=growth_medians$fructose_median-growth_sds$fructose_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$fructose_median+growth_sds$fructose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$fructose_median+growth_sds$fructose_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$fructose_median-growth_sds$fructose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$fructose_median-growth_sds$fructose_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])

lines(growth_medians$mannitol_median, type='o', col=wes_palette('FantasticFox')[3], lwd=2.5, pch=1, cex=2.5)
segments(x0=seq(1,49,1), y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, x1=seq(1,49,1), y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$mannitol_median+growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$mannitol_median-growth_sds$mannitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])

lines(growth_medians$salicin_median, type='o', col=wes_palette('FantasticFox')[5], lwd=2.5, pch=2, cex=2.5)
segments(x0=seq(1,49,1), y0=growth_medians$salicin_median+growth_sds$salicin_sd, x1=seq(1,49,1), y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$salicin_median+growth_sds$salicin_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$salicin_median+growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$salicin_median-growth_sds$salicin_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])

lines(growth_medians$acetylneuraminate_median, type='o', col='forestgreen', lwd=2.5, pch=5, cex=2.5)
segments(x0=seq(1,49,1), y0=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd, x1=seq(1,49,1), y1=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd, lwd=2.5, col='forestgreen')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetylneuraminate_median+growth_sds$acetylneuraminate_sd, lwd=2.5, col='forestgreen')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetylneuraminate_median-growth_sds$acetylneuraminate_sd, lwd=2.5, col='forestgreen')

legend('topleft', legend=c('No Carbohydrates','N-Acetyl-D-glucosamine','D-Fructose','Mannitol','Salicin','N-Acetylneuriminate'), 
       col=c('gray60','black',wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5],'forestgreen'), 
       pch=c(16,6,0,1,2,5), cex=2, pt.cex=c(0,2.9,2.9,2.9,2.9,2.9), bg='white', lwd=2.5)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_S6.pdf'
pdf(file=plot_file, width=14, height=10)
par(mar=c(7,7,1.5,2), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE, mgp=c(4,2,0))
plot(0, type='n', xaxt='n', yaxt='n', xlim=c(0,50), ylim=c(-0.03,0.83), lwd=2, pch=15, xlab='Time (hours)', ylab=expression(OD[600]), cex=2.3)
abline(h=seq(0,0.9,0.1), lty=3, col='gray68') # adding gridlines
abline(v=seq(1,50,2), lty=3, col='gray68') # adding gridlines
axis(1, at=seq(1,49,4), labels=seq(0,24,2), tck=-0.018)
axis(2, at=seq(0.0,0.8,0.2), labels=c('0.0','0.2','0.4','0.6','0.8'), tck=-0.018)

lines(growth_medians$no_aa_median, type='o', lwd=2.5, pch=15, cex=0, col='gray60')
segments(x0=seq(1,49,1), y0=growth_medians$no_aa_median+growth_sds$no_aa_sd, x1=seq(1,49,1), y1=growth_medians$no_aa_median-growth_sds$no_aa_sd, lwd=2.5, cex=2, col='gray60')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$no_aa_median+growth_sds$no_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$no_aa_median+growth_sds$no_aa_sd, lwd=2.5, col='gray60')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$no_aa_median-growth_sds$no_aa_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$no_aa_median-growth_sds$no_aa_sd, lwd=2.5, col='gray60')

lines(bhi_median, type='o', col='darkorchid3', lwd=2.5, pch=15, cex=0.8)
segments(x0=seq(1,49,1), y0=bhi_median+bhi_sd, x1=seq(1,49,1), y1=bhi_median-bhi_sd, lwd=2.5, col='darkorchid3')
segments(x0=seq(1,49,1)-0.2, y0=bhi_median+bhi_sd, x1=seq(1,49,1)+0.2, y1=bhi_median+bhi_sd, lwd=2.5, col='darkorchid3')
segments(x0=seq(1,49,1)-0.2, y0=bhi_median-bhi_sd, x1=seq(1,49,1)+0.2, y1=bhi_median-bhi_sd, lwd=2.5, col='darkorchid3')

lines(growth_medians$glucose_median, type='o', lwd=2.5, pch=15, cex=1.5, col='dodgerblue4')
segments(x0=seq(1,49,1), y0=growth_medians$glucose_median+growth_sds$glucose_sd, x1=seq(1,49,1), y1=growth_medians$glucose_median-growth_sds$glucose_sd, lwd=2.5, cex=2, col='dodgerblue4')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$glucose_median+growth_sds$glucose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$glucose_median+growth_sds$glucose_sd, lwd=2.5, col='dodgerblue4')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$glucose_median-growth_sds$glucose_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$glucose_median-growth_sds$glucose_sd, lwd=2.5, col='dodgerblue4')

#lines(growth_medians$sorbitol_median, type='o', col=wes_palette('FantasticFox')[1], lwd=2.5, pch=16, cex=1.5)
#segments(x0=seq(1,49,1), y0=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, x1=seq(1,49,1), y1=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
#segments(x0=seq(1,49,1)-0.2, y0=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
#segments(x0=seq(1,49,1)-0.2, y0=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])

lines(growth_medians$hydroxybutanoate_median, type='o', col=wes_palette('FantasticFox')[5], lwd=2.5, pch=17, cex=1.5)
segments(x0=seq(1,49,1), y0=growth_medians$hydroxybutanoate_median+growth_sds$hydroxybutanoate_sd, x1=seq(1,49,1), y1=growth_medians$hydroxybutanoate_median-growth_sds$hydroxybutanoate_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$hydroxybutanoate_median+growth_sds$hydroxybutanoate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$hydroxybutanoate_median+growth_sds$hydroxybutanoate_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$hydroxybutanoate_median-growth_sds$hydroxybutanoate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$hydroxybutanoate_median-growth_sds$hydroxybutanoate_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])

lines(growth_medians$acetate_median, type='o', col='forestgreen', lwd=2.5, pch=18, cex=1.5)
segments(x0=seq(1,49,1), y0=growth_medians$acetate_median+growth_sds$acetate_sd, x1=seq(1,49,1), y1=growth_medians$acetate_median-growth_sds$acetate_sd, lwd=2.5, col='forestgreen')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetate_median+growth_sds$acetate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetate_median+growth_sds$acetate_sd, lwd=2.5, col='forestgreen')
segments(x0=seq(1,49,1)-0.2, y0=growth_medians$acetate_median-growth_sds$acetate_sd, x1=seq(1,49,1)+0.2, y1=growth_medians$acetate_median-growth_sds$acetate_sd, lwd=2.5, col='forestgreen')

legend('topleft', legend=c('No Amino acids','BHI','Glucose','4-Hydroxybutanoic acid','Acetate'), 
       col=c('gray60','darkorchid3','dodgerblue4',wes_palette('FantasticFox')[5],'forestgreen'), 
       pch=c(16,15,15,17,18), cex=1.3, pt.cex=c(0,1,2.5,2.5,2.5), bg='white', lwd=2.5)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
rm(optimal_layout1, optimal_layout2, largest_simple_graph, network, plot_file, top_importances, 
   shared_importance, growth_medians, growth_sds, format_curve, bhi, bhi_median, bhi_sd)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()


