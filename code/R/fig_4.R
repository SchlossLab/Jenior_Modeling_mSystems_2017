
# Load dependencies
deps <- c('vegan', 'igraph', 'ggplot2', 'shape', 'wesanderson', 'matrixStats');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

# Define file variables for network plot
network_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/bipartite_graph.txt'
ko_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/cefoperazone_630.mapping.txt'
layout_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/optimal_layout.tsv'

# Read in metabolic network data
network <- read.table(network_file, header=FALSE, sep='\t')
ko <- read.table(ko_file, header=TRUE, sep='\t')
optimal_layout1 <- as.matrix(read.table(layout_file, header=TRUE, sep='\t', row.names=1))
rm(network_file, ko_file, layout_file)

# Format directed graph
raw_graph <- graph.data.frame(network, directed=TRUE)
rm(network)

# Remove loops and multiple edges to make visualzation easier
simple_graph <- simplify(raw_graph)
rm(raw_graph)

# Print a summary of nodes and edges for entire graph
summary(simple_graph)
print(length(as.vector(grep('K', V(simple_graph)$name, value=TRUE))))
print(length(as.vector(grep('C', V(simple_graph)$name, value=TRUE))))

# Decompose graph
decomp_simple_graph <- decompose.graph(simple_graph)
rm(simple_graph)

# Get largest component
largest_component <- which.max(sapply(decomp_simple_graph, vcount))
largest_simple_graph <- decomp_simple_graph[[largest_component]]
ko_simp <- as.vector(grep('K', V(largest_simple_graph)$name, value=TRUE))
substrate_simp <- as.vector(grep('C', V(largest_simple_graph)$name, value=TRUE))
nodes <- c(ko_simp, substrate_simp)
rm(decomp_simple_graph, largest_component)

# Print a summary of nodes and edges for large component
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
rm(ko, ko_simp, substrate_simp, node_size, mappings, nodes)
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

cef_importance <- cef_importance[c(1:100),]
clinda_importance <- clinda_importance[c(1:100),]
strep_importance <- strep_importance[c(1:100),]
gf_importance <- gf_importance[c(1:100),]

shared_importance <- as.data.frame(subset(cef_importance, (cef_importance[,1] %in% clinda_importance[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% clinda_importance[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% strep_importance[,1])))
shared_importance <- as.data.frame(subset(shared_importance, (shared_importance[,1] %in% gf_importance[,1])))

shared_importance$KEGG_code <- rownames(shared_importance)
write.table(shared_importance, file='~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/table_S3.tsv', quote=FALSE, sep='\t', row.names=FALSE)

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
cef_only_importance$color <- wes_palette("FantasticFox")[3]
clinda_only_importance$abx <- 'Clindamycin'
clinda_only_importance$color <- wes_palette("FantasticFox")[5]
strep_only_importance$abx <- 'Streptomycin'
strep_only_importance$color <- wes_palette("FantasticFox")[1]
gf_only_importance$abx <- 'Gnotobiotic'
gf_only_importance$color <- 'black'

top_importances <- rbind(cef_only_importance[,c(1,2,3,4,6,7)], clinda_only_importance[,c(1,2,3,4,6,7)], 
                         strep_only_importance[,c(1,2,3,4,6,7)], gf_only_importance[,c(1,2,3,4,6,7)])
top_importances$abx <- as.factor(top_importances$abx)
top_importances$abx <- ordered(top_importances$abx, levels=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Gnotobiotic'))
top_importances$Compound_name <- gsub('_',' ',top_importances$Compound_name)
top_importances$Compound_name <- gsub('beta','b',top_importances$Compound_name)
top_importances$Compound_name <- gsub('\\-phosphate','',top_importances$Compound_name)
top_importances$Compound_name <- gsub(' 6','',top_importances$Compound_name)
top_importances$Compound_name <- gsub(' 5','',top_importances$Compound_name)
top_importances$Compound_name <- gsub(' 3','',top_importances$Compound_name)
rm(cef_only_importance, clinda_only_importance, strep_only_importance, gf_only_importance)
top_importances <- top_importances[ !(rownames(top_importances) %in% c('C11436')), ]

#-------------------------------------------------------------------------------------------------------------------------------------#

# Read in growth rate data
# Define variables
growth_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/formatted_growth.tsv'

# Read in data
growth <- read.delim(growth_file, sep='\t', header=TRUE, row.names=1)
growth <- as.data.frame(t(growth))
rm(growth_file)

# Seperate to groups of each growth substrate and subset to the first 12 hours
sorbitol <- cbind(growth$B9, growth$B10, growth$B11)[1:25,]
galactitol <- cbind(growth$C9, growth$C10, growth$C11)[1:25,]
starch <- cbind(growth$D9, growth$D10, growth$D11)[1:25,]
fructose <- cbind(growth$E9, growth$E10, growth$E11)[1:25,]
combination <- cbind(growth$G3, growth$G4, growth$G5)[1:25,]
mannitol <- cbind(growth$F9, growth$F10, growth$F11)[1:25,]
salicin <- cbind(growth$G9, growth$G10, growth$G11)[1:25,]
y_glucose_y_aa <- cbind(growth$B3, growth$B4, growth$B5)[1:25,]
n_glucose_y_aa <- cbind(growth$D3, growth$D4, growth$D5)[1:25,]
y_glucose_n_aa <- cbind(growth$C3, growth$C4, growth$C5)[1:25,]
n_glucose_n_aa <- cbind(growth$E3, growth$E4, growth$E5)[1:25,]
bhi <- cbind(growth$F3, growth$F4, growth$F5)[1:25,]

#-------------------------------------------------------------------------------------------------------------------------------------#

# Prepare data for statistical tests
control_test <- c()
for (time in 1:25){
  temp <- cbind('control', time, n_glucose_y_aa[time,])
  control_test <- rbind(control_test, temp)
}

fructose_test <- c()
for (time in 1:25){
  temp <- cbind('fructose', time, fructose[time,])
  fructose_test <- rbind(fructose_test, temp)
}
fructose_test <- as.data.frame(rbind(control_test, fructose_test))
colnames(fructose_test) <- c('substrate','time','od')
fructose_test$od <- as.numeric(as.character(fructose_test$od))

sorbitol_test <- c()
for (time in 1:25){
  temp <- cbind('sorbitol', time, sorbitol[time,])
  sorbitol_test <- rbind(sorbitol_test, temp)
}
sorbitol_test <- as.data.frame(rbind(control_test, sorbitol_test))
colnames(sorbitol_test) <- c('substrate','time','od')
sorbitol_test$od <- as.numeric(as.character(sorbitol_test$od))

galactitol_test <- c()
for (time in 1:25){
  temp <- cbind('galactitol', time, galactitol[time,])
  galactitol_test <- rbind(galactitol_test, temp)
}
galactitol_test <- as.data.frame(rbind(control_test, galactitol_test))
colnames(galactitol_test) <- c('substrate','time','od')
galactitol_test$od <- as.numeric(as.character(galactitol_test$od))

starch_test <- c()
for (time in 1:25){
  temp <- cbind('starch', time, starch[time,])
  starch_test <- rbind(starch_test, temp)
}
starch_test <- as.data.frame(rbind(control_test, starch_test))
colnames(starch_test) <- c('substrate','time','od')
starch_test$od <- as.numeric(as.character(starch_test$od))

mannitol_test <- c()
for (time in 1:25){
  temp <- cbind('mannitol', time, mannitol[time,])
  mannitol_test <- rbind(mannitol_test, temp)
}
mannitol_test <- as.data.frame(rbind(control_test, mannitol_test))
colnames(mannitol_test) <- c('substrate','time','od')
mannitol_test$od <- as.numeric(as.character(mannitol_test$od))

salicin_test <- c()
for (time in 1:25){
  temp <- cbind('salicin', time, salicin[time,])
  salicin_test <- rbind(salicin_test, temp)
}
salicin_test <- as.data.frame(rbind(control_test, salicin_test))
colnames(salicin_test) <- c('substrate','time','od')
salicin_test$od <- as.numeric(as.character(salicin_test$od))
rm(temp, time, control_test)

# Calculate differences
fructose_sig <- aov(formula=od ~ substrate * time, data=fructose_test)
summary(fructose_sig) # p < 2e-16 ***, corrected = 1.2e-15 ***
sorbitol_sig <- aov(formula=od ~ substrate * time, data=sorbitol_test)
summary(sorbitol_sig) # p = 0.02474 *, corrected =  0.148 n.s.
galactitol_sig <- aov(formula=od ~ substrate * time, data=galactitol_test)
summary(galactitol_sig) # p < 2e-16 ***, corrected = 1.2e-15 ***
starch_sig <- aov(formula=od ~ substrate * time, data=starch_test)
summary(starch_sig) # n.s.
mannitol_sig <- aov(formula=od ~ substrate * time, data=mannitol_test)
summary(mannitol_sig) # p < 2e-16 ***, corrected = 1.2e-15 ***
salicin_sig <- aov(formula=od ~ substrate * time, data=salicin_test)
summary(salicin_sig) # p < 2e-16 ***, corrected = 1.2e-15 ***

p_values <- c(2e-16,0.02474,2e-16,0.663,2e-16,2e-16)
corrected_p_values <- p.adjust(p_values, method='bonferroni')

# Clean up
rm(p_values, corrected_p_values)
rm(fructose_test, sorbitol_test, galactitol_test, starch_test, mannitol_test, salicin_test)
rm(fructose_sig, sorbitol_sig, galactitol_sig, starch_sig, mannitol_sig, salicin_sig)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format growth curves

# Find medians of treatement groups and subtract blanks
sorbitol_median <- rowMedians(sorbitol, na.rm=TRUE) - growth$B8[1:25]
sorbitol_median[sorbitol_median < 0] <- 0
galactitol_median <- rowMedians(galactitol, na.rm=TRUE) - growth$C8[1:25]
galactitol_median[galactitol_median < 0] <- 0
starch_median <- rowMedians(starch, na.rm=TRUE) - growth$D8[1:25]
starch_median[starch_median < 0] <- 0
fructose_median <-  rowMedians(fructose, na.rm=TRUE) - growth$E8[1:25]
fructose_median[fructose_median < 0] <- 0
combination_median <- rowMedians(combination, na.rm=TRUE) - growth$G2[1:25]
combination_median[combination_median < 0] <- 0
mannitol_median <- rowMedians(mannitol, na.rm=TRUE) - growth$F8[1:25]
mannitol_median[mannitol_median < 0] <- 0
salicin_median <- rowMedians(salicin, na.rm=TRUE) - growth$G7[1:25]
salicin_median[salicin_median < 0] <- 0
y_glucose_y_aa_median <- rowMedians(y_glucose_y_aa, na.rm=TRUE) - growth$B2[1:25]
y_glucose_y_aa_median[y_glucose_y_aa_median < 0] <- 0
n_glucose_y_aa_median <- rowMedians(n_glucose_y_aa, na.rm=TRUE) - growth$D2[1:25]
n_glucose_y_aa_median[n_glucose_y_aa_median < 0] <- 0
y_glucose_n_aa_median <- rowMedians(y_glucose_n_aa, na.rm=TRUE) - growth$C2[1:25]
y_glucose_n_aa_median[y_glucose_n_aa_median < 0] <- 0
n_glucose_n_aa_median <- rowMedians(n_glucose_n_aa, na.rm=TRUE) - growth$E2[1:25]
n_glucose_n_aa_median[n_glucose_n_aa_median < 0] <- 0
bhi_median <- rowMedians(bhi, na.rm=TRUE) - growth$F2[1:25]
bhi_median[bhi_median < 0] <- 0
growth_medians <- as.data.frame(rbind(sorbitol_median, galactitol_median, starch_median, fructose_median, combination_median, mannitol_median, salicin_median, y_glucose_y_aa_median, n_glucose_y_aa_median, y_glucose_n_aa_median, n_glucose_n_aa_median, bhi_median))

# Determine some features of the 12 hour growth curves
# Maximum growth rate
diff(sorbitol_median)[which.max(abs(diff(sorbitol_median)))] # 0.022
diff(galactitol_median)[which.max(abs(diff(galactitol_median)))] # 0.02
diff(starch_median)[which.max(abs(diff(starch_median)))] # 0.025
diff(fructose_median)[which.max(abs(diff(fructose_median)))] # 0.089
diff(combination_median)[which.max(abs(diff(combination_median)))] # 0.095
diff(mannitol_median)[which.max(abs(diff(mannitol_median)))] # 0.044
diff(salicin_median)[which.max(abs(diff(salicin_median)))] # 0.049
diff(y_glucose_y_aa_median)[which.max(abs(diff(y_glucose_y_aa_median)))] # 0.085
diff(n_glucose_y_aa_median)[which.max(abs(diff(n_glucose_y_aa_median)))] # 0.028
diff(y_glucose_n_aa_median)[which.max(abs(diff(y_glucose_n_aa_median)))] # 0.006
diff(n_glucose_n_aa_median)[which.max(abs(diff(n_glucose_n_aa_median)))] # -0.003
diff(bhi_median)[which.max(abs(diff(bhi_median)))] # 0.057

# Time of maximum growth rate
which.max(abs(diff(sorbitol_median))) * 0.5 # 7.5 hours
which.max(abs(diff(galactitol_median))) * 0.5 # 8 hours
which.max(abs(diff(starch_median))) * 0.5 # 9 hours
which.max(abs(diff(fructose_median))) * 0.5 # 6.5 hours
which.max(abs(diff(combination_median))) * 0.5 # 7 hours
which.max(abs(diff(mannitol_median))) * 0.5 # 7.5 hours
which.max(abs(diff(salicin_median))) * 0.5 # 8.5 hours
which.max(abs(diff(y_glucose_y_aa_median))) * 0.5 # 7.5 hours
which.max(abs(diff(n_glucose_y_aa_median))) * 0.5 # 7.5 hours
which.max(abs(diff(y_glucose_n_aa_median))) * 0.5 # 0.5 hours
which.max(abs(diff(n_glucose_n_aa_median))) * 0.5 # 0.5 hours
which.max(abs(diff(bhi_median))) * 0.5 # 8.5 hours

# Maximum OD
max(sorbitol_median) # 0.199
max(galactitol_median) # 0.202
max(starch_median) # 0.204
max(fructose_median) # 0.556
max(combination_median) # 0.505
max(mannitol_median) # 0.43
max(salicin_median) # 0.549
max(salicin_median) # 0.575
max(y_glucose_y_aa_median) # 0.575
max(n_glucose_y_aa_median) # 0.211
max(y_glucose_n_aa_median) # 0.008
max(n_glucose_n_aa_median) # 0.005
max(bhi_median) # 0.662

# Growth rate at 12 hours
diff(sorbitol_median)[length(diff(sorbitol_median))] # 0.002
diff(galactitol_median)[length(diff(galactitol_median))] # 0.005
diff(starch_median)[length(diff(starch_median))] # 0.001
diff(fructose_median)[length(diff(fructose_median))] # 0.003
diff(combination_median)[length(diff(combination_median))] # 0.002
diff(mannitol_median)[length(diff(mannitol_median))] # 0.006
diff(salicin_median)[length(diff(salicin_median))] # 0.048
diff(y_glucose_y_aa_median)[length(diff(n_glucose_y_aa_median))] # 0.008
diff(n_glucose_y_aa_median)[length(diff(n_glucose_y_aa_median))] # 0.0
diff(y_glucose_n_aa_median)[length(diff(n_glucose_y_aa_median))] # 0.0
diff(n_glucose_n_aa_median)[length(diff(n_glucose_y_aa_median))] # 0.0
diff(bhi_median)[length(diff(bhi_median))] # 0.002

# Mean growth rate
mean(diff(sorbitol_median)) # 0.008291667
mean(diff(galactitol_median)) # 0.008333333
mean(diff(starch_median)) # 0.00825
mean(diff(fructose_median)) # 0.02304167
mean(diff(combination_median)) # 0.02091667
mean(diff(mannitol_median)) # 0.01775
mean(diff(salicin_median)) # 0.022125
mean(diff(n_glucose_y_aa_median)) # 0.02391667
mean(diff(y_glucose_y_aa_median)) # 0.008666667
mean(diff(y_glucose_n_aa_median)) # 0.0003333333
mean(diff(n_glucose_n_aa_median)) # -0.0002083333
mean(diff(bhi_median)) # 0.02758333

# Standard deviations
sorbitol_sd <- rowSds(sorbitol, na.rm=TRUE)
galactitol_sd <- rowSds(galactitol, na.rm=TRUE)
starch_sd <-  rowSds(starch, na.rm=TRUE)
fructose_sd <-  rowSds(fructose, na.rm=TRUE)
combination_sd <-  rowSds(combination, na.rm=TRUE)
mannitol_sd <- rowSds(mannitol, na.rm=TRUE)
salicin_sd <- rowSds(salicin, na.rm=TRUE)
y_glucose_y_aa_sd <- rowSds(y_glucose_y_aa, na.rm=TRUE)
n_glucose_y_aa_sd <- rowSds(n_glucose_y_aa, na.rm=TRUE)
y_glucose_n_aa_sd <- rowSds(y_glucose_n_aa, na.rm=TRUE)
n_glucose_n_aa_sd <- rowSds(n_glucose_n_aa, na.rm=TRUE)
bhi_sd <- rowSds(bhi, na.rm=TRUE)
growth_sds <- as.data.frame(rbind(sorbitol_sd, galactitol_sd, starch_sd, fructose_sd, combination_sd, mannitol_sd, salicin_sd, y_glucose_y_aa_sd, n_glucose_y_aa_sd, y_glucose_n_aa_sd, n_glucose_n_aa_sd, bhi_sd))
rm(sorbitol_median, galactitol_median, starch_median, fructose_median, combination_median, mannitol_median, salicin_median, y_glucose_y_aa_median, n_glucose_y_aa_median, y_glucose_n_aa_median, n_glucose_n_aa_median, bhi_median)
rm(sorbitol, galactitol, starch, fructose, combination, mannitol, salicin, y_glucose_y_aa, n_glucose_y_aa, y_glucose_n_aa, n_glucose_n_aa, bhi)
rm(sorbitol_sd, galactitol_sd, starch_sd, fructose_sd, combination_sd, mannitol_sd, salicin_sd, y_glucose_y_aa_sd, n_glucose_y_aa_sd, y_glucose_n_aa_sd, n_glucose_n_aa_sd, bhi_sd)
rm(growth)

# Subset for first 12 hours of assay
growth_medians <- as.data.frame(t(growth_medians[,1:25]))
growth_sds <- as.data.frame(t(growth_sds[,1:25]))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up plotting environment
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_4.pdf'
pdf(file=plot_file, width=17, height=16)
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
text(-0.9, 0, expression(Importance == paste(log[2],'( ',frac(Sigma * t[i], e[o]),' ','-',' ',frac(Sigma * t[o], e[i]),' )')), cex = 1.9) # Importance algorithm
text(x=-1, y=1.13, labels='Tetrathionate reductase', font=2, cex=1.6) # Enzyme 1 name
text(x=-1, y=1, labels='7', col='white', cex=1.25) # Enzyme 1 transcription
text(x=-0.5, y=-1.3, labels='Sulfate reductase', font=2, cex=1.6) # Enzyme 2 name
text(x=-0.5, y=-1, labels='94', col='white', cex=2.3) # Enzyme 2 transcription
text(x=1, y=0.75, labels='Thiosulfate Oxidase', font=2, cex=1.6) # Enzyme 3
text(x=0.99, y=0.44, labels='115', col='white', cex=2.5) # Enzyme 3 transcription
text(x=-0.17, y=0.14, 's', col='white', cex=2) # Substrate node label
text(x=c(0.1,0.1), y=c(-0.02,-0.12), labels=c('Thiosulfate','= 6.554'), cex=1.7, font=c(2,1)) # Compound & calculated importance
segments(x0=-0.05, y0=-0.17, x1=0.25, y1=-0.17, lwd=2)
legend(x=0.4, y=1.2, legend=c('KEGG ortholog', 'Reaction substrate'), 
       pt.bg=c('firebrick3', 'blue3'), col='black', pch=21, pt.cex=3, cex=2, bty='n')
text(x=-0.5, y=-1.1, expression(t[i]), col='white', cex=1.9) # labeled transcription for input reactions
text(x=0.99, y=0.34, expression(t[i]), col='white', cex=1.9)
text(x=-1, y=0.9, expression(t[o]), col='black', cex=1.9) # labeled transcription for output reactions
text(x=-0.6, y=0.7, expression(e[i]), col='black', cex=1.9) # labeled indegree
text(x=-0.4, y=-0.3, expression(e[o]), col='black', cex=1.9) # labeled outdegree
text(x=0.3, y=0.33, expression(e[o]), col='black', cex=1.9)
#mtext('B', side=2, line=2, las=2, adj=-2, padj=-16.5, cex=1.8)
Arrows(x0=0.63, y0=-0.7, x1=0.12, y1=-0.7, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.3) # Score explanation line
Arrows(x0=0.63, y0=-0.7, x1=1.14, y1=-0.7, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.3)
segments(x0=0.63, y0=-0.65, x1=0.63, y1=-0.75, lwd=3)
text(x=0.63, y=-0.83, '0', cex=2.1) 
#text(x=0.12, y=-0.8, expression(- infinity), cex=2.3)
#text(x=1.14, y=-0.8, expression(+ infinity), cex=2.3)
text(x=0.12, y=-0.81, '-', cex=2.6)
text(x=1.14, y=-0.81, '+', cex=2.6)
text(x=1.15, y=-0.6, 'More likely consumed', cex=1.6)
text(x=0.14, y=-0.6, 'More likely released', cex=1.6)
text(x=0.63, y=-0.95, 'Importance Score', cex=2.3, font=2) 

plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot

# Large component of C. difficile 630 graph
par(mar=c(1,3,1,1))
plot(largest_simple_graph, vertex.label=NA, layout=optimal_layout1,
     edge.arrow.size=0.2, edge.arrow.width=0.4, vertex.frame.color='black')
mtext('A', side=2, line=2, las=2, adj=-2, padj=-10, cex=1.8)

plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot
plot(1, type='n', axes=F, xlab='', ylab='') # Empty plot

# Boxes in A are drawn in seperate software (Gimp)

#-------------------------------------------------------------------------------------------------------------------------------------#

# B - Top compound importances
par(mar=c(4,3,1,1), xaxs='i')
dotchart(top_importances$Metabolite_score, labels=top_importances$Compound_name, 
         lcolor=NA, cex=1.5, groups=top_importances$abx, color='black', 
         xlab='Metabolite Importance Score', xlim=c(-2,10), 
         gcolor=c(wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5],'black'), pch=19)
segments(x0=rep(-2, 14), y0=c(1:8, 11, 14, 17:20), x1=rep(10, 14), y1=c(1:8, 11, 14, 17:20), lty=2)

# Add simulated means
points(x=top_importances[c(14:7),3], y=c(1:8), cex=2, col='black', pch='|') # Gnotobiotic
points(x=top_importances[2,3], y=11, cex=2, col='black', pch='|') # Clindamycin
points(x=top_importances[1,3], y=14, cex=2, col='black', pch='|') # Cefoperazone
points(x=top_importances[c(3:6),3], y=c(17:20), cex=2, col='black', pch='|') # Streptomycin

mtext('B', side=2, line=2, las=2, adj=0.5, padj=-15.5, cex=1.8)

#-------------------------------------------------------------------------------------------------------------------------------------#

# C - Growth on important compounds
par(mar=c(5,5,1,1), las=1, cex.lab=2, cex.axis=1.8, xpd=FALSE)

plot(growth_medians$y_glucose_y_aa_median, type='o', xaxt='n', xlim=c(0,29), ylim=c(-0.03,0.65), lwd=2, pch=15, xlab='Hours Postinoculation', ylab=expression(OD[600]), cex=2.3)
segments(x0=seq(1,25,1), y0=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, x1=seq(1,25,1), y1=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, lwd=2.5, cex=2)
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$y_glucose_y_aa_median+growth_sds$y_glucose_y_aa_sd, lwd=2)
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$y_glucose_y_aa_median-growth_sds$y_glucose_y_aa_sd, lwd=2)
abline(h=seq(0,0.6,0.1), lty=3, col='gray48') # adding gridlines
abline(v=seq(1,26,2), lty=3, col='gray48') # adding gridlines
lines(growth_medians$n_glucose_y_aa_median, type='o', lwd=2, pch=16, cex=2.3)
segments(x0=seq(1,25,1), y0=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, x1=seq(1,25,1), y1=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, lwd=2.5, cex=2)
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$n_glucose_y_aa_median+growth_sds$n_glucose_y_aa_sd, lwd=2)
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$n_glucose_y_aa_median-growth_sds$n_glucose_y_aa_sd, lwd=2)
lines(growth_medians$y_glucose_n_aa_median, type='o', lwd=2, pch=18, cex=2.3)
segments(x0=seq(1,25,1), y0=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, x1=seq(1,25,1), y1=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, lwd=2.5, cex=2)
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$y_glucose_n_aa_median+growth_sds$y_glucose_n_aa_sd, lwd=2)
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$y_glucose_n_aa_median-growth_sds$y_glucose_n_aa_sd, lwd=2)
lines(growth_medians$n_glucose_n_aa_median, type='o', lwd=2, pch=17, cex=2.7)
segments(x0=seq(1,25,1), y0=growth_medians$n_glucose_n_aa_median+growth_sds$n_glucose_n_aa_sd, x1=seq(1,25,1), y1=growth_medians$n_glucose_n_aa_median-growth_sds$n_glucose_n_aa_sd, lwd=2.5, cex=2)
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$n_glucose_n_aa_median+growth_sds$n_glucose_n_aa_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$n_glucose_n_aa_median+growth_sds$n_glucose_n_aa_sd, lwd=2)
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$n_glucose_n_aa_median-growth_sds$n_glucose_n_aa_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$n_glucose_n_aa_median-growth_sds$n_glucose_n_aa_sd, lwd=2)

lines(growth_medians$sorbitol_median, type='o', col=wes_palette('FantasticFox')[1], lwd=2, pch=15, cex=2.3)
segments(x0=seq(1,25,1), y0=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, x1=seq(1,25,1), y1=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$sorbitol_median+growth_sds$sorbitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$sorbitol_median-growth_sds$sorbitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
lines(growth_medians$galactitol_median, type='o', col=wes_palette('FantasticFox')[1], lwd=2, pch=16, cex=2.3)
segments(x0=seq(1,25,1), y0=growth_medians$galactitol_median+growth_sds$galactitol_sd, x1=seq(1,25,1), y1=growth_medians$galactitol_median-growth_sds$galactitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$galactitol_median+growth_sds$galactitol_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$galactitol_median+growth_sds$galactitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$galactitol_median-growth_sds$galactitol_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$galactitol_median-growth_sds$galactitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
lines(growth_medians$starch_median, type='o', col=wes_palette('FantasticFox')[1], lwd=2, pch=18, cex=2.3)
segments(x0=seq(1,25,1), y0=growth_medians$starch_median+growth_sds$starch_sd, x1=seq(1,25,1), y1=growth_medians$starch_median-growth_sds$starch_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$starch_median+growth_sds$starch_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$starch_median+growth_sds$starch_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$starch_median-growth_sds$starch_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$starch_median-growth_sds$starch_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
lines(growth_medians$fructose_median, type='o', col=wes_palette('FantasticFox')[1], lwd=2, pch=17, cex=3)
segments(x0=seq(1,25,1), y0=growth_medians$fructose_median+growth_sds$fructose_sd, x1=seq(1,25,1), y1=growth_medians$fructose_median-growth_sds$fructose_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$fructose_median+growth_sds$fructose_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$fructose_median+growth_sds$fructose_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$fructose_median-growth_sds$fructose_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$fructose_median-growth_sds$fructose_sd, lwd=2.5, col=wes_palette('FantasticFox')[1])
#lines(growth_medians$combination_median, type='o', col=color_palette[5], lwd=2, pch=19, cex=1.5)
#segments(x0=seq(1,25,1), y0=growth_medians$combination_median+growth_sds$combination_sd, x1=seq(1,25,1), y1=growth_medians$combination_median-growth_sds$combination_sd, lwd=2.5, col=color_palette[5])
#segments(x0=seq(1,25,1)-0.2, y0=growth_medians$combination_median+growth_sds$combination_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$combination_median+growth_sds$combination_sd, lwd=2.5, col=color_palette[5])
#segments(x0=seq(1,25,1)-0.2, y0=growth_medians$combination_median-growth_sds$combination_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$combination_median-growth_sds$combination_sd, lwd=2.5, col=color_palette[5])
lines(growth_medians$mannitol_median, type='o', col=wes_palette('FantasticFox')[3], lwd=2, pch=19, cex=2)
segments(x0=seq(1,25,1), y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, x1=seq(1,25,1), y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$mannitol_median+growth_sds$mannitol_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$mannitol_median+growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$mannitol_median-growth_sds$mannitol_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$mannitol_median-growth_sds$mannitol_sd, lwd=2.5, col=wes_palette('FantasticFox')[3])
lines(growth_medians$salicin_median, type='o', col=wes_palette('FantasticFox')[5], lwd=2.5, pch=19, cex=2)
segments(x0=seq(1,25,1), y0=growth_medians$salicin_median+growth_sds$salicin_sd, x1=seq(1,25,1), y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$salicin_median+growth_sds$salicin_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$salicin_median+growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])
segments(x0=seq(1,25,1)-0.2, y0=growth_medians$salicin_median-growth_sds$salicin_sd, x1=seq(1,25,1)+0.2, y1=growth_medians$salicin_median-growth_sds$salicin_sd, lwd=2.5, col=wes_palette('FantasticFox')[5])

axis(1, at=seq(1,26,4), labels=seq(0,12,2))

legend('topleft', legend=c('+Glucose +AA','-Glucose +AA','+Glucose -AA','-Glucose -AA','D-Sorbitol','Galactitol','Starch','D-Fructose','Mannitol','Salicin'), 
       col=c('black','black','black','black', wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[3],wes_palette('FantasticFox')[5]), 
       pch=c(15,16,18,17,15,16,18,17,19,19), cex=2, pt.cex=3.1, bg='white')

segments(x0=c(26,27,28), y0=c(0.556,0.549,0.430), x1=c(26,27,28), y1=c(0.211,0.211,0.211), lwd=2.5)

segments(x0=c(25.7,26.7,27.7), y0=c(0.556,0.549,0.430), x1=c(26.3,27.3,28.3), y1=c(0.556,0.549,0.430), lwd=3.5, col=c(wes_palette('FantasticFox')[1],wes_palette('FantasticFox')[5],wes_palette('FantasticFox')[3]))
segments(x0=c(25.7,26.7,27.7), y0=c(0.211,0.211,0.211), x1=c(26.3,27.3,28.3), y1=c(0.211,0.211,0.211), lwd=3.5, col='black')
text(x=c(26.3,27.3,28.3), y=c(0.384,0.38,0.321), labels=c('***','***','***'), cex=2.3, srt = 90)

mtext('C', side=2, line=2, las=2, adj=1, padj=-16, cex=1.8)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
rm(optimal_layout1, optimal_layout2, largest_simple_graph, network, plot_file, growth_sds, growth_medians, top_importances)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()

