
# Read in Monte Carlo ranges
cef_range_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/cefoperazone_630.bipartite.files/monte_carlo.score_range.tsv'
clinda_range_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/clindamycin_630.bipartite.files/monte_carlo.score_range.tsv'
strep_range_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/streptomycin_630.bipartite.files/monte_carlo.score_range.tsv'
gf_range_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/germfree_630.bipartite.files/monte_carlo.score_range.tsv'
cef_range <- as.matrix(read.delim(cef_range_file, header=TRUE, sep='\t', row.names=1))
colnames(cef_range) <- NULL
clinda_range <- as.matrix(read.delim(clinda_range_file, header=TRUE, sep='\t', row.names=1))
colnames(clinda_range) <- NULL
strep_range <- as.matrix(read.delim(strep_range_file, header=TRUE, sep='\t', row.names=1))
colnames(strep_range) <- NULL
gf_range <- as.matrix(read.delim(gf_range_file, header=TRUE, sep='\t', row.names=1))
colnames(gf_range) <- NULL
rm(cef_range_file, clinda_range_file, strep_range_file, gf_range_file)

# Combine simulated ranges by metabolite ID
sim_ranges <- merge(cef_range, clinda_range, by='row.names')
rownames(sim_ranges) <- sim_ranges$Row.names
sim_ranges$Row.names <- NULL
sim_ranges <- merge(sim_ranges, strep_range, by='row.names')
rownames(sim_ranges) <- sim_ranges$Row.names
sim_ranges$Row.names <- NULL
sim_ranges <- merge(sim_ranges, gf_range, by='row.names')
rownames(sim_ranges) <- sim_ranges$Row.names
sim_ranges$Row.names <- NULL
rm(cef_range, clinda_range, strep_range, gf_range)

# Calculate new monte carlo ranges
temp_summary <- as.vector(quantile(sim_ranges[1,], probs=c(0.5,0.25,0.5,0.75,0.95)))
temp_lower <- temp_summary[3] - abs(1.58 * (temp_summary[2] / sqrt(ncol(sim_ranges))))
temp_upper <- temp_summary[3] + abs(1.58 * (temp_summary[4] / sqrt(ncol(sim_ranges))))
sim_summary <- c(temp_summary[3],temp_summary[2],temp_summary[4],temp_summary[1],temp_summary[5],temp_lower,temp_upper)
for (x in seq(2, nrow(sim_ranges), 1)){
  temp_summary <- as.vector(quantile(sim_ranges[x,], probs=c(0.5,0.25,0.5,0.75,0.95)))
  temp_lower <- temp_summary[3] - abs(1.58 * (temp_summary[2] / sqrt(ncol(sim_ranges))))
  temp_upper <- temp_summary[3] + abs(1.58 * (temp_summary[4] / sqrt(ncol(sim_ranges))))
  sim_summary <- rbind(sim_summary, c(temp_summary[3],temp_summary[2],temp_summary[4],temp_summary[1],temp_summary[5],temp_lower,temp_upper))
}
sim_summary <- as.matrix(sim_summary)
colnames(sim_summary) <- c('Sim_Median', 'Sim_iqr_25', 'Sim_iqr_75', 'Sim_Lower_sig', 'Sim_Upper_sig', 'Sim_Lower95', 'Sim_Upper95')
rownames(sim_summary) <- rownames(sim_ranges)
write.table(sim_summary, file="~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metabolic_models/combined_sim.tsv", sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)

