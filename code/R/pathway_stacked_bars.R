
# Install packages
#install.packages("RColorBrewer")
library(RColorBrewer)

#--------------------------------------------------------------------------------------------------------------#

# Define variables
cef <- '/Users/schloss/Desktop/cefoperazone_630.mapped2cdf630.annotated.txt'
strep <- '/Users/schloss/Desktop/clindamycin_630.mapped2cdf630.annotated.txt'
clinda <- '/Users/schloss/Desktop/streptomycin_630.mapped2cdf630.annotated.txt'
gf <- '/Users/schloss/Desktop/germfree.mapped2cdf630.annotated.txt'
figure_file <- '/Users/schloss/Desktop/pathway_relabund.pdf'

#--------------------------------------------------------------------------------------------------------------#

# Define formatting functions
format_expression_aggregate <- function(file_name){

  # Load in and format data
  expression <- read.delim(file_name, sep='\t', header=F)
  colnames(expression) <- c('Normalized_abundance', 'Gene_code', 'Gene_name', 'KO', 'Pathway_code', 'Pathway_name', 'Pathway_family1', 'Pathway_family2')
  path_expression <- subset(expression, Pathway_code != 'path_key_error')
  #path_expression <- subset(expression, Pathway_code != 'pathway_unannotated')
  path_expression$Relabund <- (path_expression$Normalized_abundance  / sum(path_expression$Normalized_abundance)) * 100

  # Pool data in same pathway category and sort
  pool_path_expression <- aggregate(Relabund ~ Pathway_name, FUN = sum, data = path_expression)
  pool_path_expression <- pool_path_expression[order(-pool_path_expression$Relabund), ] 
  rownames(pool_path_expression) <- pool_path_expression$Pathway_name
  pool_path_expression$Pathway_name <- NULL
  rows <- rownames(pool_path_expression)
  
  # Bin low abundant catagories into 'Other'
  other <- sum(pool_path_expression$Relabund[which(pool_path_expression$Relabund < 0.3)])
  pool_path_expression[pool_path_expression < 0.3] <- 0
  pool_path_expression <- rbind(pool_path_expression, other)
  rows <- append(rows, 'Other')
  rownames(pool_path_expression) <- rows
  
  return(pool_path_expression)

}


format_expression<- function(file_name){
  expression <- read.delim(file_name, sep='\t', header=F)
  colnames(expression) <- c('Normalized_abundance', 'Gene_code', 'Gene_name', 'KO', 'Pathway_code', 'Pathway_name', 'Pathway_family1', 'Pathway_family2')
  path_expression <- subset(expression, Pathway_code != 'path_key_error')
  #path_expression <- subset(expression, Pathway_code != 'pathway_unannotated')
  path_expression$Relabund <- (path_expression$Normalized_abundance  / sum(path_expression$Normalized_abundance)) * 100
  rownames(path_expression) <- path_expression$Gene_code
  #path_expression$Gene_code <- NULL
  path_expression <- path_expression[order(-path_expression$Normalized_abundance), ]
  #path_expression$Gene_name <- NULL
  #path_expression$KO <- NULL
  #path_expression$Pathway_code <- NULL
  #path_expression$Pathway_name <- NULL
  #path_expression$Pathway_family1 <- NULL
  #path_expression$Pathway_family2 <- NULL
  #path_expression$Relabund <- NULL
  return(path_expression)
}

#--------------------------------------------------------------------------------------------------------------#

# Prep all of the data - aggregating by pathways
cef_expression <- format_expression_aggregate(cef)
clinda_expression <- format_expression_aggregate(clinda)
strep_expression <- format_expression_aggregate(strep)
gf_expression <- format_expression_aggregate(gf)

# Combine with other expression tables
all_expression <- merge(cef_expression , clinda_expression, by = 'row.names')
rownames(all_expression) <- all_expression$Row.names
all_expression$Row.names <- NULL
colnames(all_expression) <- c('cefoperazone', 'clindamycin')
all_expression <- merge(all_expression , strep_expression, by = 'row.names')
rownames(all_expression) <- all_expression$Row.names
all_expression$Row.names <- NULL
colnames(all_expression) <- c('cefoperazone', 'clindamycin', 'streptomycin')
all_expression <- merge(all_expression , gf_expression, by = 'row.names')
rownames(all_expression) <- all_expression$Row.names
all_expression$Row.names <- NULL
colnames(all_expression) <- c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Germfree')

# Remove categories where all of the values are 0 and sort
filter_all_expression <- all_expression[rowSums(all_expression[, -1]) > 0, ]
filter_all_expression <- filter_all_expression[order(rowSums(filter_all_expression), decreasing = T), ]
filter_all_expression <- as.matrix(filter_all_expression)

# Define color palette
color_pal <- rainbow(nrow(filter_all_expression))

# Plot the final aggregated data
pdf(file=figure_file, width=14, height=9)
par(las=1, mar=c(2.2, 3, 1, 20), mgp=c(2, 0.8, 0), xpd=TRUE)
barplot(filter_all_expression, col = color_pal, yaxt = 'n', ylim = c(-0.5, 100.5), ylab = '% Relative Abundance')
axis(side = 2, at = seq(0, 100, 20), tick = TRUE)
segments(x0 = c(0, 0, 5), y0 = c(100, 0, 100), x1 = c(5, 5, 5), y1 = c(100, 0, 0))
segments(x0 = c(0, 0, 0, 0), y0 = c(20, 40, 60, 80), x1 = c(5, 5, 5, 5), y1 = c(20, 40, 60, 80), lty = 2)
legend(5, 85, legend = rev(rownames(filter_all_expression)), pt.bg = rev(color_pal), pch = 22, pt.cex = 2, cex = 1)
dev.off()

#--------------------------------------------------------------------------------------------------------------#


# Prep all of the data - no aggregation
cef_expression_all <- format_expression(cef)
clinda_expression_all <- format_expression(clinda)
strep_expression_all <- format_expression(strep)
gf_expression_all <- format_expression(gf)

# Filter out the genes you aren't interested in
#cef_expression_all_filter <- subset(cef_expression_all, Pathway_name == 'Microbial metabolism in diverse environments')
#clinda_expression_all_filter <- subset(clinda_expression_all, Pathway_name == 'Microbial metabolism in diverse environments')
#strep_expression_all_filter <- subset(strep_expression_all, Pathway_name == 'Microbial metabolism in diverse environments')
#gf_expression_all_filter <- subset(gf_expression_all, Pathway_name == 'Microbial metabolism in diverse environments')

# Remove useless columns
cef_expression_all$Normalized_abundance <- NULL
cef_expression_all$KO <- NULL
cef_expression_all$Pathway_code <- NULL
cef_expression_all$Pathway_name <- NULL
cef_expression_all$Pathway_family1 <- NULL
cef_expression_all$Pathway_family2 <- NULL
cef_expression_all$Gene_code <- NULL
cef_expression_all$Gene_name <- NULL

clinda_expression_all$Normalized_abundance <- NULL
clinda_expression_all$KO <- NULL
clinda_expression_all$Pathway_code <- NULL
clinda_expression_all$Pathway_name <- NULL
clinda_expression_all$Pathway_family1 <- NULL
clinda_expression_all$Pathway_family2 <- NULL
clinda_expression_all$Gene_code <- NULL
clinda_expression_all$Gene_name <- NULL

strep_expression_all$Normalized_abundance <- NULL
strep_expression_all$KO <- NULL
strep_expression_all$Pathway_code <- NULL
strep_expression_all$Pathway_name <- NULL
strep_expression_all$Pathway_family1 <- NULL
strep_expression_all$Pathway_family2 <- NULL
strep_expression_all$Gene_code <- NULL
strep_expression_all$Gene_name <- NULL

gf_expression_all$Normalized_abundance <- NULL
gf_expression_all$KO <- NULL
gf_expression_all$Pathway_code <- NULL
gf_expression_all$Pathway_name <- NULL
gf_expression_all$Pathway_family1 <- NULL
gf_expression_all$Pathway_family2 <- NULL
gf_expression_all$Gene_code <- NULL
gf_expression_all$Gene_name <- NULL

