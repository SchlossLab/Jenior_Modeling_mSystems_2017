deps <- c('wesanderson','vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

#--------------------------------------------------------------------------------------------------------------#

# Define variables
growth_file <- '/home/mjenior/Desktop/growth_curve/formatted_wells.tsv'

# Read in data
growth <- read.delim(growth_file, sep='\t', header=TRUE, row.names=1)
growth <- as.data.frame(t(growth))
rm(growth_file)

#--------------------------------------------------------------------------------------------------------------#

# Average treatement groups
galactitol <- rowMeans(cbind(growth$E2, growth$F2, growth$G2), na.rm=TRUE)
starch <-  rowMeans(cbind(growth$E3, growth$F3, growth$G3), na.rm=TRUE)
fructose <-  rowMeans(cbind(growth$E4, growth$F4, growth$G4), na.rm=TRUE)
mannitol <- rowMeans(cbind(growth$E6, growth$F6, growth$G6), na.rm=TRUE)
salicin <- rowMeans(cbind(growth$E7, growth$F7, growth$G7), na.rm=TRUE)
sorbitol <- rowMeans(cbind(growth$E8, growth$F8, growth$G8), na.rm=TRUE)

y_glucose_y_aa <- rowMeans(cbind(growth$B2, growth$C2, growth$D2), na.rm=TRUE)
n_glucose_y_aa <- rowMeans(cbind(growth$B3, growth$C3, growth$D3), na.rm=TRUE)
y_glucose_n_aa <- rowMeans(cbind(growth$B4, growth$C4, growth$D4), na.rm=TRUE)
n_glucose_n_aa <- rowMeans(cbind(growth$B5, growth$C5, growth$D5), na.rm=TRUE)

y_glucose_n_MTV <- rowMeans(cbind(growth$B6, growth$C6, growth$D6), na.rm=TRUE) 
n_glucose_n_MTV <- rowMeans(cbind(growth$B7, growth$C7, growth$D7), na.rm=TRUE) 

blank <- rowMeans(cbind(growth$B9, growth$C9, growth$D9, growth$B10, growth$C10, growth$D10, growth$B11, growth$C11, growth$D11), na.rm=TRUE)
rm(growth)

# Subtract blank means
galactitol <- galactitol - blank
starch <- starch - blank
fructose <- fructose - blank
mannitol <- mannitol - blank
salicin <- salicin - blank
sorbitol <- sorbitol - blank
  
y_glucose_y_aa <- y_glucose_y_aa - blank
n_glucose_y_aa <- n_glucose_y_aa - blank
y_glucose_n_aa <- y_glucose_n_aa - blank
n_glucose_n_aa <- n_glucose_n_aa - blank
  
y_glucose_n_MTV <- y_glucose_n_MTV - blank
n_glucose_n_MTV <- n_glucose_n_MTV - blank
rm(blank)

#--------------------------------------------------------------------------------------------------------------#

# Plot growth rates









#--------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
rm(galactitol, starch, fructose, mannitol, salicin, sorbitol, y_glucose_y_aa, n_glucose_y_aa, y_glucose_n_aa, n_glucose_n_aa, y_glucose_n_MTV, n_glucose_n_MTV)

