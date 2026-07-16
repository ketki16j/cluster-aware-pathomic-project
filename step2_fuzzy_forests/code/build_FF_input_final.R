.libPaths('/home/kxj190026/R_libs')

renal <- read.csv('data/Renal_Data.csv')
df <- read.csv('data/Renal_Data_cluster_averaged_all.csv')

# Keep all clinical columns (1 to COMPUTATIONAL_MORPHOMETRIC inclusive)
comp_idx <- which(names(renal) == 'COMPUTATIONAL_MORPHOMETRIC')
clinical_part <- renal[, 1:comp_idx]
cat('Clinical columns:', ncol(clinical_part), '\n')
cat('KDPI_2024 in clinical:', 'KDPI_2024' %in% names(clinical_part), '\n')

# Function to build final file
build_file <- function(top_features, outfile) {
  # Remove any features that are already in clinical part
  top_pathomic <- top_features[top_features %in% names(df) & 
                                !top_features %in% names(clinical_part)]
  cat('Pathomic features added:', length(top_pathomic), '\n')
  
  pathomic_df <- df[, c('Slide_number', top_pathomic)]
  
  # Merge - clinical already has Slide_number
  final <- merge(clinical_part, pathomic_df, by='Slide_number', all.x=TRUE)
  final[is.na(final)] <- 0
  
  cat('Final shape:', nrow(final), 'x', ncol(final), '\n')
  cat('KDPI_2024:', 'KDPI_2024' %in% names(final), '\n')
  cat('KDPI_SUP_85_2024:', 'KDPI_SUP_85_2024' %in% names(final), '\n')
  cat('DONOR_VARIABLES:', 'DONOR_VARIABLES' %in% names(final), '\n')
  
  write.csv(final, outfile, row.names=FALSE)
  cat('Saved:', outfile, '\n\n')
}

# eGFR top 100
ff_egfr <- read.csv('results/fuzzyforest/FF_eGFR_top100_high.csv')
cat('=== eGFR features ===\n')
build_file(ff_egfr$feature_name, 'data/Renal_Data_FF_eGFR_final.csv')

# DGF top 100
ff_dgf <- read.csv('results/fuzzyforest/FF_DGF_top100_high.csv')
cat('=== DGF features ===\n')
build_file(ff_dgf$feature_name, 'data/Renal_Data_FF_DGF_final.csv')
