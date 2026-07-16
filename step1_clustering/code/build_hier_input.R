.libPaths('/home/kxj190026/R_libs')

renal <- read.csv('data/Renal_Data.csv')
df <- read.csv('data/Renal_Data_hierarchical_k4.csv')

comp_idx <- which(names(renal) == 'COMPUTATIONAL_MORPHOMETRIC')
clinical_part <- renal[, 1:comp_idx]

build_file <- function(ff_file, outfile) {
  ff <- read.csv(ff_file)
  top <- ff$feature_name[ff$feature_name %in% names(df) &
                         !ff$feature_name %in% names(clinical_part)]
  cat('Pathomic features added:', length(top), '\n')
  pathomic_df <- df[, c('Slide_number', top)]
  final <- merge(clinical_part, pathomic_df, by='Slide_number', all.x=TRUE)
  final[is.na(final)] <- 0
  cat('Shape:', nrow(final), 'x', ncol(final), '\n')
  cat('KDPI_2024:', 'KDPI_2024' %in% names(final), '\n')
  write.csv(final, outfile, row.names=FALSE)
  cat('Saved:', outfile, '\n\n')
}

cat('=== eGFR hier ===\n')
build_file('results/fuzzyforest/FF_eGFR_hier_top100.csv',
           'data/Renal_Data_FF_eGFR_hier_final.csv')

cat('=== DGF hier ===\n')
build_file('results/fuzzyforest/FF_DGF_hier_top100.csv',
           'data/Renal_Data_FF_DGF_hier_final.csv')
