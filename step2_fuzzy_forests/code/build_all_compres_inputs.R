.libPaths('/home/kxj190026/R_libs')

renal <- read.csv('data/Renal_Data.csv')
comp_idx <- which(names(renal) == 'COMPUTATIONAL_MORPHOMETRIC')
clinical_part <- renal[, 1:comp_idx]

# All combinations to build
combos <- c('G','T','A','V','GT','GA','GV','TA','TV','AV')
approaches <- c('global','hier')
outcomes <- c('eGFR','DGF')

built <- 0
for(combo in combos){
  for(approach in approaches){
    # Load cluster file
    if(approach == 'global'){
      cluster_file <- paste0('data/Renal_Data_cluster_averaged_',combo,'.csv')
    } else {
      cluster_file <- paste0('data/Renal_Data_hier_',combo,'.csv')
    }
    df <- read.csv(cluster_file)
    
    for(outcome in outcomes){
      ff_file <- paste0('results/fuzzyforest/FF_',outcome,'_',approach,'_',combo,'_top100.csv')
      out_file <- paste0('data/Renal_Data_FF_',outcome,'_',approach,'_',combo,'_final.csv')
      
      if(!file.exists(ff_file)) next
      
      ff <- read.csv(ff_file)
      top <- ff$feature_name[ff$feature_name %in% names(df) &
                             !ff$feature_name %in% names(clinical_part)]
      
      pathomic_df <- df[, c('Slide_number', top)]
      final <- merge(clinical_part, pathomic_df, by='Slide_number', all.x=TRUE)
      final[is.na(final)] <- 0
      
      write.csv(final, out_file, row.names=FALSE)
      built <- built + 1
    }
  }
}
cat('Built', built, 'ComPRePS input files\n')
