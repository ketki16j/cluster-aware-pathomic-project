.libPaths('/home/kxj190026/R_libs')
renal    <- read.csv('data/Renal_Data.csv')
comp_idx <- which(names(renal) == 'COMPUTATIONAL_MORPHOMETRIC')
glom_idx <- which(names(renal) == 'GLOMERULAR_GRANULAR_DATA_REVIEWED_MF01')
clinical_part <- renal[, 1:(glom_idx-1)]
cat('Clinical columns:', ncol(clinical_part), '\n')

build_file <- function(ff_csv, cluster_csv, outfile) {
  if(!file.exists(ff_csv))      { cat('MISSING FF:', ff_csv, '\n');      return() }
  if(!file.exists(cluster_csv)) { cat('MISSING cluster:', cluster_csv, '\n'); return() }
  df  <- read.csv(cluster_csv)
  top <- read.csv(ff_csv)$feature_name
  top_pathomic <- top[top %in% names(df) & !top %in% names(clinical_part)]
  pathomic_df  <- df[, c('Slide_number', top_pathomic), drop=FALSE]
  final <- merge(clinical_part, pathomic_df, by='Slide_number', all.x=TRUE)
  final[is.na(final)] <- 0
  write.csv(final, outfile, row.names=FALSE)
  cat('OK:', outfile, '| rows:', nrow(final), '| cols:', ncol(final), '\n')
}

combos  <- c('G','T','A','V','GT','GA','GV','TA','TV','AV',
             'GTA','GTV','GAV','TAV','GTAV')
methods <- c('global','hier')
outcomes <- c('DGF','eGFR')

for(combo in combos) {
  for(method in methods) {
    clust <- if(method=='global') paste0('data/Renal_Data_cluster_averaged_',combo,'.csv')              else paste0('data/Renal_Data_hier_',combo,'.csv')
    for(out in outcomes) {
      ff_csv <- paste0('results/fuzzyforest/FF_',out,'_',method,'_',combo,'_top100.csv')
      outfile <- paste0('data/Renal_Data_FF_',method,'_',combo,'_FF',out,'.csv')
      build_file(ff_csv, clust, outfile)
    }
  }
}
