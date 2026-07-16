.libPaths('/home/kxj190026/R_libs')

renal <- read.csv('data/Renal_Data.csv')
df_all <- read.csv('data/Renal_Data_cluster_averaged_all.csv')
df_hier <- read.csv('data/Renal_Data_hierarchical_k4.csv')
color_pat <- paste0('_feat', 28:45, '$', collapse='|')
df_all  <- df_all[,  !grepl(color_pat, names(df_all))]
df_hier <- df_hier[, !grepl(color_pat, names(df_hier))]

# Global k-means single tissue files
tissues_global <- list(
  G = list(prefix='G_', file='data/Renal_Data_cluster_averaged_G.csv'),
  T = list(prefix='T_', file='data/Renal_Data_cluster_averaged_T.csv'),
  A = list(prefix='A_', file='data/Renal_Data_cluster_averaged_A.csv'),
  V = list(prefix='V_', file='data/Renal_Data_cluster_averaged_V.csv')
)

for(t in names(tissues_global)){
  cols <- c('Slide_number', names(df_all)[grepl(paste0('^',t,'_'), names(df_all))])
  df_t <- df_all[, cols]
  write.csv(df_t, tissues_global[[t]]$file, row.names=FALSE)
  pathomic <- df_t[,-1]
  zero_pct <- round(sum(pathomic==0)/(nrow(pathomic)*ncol(pathomic))*100,1)
  cat('Global', t, ':', ncol(pathomic), 'features, zero imputation:', zero_pct, '%\n')
}

# Hierarchical single tissue files
tissues_hier <- list(
  G = list(prefix='glomeruli', file='data/Renal_Data_hier_G.csv'),
  T = list(prefix='tubules',   file='data/Renal_Data_hier_T.csv'),
  A = list(prefix='arteries',  file='data/Renal_Data_hier_A.csv'),
  V = list(prefix='scl_gloms', file='data/Renal_Data_hier_V.csv')
)

for(t in names(tissues_hier)){
  cols <- c('Slide_number', names(df_hier)[grepl(paste0('^',tissues_hier[[t]]$prefix), names(df_hier))])
  df_t <- df_hier[, cols]
  write.csv(df_t, tissues_hier[[t]]$file, row.names=FALSE)
  pathomic <- df_t[,-1]
  zero_pct <- round(sum(pathomic==0)/(nrow(pathomic)*ncol(pathomic))*100,1)
  cat('Hier', t, ':', ncol(pathomic), 'features, zero imputation:', zero_pct, '%\n')
}
