.libPaths('/home/kxj190026/R_libs')

renal <- read.csv('data/Renal_Data.csv')
df_all <- read.csv('data/Renal_Data_cluster_averaged_all.csv')
df_hier <- read.csv('data/Renal_Data_hierarchical_k4.csv')
color_pat <- paste0('_feat', 28:45, '$', collapse='|')
df_all  <- df_all[,  !grepl(color_pat, names(df_all))]
df_hier <- df_hier[, !grepl(color_pat, names(df_hier))]

# 2-tissue combinations
combos <- list(
  GT = c('G_', 'T_'),
  GA = c('G_', 'A_'),
  GV = c('G_', 'V_'),
  TA = c('T_', 'A_'),
  TV = c('T_', 'V_'),
  AV = c('A_', 'V_')
)

hier_map <- c(G='glomeruli', T='tubules', A='arteries', V='scl_gloms')

for(combo_name in names(combos)){
  prefixes <- combos[[combo_name]]
  
  # Global
  cols <- c('Slide_number', names(df_all)[grepl(paste(paste0('^',prefixes), collapse='|'), names(df_all))])
  df_t <- df_all[, cols]
  write.csv(df_t, paste0('data/Renal_Data_cluster_averaged_',combo_name,'.csv'), row.names=FALSE)
  pathomic <- df_t[,-1]
  zero_pct <- round(sum(pathomic==0)/(nrow(pathomic)*ncol(pathomic))*100,1)
  cat('Global', combo_name, ':', ncol(pathomic), 'features, zeros:', zero_pct, '%\n')
  
  # Hierarchical
  tissues <- strsplit(combo_name, '')[[1]]
  hier_prefixes <- hier_map[tissues]
  cols_h <- c('Slide_number', names(df_hier)[grepl(paste(paste0('^',hier_prefixes), collapse='|'), names(df_hier))])
  df_h <- df_hier[, cols_h]
  write.csv(df_h, paste0('data/Renal_Data_hier_',combo_name,'.csv'), row.names=FALSE)
  pathomic_h <- df_h[,-1]
  zero_pct_h <- round(sum(pathomic_h==0)/(nrow(pathomic_h)*ncol(pathomic_h))*100,1)
  cat('Hier', combo_name, ':', ncol(pathomic_h), 'features, zeros:', zero_pct_h, '%\n')
}
