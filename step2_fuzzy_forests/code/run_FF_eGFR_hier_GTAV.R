.libPaths('/home/kxj190026/R_libs')
library(fuzzyforest); library(WGCNA); library(randomForest)
set.seed(382025)
WORKDIR <- '/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml'
df   <- read.csv(file.path(WORKDIR, 'data/Renal_Data_hier_GTAV.csv'))
full <- read.csv(file.path(WORKDIR, 'data/Renal_Data.csv'))
has_outcome <- !is.na(full$eGFR_CKD_EPI_12M)
full <- full[has_outcome, ]
df   <- df[df$Slide_number %in% full$Slide_number, ]
full <- full[match(df$Slide_number, full$Slide_number), ]
train_idx    <- as.vector(unlist(read.csv(file.path(WORKDIR,'train_indices.csv'),row.names=1)))
train_slides <- df$Slide_number[train_idx]
df_t   <- df[df$Slide_number %in% train_slides, ]
full_t <- full[full$Slide_number %in% train_slides, ]
y <- as.numeric(full_t$eGFR_CKD_EPI_12M)
donor_clinical <- c('Donor_age','Expanded_criteria_donor','Donor_Sex_0_Male_1_Female',
  'Donor_height','Donor_weight_Kg','Donor_BMI','Donor_Hypertension_History',
  'Donor_Diabetes_History','Donor_final_creatinine','Donor_eGFR_CKD_EPI_final_creatinine',
  'REMUZZI_GT1')
clin       <- as.data.frame(lapply(full_t[, donor_clinical], as.numeric))
pathomic_cols <- names(df_t)[grepl('^glomeruli|^tubules|^arteries|^scl_gloms', names(df_t))]
cat('Pathomic features:', length(pathomic_cols), '\n')
X <- cbind(clin, df_t[, pathomic_cols])
X <- as.data.frame(lapply(X, as.numeric))
X <- X[, colSums(is.na(X))==0]
X <- X[, sapply(X, var, na.rm=TRUE)>0]
cat('Samples:', nrow(X), '| Total features:', ncol(X), '\n')
sp <- screen_control(min_ntree=500,keep_fraction=0.25,mtry_factor=1,ntree_factor=1,drop_fraction=0.25)
lp <- select_control(min_ntree=500,number_selected=100,mtry_factor=1,ntree_factor=1,drop_fraction=0.25)
wp <- WGCNA_control(power=6,TOMType='unsigned',minModuleSize=10)
cat('Running FF eGFR hier_GTAV...\n')
ff  <- wff(X, y, WGCNA_params=wp, screen_params=sp, select_params=lp, final_ntree=1000)
out <- file.path(WORKDIR,'results/fuzzyforest',paste0('FF_eGFR_hier_GTAV_top100.csv'))
write.csv(ff$feature_list, out, row.names=FALSE)
cat('Saved:', out, '\n')
