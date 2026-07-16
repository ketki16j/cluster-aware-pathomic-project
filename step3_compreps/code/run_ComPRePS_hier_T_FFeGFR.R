.libPaths("/home/kxj190026/R_libs")
# Author: Jeremy Rubin
# Date: 03/16/26
# Code to get internal validation metrics for models predicting 1-year eGFR 

set.seed(382025)

### Script to help with machine learning model cross-validation
source(file.path("/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml", "Helper_CV_Functions.R"))
source(file.path("/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml", "Create_Train_Test_Exclusion_Data.R"))
source(file.path("/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml", "CKD_Helper_Functions.R"))

library(glmnet)
library(randomForest)
library(dplyr)
library(mRMRe)
library(stringr)
library(caret)
library(pROC)
library(tidyr)
library(ggplot2)
library(grid)
library(foreach)
library(doParallel)

# Number of bootstraps for internal validation procedure
nboot <- 100

### Outcome is a character string "B" for Delayed Graft Function or 
### "C" for 1-year eGFR outcome
outcome <- "C"

## Boolean variable to determine if you want to save the MRMR-selected features 
## and outcome vector from this iteration
## good for when you have done full sweep over all numbers of MRMR-selected features
## and ready to go forward and pick one trained model with the optimal number of 
## MRMR-selected features
save_optimal_training_features <- T

### Whether or not you want to make a train/test split 
make.split <- F

## Name of file containing features and outcome
input.file.name <- file.path("/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml", "data/Renal_Data_FF_hier_T_FFeGFR.csv")

## Name of files containing training subject and testing subject indices 
## If the training/testing split has been made already
train.sub.name <- file.path("/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml", "pathomic_graft_train_subjects.csv")
test.sub.name <- file.path("/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml", "pathomic_graft_test_subjects.csv")

# Set up parallel backend once for all iterations
# Use detectCores() - 1 to leave one core free for other tasks
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")); if(is.na(n_cores) || n_cores < 1) n_cores <- 4
cl <- makeCluster(n_cores, type="FORK")
registerDoParallel(cl)

WORKDIR_ABS <- "/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml"
folder_c   <- file.path(WORKDIR_ABS, "results/FF_compres/eGFR_FF_hier_T_C")
folder_mse <- file.path(WORKDIR_ABS, "results/FF_compres/eGFR_FF_hier_T_MSE")
dir.create(folder_c,   recursive=TRUE, showWarnings=FALSE)
dir.create(folder_mse, recursive=TRUE, showWarnings=FALSE)
setwd(folder_c)


# Suppress glmnet connection warnings for cleaner output
options(warn = -1)

# Loop through num.MRMR.features from 1 to 100
for(num.MRMR.features in 1:100) {
  
  cat("\n========================================\n")
  cat("Processing num.MRMR.features =", num.MRMR.features, "\n")
  cat("========================================\n")
  
  ## First entry is the number of top MRMR features you want
  # Second entry is the name of the outcome variable, 
  # Third entry is the number of training subjects and 
  ## Last five columns are performance metrics for each of the five models (four machine learning and one for KDPI)
  results.vector <- rep(0,8)
  results.vector[1] <- num.MRMR.features
  
  ## First entry is the number of top MRMR features you want
  # Second entry is the name of the outcome variable, 
  # Third entry is the number of training subjects and 
  ## Last five columns are performance metrics for each of the five models (four machine learning and one for KDPI)
  results.vector[2] <- "1_yr_eGFR"
  
  names(results.vector) <- c("Number of MRMR features","Outcome","Number of training subjects",
                             "Lasso","Ridge","Elastic net","Random forest","KDPI")
  
  training_data_list <- generate_training_testing_exclusion_data(input.file.name,
                                                                 num.MRMR.features, outcome,
                                                                 train.sub.name,
                                                                 test.sub.name,
                                                                 make.split)
  
  final_features <- training_data_list[[1]]
  cur.outcome <- training_data_list[[2]]
  KDPI <- training_data_list[[3]]
  
  if(save_optimal_training_features == T)
  {
    write.csv(cbind(final_features, KDPI, cur.outcome), file=paste0(results.vector[2], "_features_outcome_KDPI_train.csv"))
  }
  
  ## Number of training subjects
  results.vector[3] <- length(cur.outcome)
  
  ## Grid of node sizes to search over for the random forest
  nodesize.try <- seq(from=1, to=floor(sqrt(length(cur.outcome))), by=floor(log(length(cur.outcome))))
  
  ###### Model training
  ### Code adapted to handle when you only have a single MRMR-selected feature vs. multiple MRMR-selected features  
  ### One vs multiple MRMR-selected features
  if(is.vector(final_features) && !is.matrix(final_features) && !is.data.frame(final_features))
  {
    training_data_matrix <- data.matrix(cbind(0, final_features))
  } else{
    training_data_matrix <- data.matrix(final_features)
  }
  
  ## Lasso model fit
  lasso.mod <- cv.glmnet(x=training_data_matrix,
                         y=cur.outcome,
                         type.measure = "mse",
                         alpha=1,
                         nfolds = 5)  
  
  ## Ridge regression model fit 
  ridge.mod <- cv.glmnet(x=training_data_matrix,
                         y=cur.outcome,
                         type.measure = "mse",
                         alpha=0,
                         nfolds = 5)
  
  ## Elastic net model fit 
  enet.mod <- cv.glmnet(x=training_data_matrix,
                        y=cur.outcome,
                        type.measure = "mse",
                        alpha=0.5,
                        nfolds = 5)
  
  ## Use cross-validation to pick nodesize that leads to smallest average MSE for the random forest
  set.seed(382025)
  rf.MSEs <- sapply(nodesize.try, FUN=CV.avg.all.folds,
                    X.train=final_features,
                    Y.train=cur.outcome,
                    algorithm="rf",
                    outcome="continuous")
  
  nodesize.best <- nodesize.try[which(rf.MSEs == min(rf.MSEs))[1]]  
  
  #### Saving optimal nodesize for predicting 1-year eGFR with random forest
  write.csv(nodesize.best, file=paste0("nodesize_best_rf_eGFR_", num.MRMR.features, ".csv"))
  
  ############################# Get the baseline MSEs for each model trained and tested on the original data
  # Case of one MRMR-selected feature
  if(is.vector(final_features) && !is.matrix(final_features) && !is.data.frame(final_features))
  {
    training_data_matrix <- as.matrix(cbind(0, final_features))
    train.df.baseline <- cbind(0, final_features)
    colnames(train.df.baseline) <- c("zero", "X.train")
    
  } else{
    training_data_matrix <- as.matrix(final_features)
    train.df.baseline <- final_features
  }
  
  ## Lasso MSE
  baseline_mse_lasso <- 
    unlist(assess.glmnet(lasso.mod,
                         newx=training_data_matrix,
                         newy=as.vector(cur.outcome), s="lambda.min"))[[1]][1]
  
  ## Ridge regression MSE
  baseline_mse_ridge <- 
    unlist(assess.glmnet(ridge.mod,
                         newx=training_data_matrix,
                         newy=as.vector(cur.outcome), s="lambda.min"))[[1]][1]
  
  ## Elastic net MSE
  baseline_mse_enet <- 
    unlist(assess.glmnet(enet.mod,
                         newx=training_data_matrix,
                         newy=as.vector(cur.outcome), s="lambda.min"))[[1]][1]
  
  ## Random forest MSE 
  rf.mod.baseline <- randomForest(
    x = train.df.baseline,  
    y = cur.outcome,                                              
    ntree = 500,
    nodesize=nodesize.best,
    keep.forest=TRUE,
    keep.inbag=TRUE
  )  
  
  predictions.rf.baseline <- as.vector(predict(rf.mod.baseline, newdata=train.df.baseline))
  baseline_mse_rf <- mean((cur.outcome - predictions.rf.baseline)^2)
  
  ## Baseline MSE from KDPI as the only predictor in linear regression
  train.df.KDPI.baseline <- data.frame(X_train = KDPI, Y_train = cur.outcome)
  KDPI.mod.baseline <- glm(Y_train ~ X_train, data=train.df.KDPI.baseline, family=gaussian(link="identity"))
  predictions.KDPI.baseline <- predict(KDPI.mod.baseline, newdata=train.df.KDPI.baseline, type="response")
  baseline_mse_KDPI <- mean((cur.outcome - predictions.KDPI.baseline)^2)
  
  ## Concatenate all baseline MSEs in one vector
  baseline_MSEs <- c(baseline_mse_lasso, baseline_mse_ridge, baseline_mse_enet, baseline_mse_rf, baseline_mse_KDPI)
  
  #### Saving predicted 1-year eGFR values from optimal random forest model and KDPI model for training cohort
  write.csv(cbind(predictions.rf.baseline, predictions.KDPI.baseline), file=paste0("training_eGFR_rf_and_KDPI_", num.MRMR.features, ".csv"))
  
  # Export necessary objects to parallel workers
  clusterExport(cl, c("final_features", "cur.outcome", "KDPI", "nodesize.try"))
  
  # Parallelized bootstrap loop using foreach
  bootstrap_results <- foreach(p = 1:nboot,
                               .packages = c("glmnet", "randomForest"),
                               .combine = 'rbind') %dopar% {
                                 
                                 ### Bootstrap sampling
                                 bootRows <- sort(sample(1:length(cur.outcome), size=length(cur.outcome), replace=T))
                                 Y.boot <- cur.outcome[bootRows]
                                 KDPI.boot <- KDPI[bootRows]
                                 
                                 ### One vs multiple MRMR-selected features 
                                 if(is.vector(final_features) && !is.matrix(final_features) && !is.data.frame(final_features))
                                 {
                                   X.train.boot <- final_features[bootRows]
                                   
                                   training_matrix_boot <- data.matrix(cbind(0, X.train.boot))
                                   boot.matrix <- as.matrix(cbind(0, X.train.boot))
                                   boot.matrix.original <- as.matrix(cbind(0, final_features))
                                   
                                   train.df.boot <- cbind(0, X.train.boot)
                                   test.df.mse <- cbind(0, final_features)
                                   
                                   colnames(train.df.boot) <- c("zero", "X.train.boot")
                                   colnames(test.df.mse) <- c("zero", "X.train.boot")
                                   
                                 } else{
                                   
                                   X.train.boot <- final_features[bootRows,]
                                   
                                   training_matrix_boot <- data.matrix(X.train.boot)
                                   boot.matrix <- as.matrix(X.train.boot)
                                   boot.matrix.original <- as.matrix(final_features)
                                   
                                   train.df.boot <- X.train.boot
                                   test.df.mse <- final_features
                                   
                                   colnames(train.df.boot) <- colnames(X.train.boot)    
                                   colnames(test.df.mse) <- colnames(final_features)
                                   
                                 }
                                 
                                 ### Bootstrap subject-trained lasso model
                                 lasso.mod.boot <- cv.glmnet(x=training_matrix_boot,
                                                             y=Y.boot,
                                                             type.measure = "mse",
                                                             alpha=1,
                                                             nfolds = 5)  
                                 
                                 ### Bootstrap subject-trained ridge regression model
                                 ridge.mod.boot <- cv.glmnet(x=training_matrix_boot,
                                                             y=Y.boot,
                                                             type.measure = "mse",
                                                             alpha=0,
                                                             nfolds = 5)
                                 
                                 ### Bootstrap subject-trained elastic net model 
                                 enet.mod.boot <- cv.glmnet(x=training_matrix_boot,
                                                            y=Y.boot,
                                                            type.measure = "mse",
                                                            alpha=0.5,
                                                            nfolds = 5)
                                 
                                 ## Use cross-validation to pick nodesize that leads to smallest average MSE
                                 ## for random forest using bootstrapped data
                                 rf.MSEs.boot <- sapply(nodesize.try, FUN=CV.avg.all.folds,
                                                        X.train=X.train.boot,
                                                        Y.train=Y.boot,
                                                        algorithm="rf",
                                                        outcome="continuous")
                                 
                                 nodesize.best.boot <- nodesize.try[which(rf.MSEs.boot == min(rf.MSEs.boot))[1]]
                                 
                                 ############################ Bootstrap trained models tested on bootstrapped data
                                 #### lasso
                                 boot_mse_lasso <- 
                                   unlist(assess.glmnet(lasso.mod.boot,
                                                        newx=boot.matrix,
                                                        newy=as.vector(Y.boot), s="lambda.min"))[[1]][1]
                                 
                                 #### ridge
                                 boot_mse_ridge <- 
                                   unlist(assess.glmnet(ridge.mod.boot,
                                                        newx=boot.matrix,
                                                        newy=as.vector(Y.boot), s="lambda.min"))[[1]][1]
                                 
                                 #### elastic net
                                 boot_mse_enet <- 
                                   unlist(assess.glmnet(enet.mod.boot,
                                                        newx=boot.matrix,
                                                        newy=as.vector(Y.boot), s="lambda.min"))[[1]][1]
                                 
                                 ### random forest
                                 rf.mod.boot <- randomForest(
                                   x = train.df.boot,  
                                   y = Y.boot,                                              
                                   ntree = 500,
                                   nodesize=nodesize.best.boot,
                                   keep.forest=TRUE,
                                   keep.inbag=TRUE
                                 )  
                                 
                                 predictions.rf.boot <- as.vector(predict(rf.mod.boot, newdata=train.df.boot))
                                 boot_mse_rf <- mean((Y.boot - predictions.rf.boot)^2)
                                 
                                 ### bootstrapped KDPI model tested on bootstrapped data
                                 train.df.KDPI.boot <- data.frame(X_train = KDPI.boot, Y.boot = Y.boot)
                                 KDPI.mod.boot <- glm(Y.boot ~ X_train, data=train.df.KDPI.boot, family=gaussian(link="identity"))
                                 predictions.KDPI.boot <- predict(KDPI.mod.boot, newdata=train.df.KDPI.boot, type="response")
                                 boot_mse_KDPI <- mean((Y.boot - predictions.KDPI.boot)^2)
                                 
                                 ####### bootstrapped trained models tested on original data
                                 ### lasso
                                 test_mse_lasso <- 
                                   unlist(assess.glmnet(lasso.mod.boot,
                                                        newx=boot.matrix.original,
                                                        newy=as.vector(cur.outcome), s="lambda.min"))[[1]][1]
                                 
                                 ### ridge regression
                                 test_mse_ridge <- 
                                   unlist(assess.glmnet(ridge.mod.boot,
                                                        newx=boot.matrix.original,
                                                        newy=as.vector(cur.outcome), s="lambda.min"))[[1]][1]
                                 
                                 ### elastic net
                                 test_mse_enet <- 
                                   unlist(assess.glmnet(enet.mod.boot,
                                                        newx=boot.matrix.original,
                                                        newy=as.vector(cur.outcome), s="lambda.min"))[[1]][1]
                                 
                                 ### random forest
                                 predictions.rf.test <- as.vector(predict(rf.mod.boot, newdata=test.df.mse))
                                 test_mse_rf <- mean((cur.outcome - predictions.rf.test)^2)
                                 
                                 ########## bootstrapped KDPI model on original data
                                 test.df.KDPI <- data.frame(X_train = KDPI)
                                 predictions.KDPI.test <- predict(KDPI.mod.boot, newdata=test.df.KDPI, type="response")
                                 test_mse_KDPI <- mean((cur.outcome - predictions.KDPI.test)^2)
                                 
                                 ### CKD stage frequencies for true bootstrapped eGFRs, RF predictions, and KDPI predictions
                                 ckd_true <- get_ckd_frequencies(Y.boot)
                                 ckd_rf   <- get_ckd_frequencies(predictions.rf.boot)
                                 ckd_kdpi <- get_ckd_frequencies(predictions.KDPI.boot)
                                 
                                 # Return boot MSEs (cols 1-5), test MSEs (cols 6-10), and CKD frequencies (cols 11-28)
                                 c(boot_mse_lasso, boot_mse_ridge, boot_mse_enet, boot_mse_rf, boot_mse_KDPI,
                                   test_mse_lasso, test_mse_ridge, test_mse_enet, test_mse_rf, test_mse_KDPI,
                                   ckd_true, ckd_rf, ckd_kdpi)
                               }
  
  # Extract results from the combined matrix
  boot_mse_all_models <- bootstrap_results[, 1:5]
  test_mse_all_models <- bootstrap_results[, 6:10]
  
  ### Reconstruct CKD stage frequencies across bootstrap iterations and compute median
  CKD_stages_boot <- array(NA, dim = c(3, 6, nboot))
  CKD_stages_boot[1,,] <- t(bootstrap_results[, 11:16])
  CKD_stages_boot[2,,] <- t(bootstrap_results[, 17:22])
  CKD_stages_boot[3,,] <- t(bootstrap_results[, 23:28])
  
  ### Want to get the median bootstrap-averaged CKD stages for the true eGFRs, 
  ### optimal random forest-predicted eGFRs and the KDPI-predicted eGFRs
  CKD_stages_median <- apply(CKD_stages_boot, MARGIN = c(1, 2), FUN = median)
  write.csv(CKD_stages_median, file=paste0("CKD_staging_median_bootstrap_", num.MRMR.features, ".csv"))
  
  ### Applying the internal validation formula to compute final MSEs for all methods 
  results.vector[4:8] <- baseline_MSEs - apply(FUN=median, MARGIN=2, boot_mse_all_models - test_mse_all_models)
  
  write.csv(results.vector, file=paste0(num.MRMR.features, "_MRMR_features_training_MSE_results.csv"))
  
  cat("Completed num.MRMR.features =", num.MRMR.features, "\n")
  cat("Time:", format(Sys.time(), "%H:%M:%S"), "\n")
}

# Stop the cluster after all iterations are complete
stopCluster(cl)

# Turn warnings back on
options(warn = 0)

cat("\n========================================\n")
cat("ALL ANALYSES COMPLETE!\n")
cat("Processed num.MRMR.features from 1 to 100\n")
cat("Final time:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("========================================\n")