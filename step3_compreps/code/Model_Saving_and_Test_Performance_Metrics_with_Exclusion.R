# Author: Jeremy Rubin
# Date: 1/30/26
# Compute optimal random forest model and KDPI model performance in test set 
# Visualize feature importance of optimal random forest model 
# Visualize CKD distributions of true eGFR, predicted eGFR by random forest 
# and predicted eGFR by KDPI models within training set and in test set
# Visualize ROC curve on test data for random forest model and KDPI
# and calculate Youden's index. 
# Additionally, get predicted 1-year eGFRs and probabilities of DGF for 
# Excluded kidneys using optimal random forest model

set.seed(382025)

source("Helper_CV_Functions.R")
source("CKD_Helper_Functions.R")
library(randomForest)
library(permimp)
library(ggplot2)
library(pROC)
library(dplyr)
library(tidyr)
library(gridExtra)
library(stringr)

#### Character string to indicate which outcome you're working with 
#### "B" for delayed graft function and "C" for eGFR 
#### representing binary vs. continuous outcomes
outcome <- "C"

### Reading in files for eGFR prediction
if(outcome=="C")
{
  # Predicted eGFRs by optimal chosen random forest model and KDPI linear regression model 
  # on training data
  training_eGFR_pred_rf_KDPI_file <- "training_eGFR_rf_and_KDPI.csv"
  
  ## File to read in distributions of CKD stages for the true eGFR, optimal random forest model
  ## and KDPI averaged across bootstraps in the training data
  CKD_stages_training_file <- "CKD_staging_median_bootstrap.csv"  
  
  # String to add to files being written out
  write_file_tag <- "eGFR"
  
  ### Reading in files for DGF prediction
} else {
  
  # File with predicted DGF probabilities from optimal random forest model 
  training_DGF_pred_file <- "training_DGF_prob.csv"
  
  # String to add to DGF files being written out
  write_file_tag <- "DGF"
}

## Read in saved features from chosen optimal number of MRMR-selected features
## as well as the vector of KDPI values, 1-year eGFR or DGF outcome vector and 
## optimal random forest node size from the training data
training_features_outcomes_file <- paste0(write_file_tag,"_features_outcome_KDPI_train.csv") 
testing_features_outcome_file <- paste0(write_file_tag,"_features_outcome_KDPI_test.csv")
rf_nodesize_file <- paste0("nodesize_best_rf_",write_file_tag,".csv")

## Read in recipient codes for train/test data for predicting 1-year eGFR or DGF outcome
rec_code_train_file <- paste0(write_file_tag,"_rec_train_code.csv")
rec_code_test_file <- paste0(write_file_tag,"_rec_test_code.csv")

# File to read in MRMR-selected features for optimal random forest model 
# with slide number 
excluded_with_slide_file <- paste0(write_file_tag,"_features_excluded_with_slide_num.csv")

## Read in chosen features for the training set, 
# the outcome vector on the training set, and 
## KDPI on the training set
features_outcome_KDPI <- read.csv(training_features_outcomes_file, row.names=1)
nodesize.best <- read.csv(rf_nodesize_file, row.names=1)

# Converting from a data frame to a scalar
nodesize.best <- as.numeric(nodesize.best[1,1])

cur.outcome <- as.numeric(unlist(subset(features_outcome_KDPI,select = cur.outcome)))
KDPI <- as.numeric(unlist(subset(features_outcome_KDPI,select = KDPI)))
final_features <- features_outcome_KDPI[,-c((ncol(features_outcome_KDPI)-1):ncol(features_outcome_KDPI))]

## Read in chosen features for the test set, the outcome vector
# on the test set, and 
## KDPI on the test set
features_outcome_KDPI_test <- read.csv(testing_features_outcome_file, row.names=1)
cur.outcome.test <- as.numeric(unlist(subset(features_outcome_KDPI_test,select = cur.outcome)))
KDPI.test <- as.numeric(unlist(subset(features_outcome_KDPI_test,select = KDPI.test)))

# Read in features and slide number for excluded kidneys
## NOTE: Feature exclusion based on missingness was performed using training subjects
## so subjects with excluded kidneys are not guaranteed to have all observed feature
## values of interest. Must take out excluded subjects who are missing features on which
## optimal random forest model was trained on
features_slidenum_excluded <- read.csv(excluded_with_slide_file, row.names=1)
colnames(features_slidenum_excluded)[1] <- "Slide_num"
features_slidenum_excluded <- na.omit(features_slidenum_excluded)
Slide_num <- features_slidenum_excluded$Slide_num
features_excluded_only <- features_slidenum_excluded[,-1]

# Making outcome a factor if we're looking at DGF
if(outcome=="B"){
  cur.outcome <- as.factor(cur.outcome)
  cur.outcome.test <- as.factor(cur.outcome.test)
}

# Make sure you just keep the features for the test subjects which correspond to those which
# the model were trained on
test_features <- features_outcome_KDPI_test[, colnames(final_features)]

### Fit random forest model with optimal/chosen number of MRMR-selected features and 
### optimal nodesize
rf.mod <- randomForest(
      x = final_features,  
      y = cur.outcome,                                              
      ntree = 500,
      nodesize=nodesize.best,
      keep.forest=TRUE,
      keep.inbag=TRUE
    )

## Save optimal random forest model
saveRDS(rf.mod, file = paste0("random_forest_model_optimal_",write_file_tag,".rds"))

## Random forest predictions on test set    
if(outcome=="C")
{
  ## RF predictions and MSE on test transplant subjects
  predictions.rf.test <- as.vector(predict(rf.mod,newdata=test_features))
  mse_rf_test <- mean((cur.outcome.test - predictions.rf.test)^2)
  
  ## RF predictions on excluded kidneys 
  predictions.rf.excluded <- as.vector(predict(rf.mod,newdata=features_excluded_only))
  write.csv(cbind(Slide_num,predictions.rf.excluded),file=paste0("rf_predictions_excluded_",write_file_tag,".csv"))
  
  # baseline KDPI mse
  train.df.KDPI.baseline <- data.frame(X_train = KDPI, Y_train = cur.outcome)
  
  # Fit KDPI linear regression model
  KDPI.mod <- glm(Y_train ~ X_train, data=train.df.KDPI.baseline, family=gaussian(link="identity"))
  
  ## Save KDPI linear regression model
  saveRDS(KDPI.mod, file = paste0("KDPI_model_",write_file_tag,".rds"))
  
  # Make predictions on test set with KDPI linear regression model
  test_KDPI_df <- data.frame(X_train = KDPI.test)
  predictions.KDPI.test <- predict(KDPI.mod, newdata=test_KDPI_df, type="response")
  mse_KDPI_test <- mean((cur.outcome.test - predictions.KDPI.test)^2)
  
  ## Save MSEs of random forest and KDPI models
  mse_test_save <- c(mse_rf_test,mse_KDPI_test)
  names(mse_test_save) <- c("Random forest test MSE","KDPI test MSE")
  write.csv(mse_test_save,file="test_MSEs_rf_KDPI.csv")
  
  ### Report out predicted 1-year eGFRs for transplanted subjects in training and test sets 
  ### with recipient codes for optimal random forest model
  eGFR_rec_test_code <- read.csv(rec_code_test_file, row.names=1)
  training_eGFR_code <- read.csv(rec_code_train_file, row.names=1)
  training_eGFR_pred <- read.csv(training_eGFR_pred_rf_KDPI_file, row.names=1)
  
  testing_eGFR_and_rec_code <- cbind(eGFR_rec_test_code,predictions.rf.test)
  colnames(testing_eGFR_and_rec_code) <- c("Recipient code","Predicted eGFR on test data from optimal random forest model")
  
  training_eGFR_and_rec_code <- cbind(training_eGFR_code,training_eGFR_pred[,1])
  colnames(training_eGFR_and_rec_code) <- c("Recipient code","Predicted eGFR on training data from optimal random forest model")
  
  write.csv(training_eGFR_and_rec_code,file="Predicted_eGFR_rf_train_rec_code.csv")
  write.csv(testing_eGFR_and_rec_code,file="Predicted_eGFR_rf_test_rec_code.csv")
  
  #### Visualize CKD Stages from KDPI linear regression, optimal random forest model 
  #### and true eGFR on the training data 
  CKD_internal_validation_stages <- read.csv(CKD_stages_training_file,row.names = 1)
  plot_CKD_distributions(CKD_internal_validation_stages,"train")
  
  ## Now do the same CKD stage visualization but based on true/predicted CKD stages in the test set
  plot_CKD_distributions(rbind(get_ckd_frequencies(cur.outcome.test),
                               get_ckd_frequencies(predictions.rf.test),
                               get_ckd_frequencies(predictions.KDPI.test)),"test")
  
  ### Scatterplot of predicted eGFR by KDPI vs. true eGFR and predicted eGFR by RF vs. true eGFR 
  ### for test data
  
  # KDPI limits
  lims_kdpi <- range(c(cur.outcome.test, predictions.KDPI.test), na.rm = TRUE)
  
  # RF limits
  lims_rf <- range(c(cur.outcome.test, predictions.rf.test), na.rm = TRUE)
  
  p_kdpi <- ggplot(
    data.frame(
      true = cur.outcome.test,
      pred = predictions.KDPI.test
    ),
    aes(x = true, y = pred)
  ) +
    geom_point(alpha = 0.5, color = "#D55E00") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal(xlim = lims_kdpi, ylim = lims_kdpi) +
    labs(
      title = "KDPI Predicted vs True eGFR",
      x = "True eGFR",
      y = "Predicted eGFR (KDPI)"
    ) +
    theme_bw()
  
  p_rf <- ggplot(
    data.frame(
      true = cur.outcome.test,
      pred = predictions.rf.test
    ),
    aes(x = true, y = pred)
  ) +
    geom_point(alpha = 0.5, color = "#0072B2") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal(xlim = lims_rf, ylim = lims_rf) +
    labs(
      title = "Random Forest Predicted vs True eGFR",
      x = "True eGFR",
      y = "Predicted eGFR (RF)"
    ) +
    theme_bw()
  
  grid.arrange(p_kdpi, p_rf, ncol = 2)
  
  ### DGF outcome 
} else{
  
  ## Compute AUCs on test data for optimal random forest model and using KDPI alone (no model)
  predictions.rf.test <- as.data.frame(unlist(predict(rf.mod,newdata=test_features,type="prob")))[,2]
  auc_rf_test <- auc(cur.outcome.test,predictions.rf.test)
  auc_KDPI_test <- auc(cur.outcome.test,KDPI.test)
  
  ## Save AUCs of random forest model and using KDPI score
  auc_test_save <- c(auc_rf_test,auc_KDPI_test)
  names(auc_test_save) <- c("Random forest test AUC","KDPI test AUC")
  write.csv(auc_test_save,file="test_AUCs_rf_KDPI.csv")
  
  ## RF predictions on excluded kidneys 
  predictions.rf.excluded <- as.data.frame(unlist(predict(rf.mod,newdata=features_excluded_only,type="prob")))[,2]
  write.csv(cbind(Slide_num,predictions.rf.excluded),file=paste0("rf_predictions_excluded_",write_file_tag,".csv"))
  
  ### Code to plot ROC curve for test data
  # Create ROC objects
  roc_rf <- roc(cur.outcome.test, predictions.rf.test)
  roc_kdpi <- roc(cur.outcome.test, KDPI.test)
  
  # Set up the plot area
  plot(
    x = 1 - roc_rf$specificities, 
    y = roc_rf$sensitivities,
    type = "l",
    col = "blue",
    lwd = 2,
    xlab = "False Positive Rate (1 - Specificity)",
    ylab = "True Positive Rate (Sensitivity)",
    main = "ROC Curve Comparison for Delayed Graft Function Prediction in Hold-out Transplant Set",
    xlim=c(0,1),
    ylim=c(0,1)
  )
  
  # Add the KDPI ROC curve
  lines(
    x = 1 - roc_kdpi$specificities, 
    y = roc_kdpi$sensitivities,
    col = "brown",
    lwd = 2
  )
  
  # Add the diagonal reference line
  abline(a = 0, b = 1, lty = 2, col = "gray")
  
  # Add the subtitle with AUC values
  mtext(paste("Random Forest AUC =", round(auc(cur.outcome.test,predictions.rf.test),2), 
              ", KDPI AUC =", round(auc(cur.outcome.test,KDPI.test),2)),
        side = 3, line = 0.5, cex = 0.8)
  
  legend(
    "topleft",
    legend = c("Random Forest", "KDPI"),
    col = c("blue", "brown"),
    lwd = 2,
    bty = "n",
    cex = 0.8,
    seg.len = 1,    # Makes the lines shorter (default is usually 2)
    x.intersp = 0.5 # Reduces space between lines and text (default is 1)
  )
  
  #### Picking Youden's J statistic based on ROC curve for test data
  # maximizes sensitivity + specificity - 1
  coords_youden <- coords(roc_rf, "best", best.method = "youden")
  threshold_youden <- coords_youden$threshold
  
  # Saving Youden's index
  # Can be useful for classification out-of-sample by using the Youden's Index
  # as a threshold for classification
  write.csv(threshold_youden,file="DGF_Youden_Index.csv")
  
  ### Report out predicted probabilities of DGF in training and test sets with recipient codes
  DGF_rec_test_code <- read.csv(rec_code_test_file, row.names=1)
  training_DGF_code <- read.csv(rec_code_train_file, row.names=1)
  training_DGF_pred <- read.csv(training_DGF_pred_file, row.names=1)
  
  testing_DGF_and_rec_code <- cbind(DGF_rec_test_code,predictions.rf.test)
  colnames(testing_DGF_and_rec_code) <- c("Recipient code","Predicted probability of DGF from optimal random forest model")
  
  training_DGF_and_rec_code <- cbind(training_DGF_code,training_DGF_pred)
  colnames(training_DGF_and_rec_code) <- c("Recipient code","Predicted probability of DGF from optimal random forest model")
  
  write.csv(training_DGF_and_rec_code,file="Predicted_DGF_rf_train_rec_code.csv")
  write.csv(testing_DGF_and_rec_code,file="Predicted_DGF_rf_test_rec_code.csv")
  
}

### Visualize feature importance of random forest model 
## Select 1 when prompted by permimp method
set.seed(382025)
rf.coeff <- permimp(rf.mod)
rf.values <- rf.coeff$values
  
# for rf.values, negative values indicate predictor is less informative,
# These negative values can still occur even if the features were selected by MRMR 
# due to correlations and interactions between features 
# We set the negative feature importance values to zero 
# so that we can do L1 normalization only for features 
# with positive permutation importances and have consistent directions of association
# of informativeness with outcome 
rf.values[rf.values<0] <- 0
rf.values <- rf.values/sum(abs(rf.values))
  
# Create data frame from your feature importance vector
  rf_importance <- data.frame(
    Feature = names(rf.values),
    Importance = as.numeric(rf.values)
  )
 
# Ensure correct ordering for plotting (biggest on top)
  rf_importance <- rf_importance %>%
      arrange(desc(Importance)) %>%
      slice_head(n = 15) %>%          # keep top 15
      mutate(plot_order = factor(Feature, levels = rev(Feature)))
    
# Create the horizontal bar plot
  ggplot(rf_importance, aes(x = Importance, y = plot_order)) +
      geom_bar(stat = "identity", fill = "#4682B4") +
      geom_text(aes(label = sprintf("%.3f", Importance)), 
                hjust = -0.2, size = 3.5) +
      labs(
        title = paste0("Random Forest Feature Importance for Predicting ", write_file_tag),
        subtitle = "L1-Normalized Feature Importance Values",
        x = "Importance",
        y = ""
      ) +
      theme_minimal() +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size = 9)
      ) +
      scale_x_continuous(limits = c(0, max(rf_importance$Importance) * 1.15))