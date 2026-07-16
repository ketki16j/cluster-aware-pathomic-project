### Author: Jeremy Rubin
### Date: 1/30/26
### Code to make train/test/excluded kidneys data for downstream model 
# fitting and testing 

### FUNCTION INPUTS: 
### Name of input file
### Number of MRMR-selected/top ranking features that you want
### outcome you want to predict: Either eGFR_CKD_EPI_12M ("C") or Delayed_Graft_Function ("B")
### If you are inputting the training and testing indices, you need to provide the file 
### train.sub.name and test.sub.name are the file names of the training/testing split indices
### either to write to or read from 
### make.split is a boolean indicating whether you are making the training/testing split or not
### make.split should only be true for your first time predicting the Delayed Graft Function (DGF) outcome

### Important note: Start analysis by running this function once with (DGF)
### outcome so that you can create the train/test split based on the event rate of DGF
### Then, you can take the generated training/testing split and read it in for predicting 1-year eGFR

### Input data file to read in
## File should contain donor clinical factors, pathomic features and 1-year eGFR
### Assumes the variables are arranged in the following order with the following column separators 
### in the following order: RECIPIENT_DATA, DONOR_VARIABLES, REMUZZI_CALCULATIONS, AI_FEATURES, 
### COMPUTATIONAL_MORPHOMETRIC, GLOMERULAR_GRANULAR_DATA_REVIEWED_MF01, CT_GRANULAR_FEATURES_REV, 
### GRANULAR_FEATURES_VASCULAR_MEAN_VALUES_ALL_SIZES, GRANULAR_FEATURES_VASCULAR_ARTERY_WITH_MAX_AREA_R
### GRANULAR_FEATURES_VASCULAR_ARTERIOLAR_AREA_INF_35000_PX_R, GRANULAR_FEATURES_VASCULAR_ARTERY_AREA_SUP_35000_PX_R

generate_training_testing_exclusion_data = function(input.file.name,num.MRMR.features,outcome,
                                  train.sub.name,test.sub.name,make.split)
{
  Renal_Data <- read.csv(input.file.name)

  ### If one-year eGFR is our outcome #####
  if(outcome=="C")
  {
    ### Assumes that outcome variable (1-year eGFR) is called eGFR_CKD_EPI_12M
    ### Want to keep only subjects for whom we have a 1-year eGFR recorded
    has.outcome <- !is.na(Renal_Data$eGFR_CKD_EPI_12M)
    cur.outcome <- na.omit(as.numeric(unlist(subset(Renal_Data,select=eGFR_CKD_EPI_12M))))
    filename.tag <- "eGFR"
  } 
  
    ## Otherwise, we want to regress on delayed graft function as the outcome
    else{
      has.outcome <- !is.na(Renal_Data$Delayed_Graft_Function)
      cur.outcome <- na.omit(as.numeric(unlist(subset(Renal_Data,select=Delayed_Graft_Function))))
      filename.tag <- "DGF"
    }
  
  ### Reading in training/testing subject indices if you've done the split already
  if(outcome=="B" & make.split==T)
  {
    set.seed(382025)
    ## Doing an 80% training/testing split
    ## Where you're 
    train_ind <- createDataPartition(y = cur.outcome, 
                                     p = 0.8,  # 80% for training
                                     list = FALSE)
    
    # Convert to vector if it's not already
    train_ind <- as.vector(train_ind)
    # Create test indices
    test_ind <- setdiff(1:length(cur.outcome), train_ind)
    
    # Sort indices
    train_ind <- sort(train_ind)
    test_ind <- sort(test_ind)
    
    write.csv(train_ind,file=train.sub.name)
    write.csv(test_ind,file=test.sub.name)
  } 
  
  ### Otherwise, making training/testing split and saving the results 
  else {
    
    ### Pre-specified vectors of indices for transplanted subjects to use in train/validation sets
    train_ind <- as.vector(unlist(read.csv(train.sub.name, row.names=1)))
    test_ind <- as.vector(unlist(read.csv(test.sub.name, row.names=1)))     
  }
  
  ## Get rid of variables that are dates
  date_columns <- grep("date", colnames(Renal_Data))
  Renal_Data <- Renal_Data[,-c(date_columns)]
  
  # Save slide number to save later with excluded kidney features 
  Slide_num <- Renal_Data$Slide_number
  
  # Save KDPI for transplanted subjects 
  # to use for later as a competitor model to machine learning models
  KDPI <- subset(Renal_Data, select = KDPI_2024)
  KDPI <- as.numeric(unlist(KDPI))[has.outcome]
  
  ## Additional variables that shouldn't be regressed on
  Renal_Data <- subset(Renal_Data, select = -c(Slide_number,
                                               Allocated_Discarded,
                                               Number_of_Sections_Per_Slide,
                                               Kidney_used,
                                               Donation_after_Cardiac_Death_donor,
                                               KDRI_2024,
                                               KDPI_2024,
                                               KDPI_SUP_85_2024))
  
  ### Save transplant recipient codes
  Rec_code <- Renal_Data$Recipient_code
  Rec_code <- Renal_Data$Recipient_code[has.outcome]
  
  # Get rid of recipient features that aren't 1-year eGFR
  Renal_Data <- Renal_Data[, !(seq_along(Renal_Data) >= 
                                 which(names(Renal_Data) == "Recipient_code") & 
                                 seq_along(Renal_Data) < which(names(Renal_Data) == "DONOR_VARIABLES"))]
  
  # Donor variables to remove before model fitting because they are either empty columns 
  # or have no variance
  Renal_Data <- subset(Renal_Data, select = -c(Donor_cause_of_death_code,
                                               REMUZZI_CALCULATIONS,
                                               Donor_HCV_Status,
                                               Donor_race,
                                               Glomerular_GT1,
                                               Vascular_GT1))
  
  # Keep relevant donor variables for machine learning model fitting  
  Donor_Features <- Renal_Data[, (seq_along(Renal_Data) > 
                                    which(names(Renal_Data) == "DONOR_VARIABLES") & 
                                    seq_along(Renal_Data) < which(names(Renal_Data) == "AI_FEATURES"))] 
  
  # Making REMUZZI_GT1 a multicategorical variable
  Donor_Features$REMUZZI_GT1_cat <- cut(
    Donor_Features$REMUZZI_GT1, 
    breaks = c(-1, 4, 5, 12),
    labels = c("Mild", "Moderate", "Severe"),
    right = TRUE
  )
  
  ## Collapsing donor cause of death to multicategorical variable with levels stroke, thrauma, 
  ## and other 
  Donor_Features$Donor_cause_of_death <- dplyr::case_when(
    grepl("Stroke", Donor_Features$Donor_cause_of_death, ignore.case = TRUE) ~ 1,
    grepl("Thrauma", Donor_Features$Donor_cause_of_death, ignore.case = TRUE) ~ 2,
    TRUE ~ 3  # All other cases
  )
  
  # Converting multilevel categorical variables to factors   
  categorical.don.var <- subset(Donor_Features,select = c(Donor_cause_of_death,REMUZZI_GT1_cat))
  for(i in 1:ncol(categorical.don.var)){categorical.don.var[,i] <- as.factor(categorical.don.var[,i])}
  
  ## Dummy coding multicategorical variables
  categorical.dummy.variables <- model.matrix(~ 0 + Donor_cause_of_death + REMUZZI_GT1_cat 
                                              , data = categorical.don.var)[,-1]
  
  # Selecting donor continuous variables or binary variables that have been proprely dummy coded
  cont.don.var <- subset(Donor_Features, select = c(Donor_age,Expanded_criteria_donor,
                                                    Donor_Sex_0_Male_1_Female,
                                                    Donor_height,
                                                    Donor_weight_Kg,
                                                    Donor_BMI,
                                                    Donor_Hypertension_History,
                                                    Donor_Diabetes_History,
                                                    Donor_final_creatinine,
                                                    Donor_eGFR_CKD_EPI_final_creatinine))
  
  Donor_modified_var <- cbind(categorical.dummy.variables,cont.don.var)
  
  ## Removing all columns strictly after DONOR_VARIABLES up to and including COMPUTATIONAL_MORPHOMETRIC   
  Renal_Data <- Renal_Data[, !(seq_along(Renal_Data) > 
                                 which(names(Renal_Data) == "DONOR_VARIABLES") & 
                                 seq_along(Renal_Data) <= which(names(Renal_Data) == "COMPUTATIONAL_MORPHOMETRIC"))]
  
  # Get rid of columns that are all NAs
  all.missing = function(x){sum(is.na(x))}
  column.spacers <- (apply(FUN=all.missing,MARGIN=2,Renal_Data) == nrow(Renal_Data))
  Renal_Data <- Renal_Data[,!column.spacers]
  
  ## Combine the corrected donor clinical features with the image features we want
  Renal_Data <- cbind(Donor_modified_var,Renal_Data)
  
  # Number of subjects
  n_subjects <- length(cur.outcome)
  
  # Get the missingness percentages for each variable
  # among subjects who have the outcome of interest
  Renal_Data_with_outcome <- Renal_Data[has.outcome,]
  missingness_rates <- colSums(is.na(Renal_Data_with_outcome)) / n_subjects * 100
  
  # Keep features with no missingness 
  no_missingness_columns <- missingness_rates == 0
  Renal_Data <- Renal_Data[, no_missingness_columns]
  Renal_Data_with_outcome <- Renal_Data[has.outcome,]
  
  # The target variable for MRMR needs to be the outcome column 
  Renal_Data_with_outcome <- cbind(Renal_Data_with_outcome,cur.outcome)
  for(i in 1:ncol(Renal_Data_with_outcome))
    {Renal_Data_with_outcome[,i] <- as.numeric(Renal_Data_with_outcome[,i])}
  
  # Saving recipient codes for train/test sets
  recipient.train.code.filename <- paste0(filename.tag,"_rec_train_code.csv")
  recipient.test.code.filename <- paste0(filename.tag,"_rec_test_code.csv")
    
  write.csv(Rec_code[train_ind],file=recipient.train.code.filename)
  write.csv(Rec_code[test_ind],file=recipient.test.code.filename)  
  
  # Save test features and KDPI values
  Renal_Data_test <- Renal_Data_with_outcome[test_ind,]
  KDPI.test <- as.numeric(unlist(KDPI))
  KDPI.test <- KDPI.test[test_ind]
  
  write.csv(file = paste0(filename.tag,"_features_outcome_KDPI_test.csv"),
              cbind(Renal_Data_test,KDPI.test))  
  
  ### Set up model training to happen only on the training subjects
  Renal_Data_train <- Renal_Data_with_outcome[train_ind,]
  cur.outcome <- as.numeric(unlist(cur.outcome))[train_ind]
  KDPI <- as.numeric(unlist(KDPI))[train_ind]
  
  ### Must ensure all variables are numeric before running MRMR feature selection
  ### Assumes the last column is your target variable (in our case, 1-year eGFR)
  feature_data <- mRMR.data(data = Renal_Data_train)
  
  num.MRMR.features <- as.numeric(num.MRMR.features)
  
  #### Run mRMR
  # Adjust the feature_count parameter to select the desired number of features
  set.seed(382025)
  result <- mRMR.classic(data = feature_data, 
                         target_indices = ncol(Renal_Data_train), 
                         feature_count = num.MRMR.features)  
  
  # Features selected by MRMR by column name
  selected_features <- solutions(result)[[1]]
  
  # Get the actual feature names
  feature_names <- featureNames(feature_data)[selected_features]
  
  ### Check for any multicategorical variables present - if only part of one is selected, then want to make sure
  ### you include all binary indicators that make up that predictor after running MRMR and getting feature_names
  search_strings <- c("REMUZZI_GT1", "Donor_cause_of_death")
  
  # Create a single pattern with all search strings joined by '|' (OR operator)
  pattern <- paste(search_strings, collapse="|")
  
  # Find all matches in one operation
  matches <- grep(pattern, feature_names, value=TRUE)
  
  # Find which search strings match any feature names
  # matched_strings should contain the predictors for the multicategorical variables that you want to keep all levels of
  matched_strings <- character(0)
  for(s in search_strings) {
    if(any(grepl(s, feature_names))) {
      matched_strings <- c(matched_strings, s)
    }
  }
  
  # Get all column names
  all_column_names <- colnames(Renal_Data)
  
  if(length(matched_strings)>0){
    
    matching_columns <- which(str_detect(all_column_names, str_c(matched_strings, collapse = "|")))
    
    ##### KEEP FEATURES SELECTED BY MRMR AND ANY MISSING LEVELS OF SELECTED MULTICATEGORICAL VARIABLES
    ##### Doing so for training data and for kidneys which were excluded
    final_features_train <- Renal_Data_train[,unique(sort(c(matching_columns,selected_features)))]
    final_features_excluded <- Renal_Data[!has.outcome,unique(sort(c(matching_columns,selected_features)))]
    
  } else{
    ##### If no multicategorical predictors were selected, 
    #### then MRMR feature selection gave you all of the features you should need 
    #### Also want to save the corresponding selected features for the kidneys which were excluded
    final_features_train <- Renal_Data_train[,sort(selected_features)]
    final_features_excluded <- Renal_Data[!has.outcome,sort(selected_features)]
  }
  
  ## Save the MRMR-selected features, outcome vector and KDPI values for the training data  
  write.csv(file = paste0(filename.tag,"_features_outcome_KDPI_train.csv"),
              cbind(final_features_train,cur.outcome,KDPI))  
  
  ## Save the MRMR-selected features and slide number for excluded kidneys 
  write.csv(file = paste0(filename.tag,"_features_excluded_with_slide_num.csv"),
            cbind(Slide_num[!has.outcome],final_features_excluded))
  
  ### Returning as a list the final design matrix of features, outcome vector, and 
  ### KDPI for the training subjects
  return(list(final_features_train,cur.outcome,KDPI))
}
