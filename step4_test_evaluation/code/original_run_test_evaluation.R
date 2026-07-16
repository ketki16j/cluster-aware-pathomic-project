.libPaths('/home/kxj190026/R_libs')
library(randomForest)
library(glmnet)
library(pROC)

set.seed(382025)
WORKDIR <- '/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml'

configs <- list(
  list(outcome='DGF',  method='hier',   combo='T',    n=12),
  list(outcome='DGF',  method='hier',   combo='TV',   n=8),
  list(outcome='DGF',  method='hier',   combo='GT',   n=13),
  list(outcome='DGF',  method='global', combo='V',    n=22),
  list(outcome='DGF',  method='global', combo='GTAV', n=32),
  list(outcome='eGFR', method='hier',   combo='T',    n=12),
  list(outcome='eGFR', method='hier',   combo='GT',   n=29),
  list(outcome='eGFR', method='global', combo='V',    n=11),
  list(outcome='eGFR', method='global', combo='GTAV', n=15)
)

results <- data.frame(
  Outcome=character(), Method=character(), Combo=character(),
  Optimal_N=integer(), Train_Perf=numeric(), Test_RF=numeric(),
  Test_Lasso=numeric(), Test_Ridge=numeric(), Test_Enet=numeric(),
  Test_KDPI=numeric(), stringsAsFactors=FALSE
)

for(cfg in configs) {
  outcome <- cfg$outcome
  method  <- cfg$method
  combo   <- cfg$combo
  n       <- cfg$n

  cat('\n========================================\n')
  cat(sprintf('Testing: %s %s %s N=%d\n', outcome, method, combo, n))

  folder <- file.path(WORKDIR, sprintf('results/FF_compres/%s_FF_%s_%s_C',
                                        outcome, method, combo))

  if(outcome == 'DGF') {
    train_file  <- file.path(folder, 'DGF_features_outcome_KDPI_train.csv')
    test_file   <- file.path(folder, 'DGF_features_outcome_KDPI_test.csv')
    outcome_col <- 'cur.outcome'
    kdpi_col    <- 'KDPI.test'
  } else {
    train_file  <- file.path(folder, 'eGFR_features_outcome_KDPI_train.csv')
    test_file   <- file.path(folder, 'eGFR_features_outcome_KDPI_test.csv')
    outcome_col <- 'cur.outcome'
    kdpi_col    <- 'KDPI.test'
  }

  train_df <- read.csv(train_file, row.names=1, check.names=FALSE)
  test_df  <- read.csv(test_file,  row.names=1, check.names=FALSE)

  cat(sprintf('  Train: %d x %d | Test: %d x %d\n',
              nrow(train_df), ncol(train_df),
              nrow(test_df),  ncol(test_df)))

  # Extract outcomes and KDPI
  train_y    <- train_df[[outcome_col]]
  test_y     <- test_df[[outcome_col]]
  train_KDPI <- train_df[['KDPI']]
  test_KDPI  <- test_df[['KDPI.test']]

  # Feature columns — exclude non-feature cols
  exclude <- c(outcome_col, kdpi_col, 'cur.outcome', 'KDPI', 'KDPI.test',
               'Slide_number', 'DONOR_VARIABLES', 'X')
  feat_cols <- setdiff(names(train_df), exclude)
  feat_cols <- feat_cols[1:min(n, length(feat_cols))]

  cat(sprintf('  Using %d features\n', length(feat_cols)))

  X_train <- as.matrix(train_df[, feat_cols])
  X_test  <- as.matrix(test_df[,  feat_cols])

  # Get best nodesize for RF
  nodesize <- 5
  if(outcome == 'DGF') {
    ns_file <- file.path(folder, sprintf('nodesize_best_rf_DGF_%d.csv', n))
    if(file.exists(ns_file)) {
      ns_df <- read.csv(ns_file)
      nodesize <- max(1, as.integer(ns_df[1,1]))
      cat(sprintf('  Nodesize from file: %d\n', nodesize))
    }
  }

  if(outcome == 'DGF') {
    # ── Classification ──────────────────────────────────
    train_y_fac <- factor(train_y)

    lasso_mod <- cv.glmnet(X_train, train_y, family='binomial',
                           alpha=1, type.measure='auc', nfolds=5)
    ridge_mod <- cv.glmnet(X_train, train_y, family='binomial',
                           alpha=0, type.measure='auc', nfolds=5)
    enet_mod  <- cv.glmnet(X_train, train_y, family='binomial',
                           alpha=0.5, type.measure='auc', nfolds=5)
    train_df2      <- as.data.frame(X_train)
    train_df2$y    <- train_y_fac
    rf_mod <- randomForest(y ~ ., data=train_df2, ntree=500,
                           nodesize=nodesize)

    # Test predictions
    test_df2    <- as.data.frame(X_test)
    pred_lasso  <- as.numeric(predict(lasso_mod, newx=X_test,
                                      s='lambda.min', type='response'))
    pred_ridge  <- as.numeric(predict(ridge_mod, newx=X_test,
                                      s='lambda.min', type='response'))
    pred_enet   <- as.numeric(predict(enet_mod,  newx=X_test,
                                      s='lambda.min', type='response'))
    pred_rf     <- as.data.frame(predict(rf_mod, newdata=test_df2,
                                         type='prob'))[,2]

    c_lasso <- as.numeric(auc(test_y, pred_lasso))
    c_ridge <- as.numeric(auc(test_y, pred_ridge))
    c_enet  <- as.numeric(auc(test_y, pred_enet))
    c_rf    <- as.numeric(auc(test_y, pred_rf))
    c_kdpi  <- as.numeric(auc(test_y, test_KDPI))

    # Training C for reference
    pred_rf_train <- as.data.frame(predict(rf_mod,
                      newdata=as.data.frame(X_train), type='prob'))[,2]
    train_perf <- as.numeric(auc(train_y, pred_rf_train))

    cat(sprintf('  Train C (RF): %.3f\n', train_perf))
    cat(sprintf('  Test  C: RF=%.3f  Lasso=%.3f  Ridge=%.3f  Enet=%.3f  KDPI=%.3f\n',
                c_rf, c_lasso, c_ridge, c_enet, c_kdpi))

    results <- rbind(results, data.frame(
      Outcome=outcome, Method=method, Combo=combo, Optimal_N=n,
      Train_Perf=round(train_perf,3), Test_RF=round(c_rf,3),
      Test_Lasso=round(c_lasso,3), Test_Ridge=round(c_ridge,3),
      Test_Enet=round(c_enet,3), Test_KDPI=round(c_kdpi,3),
      stringsAsFactors=FALSE))

    # Save predictions
    write.csv(data.frame(
      Slide=rownames(test_df), True_DGF=test_y,
      Pred_RF=round(pred_rf,4), Pred_Lasso=round(pred_lasso,4),
      Pred_Ridge=round(pred_ridge,4), Pred_Enet=round(pred_enet,4),
      KDPI=test_KDPI),
      file.path(folder, sprintf('test_predictions_DGF_N%d.csv', n)),
      row.names=FALSE)

  } else {
    # ── Regression ──────────────────────────────────────
    lasso_mod <- cv.glmnet(X_train, train_y, family='gaussian',
                           alpha=1, type.measure='mse', nfolds=5)
    ridge_mod <- cv.glmnet(X_train, train_y, family='gaussian',
                           alpha=0, type.measure='mse', nfolds=5)
    enet_mod  <- cv.glmnet(X_train, train_y, family='gaussian',
                           alpha=0.5, type.measure='mse', nfolds=5)
    train_df2   <- as.data.frame(X_train)
    train_df2$y <- train_y
    rf_mod <- randomForest(y ~ ., data=train_df2, ntree=500,
                           nodesize=nodesize)

    test_df2   <- as.data.frame(X_test)
    pred_lasso <- as.numeric(predict(lasso_mod, newx=X_test, s='lambda.min'))
    pred_ridge <- as.numeric(predict(ridge_mod, newx=X_test, s='lambda.min'))
    pred_enet  <- as.numeric(predict(enet_mod,  newx=X_test, s='lambda.min'))
    pred_rf    <- as.numeric(predict(rf_mod, newdata=test_df2))

    # KDPI linear model benchmark
    kdpi_lm   <- lm(train_y ~ train_KDPI)
    pred_kdpi <- predict(kdpi_lm, newdata=data.frame(train_KDPI=test_KDPI))

    mse_lasso <- round(mean((test_y - pred_lasso)^2), 1)
    mse_ridge <- round(mean((test_y - pred_ridge)^2), 1)
    mse_enet  <- round(mean((test_y - pred_enet)^2),  1)
    mse_rf    <- round(mean((test_y - pred_rf)^2),    1)
    mse_kdpi  <- round(mean((test_y - pred_kdpi)^2),  1)

    train_pred_rf <- as.numeric(predict(rf_mod,
                                        newdata=as.data.frame(X_train)))
    train_perf <- round(mean((train_y - train_pred_rf)^2), 1)

    cat(sprintf('  Train MSE (RF): %.1f\n', train_perf))
    cat(sprintf('  Test MSE: RF=%.1f  Lasso=%.1f  Ridge=%.1f  Enet=%.1f  KDPI=%.1f\n',
                mse_rf, mse_lasso, mse_ridge, mse_enet, mse_kdpi))

    results <- rbind(results, data.frame(
      Outcome=outcome, Method=method, Combo=combo, Optimal_N=n,
      Train_Perf=train_perf, Test_RF=mse_rf, Test_Lasso=mse_lasso,
      Test_Ridge=mse_ridge, Test_Enet=mse_enet, Test_KDPI=mse_kdpi,
      stringsAsFactors=FALSE))

    write.csv(data.frame(
      Slide=rownames(test_df), True_eGFR=test_y,
      Pred_RF=round(pred_rf,2), Pred_Lasso=round(pred_lasso,2),
      Pred_Ridge=round(pred_ridge,2), Pred_Enet=round(pred_enet,2),
      KDPI_pred=round(pred_kdpi,2), KDPI=test_KDPI),
      file.path(folder, sprintf('test_predictions_eGFR_N%d.csv', n)),
      row.names=FALSE)
  }
}

# Save summary
dir.create(file.path(WORKDIR, 'results'), showWarnings=FALSE)
out_file <- file.path(WORKDIR, 'results/test_evaluation_summary.csv')
write.csv(results, out_file, row.names=FALSE)

cat('\n========================================\n')
cat('TEST EVALUATION COMPLETE\n')
cat('========================================\n')
print(results)

# ── NAIVE TEST EVALUATION ────────────────────────────────────────────────
cat('\n========================================\n')
cat('Running NAIVE test evaluation\n')
cat('========================================\n')

naive_configs <- list(
  list(outcome='DGF',  combo='T',    n=7),
  list(outcome='DGF',  combo='GTAV', n=27),
  list(outcome='eGFR', combo='T',    n=5),
  list(outcome='eGFR', combo='GTAV', n=6)
)

naive_results <- data.frame(
  Outcome=character(), Method=character(), Combo=character(),
  Optimal_N=integer(), Train_Perf=numeric(), Test_RF=numeric(),
  Test_Lasso=numeric(), Test_Ridge=numeric(), Test_Enet=numeric(),
  Test_KDPI=numeric(), stringsAsFactors=FALSE
)

for(cfg in naive_configs) {
  outcome <- cfg$outcome
  combo   <- cfg$combo
  n       <- cfg$n

  cat(sprintf('\nTesting: %s naive %s N=%d\n', outcome, combo, n))

  folder <- file.path(WORKDIR, sprintf('results/naive_compres/%s_naive_%s_C',
                                        outcome, combo))

  if(outcome == 'DGF') {
    train_file <- file.path(folder, 'DGF_features_outcome_KDPI_train.csv')
    test_file  <- file.path(folder, 'DGF_features_outcome_KDPI_test.csv')
  } else {
    # Try eGFR first, fall back to 1_yr_eGFR
    train_file <- file.path(folder, 'eGFR_features_outcome_KDPI_train.csv')
    if(!file.exists(train_file))
      train_file <- file.path(folder, '1_yr_eGFR_features_outcome_KDPI_train.csv')
    test_file <- file.path(folder, 'eGFR_features_outcome_KDPI_test.csv')
  }

  if(!file.exists(train_file)) { cat('  MISSING train\n'); next }
  if(!file.exists(test_file))  { cat('  MISSING test\n');  next }

  train_df <- read.csv(train_file, row.names=1, check.names=FALSE)
  test_df  <- read.csv(test_file,  row.names=1, check.names=FALSE)

  train_y    <- train_df[['cur.outcome']]
  test_y     <- test_df[['cur.outcome']]
  train_KDPI <- train_df[['KDPI']]
  test_KDPI  <- test_df[['KDPI.test']]

  exclude <- c('cur.outcome','KDPI','KDPI.test','Slide_number','DONOR_VARIABLES','X','')
  feat_cols <- setdiff(names(train_df), exclude)
  feat_cols <- feat_cols[1:min(n, length(feat_cols))]

  X_train <- as.matrix(train_df[, feat_cols])
  X_test  <- as.matrix(test_df[,  feat_cols])

  # Nodesize
  nodesize <- 5
  ns_file <- file.path(folder, sprintf('nodesize_best_rf_%s_%d.csv',
                                        if(outcome=='DGF') 'DGF' else 'eGFR', n))
  if(file.exists(ns_file)) {
    ns_df <- read.csv(ns_file)
    nodesize <- max(1, as.integer(ns_df[1,1]))
  }

  if(outcome == 'DGF') {
    lasso_mod <- cv.glmnet(X_train, train_y, family='binomial', alpha=1,
                           type.measure='auc', nfolds=5)
    ridge_mod <- cv.glmnet(X_train, train_y, family='binomial', alpha=0,
                           type.measure='auc', nfolds=5)
    enet_mod  <- cv.glmnet(X_train, train_y, family='binomial', alpha=0.5,
                           type.measure='auc', nfolds=5)
    train_df2 <- as.data.frame(X_train); train_df2$y <- factor(train_y)
    rf_mod <- randomForest(y ~ ., data=train_df2, ntree=500, nodesize=nodesize)

    test_df2   <- as.data.frame(X_test)
    pred_lasso <- as.numeric(predict(lasso_mod, newx=X_test, s='lambda.min', type='response'))
    pred_ridge <- as.numeric(predict(ridge_mod, newx=X_test, s='lambda.min', type='response'))
    pred_enet  <- as.numeric(predict(enet_mod,  newx=X_test, s='lambda.min', type='response'))
    pred_rf    <- as.data.frame(predict(rf_mod, newdata=test_df2, type='prob'))[,2]

    c_lasso <- round(as.numeric(auc(test_y, pred_lasso)), 3)
    c_ridge <- round(as.numeric(auc(test_y, pred_ridge)), 3)
    c_enet  <- round(as.numeric(auc(test_y, pred_enet)),  3)
    c_rf    <- round(as.numeric(auc(test_y, pred_rf)),    3)
    c_kdpi  <- round(as.numeric(auc(test_y, test_KDPI)), 3)
    train_perf <- round(as.numeric(auc(train_y,
      as.data.frame(predict(rf_mod, newdata=as.data.frame(X_train),
                            type='prob'))[,2])), 3)

    cat(sprintf('  Train C: %.3f | Test RF=%.3f Lasso=%.3f Ridge=%.3f Enet=%.3f KDPI=%.3f\n',
                train_perf, c_rf, c_lasso, c_ridge, c_enet, c_kdpi))

    naive_results <- rbind(naive_results, data.frame(
      Outcome=outcome, Method='naive', Combo=combo, Optimal_N=n,
      Train_Perf=train_perf, Test_RF=c_rf, Test_Lasso=c_lasso,
      Test_Ridge=c_ridge, Test_Enet=c_enet, Test_KDPI=c_kdpi,
      stringsAsFactors=FALSE))

  } else {
    lasso_mod <- cv.glmnet(X_train, train_y, family='gaussian', alpha=1,
                           type.measure='mse', nfolds=5)
    ridge_mod <- cv.glmnet(X_train, train_y, family='gaussian', alpha=0,
                           type.measure='mse', nfolds=5)
    enet_mod  <- cv.glmnet(X_train, train_y, family='gaussian', alpha=0.5,
                           type.measure='mse', nfolds=5)
    train_df2 <- as.data.frame(X_train); train_df2$y <- train_y
    rf_mod <- randomForest(y ~ ., data=train_df2, ntree=500, nodesize=nodesize)

    test_df2   <- as.data.frame(X_test)
    pred_lasso <- as.numeric(predict(lasso_mod, newx=X_test, s='lambda.min'))
    pred_ridge <- as.numeric(predict(ridge_mod, newx=X_test, s='lambda.min'))
    pred_enet  <- as.numeric(predict(enet_mod,  newx=X_test, s='lambda.min'))
    pred_rf    <- as.numeric(predict(rf_mod, newdata=test_df2))
    kdpi_lm    <- lm(train_y ~ train_KDPI)
    pred_kdpi  <- predict(kdpi_lm, newdata=data.frame(train_KDPI=test_KDPI))

    mse_rf    <- round(mean((test_y - pred_rf)^2),    1)
    mse_lasso <- round(mean((test_y - pred_lasso)^2), 1)
    mse_ridge <- round(mean((test_y - pred_ridge)^2), 1)
    mse_enet  <- round(mean((test_y - pred_enet)^2),  1)
    mse_kdpi  <- round(mean((test_y - pred_kdpi)^2),  1)
    train_perf <- round(mean((train_y -
      as.numeric(predict(rf_mod, newdata=as.data.frame(X_train))))^2), 1)

    cat(sprintf('  Train MSE: %.1f | Test RF=%.1f Lasso=%.1f Ridge=%.1f Enet=%.1f KDPI=%.1f\n',
                train_perf, mse_rf, mse_lasso, mse_ridge, mse_enet, mse_kdpi))

    naive_results <- rbind(naive_results, data.frame(
      Outcome=outcome, Method='naive', Combo=combo, Optimal_N=n,
      Train_Perf=train_perf, Test_RF=mse_rf, Test_Lasso=mse_lasso,
      Test_Ridge=mse_ridge, Test_Enet=mse_enet, Test_KDPI=mse_kdpi,
      stringsAsFactors=FALSE))
  }
}

# Combine and save all results
all_results <- rbind(results, naive_results)
write.csv(all_results,
          file.path(WORKDIR, 'results/test_evaluation_summary.csv'),
          row.names=FALSE)

cat('\n========================================\n')
cat('NAIVE EVALUATION COMPLETE\n')
cat('========================================\n')
print(naive_results)
