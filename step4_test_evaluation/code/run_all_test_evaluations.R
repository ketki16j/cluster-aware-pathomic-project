.libPaths('/home/kxj190026/R_libs')
library(randomForest)
library(glmnet)
library(pROC)

# ── DESIGN ─────────────────────────────────────────────────────────────────
# Single set.seed(382025) at top — NO per-config reset.
# Original 9 FF configs run FIRST in original order (matches run_test_evaluation.R)
# Original 4 naive configs run SECOND in original order (matches run_test_evaluation_allT.R)
# Remaining new configs run after — N selected by best RF bootstrap from _C folders.
# This guarantees the 9+4 original configs reproduce exact original results.
# ──────────────────────────────────────────────────────────────────────────

WORKDIR <- '/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml'
set.seed(382025)

find_optimal_N <- function(folder, file_metric, eval_metric=file_metric) {
  csvs <- Sys.glob(file.path(folder,
            paste0("*_MRMR_features_training_", file_metric, "_results.csv")))
  if (length(csvs) == 0) return(NA)
  best_n <- NA; best_val <- if (eval_metric=="C") -Inf else Inf
  for (f in csvs) {
    tryCatch({
      df <- read.csv(f, header=FALSE, stringsAsFactors=FALSE)
      rd <- setNames(trimws(df[,2]), trimws(df[,1]))
      n  <- as.integer(as.numeric(rd["Number of MRMR features"]))
      rf <- suppressWarnings(as.numeric(rd["Random forest"]))
      if (is.na(n) || is.na(rf) || rf <= 0) next
      if (eval_metric=="C"   && rf > best_val) { best_val <- rf; best_n <- n }
      if (eval_metric=="MSE" && rf < best_val) { best_val <- rf; best_n <- n }
    }, error=function(e){})
  }
  best_n
}

run_config <- function(outcome, folder, n, nodesize=5) {
  if (outcome == "DGF") {
    train_file <- file.path(folder, 'DGF_features_outcome_KDPI_train.csv')
    test_file  <- file.path(folder, 'DGF_features_outcome_KDPI_test.csv')
  } else {
    train_file <- file.path(folder, 'eGFR_features_outcome_KDPI_train.csv')
    if (!file.exists(train_file))
      train_file <- file.path(folder, '1_yr_eGFR_features_outcome_KDPI_train.csv')
    test_file  <- file.path(folder, 'eGFR_features_outcome_KDPI_test.csv')
  }
  if (!file.exists(train_file)) { cat('  SKIP: missing train CSV\n'); return(NULL) }
  if (!file.exists(test_file))  { cat('  SKIP: missing test CSV\n');  return(NULL) }

  train_df <- read.csv(train_file, row.names=1, check.names=FALSE)
  test_df  <- read.csv(test_file,  row.names=1, check.names=FALSE)

  train_y    <- train_df[['cur.outcome']]
  test_y     <- test_df[['cur.outcome']]
  train_KDPI <- train_df[['KDPI']]
  test_KDPI  <- test_df[['KDPI.test']]

  exclude   <- c('cur.outcome','KDPI','KDPI.test','Slide_number',
                 'DONOR_VARIABLES','X','')
  feat_cols <- setdiff(names(train_df), exclude)
  feat_cols <- feat_cols[seq_len(min(n, length(feat_cols)))]
  cat(sprintf('  Train: %d x %d | Test: %d x %d | N=%d | nodesize=%d\n',
              nrow(train_df), ncol(train_df), nrow(test_df), ncol(test_df),
              length(feat_cols), nodesize))

  X_train <- as.matrix(train_df[, feat_cols])
  X_test  <- as.matrix(test_df[,  feat_cols])

  tryCatch({
    if (outcome == "DGF") {
      lasso_mod <- cv.glmnet(X_train, train_y, family='binomial',
                             alpha=1,   type.measure='auc', nfolds=5)
      ridge_mod <- cv.glmnet(X_train, train_y, family='binomial',
                             alpha=0,   type.measure='auc', nfolds=5)
      enet_mod  <- cv.glmnet(X_train, train_y, family='binomial',
                             alpha=0.5, type.measure='auc', nfolds=5)
      train_df2 <- as.data.frame(X_train); train_df2$y <- factor(train_y)
      rf_mod    <- randomForest(y ~ ., data=train_df2, ntree=500, nodesize=nodesize)
      test_df2  <- as.data.frame(X_test)
      pred_lasso <- as.numeric(predict(lasso_mod, newx=X_test, s='lambda.min', type='response'))
      pred_ridge <- as.numeric(predict(ridge_mod, newx=X_test, s='lambda.min', type='response'))
      pred_enet  <- as.numeric(predict(enet_mod,  newx=X_test, s='lambda.min', type='response'))
      pred_rf    <- as.data.frame(predict(rf_mod, newdata=test_df2, type='prob'))[,2]
      c_rf    <- round(as.numeric(auc(test_y, pred_rf)),    3)
      c_lasso <- round(as.numeric(auc(test_y, pred_lasso)), 3)
      c_ridge <- round(as.numeric(auc(test_y, pred_ridge)), 3)
      c_enet  <- round(as.numeric(auc(test_y, pred_enet)),  3)
      c_kdpi  <- round(as.numeric(auc(test_y, test_KDPI)),  3)
      cat(sprintf('  Test C: RF=%.3f Lasso=%.3f Ridge=%.3f Enet=%.3f KDPI=%.3f\n',
                  c_rf, c_lasso, c_ridge, c_enet, c_kdpi))
      pred_label <- 'DGF'
      write.csv(data.frame(
        Slide=rownames(test_df), True_DGF=test_y,
        Pred_RF=round(pred_rf,4), Pred_Lasso=round(pred_lasso,4),
        Pred_Ridge=round(pred_ridge,4), Pred_Enet=round(pred_enet,4),
        KDPI=test_KDPI),
        file.path(folder, sprintf('test_predictions_%s_N%d.csv', pred_label, n)),
        row.names=FALSE)
      return(list(rf=c_rf, lasso=c_lasso, ridge=c_ridge, enet=c_enet, kdpi=c_kdpi))
    } else {
      lasso_mod <- cv.glmnet(X_train, train_y, family='gaussian',
                             alpha=1,   type.measure='mse', nfolds=5)
      ridge_mod <- cv.glmnet(X_train, train_y, family='gaussian',
                             alpha=0,   type.measure='mse', nfolds=5)
      enet_mod  <- cv.glmnet(X_train, train_y, family='gaussian',
                             alpha=0.5, type.measure='mse', nfolds=5)
      train_df2 <- as.data.frame(X_train); train_df2$y <- train_y
      rf_mod    <- randomForest(y ~ ., data=train_df2, ntree=500, nodesize=nodesize)
      test_df2  <- as.data.frame(X_test)
      pred_lasso <- as.numeric(predict(lasso_mod, newx=X_test, s='lambda.min'))
      pred_ridge <- as.numeric(predict(ridge_mod, newx=X_test, s='lambda.min'))
      pred_enet  <- as.numeric(predict(enet_mod,  newx=X_test, s='lambda.min'))
      pred_rf    <- as.numeric(predict(rf_mod, newdata=test_df2))
      kdpi_lm   <- lm(train_y ~ train_KDPI)
      pred_kdpi <- predict(kdpi_lm, newdata=data.frame(train_KDPI=test_KDPI))
      mse_rf    <- round(mean((test_y - pred_rf)^2),    1)
      mse_lasso <- round(mean((test_y - pred_lasso)^2), 1)
      mse_ridge <- round(mean((test_y - pred_ridge)^2), 1)
      mse_enet  <- round(mean((test_y - pred_enet)^2),  1)
      mse_kdpi  <- round(mean((test_y - pred_kdpi)^2),  1)
      cat(sprintf('  Test MSE: RF=%.1f Lasso=%.1f Ridge=%.1f Enet=%.1f KDPI=%.1f\n',
                  mse_rf, mse_lasso, mse_ridge, mse_enet, mse_kdpi))
      write.csv(data.frame(
        Slide=rownames(test_df), True_eGFR=test_y,
        Pred_RF=round(pred_rf,2), Pred_Lasso=round(pred_lasso,2),
        Pred_Ridge=round(pred_ridge,2), Pred_Enet=round(pred_enet,2),
        KDPI_pred=round(pred_kdpi,2), KDPI=test_KDPI),
        file.path(folder, sprintf('test_predictions_eGFR_N%d.csv', n)),
        row.names=FALSE)
      return(list(rf=mse_rf, lasso=mse_lasso, ridge=mse_ridge,
                  enet=mse_enet, kdpi=mse_kdpi))
    }
  }, error=function(e) { cat(sprintf('  ERROR: %s\n', e$message)); NULL })
}

get_folder <- function(outcome, method, combo) {
  if (method == "naive") {
    file.path(WORKDIR, sprintf('results/naive_compres/%s_naive_%s_C', outcome, combo))
  } else {
    file.path(WORKDIR, sprintf('results/FF_compres/%s_FF_%s_%s_C', outcome, method, combo))
  }
}

get_nodesize <- function(folder, outcome, n) {
  ns_label <- if (outcome=="DGF") "DGF" else "eGFR"
  ns_file  <- file.path(folder, sprintf('nodesize_best_rf_%s_%d.csv', ns_label, n))
  if (file.exists(ns_file)) max(1L, as.integer(read.csv(ns_file)[1,1])) else 5L
}

results <- data.frame(
  Outcome=character(), Method=character(), Combo=character(),
  Optimal_N=integer(), Test_RF=numeric(), Test_Lasso=numeric(),
  Test_Ridge=numeric(), Test_Enet=numeric(), Test_KDPI=numeric(),
  stringsAsFactors=FALSE
)

add_result <- function(outcome, method, combo, n, r) {
  if (is.null(r)) return()
  results <<- rbind(results, data.frame(
    Outcome=outcome, Method=method, Combo=combo, Optimal_N=n,
    Test_RF=r$rf, Test_Lasso=r$lasso, Test_Ridge=r$ridge,
    Test_Enet=r$enet, Test_KDPI=r$kdpi, stringsAsFactors=FALSE))
}

# ── PHASE 1: Original FF configs in original order (run_test_evaluation.R) ─
cat('\n======== PHASE 1: Original FF configs ========\n')
ff_orig <- list(
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
for (cfg in ff_orig) {
  cat(sprintf('\n=== %s | %s | %s (orig N=%d) ===\n',
              cfg$outcome, cfg$method, cfg$combo, cfg$n))
  folder <- get_folder(cfg$outcome, cfg$method, cfg$combo)
  ns     <- get_nodesize(folder, cfg$outcome, cfg$n)
  r      <- run_config(cfg$outcome, folder, cfg$n, ns)
  add_result(cfg$outcome, cfg$method, cfg$combo, cfg$n, r)
}

# ── PHASE 2: Original naive configs in original order (run_test_evaluation_allT.R) ─
cat('\n======== PHASE 2: Original naive configs ========\n')
naive_orig <- list(
  list(outcome='DGF',  combo='T',    n=7),
  list(outcome='DGF',  combo='GTAV', n=27),
  list(outcome='eGFR', combo='T',    n=5),
  list(outcome='eGFR', combo='GTAV', n=6)
)
for (cfg in naive_orig) {
  cat(sprintf('\n=== %s | naive | %s (orig N=%d) ===\n',
              cfg$outcome, cfg$combo, cfg$n))
  folder <- get_folder(cfg$outcome, 'naive', cfg$combo)
  ns     <- get_nodesize(folder, cfg$outcome, cfg$n)
  r      <- run_config(cfg$outcome, folder, cfg$n, ns)
  add_result(cfg$outcome, 'naive', cfg$combo, cfg$n, r)
}

# ── PHASE 3: All remaining configs ─────────────────────────────────────────
cat('\n======== PHASE 3: Remaining configs ========\n')

# Track which configs are already done
done <- paste(results$Outcome, results$Method, results$Combo)

COMBOS  <- c("G","T","A","V","GT","GA","GV","TA","TV","AV",
             "GTA","GTV","GAV","TAV","GTAV")
METHODS <- c("global","hier","naive")
OUTCOMES <- c("DGF","eGFR")

for (method in METHODS) {
  for (combo in COMBOS) {
    for (outcome in OUTCOMES) {
      key <- paste(outcome, method, combo)
      if (key %in% done) {
        cat(sprintf('\n=== %s | %s | %s — SKIP (already done)\n',
                    outcome, method, combo)); next
      }

      metric      <- if (outcome=="DGF") "C" else "MSE"
      file_metric <- if (outcome=="DGF") "C" else "MSE"
      folder      <- get_folder(outcome, method, combo)

      cat(sprintf('\n=== %s | %s | %s ===\n', outcome, method, combo))
      if (!dir.exists(folder)) { cat('  SKIP: folder not found\n'); next }

      n <- find_optimal_N(folder, file_metric, metric)
      if (is.na(n)) { cat('  SKIP: no bootstrap CSVs\n'); next }
      cat(sprintf('  Optimal N=%d\n', n))

      ns <- get_nodesize(folder, outcome, n)
      r  <- run_config(outcome, folder, n, ns)
      add_result(outcome, method, combo, n, r)
    }
  }
}

# ── Save & print ───────────────────────────────────────────────────────────
write.csv(results,
  file.path(WORKDIR, 'results/test_evaluation_summary.csv'),
  row.names=FALSE)

cat('\n\n========================================\n')
cat('ALL TEST EVALUATIONS COMPLETE\n')
cat('KDPI benchmark:      DGF C=0.410  eGFR MSE=328.5\n')
cat('Naive GTAV baseline: DGF C=0.847  eGFR RF=324.7\n')
cat('Hier T baseline:     DGF C=0.792  eGFR RF=337.9\n')
cat('========================================\n\n')
print(results)
cat('\nSaved: results/test_evaluation_summary.csv\n')
cat('Paste the printed table back to Claude to fill slide 31.\n')

