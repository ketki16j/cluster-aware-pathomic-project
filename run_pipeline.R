# ============================================================
# run_pipeline.R
# Master pipeline script — cluster-aware pathomic ML framework
#
# Runs end-to-end: FF feature selection → ComPRePS training → test evaluation
#
# Usage:
#   Rscript run_pipeline.R [options]
#
# Options:
#   --pipeline  : hier | global | naive | all  (default: all)
#   --outcome   : DGF  | eGFR   | all          (default: all)
#   --tissues   : T | G | A | V | GT | GV | AV | GTAV | all  (default: all)
#   --step      : FF | compreps | test | figures | all  (default: all)
#
# Examples:
#   Rscript run_pipeline.R --pipeline hier --outcome DGF --tissues T
#   Rscript run_pipeline.R --pipeline naive --outcome all --tissues GTAV
#   Rscript run_pipeline.R --pipeline all --outcome all --tissues all --step test
#
# Author: Ketki Joshi
# Date:   June 2026
# ============================================================

.libPaths('/home/kxj190026/R_libs')
library(optparse)

# ── 1. Parse arguments ────────────────────────────────────────────────────
option_list <- list(
  make_option('--pipeline', type='character', default='all',
    help='Pipeline type: hier | global | naive | all [default: all]'),
  make_option('--outcome', type='character', default='all',
    help='Outcome: DGF | eGFR | all [default: all]'),
  make_option('--tissues', type='character', default='all',
    help='Tissue combo: T | G | A | V | GT | GV | AV | GTAV | all [default: all]'),
  make_option('--step', type='character', default='all',
    help='Pipeline step: FF | compreps | test | figures | all [default: all]')
)

opt <- parse_args(OptionParser(option_list=option_list))

cat('\n========================================\n')
cat('Cluster-Aware Pathomic ML Pipeline\n')
cat('========================================\n')
cat(sprintf('Pipeline : %s\n', opt$pipeline))
cat(sprintf('Outcome  : %s\n', opt$outcome))
cat(sprintf('Tissues  : %s\n', opt$tissues))
cat(sprintf('Step     : %s\n', opt$step))
cat('========================================\n\n')

# ── 2. Setup ──────────────────────────────────────────────────────────────
WORKDIR <- '/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml'
set.seed(382025)

# All available options
ALL_PIPELINES <- c('hier', 'global', 'naive')
ALL_OUTCOMES  <- c('DGF', 'eGFR')
ALL_TISSUES   <- c('G','T','A','V',
                   'GT','GA','GV','TA','TV','AV',
                   'GTA','GTV','GAV','TAV','GTAV')

# Resolve "all" selections
pipelines <- if (opt$pipeline == 'all') ALL_PIPELINES else strsplit(opt$pipeline, ',')[[1]]
outcomes  <- if (opt$outcome  == 'all') ALL_OUTCOMES  else strsplit(opt$outcome,  ',')[[1]]
tissues   <- if (opt$tissues  == 'all') ALL_TISSUES   else strsplit(opt$tissues,  ',')[[1]]
steps     <- if (opt$step     == 'all') c('FF','compreps','test','figures')
             else strsplit(opt$step, ',')[[1]]

# Validate inputs
valid_pipelines <- c('hier','global','naive')
valid_outcomes  <- c('DGF','eGFR')
valid_steps     <- c('FF','compreps','test','figures')

for (p in pipelines) if (!p %in% valid_pipelines) stop(sprintf('Unknown pipeline: %s', p))
for (o in outcomes)  if (!o %in% valid_outcomes)  stop(sprintf('Unknown outcome: %s', o))
for (s in steps)     if (!s %in% valid_steps)     stop(sprintf('Unknown step: %s', s))

# Tissue prefix map (for FF feature column detection)
tissue_prefix_map <- list(
  G='glomeruli', T='tubules', A='arteries', V='arterioles',
  GT='(glomeruli|tubules)', GA='(glomeruli|arteries)',
  GV='(glomeruli|arterioles)', TA='(tubules|arteries)',
  TV='(tubules|arterioles)', AV='(arteries|arterioles)',
  GTA='(glomeruli|tubules|arteries)', GTV='(glomeruli|tubules|arterioles)',
  GAV='(glomeruli|arteries|arterioles)', TAV='(tubules|arteries|arterioles)',
  GTAV='(glomeruli|tubules|arteries|arterioles)'
)
global_prefix_map <- list(
  G='^G_cluster', T='^T_cluster', A='^A_cluster', V='^V_cluster',
  GT='^(G|T)_cluster', GA='^(G|A)_cluster', GV='^(G|V)_cluster',
  TA='^(T|A)_cluster', TV='^(T|V)_cluster', AV='^(A|V)_cluster',
  GTA='^(G|T|A)_cluster', GTV='^(G|T|V)_cluster', GAV='^(G|A|V)_cluster',
  TAV='^(T|A|V)_cluster', GTAV='^(G|T|A|V)_cluster'
)

# Outcome column map
outcome_col_map <- list(
  DGF  = 'Delayed_Graft_Function',
  eGFR = 'eGFR_CKD_EPI_12M'
)

# Clinical variables
donor_clinical <- c('Donor_age','Expanded_criteria_donor',
  'Donor_Sex_0_Male_1_Female','Donor_height','Donor_weight_Kg','Donor_BMI',
  'Donor_Hypertension_History','Donor_Diabetes_History',
  'Donor_final_creatinine','Donor_eGFR_CKD_EPI_final_creatinine','REMUZZI_GT1')

# ── 3. Build config list ──────────────────────────────────────────────────
configs <- list()

for (pipeline in pipelines) {
  for (outcome in outcomes) {
    for (tissue in tissues) {

      if (pipeline == 'naive') {
        cfg <- list(
          pipeline    = pipeline,
          outcome     = outcome,
          tissue      = tissue,
          label       = sprintf('%s naive_%s', outcome, tissue),
          data_file   = sprintf('data/Renal_Data_naive_%s.csv', tissue),
          ff_out      = NULL,  # naive has no FF step
          input_file  = sprintf('data/Renal_Data_naive_%s.csv', tissue),
          folder_c    = sprintf('results/naive_compres/%s_naive_%s_C', outcome, tissue),
          folder_mse  = sprintf('results/naive_compres/%s_naive_%s_MSE', outcome, tissue),
          outcome_B_C = if (outcome=='DGF') 'B' else 'C',
          tissue_prefix = NULL
        )
      } else {
        # hier or global
        data_file <- if (pipeline == 'hier')
          sprintf('data/Renal_Data_hier_%s.csv', tissue)
        else
          sprintf('data/Renal_Data_cluster_averaged_%s.csv', tissue)

        ff_suffix <- if (outcome == 'DGF') 'FFDGF' else 'FFeGFR'
        input_file <- sprintf('data/Renal_Data_FF_%s_%s_%s_%s.csv',
                              outcome, pipeline, tissue, ff_suffix)
        tprefix <- if (pipeline == 'hier')
          paste0('^', tissue_prefix_map[[tissue]])
        else
          global_prefix_map[[tissue]]

        cfg <- list(
          pipeline    = pipeline,
          outcome     = outcome,
          tissue      = tissue,
          label       = sprintf('%s %s_%s', outcome, pipeline, tissue),
          data_file   = data_file,
          ff_out      = sprintf('results/fuzzyforest/FF_%s_%s_%s_top100.csv',
                                outcome, pipeline, tissue),
          input_file  = input_file,
          folder_c    = sprintf('results/FF_compres/%s_FF_%s_%s_C', outcome, pipeline, tissue),
          folder_mse  = sprintf('results/FF_compres/%s_FF_%s_%s_MSE', outcome, pipeline, tissue),
          outcome_B_C = if (outcome == 'DGF') 'B' else 'C',
          tissue_prefix = tprefix,
          outcome_col   = outcome_col_map[[outcome]]
        )
      }
      configs[[length(configs)+1]] <- cfg
    }
  }
}

cat(sprintf('Total configs to run: %d\n\n', length(configs)))

# ── 4. Step FF: Fuzzy Forests feature selection ───────────────────────────
if ('FF' %in% steps) {
  cat('======== STEP 2: Fuzzy Forests Feature Selection ========\n\n')

  library(fuzzyforest); library(WGCNA)

  sp <- screen_control(min_ntree=500, keep_fraction=0.25,
                       mtry_factor=1, ntree_factor=1, drop_fraction=0.25)
  lp <- select_control(min_ntree=500, number_selected=100,
                       mtry_factor=1, ntree_factor=1, drop_fraction=0.25)
  wp <- WGCNA_control(power=6, TOMType='unsigned', minModuleSize=10)

  full      <- read.csv(file.path(WORKDIR, 'data/Renal_Data.csv'))
  train_idx <- as.vector(unlist(read.csv(
    file.path(WORKDIR,'train_indices.csv'), row.names=1)))

  dir.create(file.path(WORKDIR,'results/fuzzyforest'),
             showWarnings=FALSE, recursive=TRUE)

  ff_configs <- configs[sapply(configs, function(c) !is.null(c$ff_out))]
  cat(sprintf('FF configs: %d\n\n', length(ff_configs)))

  for (i in seq_along(ff_configs)) {
    cfg <- ff_configs[[i]]
    out_file <- file.path(WORKDIR, cfg$ff_out)

    if (file.exists(out_file)) {
      cat(sprintf('[%2d/%d] SKIP: %s\n', i, length(ff_configs), cfg$label))
      next
    }

    cat(sprintf('[%2d/%d] FF: %s\n', i, length(ff_configs), cfg$label))

    tryCatch({
      df       <- read.csv(file.path(WORKDIR, cfg$data_file))
      full_sub <- full[!is.na(full[[cfg$outcome_col]]), ]
      df       <- df[df$Slide_number %in% full_sub$Slide_number, ]
      full_sub <- full_sub[match(df$Slide_number, full_sub$Slide_number), ]

      train_slides <- df$Slide_number[train_idx]
      df_t   <- df[df$Slide_number %in% train_slides, ]
      full_t <- full_sub[full_sub$Slide_number %in% train_slides, ]

      y             <- as.numeric(full_t[[cfg$outcome_col]])
      clin          <- as.data.frame(lapply(full_t[, donor_clinical], as.numeric))
      pathomic_cols <- names(df_t)[grepl(cfg$tissue_prefix, names(df_t))]

      X <- cbind(clin, df_t[, pathomic_cols])
      X <- as.data.frame(lapply(X, as.numeric))
      X <- X[, colSums(is.na(X)) == 0]
      X <- X[, sapply(X, var, na.rm=TRUE) > 0]

      module_membership <- rep('clinical', ncol(X))
      module_membership[grepl(cfg$tissue_prefix, names(X))] <- 'pathomic'

      set.seed(382025)
      ff_mod <- wff(X, y,
                    module_membership = module_membership,
                    screen_params     = sp,
                    select_params     = lp,
                    WGCNA_params      = wp)

      write.csv(ff_mod$feature_list, out_file, row.names=FALSE)
      cat(sprintf('  Saved: %s\n', basename(out_file)))

    }, error=function(e) cat(sprintf('  ERROR: %s\n', conditionMessage(e))))
  }
  cat('\nFF complete.\n\n')
}

# ── 5. Step ComPRePS: Bootstrap training ──────────────────────────────────
if ('compreps' %in% steps) {
  cat('======== STEP 3: ComPRePS Bootstrap Training ========\n\n')

  source(file.path(WORKDIR, 'Helper_CV_Functions.R'))
  source(file.path(WORKDIR, 'Create_Train_Test_Exclusion_Data.R'))
  library(pROC); library(glmnet); library(mRMRe)
  library(stringr); library(caret); library(foreach); library(doParallel)

  cat(sprintf('ComPRePS configs: %d (%d bootstraps each)\n\n',
              length(configs), 100))

  for (i in seq_along(configs)) {
    cfg <- configs[[i]]

    # Skip if output already exists
    folder_c <- file.path(WORKDIR, cfg$folder_c)
    if (dir.exists(folder_c) &&
        length(list.files(folder_c, pattern='results\\.csv$')) >= 10) {
      cat(sprintf('[%3d/%d] SKIP: %s\n', i, length(configs), cfg$label))
      next
    }

    cat(sprintf('[%3d/%d] ComPRePS: %s\n', i, length(configs), cfg$label))

    tryCatch({
      dir.create(file.path(WORKDIR, cfg$folder_c),
                 showWarnings=FALSE, recursive=TRUE)
      dir.create(file.path(WORKDIR, cfg$folder_mse),
                 showWarnings=FALSE, recursive=TRUE)

      input_file <- file.path(WORKDIR, cfg$input_file)
      if (!file.exists(input_file)) {
        cat(sprintf('  SKIP: input file missing: %s\n', cfg$input_file))
        next
      }

      # Source the appropriate individual script for this config
      # (uses existing per-config scripts as backend)
      script <- if (cfg$pipeline == 'naive')
        file.path(WORKDIR, sprintf('run_naive_%s_%s.R', cfg$outcome, cfg$tissue))
      else
        file.path(WORKDIR, sprintf('run_ComPRePS_%s_%s_FF%s.R',
                                   cfg$pipeline, cfg$tissue, cfg$outcome))

      if (file.exists(script)) {
        source(script)
        cat(sprintf('  Done: %s\n', cfg$label))
      } else {
        cat(sprintf('  SKIP: script not found: %s\n', basename(script)))
      }

    }, error=function(e) cat(sprintf('  ERROR: %s\n', conditionMessage(e))))
  }
  cat('\nComPRePS complete.\n\n')
}

# ── 6. Step Test: Held-out test evaluation ────────────────────────────────
if ('test' %in% steps) {
  cat('======== STEP 4: Test Set Evaluation ========\n\n')
  source(file.path(WORKDIR, 'run_all_test_evaluations.R'))
  cat('Test evaluation complete.\n\n')
}

# ── 7. Step Figures: Generate publication figures ─────────────────────────
if ('figures' %in% steps) {
  cat('======== STEP 5: Generate Figures ========\n\n')
  setwd(file.path(WORKDIR))
  for (fig_script in c('plots/fig4a.R','plots/fig5a.R',
                        'plots/FF_impt.R',
                        'plots/roc_curvesonesinglepanelnew.R')) {
    if (file.exists(fig_script)) {
      cat(sprintf('Running: %s\n', fig_script))
      source(fig_script)
    }
  }
  cat('Figures complete.\n\n')
}

cat('========================================\n')
cat('Pipeline complete!\n')
cat('========================================\n')
