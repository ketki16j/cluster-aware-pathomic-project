# Cluster-Aware Pathomic Pipeline — Full Reproducible Directory

## Overview
End-to-end pipeline for cluster-aware procurement biopsy ML analysis.
Predicts DGF (C-statistic) and 12-month eGFR (MSE) from n=139 transplanted kidneys
at Coimbra University Hospital (2011-2023).
Seed: 382025 | Train: n=112 | Test (held-out): n=27

## Pipeline Steps

### Step 1: Clustering (step1_clustering/)
Build global K-means and hierarchical cluster-averaged tissue data files.
- Input:  Raw CAMELOMICS features + clinical data + train/test subject IDs
- Code:
  * run_hierarchical_clustering.py   → hierarchical within-patient clustering
  * step6B_prepare_camelomics_data.R → prepare CAMELOMICS feature matrix
  * build_single_tissue_files.R      → build per-tissue and multi-tissue CSVs
  * build_2tissue_files.R            → build 2-tissue combination files
  * build_all_FF_inputs.R            → master script for all FF input files
  * build_hier_input.R               → hierarchical clustering input builder
- Output: Renal_Data_hier_*.csv, Renal_Data_naive_*.csv,
          Renal_Data_cluster_averaged_*.csv (47 files total)

### Step 2: Fuzzy Forests Feature Selection (step2_fuzzy_forests/)
Run Fuzzy Forests on clustered data to select top 100 features per config.
- Input:  Clustered tissue data files + FF input files (60 files)
- Code:
  * run_FF_DGF_{method}_{combo}.R  → FF for DGF (30 configs)
  * run_FF_eGFR_{method}_{combo}.R → FF for eGFR (30 configs)
  * build_FF_input_final.R         → build final ComPRePS input files
  * build_all_FF_inputs.R          → master FF input builder
  * build_hier_input.R             → hierarchical FF input builder
- Output: FF_{DGF|eGFR}_{method}_{combo}_top100.csv (69 files)

### Step 3: ComPRePS Bootstrap Training (step3_compreps/)
100-bootstrap training of RF, Lasso, Ridge, Elastic Net per config.
- Input:  FF-selected train/test CSVs (90 config subfolders)
- Code:
  * run_ComPRePS_{method}_{combo}_FF{outcome}.R → FF pipeline (60 configs)
  * run_naive_{outcome}_{combo}.R               → naive pipeline (30 configs)
  * Helper_CV_Functions.R                       → core CV functions
  * Binary/Continuous_Internal_Validation_Metrics.R → metric computation
  * Create_Train_Test_Exclusion_Data.R          → data splitting
  * Model_Saving_and_Test_Performance_Metrics_with_Exclusion.R
- Output: Bootstrap performance CSVs per config (92 folders)
- Key results: Hier T DGF C=0.938 (train), Naive GTAV DGF C=0.949 (train)

### Step 4: Test Set Evaluation (step4_test_evaluation/)
Apply trained models to locked held-out test set (n=27).
- Input:  Train/test CSVs (90 config subfolders)
- Code:
  * run_all_test_evaluations.R          → final reproducible evaluation script
  * original_run_test_evaluation.R      → Jeremy original (9 FF configs)
  * original_run_test_evaluation_allT.R → Jeremy original (4 naive configs)
- Output: test_evaluation_summary.csv + 90 per-config prediction files
- Key results:
    Hier T DGF C=0.792 vs KDPI C=0.410 (test set)
    Naive GTAV DGF C=0.847 (test set)
- Seed: 382025 (single global seed, original configs run first)

### Step 5: Figures (step5_figures/)
Publication-ready figures for manuscript.
- Code:
  * fig4a.R                        → eGFR MSE elbow plots (Fig 4A)
  * fig5a.R / make_fig5a.R         → DGF C-stat elbow plots (Fig 5A)
  * FF_impt.R                      → FF importance bar plots (Fig 5B/6B)
  * roc_curvesonesinglepanelnew.R  → 6-config ROC curves
  * make_roc_hierT.R               → Hier T single ROC
- Output: 8 final publication figures (PNG, 200 DPI)

## To reproduce test evaluation from scratch
  cd step4_test_evaluation/
  module load R/4.5.0
  Rscript code/run_all_test_evaluations.R

## Validated benchmarks
  KDPI:            DGF C=0.410  eGFR MSE=328.5
  Naive GTAV test: DGF C=0.847  eGFR RF MSE=324.0
  Hier T test:     DGF C=0.792  eGFR RF MSE=347.9 (N=12)
