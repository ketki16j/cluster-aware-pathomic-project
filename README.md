# procurement-biopsy-pathomics-ml

Code used for training and testing cluster-aware machine learning models with donor clinical and pathomic biopsy features to predict kidney transplant recipient outcomes.

**Manuscript in preparation:** *"Cluster-Aware Pathomic Modeling of Procurement Biopsies for Prediction of Post-Transplant Allograft Outcomes"*

**In collaboration with:** Jeremy Rubin's Lab, University of Maryland

---

## Overview

This repository implements a cluster-aware pathomic ML framework for analyzing deceased donor procurement biopsies. Rather than averaging morphometric features across all tissue objects (naive approach), we apply unsupervised clustering (global K-means and hierarchical within-patient) per tissue type to capture morphological subpopulation structure. Cluster-specific features are combined with donor clinical variables and passed through Fuzzy Forests feature selection prior to multimodal ML training using the [ComPRePS](https://github.com/jeremysrubin/procurement-biopsy-pathomics-ml) framework.

### Key Results (held-out test set, n=27)
| Configuration | DGF C-statistic | vs KDPI (C=0.410) |
|---|---|---|
| Hierarchical T (cluster-aware) | 0.792 | +93% |
| Naive GTAV (all-tissue average) | 0.847 | +106% |

---

## Pipeline Structure

```
pipeline/
├── step1_clustering/        ← Global K-means + hierarchical clustering
├── step2_fuzzy_forests/     ← Fuzzy Forests feature selection (top 100)
├── step3_compreps/          ← ComPRePS 100-bootstrap ML training
├── step4_test_evaluation/   ← Held-out test set evaluation (n=27)
└── step5_figures/           ← Publication figures (elbow, ROC, importance)
```

### Step 1: Clustering
Builds cluster-averaged tissue data files for 15 tissue combinations (G, T, A, V and all multi-tissue combos) across 3 pipelines:
- **Global K-means**: cluster across all subjects, average within cluster
- **Hierarchical**: cluster within each patient, average within cluster  
- **Naive**: slide-level average across all objects (baseline)

Key scripts:
- `run_hierarchical_clustering.py` — hierarchical within-patient clustering
- `build_single_tissue_files.R` — build per-tissue combination CSVs
- `build_all_FF_inputs.R` — master builder for all FF input files

### Step 2: Fuzzy Forests Feature Selection
Runs Fuzzy Forests on clustered data to select top 100 features per config (15 combos × 2 methods × 2 outcomes = 60 configs).

Key scripts:
- `run_FF_{DGF|eGFR}_{method}_{combo}.R` — one script per config

### Step 3: ComPRePS Bootstrap Training
100-bootstrap training of RF, Lasso, Ridge, and Elastic Net for all 90 configs (45 FF + 45 naive × 2 outcomes).

Key scripts:
- `run_ComPRePS_{method}_{combo}_FF{outcome}.R` — FF pipeline configs
- `run_naive_{outcome}_{combo}.R` — naive pipeline configs
- `Helper_CV_Functions.R` — core cross-validation functions

### Step 4: Test Set Evaluation
Applies trained models to locked held-out test set (n=27, seed=382025).

Key scripts:
- `run_all_test_evaluations.R` — **main reproducible evaluation script**
- `original_run_test_evaluation.R` — Jeremy's original 9-config script (reference)

To reproduce:
```bash
module load R/4.5.0
Rscript pipeline/step4_test_evaluation/code/run_all_test_evaluations.R
```

### Step 5: Figures
Publication-ready figures.

| Figure | Script | Description |
|---|---|---|
| Fig 4A | `fig4a.R` | eGFR MSE vs N MRMR features (elbow plots) |
| Fig 5A | `fig5a.R` | DGF C-statistic vs N MRMR features (elbow plots) |
| Fig 5B/6B | `FF_impt.R` | FF variable importance bar plots |
| ROC curves | `roc_curvesonesinglepanelnew.R` | DGF ROC curves, held-out test set |

---

## Dependencies

```r
library(randomForest)
library(glmnet)
library(pROC)
library(fuzzyforest)
library(ggplot2)
library(dplyr)
library(tidyr)
```

Python:
```
numpy, pandas, scikit-learn, scipy
```

---

## Data Availability

Raw data (procurement biopsy WSIs and clinical variables) from Coimbra University Hospital (2011–2023) are not publicly available due to patient privacy restrictions. Processed feature matrices and model outputs are available upon reasonable request.

- **Cohort**: n=139 transplanted kidneys (112 train, 27 held-out test)
- **Features**: 59 morphometric, textural, and spatial features per tissue object
- **Tissue types**: Glomeruli (G), Tubules (T), Arteries (A), Arterioles (V)
- **Outcomes**: Delayed Graft Function (DGF, binary), 12-month eGFR (continuous)

---

## Validated Benchmarks

| Metric | Value |
|---|---|
| Random seed | 382025 |
| KDPI DGF C-statistic (test) | 0.410 |
| KDPI eGFR MSE (test) | 328.5 |
| Best cluster-aware DGF (Hier T, test) | 0.792 |
| Best naive DGF (Naive GTAV, test) | 0.847 |

---

## Related Work

This framework builds on the ComPRePS pipeline:

> Rodrigues L, Paul AS, Rubin J et al. Multimodal ComPRePS: Integrating High-dimensional Procurement Biopsy Pathomics and Clinical Data for Prediction of Post Transplant Allograft Outcomes. *CJASN*, 2026.
> [GitHub](https://github.com/jeremysrubin/procurement-biopsy-pathomics-ml)

---

## Authors

**Ketki Joshi** — Postdoctoral Associate, Department of Computational Biology, Cornell University  
In collaboration with the Rubin Lab, University of Maryland
