# ============================================================
# make_ff_importance.R
# Fuzzy Forests variable importance — 2-panel figure
# Panel 1: DGF | Panel 2: eGFR
# 4 configs each: Hier A, Hier T, Global AV, Global GV
# Single blue color, clean feature names (Jeremy style)
#
# Usage:
#   module load R/4.5.0
#   Rscript make_ff_importance.R
#
# Output:
#   plots/ff_importance_panel_DGF.png
#   plots/ff_importance_panel_eGFR.png
#
# Author: Ketki Joshi
# Date:   June 2026
# ============================================================

.libPaths('/home/kxj190026/R_libs')
library(ggplot2)
library(dplyr)

WORKDIR <- '/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml'
FF_DIR  <- file.path(WORKDIR, 'results/fuzzyforest')
setwd(WORKDIR)
dir.create("plots", showWarnings=FALSE)

# ── 1. Feature name maps per tissue ───────────────────────────────────────

# Helper: clean up raw column name to short readable label
clean_name <- function(x) {
  x <- gsub("_reviewed_MF01$", "", x)
  x <- gsub("_V_Mean_rev$",    "", x)
  x <- gsub("_TF_CT_rev$",     "", x)
  x <- gsub("_TF_Cortical$",   "", x)
  x <- gsub("_TF_Rev$",        "", x)
  x <- gsub("^TUBULAR_",       "", x)
  x <- gsub("^VASCULAR_",      "", x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x
}

# Tubule features (74 total)
tubule_raw <- c(
  "Cortical Tubular Density","Avg Cortical Tubular Area",
  "Std Cortical Tubular Area","Avg Cortical Tubular Radius",
  "Avg Medullary Tubular Area","Avg Medullary Tubular Area 2",
  "Std Medullary Tubular Area","Avg Medullary Tubular Radius",
  "Std Medullary Tubular Radius","Medullary Tubular Density",
  "Medullary Tubular Density Area","Gloms/Cortex/Tubules Ratio",
  "Cortical Tubule Count","Avg TBM Thickness",
  "Avg Cell Thickness","Luminal Fraction",
  "Avg TBM Thickness (Rev)","Luminal Fraction (Rev)",
  "Luminal Fraction (Rev2)","Sum DT/ObjArea × Nuclei",
  "Sum DT/NuclArea × Nuclei","Mean DT/ObjArea × Nuclei",
  "Mean DT/NuclArea × Nuclei","Mean DT × Nuclei",
  "Max DT/ObjArea × Nuclei","Max DT/NuclArea × Nuclei",
  "Max DT × Nuclei","Sum DT/ObjArea × Eosin",
  "Sum DT/EosinArea × Eosin","Sum DT × Eosin",
  "Mean DT/ObjArea × Eosin","Mean DT/EosinArea × Eosin",
  "Mean DT × Eosin","Max DT/ObjArea × Eosin",
  "Max DT/EosinArea × Eosin","Max DT × Eosin",
  "Sum DT/ObjArea × Luminal","Sum DT/LumArea × Luminal",
  "Sum DT × Luminal","Mean DT/ObjArea × Luminal",
  "Mean DT/LumArea × Luminal","Mean DT × Luminal",
  "Max DT/ObjArea × Luminal","Max DT/LumArea × Luminal",
  "Max DT × Luminal","Homogeneity Eosinophilic",
  "Homogeneity Nuclei","Correlation Nuclei",
  "Energy Nuclei","Contrast Eosinophilic",
  "Contrast Luminal Space","Nuclei Number",
  "Avg Aspect Ratio Nuclei","Std Aspect Ratio Nuclei",
  "Mean Nuclear Area","Correlation Luminal Space",
  "Energy Luminal Space","Nuclei Area / Object Area",
  "Nuclei Area","Eosinophilic Area / Object Area",
  "Eosinophilic Area","Luminal Space Area / Object Area",
  "Luminal Space Area","Nuclei Number (Rev)",
  "Mean Aspect Ratio Nuclei","Std Aspect Ratio Nuclei (Rev)",
  "Mean Nuclear Area (Rev)","Total Object Area",
  "Total Object Perimeter","Total Object Aspect Ratio",
  "Major Axis Length","Minor Axis Length",
  "Area (pixel²)","Radius (pixel)"
)
tubule_map <- setNames(tubule_raw, paste0("feat", seq_along(tubule_raw)))

# Glomeruli features (77 total)
glom_raw <- c(
  "Avg TBM Thickness","Avg Cell Thickness","Luminal Fraction",
  "Avg TBM Thickness (Rev)","Avg Cell Thickness (Rev)","Luminal Fraction (Rev)",
  "Luminal Ratio (Mean)","Luminal Ratio (Mean Rev)",
  "Luminal Ratio (Max Area)","Luminal Ratio (MaxAr Rev)",
  "Luminal Ratio (Arteriolar)","Luminal Ratio (INF35 Rev)",
  "Luminal Ratio (Artery Mean)","Luminal Ratio (SUP35 Rev)",
  "Sum DT/ObjArea × Nuclei","Sum DT/NuclArea × Nuclei",
  "Sum DT × Nuclei","Mean DT/ObjArea × Nuclei",
  "Mean DT/NuclArea × Nuclei","Mean DT × Nuclei",
  "Max DT/ObjArea × Nuclei","Max DT/NuclArea × Nuclei",
  "Max DT × Nuclei","Sum DT/ObjArea × Eosin",
  "Sum DT/EosinArea × Eosin","Sum DT × Eosin",
  "Mean DT/ObjArea × Eosin","Mean DT/EosinArea × Eosin",
  "Mean DT × Eosin","Max DT/ObjArea × Eosin",
  "Max DT/EosinArea × Eosin","Max DT × Eosin",
  "Sum DT/ObjArea × Luminal","Sum DT/LumArea × Luminal",
  "Sum DT × Luminal","Mean DT/ObjArea × Luminal",
  "Mean DT/LumArea × Luminal","Mean DT × Luminal",
  "Max DT/ObjArea × Luminal","Max DT/LumArea × Luminal",
  "Max DT × Luminal","Sum DT/LumArea × Luminal (2)",
  "Sum DT × Luminal (2)","Mean DT/ObjArea × Luminal (2)",
  "Mean DT/LumArea × Luminal (2)","Mean DT × Luminal (2)",
  "Max DT/ObjArea × Luminal (2)","Max DT/LumArea × Luminal (2)",
  "Max DT × Luminal (2)","Contrast Nuclei",
  "Homogeneity Nuclei","Correlation Nuclei",
  "Energy Nuclei","Contrast Eosinophilic",
  "Homogeneity Eosinophilic","Correlation Eosinophilic",
  "Energy Eosinophilic","Contrast Luminal Space",
  "Homogeneity Luminal Space","Correlation Luminal Space",
  "Energy Luminal Space","Nuclei Area / Object Area",
  "Nuclei Area","Eosinophilic Area / Object Area",
  "Eosinophilic Area","Luminal Space Area / Object Area",
  "Luminal Space Area","Nuclei Number",
  "Mean Aspect Ratio Nuclei","Std Aspect Ratio Nuclei",
  "Mean Nuclear Area","Total Object Area",
  "Total Object Perimeter","Total Object Aspect Ratio",
  "Major Axis Length","Area (pixel²)","Radius (pixel)"
)
glom_map <- setNames(glom_raw, paste0("feat", seq_along(glom_raw)))

# Artery features (250 total — map key ones)
artery_raw <- c(
  "Cortical Artery Density","Cortical Artery Count",
  "Avg TBM Thickness","Avg Cell Thickness","Luminal Fraction",
  "Avg TBM Thickness (Rev)","Avg Cell Thickness (Rev)","Luminal Fraction (Rev)",
  "Luminal Ratio (Mean)","Arterial Area (Mean)",
  "Luminal Ratio (Mean Rev)","Arterial Area (Mean Rev)",
  "Luminal Ratio (Max Area)","Arterial Area (Max Area)",
  "Luminal Ratio (MaxAr Rev)","Arterial Area (MaxAr Rev)",
  "Luminal Ratio (Arteriolar)","Arterial Area (Arteriolar)",
  "Luminal Ratio (INF35 Rev)","Arterial Area (INF35 Rev)",
  "Luminal Ratio (Artery Mean)","Arterial Area (Artery Mean)",
  "Luminal Ratio (SUP35 Rev)","Arterial Area (SUP35 Rev)",
  "Sum DT/ObjArea × Nuclei","Sum DT/NuclArea × Nuclei",
  "Sum DT × Nuclei","Mean DT/ObjArea × Nuclei",
  "Mean DT/NuclArea × Nuclei","Mean DT × Nuclei",
  "Max DT/ObjArea × Nuclei","Max DT/NuclArea × Nuclei",
  "Max DT × Nuclei","Sum DT/ObjArea × Eosin",
  "Sum DT/EosinArea × Eosin","Sum DT × Eosin",
  "Mean DT/ObjArea × Eosin","Mean DT/EosinArea × Eosin",
  "Mean DT × Eosin","Max DT/ObjArea × Eosin",
  "Max DT/EosinArea × Eosin","Max DT × Eosin",
  "Sum DT/ObjArea × Luminal","Sum DT/LumArea × Luminal",
  "Sum DT × Luminal","Mean DT/ObjArea × Luminal",
  "Mean DT/LumArea × Luminal","Mean DT × Luminal",
  "Max DT/ObjArea × Luminal","Max DT/LumArea × Luminal",
  "Max DT × Luminal","Contrast Nuclei",
  "Homogeneity Nuclei","Correlation Nuclei",
  "Energy Nuclei","Contrast Eosinophilic",
  "Homogeneity Eosinophilic","Correlation Eosinophilic",
  "Energy Eosinophilic","Contrast Luminal Space",
  "Homogeneity Luminal Space","Correlation Luminal Space",
  "Energy Luminal Space","Nuclei Area / Object Area",
  "Nuclei Area","Eosinophilic Area / Object Area",
  "Eosinophilic Area","Luminal Space Area / Object Area",
  "Luminal Space Area","Nuclei Number",
  "Mean Aspect Ratio Nuclei","Std Aspect Ratio Nuclei",
  "Mean Nuclear Area","Total Object Area",
  "Total Object Perimeter","Total Object Aspect Ratio",
  "Major Axis Length","Area (pixel²)","Radius (pixel)",
  "Granular Features Artery"
)
artery_map <- setNames(artery_raw, paste0("feat", seq_along(artery_raw)))

# Arteriole features (22 total)
arteriole_raw <- c(
  "Avg TBM Thickness","Avg Cell Thickness","Luminal Fraction",
  "Avg TBM Thickness (Rev)","Avg Cell Thickness (Rev)","Luminal Fraction (Rev)",
  "Luminal Ratio (Mean)","Arterial Area (Mean)",
  "Luminal Ratio (Mean Rev)","Arterial Area (Mean Rev)",
  "Luminal Ratio (Max Area)","Arterial Area (Max Area)",
  "Luminal Ratio (MaxAr Rev)","Arterial Area (MaxAr Rev)",
  "Luminal Ratio (Arteriolar)","Arterial Area (Arteriolar)",
  "Luminal Ratio (INF35 Rev)","Arterial Area (INF35 Rev)",
  "Luminal Ratio (Artery Mean)","Arterial Area (Artery Mean)",
  "Luminal Ratio (SUP35 Rev)","Arterial Area (SUP35 Rev)"
)
arteriole_map <- setNames(arteriole_raw,
                          paste0("feat", seq_along(arteriole_raw)))

# Combined lookup by tissue
tissue_maps <- list(
  tubules    = tubule_map,
  glomeruli  = glom_map,
  arteries   = artery_map,
  arterioles = arteriole_map,
  scl_gloms  = glom_map   # sclerotic gloms use same morphometrics
)

# ── 2. Clinical feature prefixes ──────────────────────────────────────────
clinical_patterns <- c(
  "Donor_","REMUZZI","Expanded_criteria","Fibro","IFTA",
  "Glomerulosclerosis_PERC","Glomerulus_count","sGlomerulus",
  "Average_glomerular","Std_glomerular","Cortical_interstitial",
  "Total_nephron","Cortical_tubule_count","Gloms_cortex",
  "Fibrointimal","VASCULAR","Max_Distance_TF_By"
)

# ── 3. Parse feature name → clean readable label ──────────────────────────
parse_feature <- function(feat_name) {
  # Clinical
  is_clinical <- any(sapply(clinical_patterns,
                            function(p) startsWith(feat_name, p)))
  if (is_clinical) {
    label <- feat_name
    label <- gsub("^Donor_", "", label)
    label <- gsub("_0_Male_1_Female$", " (M/F)", label)
    label <- gsub("_", " ", label)
    label <- tools::toTitleCase(tolower(label))
    return(trimws(label))
  }

  # Cluster pathomic — two naming conventions:
  # Long:  tubules_hc2_feat58, glomeruli_cluster6_feat22
  # Short: T_cluster2_feat58,  G_cluster6_feat22, A_cluster2_feat71
  pat_long  <- "^(tubules|glomeruli|arteries|arterioles|scl_gloms)_(hc|cluster)(\\d+)_feat(\\d+)$"
  pat_short <- "^([TGAVS])_(hc|cluster)(\\d+)_feat(\\d+)$"

  tissue_raw <- NULL; feat_n <- NULL

  if (grepl(pat_long, feat_name, perl=TRUE)) {
    tissue_raw <- sub(pat_long, "\\1", feat_name, perl=TRUE)
    feat_n     <- sub(pat_long, "\\4", feat_name, perl=TRUE)
  } else if (grepl(pat_short, feat_name, perl=TRUE)) {
    prefix <- sub(pat_short, "\\1", feat_name, perl=TRUE)
    feat_n <- sub(pat_short, "\\4", feat_name, perl=TRUE)
    tissue_raw <- switch(prefix,
      T = "tubules", G = "glomeruli",
      A = "arteries", V = "arterioles",
      S = "scl_gloms", "tubules")
  }

  if (!is.null(tissue_raw) && !is.null(feat_n)) {
    feat_key <- paste0("feat", feat_n)
    tmap <- tissue_maps[[tissue_raw]]
    if (!is.null(tmap) && feat_key %in% names(tmap))
      return(tmap[[feat_key]])
    else
      return(feat_key)
  }

  gsub("_", " ", feat_name)
}

# ── 4. Load FF file → top N features ─────────────────────────────────────
load_ff <- function(filepath, config_label, top_n=15) {
  df <- read.csv(filepath, stringsAsFactors=FALSE)
  df <- df[order(-df$variable_importance), ]
  df <- head(df, top_n)
  df$vi_norm <- 100 * df$variable_importance / max(df$variable_importance)
  df$label   <- sapply(df$feature_name, parse_feature)
  df$Config  <- config_label
  df$label   <- ifelse(nchar(df$label) > 45,
                       paste0(substr(df$label, 1, 42), "..."),
                       df$label)
  df
}

# ── 5. Plot panel — single blue color (Jeremy style) ─────────────────────
make_panel <- function(configs, outcome_label, out_file) {
  all_data <- list()
  for (cfg in configs) {
    fp <- file.path(FF_DIR, cfg$file)
    if (!file.exists(fp)) {
      cat(sprintf("  WARNING: %s not found\n", fp)); next
    }
    df <- load_ff(fp, cfg$label)
    cat(sprintf("  %-15s  top: %s (%.1f)\n",
                cfg$label, df$label[1], df$vi_norm[1]))
    all_data[[length(all_data)+1]] <- df
  }

  plot_df <- do.call(rbind, all_data)

  plot_df <- plot_df %>%
    group_by(Config) %>%
    mutate(label=factor(label,
                        levels=unique(label[order(vi_norm)]))) %>%
    ungroup()

  plot_df$Config <- factor(plot_df$Config,
                           levels=sapply(configs, `[[`, "label"))

  p <- ggplot(plot_df, aes(x=label, y=vi_norm)) +
    geom_col(width=0.72, fill="#2166AC",
             color="white", linewidth=0.2) +
    geom_text(aes(label=sprintf("%.1f", vi_norm)),
              hjust=-0.1, size=2.6, color="gray30") +
    coord_flip() +
    facet_wrap(~Config, nrow=1, scales="free_y") +
    scale_y_continuous(limits=c(0, 118),
                       breaks=seq(0, 100, 25)) +
    labs(
      title    = sprintf(
        "Top 15 Fuzzy Forests Features — %s Prediction",
        outcome_label),
      subtitle = paste0(
        "FF-selected features input to ComPRePS multimodal ML ",
        "framework  |  Importance normalized to top feature = 100%"),
      x        = NULL,
      y        = "Variable Importance (normalized)"
    ) +
    theme_bw(base_size=10.5) +
    theme(
      strip.background   = element_rect(fill="#EEF4FB", color="gray70"),
      strip.text         = element_text(face="bold", size=10),
      legend.position    = "none",
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color="gray92"),
      plot.title         = element_text(face="bold", size=13),
      plot.subtitle      = element_text(size=8.5, color="gray40"),
      axis.text.y        = element_text(size=7.5)
    )

  ggsave(out_file, p, width=18, height=7, dpi=200, bg="white")
  cat(sprintf("  Saved: %s\n\n", out_file))
}

# ── 6. Generate panels ────────────────────────────────────────────────────
cat("\n=== Panel 1: DGF ===\n")
make_panel(
  list(
    list(file="FF_DGF_hier_A_top100.csv",    label="Hier A"),
    list(file="FF_DGF_hier_T_top100.csv",    label="Hier T"),
    list(file="FF_DGF_global_AV_top100.csv", label="Global AV"),
    list(file="FF_DGF_global_GV_top100.csv", label="Global GV")
  ),
  "DGF",
  "plots/ff_importance_panel_DGF.png"
)

cat("=== Panel 2: eGFR ===\n")
make_panel(
  list(
    list(file="FF_eGFR_hier_A_top100.csv",    label="Hier A"),
    list(file="FF_eGFR_hier_T_top100.csv",    label="Hier T"),
    list(file="FF_eGFR_global_AV_top100.csv", label="Global AV"),
    list(file="FF_eGFR_global_GV_top100.csv", label="Global GV")
  ),
  "eGFR (12-month)",
  "plots/ff_importance_panel_eGFR.png"
)

cat("Done!\n")
cat("  cp plots/ff_importance_panel_*.png testing/output/\n")
cat("  cp make_ff_importance.R testing/code/\n")

