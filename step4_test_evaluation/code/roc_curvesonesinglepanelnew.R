# ============================================================
# make_roc_single_panel.R
# Single-panel ROC curve — held-out test set (n=27)
# 7 curves: 5 config RFs + KDPI + Naive GTAV RF (reference)
#
# Usage:
#   module load R/4.5.0
#   Rscript make_roc_single_panel.R
#
# Output:
#   plots/roc_single_panel_DGF.png
#
# Author: Ketki Joshi
# Date:   June 2026
# ============================================================

.libPaths('/home/kxj190026/R_libs')
library(pROC)
library(ggplot2)
library(dplyr)

WORKDIR  <- '/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml'
PRED_DIR <- file.path(WORKDIR, 'testing/output/predictions')
setwd(WORKDIR)
dir.create("plots", showWarnings=FALSE)

# ── 1. Config definitions (5 selected configs, excluding Naive GTAV) ──────
configs <- list(
  list(label="Hier A",     file="DGF_FF_hier_A_C_test_predictions_DGF_N24.csv"),
  list(label="Hier T",     file="DGF_FF_hier_T_C_test_predictions_DGF_N12.csv"),
  list(label="Global AV",  file="DGF_FF_global_AV_C_test_predictions_DGF_N32.csv"),
  list(label="Global GV",  file="DGF_FF_global_GV_C_test_predictions_DGF_N13.csv"),
  list(label="Naive GT",   file="DGF_naive_GT_C_test_predictions_DGF_N17.csv")
)

# ── 2. Load Naive GTAV RF as reference ────────────────────────────────────
naive_gtav_df  <- read.csv(file.path(PRED_DIR,
  "DGF_naive_GTAV_C_test_predictions_DGF_N27.csv"))
naive_gtav_roc <- roc(naive_gtav_df$True_DGF, naive_gtav_df$Pred_RF, quiet=TRUE)
naive_gtav_auc <- round(as.numeric(auc(naive_gtav_roc)), 3)
cat(sprintf("Naive GTAV RF reference C = %.3f\n", naive_gtav_auc))

# ── 3. Helper: ROC object → sorted dataframe ──────────────────────────────
roc_to_df <- function(roc_obj, model, auc_val) {
  df <- data.frame(
    FPR     = 1 - roc_obj$specificities,
    TPR     = roc_obj$sensitivities,
    Model   = model,
    AUC_val = auc_val,
    stringsAsFactors = FALSE
  )
  df[order(df$FPR, df$TPR), ]
}

# ── 4. Build plot data ─────────────────────────────────────────────────────
all_roc <- list()

# 5 config RF curves
for (cfg in configs) {
  df      <- read.csv(file.path(PRED_DIR, cfg$file))
  roc_rf  <- roc(df$True_DGF, df$Pred_RF, quiet=TRUE)
  auc_rf  <- round(as.numeric(auc(roc_rf)), 3)
  cat(sprintf("  %s  RF C = %.3f\n", cfg$label, auc_rf))
  all_roc[[length(all_roc)+1]] <- roc_to_df(roc_rf, cfg$label, auc_rf)
}

# KDPI — use first config file (KDPI is same across all)
df_kdpi  <- read.csv(file.path(PRED_DIR, configs[[1]]$file))
roc_kdpi <- roc(df_kdpi$True_DGF, df_kdpi$KDPI, direction=">", quiet=TRUE)
auc_kdpi <- round(as.numeric(auc(roc_kdpi)), 3)
cat(sprintf("  KDPI  C = %.3f\n", auc_kdpi))
all_roc[[length(all_roc)+1]] <- roc_to_df(roc_kdpi, "KDPI", auc_kdpi)

# Naive GTAV RF reference
all_roc[[length(all_roc)+1]] <- roc_to_df(naive_gtav_roc, "Naive GTAV RF", naive_gtav_auc)

plot_df <- do.call(rbind, all_roc)

# ── 5. Factor ordering ─────────────────────────────────────────────────────
model_levels <- c("Hier A","Hier T","Global AV","Global GV","Naive GT",
                  "KDPI","Naive GTAV RF")
plot_df$Model <- factor(plot_df$Model, levels=model_levels)

# 5 distinct colors for configs, gray for KDPI, black for reference
model_colors <- c(
  "Hier A"       = "#E41A1C",   # red
  "Hier T"       = "#FF7F00",   # orange
  "Global AV"    = "#4DAF4A",   # green
  "Global GV"    = "#984EA3",   # purple
  "Naive GT"     = "#2166AC",   # blue
  "KDPI"         = "#808080",   # gray
  "Naive GTAV RF"= "#000000"    # black
)
model_lines <- c(
  "Hier A"       = "solid",
  "Hier T"       = "solid",
  "Global AV"    = "solid",
  "Global GV"    = "solid",
  "Naive GT"     = "solid",
  "KDPI"         = "solid",
  "Naive GTAV RF"= "dashed"
)
model_widths <- c(
  "Hier A"       = 0.8,
  "Hier T"       = 0.8,
  "Global AV"    = 0.8,
  "Global GV"    = 0.8,
  "Naive GT"     = 0.8,
  "KDPI"         = 0.7,
  "Naive GTAV RF"= 1.0
)

# Legend labels with C-statistic
auc_lookup <- plot_df %>%
  group_by(Model) %>%
  slice(1) %>%
  select(Model, AUC_val)

model_labels <- setNames(
  sprintf("%s (C = %.3f)", as.character(auc_lookup$Model), auc_lookup$AUC_val),
  as.character(auc_lookup$Model)
)

# ── 6. Plot ────────────────────────────────────────────────────────────────
# Separate reference and non-reference for layering
plot_df_ref    <- plot_df %>% filter(Model == "Naive GTAV RF")
plot_df_kdpi   <- plot_df %>% filter(Model == "KDPI")
plot_df_config <- plot_df %>% filter(!Model %in% c("Naive GTAV RF", "KDPI"))

p <- ggplot() +

  # Random classifier diagonal (bottom layer)
  geom_abline(slope=1, intercept=0,
              color="gray80", linetype="dotted", linewidth=0.4) +

  # KDPI — solid gray, semi-transparent
  geom_line(data=plot_df_kdpi,
            aes(x=FPR, y=TPR, color=Model, linetype=Model, linewidth=Model),
            alpha=0.8) +

  # 5 config RF curves — semi-transparent so overlaps are visible
  geom_line(data=plot_df_config,
            aes(x=FPR, y=TPR, color=Model, linetype=Model, linewidth=Model),
            alpha=0.75) +

  # Naive GTAV RF reference — on top, full opacity, dashed black
  geom_line(data=plot_df_ref,
            aes(x=FPR, y=TPR, color=Model, linetype=Model, linewidth=Model),
            alpha=1.0) +

  scale_color_manual(values=model_colors, labels=model_labels, name=NULL) +
  scale_linetype_manual(values=model_lines, labels=model_labels, name=NULL) +
  scale_linewidth_manual(values=model_widths, labels=model_labels, name=NULL) +

  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.25),
                     labels=c("0","0.25","0.50","0.75","1.00")) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25),
                     labels=c("0","0.25","0.50","0.75","1.00")) +

  labs(
    title   = "ROC Curves for DGF Prediction — Held-Out Test Set (n=27)",
    x       = "False Positive Rate (1 - Specificity)",
    y       = "True Positive Rate (Sensitivity)",
    caption = "Dashed black = Naive GTAV RF reference (C=0.847); dotted = random classifier"
  ) +

  theme_bw(base_size=13) +
  theme(
    legend.position   = c(0.72, 0.30),     # inside plot, bottom-right area
    legend.background = element_rect(fill=alpha("white", 0.85), color="gray80"),
    legend.key.width  = unit(1.2, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.text       = element_text(size=9.5),
    legend.spacing.y  = unit(0.15, "cm"),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color="gray92"),
    plot.title        = element_text(face="bold", size=13),
    plot.caption      = element_text(size=8, color="gray50"),
    axis.title        = element_text(size=11)
  ) +
  guides(
    color     = guide_legend(override.aes=list(linewidth=1.0, alpha=1.0)),
    linetype  = guide_legend(),
    linewidth = guide_legend()
  )

# ── 7. Save ────────────────────────────────────────────────────────────────
out_file <- "plots/roc_single_panel_DGF.png"
ggsave(out_file, p, width=9, height=7, dpi=200, bg="white")
cat(sprintf("\nSaved: %s\n", out_file))

