# ============================================================
# make_roc_hierT_naiveGTAV.R
# Side-by-side ROC panels: Hier T | Naive GTAV
# Each panel: RF (blue) + KDPI (gray) + Naive GTAV RF ref (black dashed)
#
# Usage:
#   module load R/4.5.0
#   Rscript make_roc_hierT_naiveGTAV.R
#
# Output:
#   plots/roc_hierT_naiveGTAV_DGF.png
# ============================================================
.libPaths('/home/kxj190026/R_libs')
library(pROC)
library(ggplot2)
library(dplyr)

WORKDIR  <- '/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml'
PRED_DIR <- file.path(WORKDIR, 'testing/output/predictions')
setwd(WORKDIR)
dir.create("plots", showWarnings=FALSE)

# ── 1. Config definitions ─────────────────────────────────────────────────
configs <- list(
  list(label="Hier T",
       file="DGF_FF_hier_T_C_test_predictions_DGF_N12.csv"),
  list(label="Naive GTAV",
       file="DGF_naive_GTAV_C_test_predictions_DGF_N27.csv")
)

# Naive GTAV RF as reference line
naive_gtav_df  <- read.csv(file.path(PRED_DIR,
  "DGF_naive_GTAV_C_test_predictions_DGF_N27.csv"))
naive_gtav_roc <- roc(naive_gtav_df$True_DGF, naive_gtav_df$Pred_RF,
                      quiet=TRUE)
naive_gtav_auc <- round(as.numeric(auc(naive_gtav_roc)), 3)

# ── 2. ROC to dataframe ───────────────────────────────────────────────────
roc_to_df <- function(roc_obj, model, auc_val, config) {
  df <- data.frame(
    FPR     = 1 - roc_obj$specificities,
    TPR     = roc_obj$sensitivities,
    Model   = model,
    AUC_val = auc_val,
    Config  = config,
    stringsAsFactors = FALSE
  )
  df[order(df$FPR, df$TPR), ]
}

# ── 3. Build plot data ─────────────────────────────────────────────────────
all_roc <- list()

for (cfg in configs) {
  df      <- read.csv(file.path(PRED_DIR, cfg$file))
  roc_rf  <- roc(df$True_DGF, df$Pred_RF,             quiet=TRUE)
  roc_kdpi<- roc(df$True_DGF, df$KDPI, direction=">", quiet=TRUE)
  auc_rf  <- round(as.numeric(auc(roc_rf)),   3)
  auc_kdpi<- round(as.numeric(auc(roc_kdpi)), 3)
  cat(sprintf("%-15s  RF C=%.3f  KDPI C=%.3f\n",
              cfg$label, auc_rf, auc_kdpi))

  # RF curve
  all_roc[[length(all_roc)+1]] <- roc_to_df(
    roc_rf, "RF", auc_rf, cfg$label)

  # KDPI curve
  all_roc[[length(all_roc)+1]] <- roc_to_df(
    roc_kdpi, "KDPI", auc_kdpi, cfg$label)

  # Naive GTAV RF reference (same line in both panels)
  # Skip for Naive GTAV panel itself — RF IS the reference
  if (cfg$label != "Naive GTAV") {
    all_roc[[length(all_roc)+1]] <- roc_to_df(
      naive_gtav_roc, "Naive GTAV RF", naive_gtav_auc, cfg$label)
  }
}

plot_df <- do.call(rbind, all_roc)

# ── 4. Factor ordering ─────────────────────────────────────────────────────
cfg_levels   <- c("Hier T","Naive GTAV")
model_levels <- c("RF","KDPI","Naive GTAV RF")
plot_df$Config <- factor(plot_df$Config, levels=cfg_levels)
plot_df$Model  <- factor(plot_df$Model,  levels=model_levels)

# ── 5. Legend labels with C-statistics ────────────────────────────────────
auc_lookup <- plot_df %>%
  group_by(Config, Model) %>%
  slice(1) %>%
  ungroup()

# Build per-panel legend labels
plot_df <- plot_df %>%
  left_join(auc_lookup %>% select(Config, Model, AUC_val),
            by=c("Config","Model")) %>%
  mutate(AUC_val = coalesce(AUC_val.x, AUC_val.y)) %>%
  select(-AUC_val.x, -AUC_val.y)

# Global model labels (AUC shown per panel via annotation)
model_labels <- c(
  "RF"            = "Random Forest",
  "KDPI"          = "KDPI (clinical)",
  "Naive GTAV RF" = "Naive GTAV RF (ref)"
)

# ── 6. Annotation: C values inside each panel ─────────────────────────────
annot_df <- plot_df %>%
  group_by(Config, Model) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    x     = 0.52,
    y     = case_when(
      Model == "RF"            ~ 0.22,
      Model == "KDPI"          ~ 0.14,
      Model == "Naive GTAV RF" ~ 0.06
    ),
    label = sprintf("C = %.3f", AUC_val)
  )

# ── 7. Styling ─────────────────────────────────────────────────────────────
model_colors <- c(
  "RF"            = "#2166AC",
  "KDPI"          = "#808080",
  "Naive GTAV RF" = "#000000"
)
model_lines <- c(
  "RF"            = "solid",
  "KDPI"          = "solid",
  "Naive GTAV RF" = "dashed"
)
model_widths <- c(
  "RF"            = 1.2,
  "KDPI"          = 0.8,
  "Naive GTAV RF" = 1.0
)

# ── 8. Plot ────────────────────────────────────────────────────────────────
p <- ggplot(plot_df,
            aes(x=FPR, y=TPR,
                color=Model, linetype=Model, linewidth=Model)) +
  geom_line(alpha=0.9) +

  # Random classifier diagonal
  geom_abline(slope=1, intercept=0,
              color="gray80", linetype="dotted", linewidth=0.4) +

  # C-statistic annotations
  geom_text(data=annot_df,
            aes(x=x, y=y, label=label, color=Model),
            size=3.5, fontface="bold",
            inherit.aes=FALSE, show.legend=FALSE) +

  scale_color_manual(values=model_colors,
                     labels=model_labels, name=NULL) +
  scale_linetype_manual(values=model_lines,
                        labels=model_labels, name=NULL) +
  scale_linewidth_manual(values=model_widths,
                         labels=model_labels, name=NULL) +

  facet_wrap(~Config, ncol=2) +

  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.25),
                     labels=c("0","0.25","0.50","0.75","1.00")) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25),
                     labels=c("0","0.25","0.50","0.75","1.00")) +

  labs(
    title   = "ROC Curves for DGF Prediction — Held-Out Test Set",
    x       = "False Positive Rate (1 - Specificity)",
    y       = "True Positive Rate (Sensitivity)",
    caption = "Dashed = Naive GTAV RF reference; dotted = random classifier"
  ) +

  theme_bw(base_size=12) +
  theme(
    strip.background  = element_rect(fill="#2166AC", color="#2166AC"),
    strip.text        = element_text(face="bold", size=11,
                                     color="white"),
    legend.position   = "bottom",
    legend.key.width  = unit(1.5, "cm"),
    legend.text       = element_text(size=10),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color="gray92"),
    plot.title        = element_text(face="bold", size=13),
    plot.caption      = element_text(size=8, color="gray50"),
    axis.title        = element_text(size=11)
  ) +
  guides(
    color     = guide_legend(nrow=1,
                             override.aes=list(linewidth=1.2,
                                               alpha=1.0)),
    linetype  = guide_legend(nrow=1),
    linewidth = guide_legend(nrow=1)
  )

# ── 9. Save ────────────────────────────────────────────────────────────────
ggsave("plots/roc_hierT_naiveGTAV_DGF.png", p,
       width=10, height=6, dpi=200, bg="white")
cat("Saved: plots/roc_hierT_naiveGTAV_DGF.png\n")
