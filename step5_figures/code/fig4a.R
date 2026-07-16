# ============================================================
# make_fig4a.R
# Reproduces Jeremy's Figure 4A style — bootstrap MSE vs N MRMR features
# for 6 selected eGFR configurations
#
# IMPORTANT FOLDER NOTE:
#   FF configs (hier/global): bootstrap MSE CSVs live in _MSE folders
#   Naive configs:            bootstrap MSE CSVs live in _C folders
#
# Usage:
#   module load R/4.5.0
#   Rscript make_fig4a.R
#
# Output:
#   plots/fig4a_eGFR_elbow_6configs.png
#
# Author: Ketki Joshi
# Date:   June 2026
# ============================================================

.libPaths('/home/kxj190026/R_libs')
library(ggplot2)
library(tidyr)
library(dplyr)

WORKDIR <- '/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml'
setwd(WORKDIR)
dir.create("plots", showWarnings=FALSE)

# ── 1. Config definitions ─────────────────────────────────────────────────
# bootstrap_folder: where MSE bootstrap CSVs live (for N selection + plotting)
# data_folder:      where train/test CSVs live (for model evaluation)
configs <- list(
  list(label="Hier A",
       bootstrap_folder = "results/FF_compres/eGFR_FF_hier_A_MSE",
       opt_n = 29),
  list(label="Hier T",
       bootstrap_folder = "results/FF_compres/eGFR_FF_hier_T_MSE",
       opt_n = 12),
  list(label="Global AV",
       bootstrap_folder = "results/FF_compres/eGFR_FF_global_AV_MSE",
       opt_n = 9),
  list(label="Global GV",
       bootstrap_folder = "results/FF_compres/eGFR_FF_global_GV_MSE",
       opt_n = 11),
  list(label="Naive GT",
       bootstrap_folder = "results/naive_compres/eGFR_naive_GT_C",
       opt_n = 9),
  list(label="Naive GTAV",
       bootstrap_folder = "results/naive_compres/eGFR_naive_GTAV_C",
       opt_n = 6)
)

# ── 2. Load bootstrap MSE values ──────────────────────────────────────────
all_data   <- list()
opt_points <- list()

for (cfg in configs) {
  csvs <- Sys.glob(file.path(cfg$bootstrap_folder,
            "*_MRMR_features_training_MSE_results.csv"))
  if (length(csvs) == 0) {
    cat(sprintf("WARNING: no CSVs found in %s\n", cfg$bootstrap_folder))
    next
  }

  rows <- list()
  for (f in csvs) {
    tryCatch({
      df <- read.csv(f, header=FALSE, stringsAsFactors=FALSE)
      rd <- setNames(trimws(df[,2]), trimws(df[,1]))
      rows[[length(rows)+1]] <- data.frame(
        N     = as.integer(as.numeric(rd["Number of MRMR features"])),
        RF    = as.numeric(rd["Random forest"]),
        Lasso = as.numeric(rd["Lasso"]),
        Ridge = as.numeric(rd["Ridge"]),
        Enet  = as.numeric(rd["Elastic net"]),
        KDPI  = as.numeric(rd["KDPI"])
      )
    }, error=function(e){})
  }

  df_wide <- do.call(rbind, rows)
  df_wide <- df_wide[order(df_wide$N), ]

  # Long format
  df_long <- pivot_longer(df_wide,
                          cols      = c(RF, Lasso, Ridge, Enet, KDPI),
                          names_to  = "Model",
                          values_to = "MSE")
  df_long$Config <- cfg$label
  all_data[[length(all_data)+1]] <- df_long

  # Optimal RF point (minimum MSE at opt_n)
  opt_mse <- df_wide$RF[df_wide$N == cfg$opt_n]
  opt_points[[length(opt_points)+1]] <- data.frame(
    Config = cfg$label,
    N      = cfg$opt_n,
    MSE    = opt_mse)

  cat(sprintf("%-12s  N=%2d  RF_MSE=%.1f\n",
              cfg$label, cfg$opt_n, opt_mse))
}

plot_df <- do.call(rbind, all_data)
opt_df  <- do.call(rbind, opt_points)

# ── 3. Styling ─────────────────────────────────────────────────────────────
cfg_levels <- c("Hier A","Hier T","Global AV","Global GV","Naive GT","Naive GTAV")
plot_df$Config <- factor(plot_df$Config, levels=cfg_levels)
opt_df$Config  <- factor(opt_df$Config,  levels=cfg_levels)

model_levels <- c("RF","Lasso","Ridge","Enet","KDPI")
plot_df$Model <- factor(plot_df$Model, levels=model_levels)

# Bright solid colors — same palette as fig5a for consistency
model_colors <- c(
  RF    = "#2166AC",   # strong blue
  Lasso = "#E31A1C",   # bright red
  Ridge = "#FF7F00",   # vivid orange
  Enet  = "#33A02C",   # bright green
  KDPI  = "#808080"    # gray reference
)

model_lines <- c(
  RF    = "solid",
  Lasso = "solid",
  Ridge = "solid",
  Enet  = "solid",
  KDPI  = "solid"
)

model_widths <- c(
  RF    = 1.4,
  Lasso = 0.9,
  Ridge = 0.9,
  Enet  = 0.9,
  KDPI  = 0.9
)

model_labels <- c(
  RF    = "Random Forest",
  Lasso = "Lasso",
  Ridge = "Ridge",
  Enet  = "Elastic Net",
  KDPI  = "KDPI (clinical)"
)

# ── 4. Plot ────────────────────────────────────────────────────────────────
p <- ggplot(plot_df, aes(x = N, y = MSE,
                          color     = Model,
                          linetype  = Model,
                          linewidth = Model)) +
  geom_line() +

  # Vertical dashed line at optimal N
  geom_vline(data        = opt_df,
             aes(xintercept = N),
             color       = "black",
             linetype    = "dashed",
             linewidth   = 0.5,
             inherit.aes = FALSE) +

  # Star at optimal RF MSE point
  geom_point(data        = opt_df,
             aes(x = N, y = MSE),
             shape       = 8,
             size        = 3,
             color       = "black",
             stroke      = 1.2,
             inherit.aes = FALSE) +

  scale_color_manual(values = model_colors, labels = model_labels) +
  scale_linetype_manual(values = model_lines,  labels = model_labels) +
  scale_linewidth_manual(values = model_widths, labels = model_labels) +

  facet_wrap(~ Config, ncol = 2, scales = "free_x") +

  scale_y_continuous(limits = c(100, 450),
                     breaks = seq(100, 450, 50)) +
  scale_x_continuous(breaks = seq(0, 100, 25)) +

  labs(
    title     = "eGFR Prediction: Bootstrap MSE vs. Number of MRMR Features",
    x         = "Number of MRMR-selected features",
    y         = expression("Bootstrap-validated MSE (mL/min/1.73m"^2*")²"),
    color     = NULL,
    linetype  = NULL,
    linewidth = NULL,
    caption   = "★ = optimal N (lowest RF bootstrap MSE); dashed line = selected N cutoff"
  ) +

  theme_bw(base_size = 12) +
  theme(
    strip.background  = element_rect(fill = "#EEF4FB", color = "gray70"),
    strip.text        = element_text(face = "bold", size = 11),
    legend.position   = "bottom",
    legend.key.width  = unit(1.8, "cm"),
    legend.text       = element_text(size = 10),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(color = "gray90"),
    plot.title        = element_text(face = "bold", size = 13),
    plot.caption      = element_text(size = 8, color = "gray50"),
    axis.title        = element_text(size = 11)
  ) +
  guides(
    color     = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)),
    linetype  = guide_legend(nrow = 1),
    linewidth = guide_legend(nrow = 1)
  )

# ── 5. Save ────────────────────────────────────────────────────────────────
out_file <- "plots/fig4a_eGFR_elbow_6configs.png"
ggsave(out_file, p, width = 10, height = 12, dpi = 200, bg = "white")
cat(sprintf("\nSaved: %s\n", out_file))

# ── 6. Summary table ───────────────────────────────────────────────────────
cat("\nOptimal N summary (lowest RF MSE):\n")
print(opt_df)

