# ============================================================
# make_fig5a.R
# Reproduces Jeremy's Figure 5A style — C-statistic vs N MRMR features
# for 6 selected DGF configurations
#
# Usage:
#   module load R/4.5.0
#   Rscript make_fig5a.R
#
# Output:
#   results/fig5a_DGF_elbow_6configs.png
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

# ── 1. Config definitions ─────────────────────────────────────────────────
configs <- list(
  list(label="Hier A",
       folder="results/FF_compres/DGF_FF_hier_A_C",     opt_n=24),
  list(label="Hier T",
       folder="results/FF_compres/DGF_FF_hier_T_C",     opt_n=12),
  list(label="Global AV",
       folder="results/FF_compres/DGF_FF_global_AV_C",  opt_n=32),
  list(label="Global GV",
       folder="results/FF_compres/DGF_FF_global_GV_C",  opt_n=13),
  list(label="Naive GT",
       folder="results/naive_compres/DGF_naive_GT_C",   opt_n=17),
  list(label="Naive GTAV",
       folder="results/naive_compres/DGF_naive_GTAV_C", opt_n=27)
)

# ── 2. Load bootstrap C-statistics ────────────────────────────────────────
all_data   <- list()
opt_points <- list()

for (cfg in configs) {
  csvs <- Sys.glob(file.path(cfg$folder,
            "*_MRMR_features_training_C_results.csv"))
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
                          cols = c(RF, Lasso, Ridge, Enet, KDPI),
                          names_to  = "Model",
                          values_to = "C_stat")
  df_long$Config <- cfg$label
  all_data[[length(all_data)+1]] <- df_long

  # Optimal RF point
  opt_c <- df_wide$RF[df_wide$N == cfg$opt_n]
  opt_points[[length(opt_points)+1]] <- data.frame(
    Config = cfg$label,
    N      = cfg$opt_n,
    C_stat = opt_c)
}

plot_df <- do.call(rbind, all_data)
opt_df  <- do.call(rbind, opt_points)

# ── 3. Styling ─────────────────────────────────────────────────────────────
# Factor order — controls facet layout
cfg_levels <- c("Hier A","Hier T","Global AV","Global GV","Naive GT","Naive GTAV")
plot_df$Config <- factor(plot_df$Config, levels=cfg_levels)
opt_df$Config  <- factor(opt_df$Config,  levels=cfg_levels)

# Model factor order
model_levels <- c("RF","Lasso","Ridge","Enet","KDPI")
plot_df$Model <- factor(plot_df$Model, levels=model_levels)

# Bright solid colors matching Jeremy's paper palette
model_colors <- c(
  RF    = "#2166AC",   # strong blue
  Lasso = "#E31A1C",   # bright red
  Ridge = "#FF7F00",   # vivid orange
  Enet  = "#33A02C",   # bright green
  KDPI  = "#808080"    # gray reference
)

# ALL solid lines — no dashes
model_lines <- c(
  RF    = "solid",
  Lasso = "solid",
  Ridge = "solid",
  Enet  = "solid",
  KDPI  = "solid"
)

# Line widths — RF slightly thicker
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
p <- ggplot(plot_df, aes(x = N, y = C_stat,
                          color    = Model,
                          linetype = Model,
                          linewidth = Model)) +
  geom_line() +

  # Vertical line at optimal N
  geom_vline(data    = opt_df,
             aes(xintercept = N),
             color   = "black",
             linetype = "dashed",
             linewidth = 0.5,
             inherit.aes = FALSE) +

  # Star at optimal RF value
  geom_point(data  = opt_df,
             aes(x = N, y = C_stat),
             shape = 8,          # star
             size  = 3,
             color = "black",
             stroke = 1.2,
             inherit.aes = FALSE) +

  scale_color_manual(values = model_colors, labels = model_labels) +
  scale_linetype_manual(values = model_lines,  labels = model_labels) +
  scale_linewidth_manual(values = model_widths, labels = model_labels) +

  facet_wrap(~ Config, ncol = 2, scales = "free_x") +

  scale_y_continuous(limits = c(0.40, 1.00),
                     breaks = seq(0.40, 1.00, 0.10),
                     labels = sprintf("%.2f", seq(0.40, 1.00, 0.10))) +
  scale_x_continuous(breaks = seq(0, 100, 25)) +

  labs(
    title    = "DGF Prediction: Bootstrap C-statistic vs. Number of MRMR Features",
    x        = "Number of MRMR-selected features",
    y        = "Bootstrap-validated C-statistic",
    color    = NULL,
    linetype = NULL,
    linewidth = NULL,
    caption  = "★ = optimal N (best RF bootstrap); dashed line = selected N cutoff"
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
out_file <- "plots/fig5a_DGF_elbow_6configs.png"
ggsave(out_file, p, width = 10, height = 12, dpi = 200, bg = "white")
cat(sprintf("Saved: %s\n", out_file))

# ── 6. Summary table ───────────────────────────────────────────────────────
cat("\nOptimal N summary:\n")
print(opt_df)
