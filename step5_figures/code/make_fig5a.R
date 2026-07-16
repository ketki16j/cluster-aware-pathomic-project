.libPaths('/home/kxj190026/R_libs')
library(ggplot2)
library(tidyr)
library(dplyr)

configs <- list(
  list(label="Hier A",     folder="results/FF_compres/DGF_FF_hier_A_C",     opt_n=24),
  list(label="Hier T",     folder="results/FF_compres/DGF_FF_hier_T_C",     opt_n=12),
  list(label="Global AV",  folder="results/FF_compres/DGF_FF_global_AV_C",  opt_n=32),
  list(label="Global GV",  folder="results/FF_compres/DGF_FF_global_GV_C",  opt_n=13),
  list(label="Naive GT",   folder="results/naive_compres/DGF_naive_GT_C",   opt_n=17),
  list(label="Naive GTAV", folder="results/naive_compres/DGF_naive_GTAV_C", opt_n=27)
)

all_data  <- list()
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
  df_long <- pivot_longer(df_wide, cols=c(RF,Lasso,Ridge,Enet,KDPI),
                          names_to="Model", values_to="C_stat")
  df_long$Config <- cfg$label
  all_data[[length(all_data)+1]] <- df_long

  # Optimal point (RF only)
  opt_rf <- df_wide$RF[df_wide$N == cfg$opt_n]
  opt_points[[length(opt_points)+1]] <- data.frame(
    Config=cfg$label, N=cfg$opt_n, C_stat=opt_rf, Model="RF")
}

plot_df  <- do.call(rbind, all_data)
opt_df   <- do.call(rbind, opt_points)

# Factor ordering for facets
plot_df$Config <- factor(plot_df$Config,
  levels=c("Hier A","Hier T","Global AV","Global GV","Naive GT","Naive GTAV"))
opt_df$Config  <- factor(opt_df$Config,
  levels=levels(plot_df$Config))

# Model styling — match Jeremy's paper style
model_colors <- c(RF="#2C7BB6", Lasso="#D7191C", Ridge="#F46D43",
                  Enet="#1A9641", KDPI="gray50")
model_lines  <- c(RF="solid", Lasso="dashed", Ridge="dotdash",
                  Enet="dotted", KDPI="longdash")
model_size   <- c(RF=1.2, Lasso=0.8, Ridge=0.8, Enet=0.8, KDPI=0.8)

plot_df$Model <- factor(plot_df$Model,
  levels=c("RF","Lasso","Ridge","Enet","KDPI"))
opt_df$Model  <- factor(opt_df$Model, levels=levels(plot_df$Model))

p <- ggplot(plot_df, aes(x=N, y=C_stat, color=Model,
                          linetype=Model, linewidth=Model)) +
  geom_line() +
  # Vertical dashed line at optimal N
  geom_vline(data=opt_df, aes(xintercept=N),
             color="black", linetype="dashed", linewidth=0.5) +
  # Star at optimal RF point
  geom_point(data=opt_df, aes(x=N, y=C_stat),
             color="black", size=3, shape=8,
             inherit.aes=FALSE) +
  # KDPI reference line label
  scale_color_manual(values=model_colors,
    labels=c("Random Forest","Lasso","Ridge","Elastic Net","KDPI")) +
  scale_linetype_manual(values=model_lines,
    labels=c("Random Forest","Lasso","Ridge","Elastic Net","KDPI")) +
  scale_linewidth_manual(values=model_size,
    labels=c("Random Forest","Lasso","Ridge","Elastic Net","KDPI")) +
  facet_wrap(~Config, ncol=2, scales="free_x") +
  scale_y_continuous(limits=c(0.4, 1.0), breaks=seq(0.4,1.0,0.1)) +
  scale_x_continuous(breaks=seq(0,100,20)) +
  labs(
    title   = "DGF Prediction: C-statistic vs. Number of MRMR Features",
    x       = "Number of MRMR-selected features",
    y       = "Bootstrap-validated C-statistic",
    color   = "Model", linetype = "Model", linewidth = "Model",
    caption = "★ = optimal N (best RF); dashed line = optimal N cutoff; KDPI = clinical benchmark"
  ) +
  theme_bw(base_size=11) +
  theme(
    strip.background  = element_rect(fill="#EEF4FB", color="gray70"),
    strip.text        = element_text(face="bold", size=11),
    legend.position   = "bottom",
    legend.title      = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face="bold", size=13),
    plot.caption      = element_text(size=8, color="gray40")
  ) +
  guides(color    = guide_legend(nrow=1),
         linetype = guide_legend(nrow=1),
         linewidth= guide_legend(nrow=1))

ggsave("results/fig5a_DGF_elbow_6configs.png", p,
       width=10, height=12, dpi=200, bg="white")
cat("Saved: results/fig5a_DGF_elbow_6configs.png\n")

# Print optimal point summary
cat("\nOptimal N summary:\n")
print(opt_df[,c("Config","N","C_stat")])
