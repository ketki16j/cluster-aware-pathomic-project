.libPaths('/home/kxj190026/R_libs')
library(pROC)
library(ggplot2)

WORKDIR  <- '/scratch/juno/kxj190026/others/jeremy/procurement-biopsy-pathomics-ml'
PRED_DIR <- file.path(WORKDIR, 'testing/output/predictions')
setwd(WORKDIR)
dir.create("plots", showWarnings=FALSE)

# Load predictions
df         <- read.csv(file.path(PRED_DIR,
               "DGF_FF_hier_T_C_test_predictions_DGF_N12.csv"))
naive_df   <- read.csv(file.path(PRED_DIR,
               "DGF_naive_GTAV_C_test_predictions_DGF_N27.csv"))

# ROC objects
roc_rf    <- roc(df$True_DGF, df$Pred_RF,           quiet=TRUE)
roc_kdpi  <- roc(df$True_DGF, df$KDPI, direction=">",quiet=TRUE)
roc_naive <- roc(naive_df$True_DGF, naive_df$Pred_RF,quiet=TRUE)

auc_rf    <- round(as.numeric(auc(roc_rf)),    3)
auc_kdpi  <- round(as.numeric(auc(roc_kdpi)),  3)
auc_naive <- round(as.numeric(auc(roc_naive)), 3)

cat(sprintf("Hier T RF:         C=%.3f\n", auc_rf))
cat(sprintf("KDPI:              C=%.3f\n", auc_kdpi))
cat(sprintf("Naive GTAV RF ref: C=%.3f\n", auc_naive))

# Build plot data
roc_df <- rbind(
  data.frame(FPR=1-roc_rf$specificities,
             TPR=roc_rf$sensitivities,
             Model=sprintf("Hier T RF (C=%.3f)", auc_rf)),
  data.frame(FPR=1-roc_kdpi$specificities,
             TPR=roc_kdpi$sensitivities,
             Model=sprintf("KDPI (C=%.3f)", auc_kdpi)),
  data.frame(FPR=1-roc_naive$specificities,
             TPR=roc_naive$sensitivities,
             Model=sprintf("Naive GTAV RF (C=%.3f)", auc_naive))
)

model_colors <- c(
  sprintf("Hier T RF (C=%.3f)", auc_rf)        = "#2166AC",
  sprintf("KDPI (C=%.3f)", auc_kdpi)            = "#808080",
  sprintf("Naive GTAV RF (C=%.3f)", auc_naive)  = "#000000"
)
model_lines <- c(
  sprintf("Hier T RF (C=%.3f)", auc_rf)        = "solid",
  sprintf("KDPI (C=%.3f)", auc_kdpi)            = "solid",
  sprintf("Naive GTAV RF (C=%.3f)", auc_naive)  = "dashed"
)
model_widths <- c(
  sprintf("Hier T RF (C=%.3f)", auc_rf)        = 1.4,
  sprintf("KDPI (C=%.3f)", auc_kdpi)            = 1.0,
  sprintf("Naive GTAV RF (C=%.3f)", auc_naive)  = 1.0
)

roc_df$Model <- factor(roc_df$Model,
  levels=c(sprintf("Hier T RF (C=%.3f)", auc_rf),
           sprintf("KDPI (C=%.3f)", auc_kdpi),
           sprintf("Naive GTAV RF (C=%.3f)", auc_naive)))

p <- ggplot(roc_df, aes(x=FPR, y=TPR,
                         color=Model, linetype=Model,
                         linewidth=Model)) +
  geom_line() +
  geom_abline(slope=1, intercept=0,
              color="gray80", linetype="dotted", linewidth=0.4) +
  scale_color_manual(values=model_colors, name=NULL) +
  scale_linetype_manual(values=model_lines, name=NULL) +
  scale_linewidth_manual(values=model_widths, name=NULL) +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
  labs(
    title   = "ROC Curve — DGF Prediction (Hier T, n=27 test set)",
    x       = "False Positive Rate (1 - Specificity)",
    y       = "True Positive Rate (Sensitivity)",
    caption = "Dashed = Naive GTAV RF reference; dotted = random classifier"
  ) +
  theme_bw(base_size=13) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size=11),
    legend.key.width = unit(1.8, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color="gray92"),
    plot.title       = element_text(face="bold", size=14),
    plot.caption     = element_text(size=9, color="gray50")
  ) +
  guides(
    color     = guide_legend(nrow=1,
                             override.aes=list(linewidth=1.2)),
    linetype  = guide_legend(nrow=1),
    linewidth = guide_legend(nrow=1)
  )

ggsave("plots/roc_hierT_DGF.png", p,
       width=6, height=6, dpi=200, bg="white")
cat("Saved: plots/roc_hierT_DGF.png\n")
