#===============================================================
# 5plots_all.R
# Combined plotting script for:
# AUC / Accuracy / Sensitivity / Specificity / Line plot
# Based on full held-out test set results
#===============================================================

rm(list = ls())

library(dplyr)
library(ggplot2)

setwd("D:/XiangResearch/X2VTE/2vte")

#---------------------------------------------------------------
# Read metrics
#---------------------------------------------------------------
df <- read.csv("4metrics.csv", stringsAsFactors = FALSE)

#---------------------------------------------------------------
# Relabel datasets
#---------------------------------------------------------------
df <- df %>%
  mutate(
    Dataset_label = recode(
      Dataset,
      "CAPRINI" = "Caprini",
      "TOP5"    = "Top 5",
      "TOP10"   = "Top 10",
      "TOP13"   = "Top 13"
    ),
    Dataset_label = factor(Dataset_label, levels = c("Caprini", "Top 5", "Top 10", "Top 13")),
    n_vars = case_when(
      Dataset == "CAPRINI" ~ 1,
      Dataset == "TOP5"    ~ 5,
      Dataset == "TOP10"   ~ 10,
      Dataset == "TOP13"   ~ 13,
      TRUE ~ NA_real_
    )
  )

#---------------------------------------------------------------
# Colors
#---------------------------------------------------------------
dataset_colors <- c(
  "Caprini" = "#0072B2",
  "Top 5"   = "#E69F00",
  "Top 10"  = "#009E73",
  "Top 13"  = "#D55E00"
)

alg_colors <- c(
  "GLM" = "#3B5BA5",  
  "RDF" = "#E64B35",  
  "GBA" = "#4DBBD5",  
  "SVM" = "#7E6148"   
)

#---------------------------------------------------------------
# Generic bar plot function
#---------------------------------------------------------------
make_metric_plot <- function(data, yvar, ylab, file_name, add_auc_ci = FALSE) {
  
  p <- ggplot(data, aes(x = Algorithm, y = .data[[yvar]], fill = Dataset_label)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68)
  
  if (add_auc_ci) {
    p <- p +
      geom_errorbar(
        aes(ymin = AUC_LCL, ymax = AUC_UCL),
        position = position_dodge(width = 0.75),
        width = 0.18,
        linewidth = 0.7
      )
  }
  
  p <- p +
    labs(
      x = NULL,
      y = ylab,
      fill = "Predictor set"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
    ) +
    scale_fill_manual(values = dataset_colors)
  
  pdf(file_name, height = 5, width = 6)
  print(p)
  dev.off()
  
  return(p)
}

#---------------------------------------------------------------
# 1. AUC plot
#---------------------------------------------------------------
p_auc <- make_metric_plot(
  data = df,
  yvar = "AUC",
  ylab = "AUC on full held-out test set",
  file_name = "5AUC.pdf",
  add_auc_ci = TRUE
)

#---------------------------------------------------------------
# 2. Accuracy plot
# Use Youden-threshold metrics if you want clinically actionable cutoff
# or switch to Accuracy_default if preferred
#---------------------------------------------------------------
p_acc <- make_metric_plot(
  data = df,
  yvar = "Accuracy_Youden",
  ylab = "Accuracy on full held-out test set",
  file_name = "5Accuracy.pdf",
  add_auc_ci = FALSE
)

#---------------------------------------------------------------
# 3. Sensitivity plot
#---------------------------------------------------------------
p_sens <- make_metric_plot(
  data = df,
  yvar = "Sensitivity_Youden",
  ylab = "Sensitivity on full held-out test set",
  file_name = "5Sensitivity.pdf",
  add_auc_ci = FALSE
)

#---------------------------------------------------------------
# 4. Specificity plot
#---------------------------------------------------------------
p_spec <- make_metric_plot(
  data = df,
  yvar = "Specificity_Youden",
  ylab = "Specificity on full held-out test set",
  file_name = "5Specificity.pdf",
  add_auc_ci = FALSE
)

#---------------------------------------------------------------
# 5. Line plot: AUC vs number of predictors
#---------------------------------------------------------------
df_line <- df %>%
  filter(!is.na(n_vars)) %>%
  arrange(Algorithm, n_vars)

p_line_auc <- ggplot(df_line, aes(x = n_vars, y = AUC, color = Algorithm, group = Algorithm)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = c(1, 5, 10, 13)) +
  scale_y_continuous(limits = c(0.50, 1.00), breaks = seq(0.5, 1.0, by = 0.1)) +
  scale_color_manual(values = alg_colors) +
  labs(
    x = "Number of predictors",
    y = "AUC on full held-out test set",
    color = "Algorithm"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  ) +
  annotate(
    "segment",
    x = 10, xend = 10,
    y = 0.50, yend = 1.00,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 5,
    y = 0.985,
    label = "Parsimonious choice",
    vjust = 0,
    hjust = -0.05,
    size = 4
  ) +
  coord_cartesian(clip = "off")

pdf("5Line.pdf", height = 5, width = 6)
print(p_line_auc)
dev.off()

#---------------------------------------------------------------
# Optional: export plotting table
#---------------------------------------------------------------
write.csv(df, "5plot_input.csv", row.names = FALSE)

cat("\nDone.\n")
cat("Saved files:\n")
cat("- 5AUC.pdf\n")
cat("- 5Accuracy.pdf\n")
cat("- 5Sensitivity.pdf\n")
cat("- 5Specificity.pdf\n")
cat("- 5Line.pdf\n")
cat("- 5plot_input.csv\n")
