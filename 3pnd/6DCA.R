#===============================================================
# 6DCA.R
# Calibration and DCA using full held-out test set predictions
#===============================================================

rm(list = ls())

library(dplyr)
library(ggplot2)

setwd("D:/XiangResearch/X2VTE/3PND")

pred_long <- read.csv("4pred_long.csv", stringsAsFactors = FALSE)

#---------------------------------------------------------------
# Keep only the final model of interest
#---------------------------------------------------------------
plot_models <- pred_long %>%
  filter(Dataset == "TOP10", Algorithm == "GLM") %>%
  mutate(
    Model = "Top 10 GLM",
    Dataset_label = "Top 10"
  )

#===============================================================
# Part 1. Calibration curve
#===============================================================

make_calibration_df <- function(df, n_bins = 4) {
  df %>%
    mutate(
      truth_num = ifelse(truth == "Yes", 1, 0),
      bin = ntile(prob_yes, n_bins)
    ) %>%
    group_by(Model, Dataset_label, Algorithm, bin) %>%
    summarise(
      pred_mean = mean(prob_yes, na.rm = TRUE),
      obs_rate = mean(truth_num, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
}

cal_df <- make_calibration_df(plot_models, n_bins = 4)

p_cal <- ggplot(cal_df, aes(x = pred_mean, y = obs_rate, color = Model, group = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.7, color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  labs(
    x = "Predicted probability",
    y = "Observed event rate",
    color = "Model"
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

pdf("6Calibration.pdf", width = 6, height = 5)
print(p_cal)
dev.off()

#===============================================================
# Part 2. Decision curve analysis
#===============================================================

calc_nb <- function(y_true, p, thresholds) {
  n <- length(y_true)
  y_num <- ifelse(y_true == "Yes", 1, 0)
  
  out <- lapply(thresholds, function(pt) {
    pred_pos <- p >= pt
    TP <- sum(pred_pos & y_num == 1)
    FP <- sum(pred_pos & y_num == 0)
    NB <- TP / n - FP / n * (pt / (1 - pt))
    
    data.frame(
      threshold = pt,
      net_benefit = NB
    )
  })
  
  bind_rows(out)
}

calc_nb_all <- function(y_true, thresholds) {
  n <- length(y_true)
  y_num <- ifelse(y_true == "Yes", 1, 0)
  prevalence <- mean(y_num)
  
  bind_rows(lapply(thresholds, function(pt) {
    nb_all <- prevalence - (1 - prevalence) * (pt / (1 - pt))
    data.frame(
      threshold = pt,
      net_benefit = nb_all
    )
  }))
}

calc_nb_none <- function(thresholds) {
  data.frame(
    threshold = thresholds,
    net_benefit = 0
  )
}

thresholds <- seq(0.01, 0.80, by = 0.01)

dca_model_df <- plot_models %>%
  group_by(Model, Dataset_label, Algorithm) %>%
  group_modify(~{
    calc_nb(.x$truth, .x$prob_yes, thresholds)
  }) %>%
  ungroup()

# add treat-all and treat-none separately for each dataset
dca_ref_df <- plot_models %>%
  group_by(Dataset_label) %>%
  group_modify(~{
    y_true <- .x$truth
    
    bind_rows(
      calc_nb_all(y_true, thresholds) %>%
        mutate(Model = "Treat all", Algorithm = "Reference"),
      calc_nb_none(thresholds) %>%
        mutate(Model = "Treat none", Algorithm = "Reference")
    )
  }) %>%
  ungroup()

dca_plot_df <- bind_rows(
  dca_model_df,
  dca_ref_df
)

p_dca <- ggplot(dca_plot_df, aes(x = threshold, y = net_benefit,
                                 color = Model, group = Model)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c(
    "Top 10 GLM" = "#3B5BA5",
    "Treat all"  = "#00A087",
    "Treat none" = "#3C5488"
  )) +
  labs(
    x = "Threshold probability",
    y = "Net benefit",
    color = "Model"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )+ ylim(-0.5, 0.5)

pdf("6DCA.pdf", width = 6, height = 5)
print(p_dca)
dev.off()

#===============================================================
# Part 3. Export threshold ranges with positive incremental benefit
# model net benefit > treat-all and > treat-none(=0)
#===============================================================
dca_model_only <- dca_model_df

dca_all_ref <- dca_ref_df %>%
  filter(Model == "Treat all") %>%
  select(Dataset_label, threshold, nb_all = net_benefit)

dca_sig_df <- dca_model_only %>%
  left_join(dca_all_ref, by = c("Dataset_label", "threshold")) %>%
  mutate(
    better_than_all = net_benefit > nb_all,
    better_than_none = net_benefit > 0
  ) %>%
  filter(better_than_all, better_than_none)

dca_range <- dca_sig_df %>%
  group_by(Model, Dataset_label, Algorithm) %>%
  summarise(
    threshold_min = min(threshold, na.rm = TRUE),
    threshold_max = max(threshold, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(cal_df, "6Calibration_data.csv", row.names = FALSE)
write.csv(dca_plot_df, "6DCA_data.csv", row.names = FALSE)
write.csv(dca_range, "6DCA_threshold_range.csv", row.names = FALSE)

cat("\nDone.\n")
cat("Saved files:\n")
cat("- 6Calibration.pdf\n")
cat("- 6DCA.pdf\n")
cat("- 6Calibration_data.csv\n")
cat("- 6DCA_data.csv\n")
cat("- 6DCA_threshold_range.csv\n")
