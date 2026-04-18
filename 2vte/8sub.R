#===============================================================
# 8sub.R
# Sensitivity analysis for symptomatic VTE
# VTE_clinical is used ONLY as an evaluation label, not as a predictor
#===============================================================

rm(list = ls())

library(dplyr)
library(pROC)
library(caret)
library(ggplot2)

#---------------------------------------------------------------
# 1. Read raw df_VTE to get VTE_clinical
#---------------------------------------------------------------
setwd("D:/XiangResearch/X2VTE/1data")

raw <- read.csv("df_VTE.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
raw$row_id <- rownames(raw)

cat("\nCheck VTE_clinical values:\n")
print(table(raw$VTE_clinical, useNA = "always"))

# In your file, VTE_clinical is coded as 0/1
# 1 = clinically symptomatic VTE
# 0 = not clinically symptomatic
raw$truth_symptomatic <- ifelse(raw$VTE_clinical == 1, "Yes", "No")
raw$truth_symptomatic <- factor(raw$truth_symptomatic, levels = c("No", "Yes"))

raw_symp <- raw[, c("row_id", "truth_symptomatic")]

#---------------------------------------------------------------
# 2. Read model predictions
#---------------------------------------------------------------
setwd("D:/XiangResearch/X2VTE/2vte")

pred <- read.csv("4pred_long.csv", stringsAsFactors = FALSE)

cat("\nPrediction table dimensions:\n")
print(dim(pred))

#---------------------------------------------------------------
# 3. Merge symptomatic label to prediction table
#---------------------------------------------------------------
pred2 <- merge(pred, raw_symp, by = "row_id", all.x = TRUE)

cat("\nMerged table dimensions:\n")
print(dim(pred2))

cat("\nCheck truth_symptomatic after merge:\n")
print(table(pred2$truth_symptomatic, useNA = "always"))

pred2$truth_symptomatic <- factor(pred2$truth_symptomatic, levels = c("No", "Yes"))

#---------------------------------------------------------------
# 4. Function to calculate symptomatic metrics
#---------------------------------------------------------------
calc_symptomatic_metrics <- function(df_one) {
  
  df_one <- df_one[!is.na(df_one$truth_symptomatic) & !is.na(df_one$prob_yes), , drop = FALSE]
  
  if (length(unique(df_one$truth_symptomatic)) < 2) {
    return(data.frame(
      AUC_symptomatic = NA_real_,
      AUC_LCL_symptomatic = NA_real_,
      AUC_UCL_symptomatic = NA_real_,
      Threshold_Youden_symptomatic = NA_real_,
      Sensitivity_Youden_symptomatic = NA_real_,
      Specificity_Youden_symptomatic = NA_real_,
      Accuracy_Youden_symptomatic = NA_real_
    ))
  }
  
  roc_obj <- roc(
    response = df_one$truth_symptomatic,
    predictor = df_one$prob_yes,
    levels = c("No", "Yes"),
    direction = "<",
    quiet = TRUE
  )
  
  auc_val <- as.numeric(auc(roc_obj))
  ci_vals <- as.numeric(ci.auc(roc_obj))
  
  best_cut <- coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  cutoff <- as.numeric(best_cut["threshold"])
  
  pred_youden <- ifelse(df_one$prob_yes >= cutoff, "Yes", "No")
  pred_youden <- factor(pred_youden, levels = c("No", "Yes"))
  
  cm <- confusionMatrix(pred_youden, df_one$truth_symptomatic, positive = "Yes")
  
  data.frame(
    AUC_symptomatic = auc_val,
    AUC_LCL_symptomatic = ci_vals[1],
    AUC_UCL_symptomatic = ci_vals[3],
    Threshold_Youden_symptomatic = cutoff,
    Sensitivity_Youden_symptomatic = unname(cm$byClass["Sensitivity"]),
    Specificity_Youden_symptomatic = unname(cm$byClass["Specificity"]),
    Accuracy_Youden_symptomatic = unname(cm$overall["Accuracy"])
  )
}

#---------------------------------------------------------------
# 5. Calculate symptomatic AUC and threshold-based metrics
#---------------------------------------------------------------
res_symp <- pred2 %>%
  group_by(Dataset, Algorithm) %>%
  group_modify(~ calc_symptomatic_metrics(.x)) %>%
  ungroup()

cat("\nSymptomatic VTE sensitivity-analysis results:\n")
print(res_symp)

write.csv(res_symp, "8AUC_symptomatic.csv", row.names = FALSE)

#---------------------------------------------------------------
# 6. Merge with main full-test metrics
#---------------------------------------------------------------
metrics_main <- read.csv("4metrics.csv", stringsAsFactors = FALSE)

compare_df <- left_join(metrics_main, res_symp, by = c("Dataset", "Algorithm"))

# relabel dataset
compare_df <- compare_df %>%
  mutate(
    Dataset_label = recode(
      Dataset,
      "CAPRINI" = "Caprini",
      "TOP5" = "Top 5",
      "TOP10" = "Top 10",
      "TOP13" = "Top 13"
    ),
    Dataset_label = factor(Dataset_label, levels = c("Caprini", "Top 5", "Top 10", "Top 13"))
  )

write.csv(compare_df, "8compare_main_vs_symptomatic.csv", row.names = FALSE)

#---------------------------------------------------------------
# 7. Build AUC comparison plot (point + errorbar)
#---------------------------------------------------------------
plot_main <- compare_df %>%
  transmute(
    Dataset,
    Dataset_label,
    Algorithm,
    Analysis = "All VTE",
    AUC_value = AUC,
    AUC_LCL = AUC_LCL,
    AUC_UCL = AUC_UCL
  )

plot_symp <- compare_df %>%
  transmute(
    Dataset,
    Dataset_label,
    Algorithm,
    Analysis = "Symptomatic VTE",
    AUC_value = AUC_symptomatic,
    AUC_LCL = AUC_LCL_symptomatic,
    AUC_UCL = AUC_UCL_symptomatic
  )

plot_auc <- bind_rows(plot_main, plot_symp)

head(plot_auc)
tail(plot_auc)

p_auc_compare <- ggplot(
  plot_auc,
  aes(x = Algorithm, y = AUC_value, color = Analysis, group = Analysis)
) +
  geom_point(position = position_dodge(width = 0.45), size = 2.8) +
  geom_errorbar(
    aes(ymin = AUC_LCL, ymax = AUC_UCL),
    position = position_dodge(width = 0.45),
    width = 0.18,
    linewidth = 0.7
  ) +
  facet_wrap(~ Dataset_label) +
  scale_color_manual(values = c(
    "All VTE"        = "#D55E00",
    "Symptomatic VTE" = "#0072B2"
  )) +
  labs(
    x = NULL,
    y = "AUC",
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

pdf("8AUC_compare_main_vs_symptomatic.pdf", width = 10, height = 6)
print(p_auc_compare)
dev.off()
