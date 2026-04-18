#===============================================================
# 7shap.R
# SHAP analysis for the final TOP10 GLM model using training data
#===============================================================

rm(list = ls())

library(dplyr)
library(fastshap)
library(ggplot2)

setwd("D:/XiangResearch/X2VTE/3pnd/save")

#---------------------------------------------------------------
# Fixed model / dataset
#---------------------------------------------------------------
dataset_name <- "TOP10"
algorithm_name <- "GLM"

train_file <- switch(
  dataset_name,
  "TOP10" = "df_top10_train.csv"
)

df <- read.csv(train_file, row.names = 1)
df$PND <- factor(df$PND, levels = c("No", "Yes"))

X <- df %>% select(-PND)
y <- df$PND

#---------------------------------------------------------------
# Load final model
#---------------------------------------------------------------
setwd("D:/XiangResearch/X2VTE/3pnd")
model_file <- file.path("save", paste0(algorithm_name, "_", dataset_name, "_model.rds"))
fit <- readRDS(model_file)

#---------------------------------------------------------------
# Prediction wrapper for fastshap
#---------------------------------------------------------------
pred_wrapper <- function(object, newdata) {
  as.numeric(predict(object, newdata = newdata, type = "prob")[, "Yes"])
}

#---------------------------------------------------------------
# Compute SHAP values
#---------------------------------------------------------------
set.seed(8)

shap_values <- fastshap::explain(
  object = fit,
  X = X,
  pred_wrapper = pred_wrapper,
  nsim = 100,
  adjust = TRUE
)

#---------------------------------------------------------------
# Long format for plotting
#---------------------------------------------------------------
shap_long <- as.data.frame(shap_values) %>%
  mutate(row_id = rownames(X)) %>%
  tidyr::pivot_longer(
    cols = -row_id,
    names_to = "feature",
    values_to = "shap"
  ) %>%
  left_join(
    X %>%
      mutate(row_id = rownames(X)) %>%
      tidyr::pivot_longer(
        cols = -row_id,
        names_to = "feature",
        values_to = "value"
      ),
    by = c("row_id", "feature")
  )

#---------------------------------------------------------------
# Feature labels
#---------------------------------------------------------------
feature_labels <- c(
  age                   = "Age",
  hypertension_history  = "Hypertension history",
  pre_caprini_score     = "Preop Caprini score",
  pre_lymphocyte_pct    = "Preop lymphocyte",
  pre_albumin_g_l       = "Preop albumin",
  pre_neutrophil_pct    = "Preop neutrophil",
  pre_total_protein_g_l = "Preop total protein",
  tumor_spread          = "Tumor spread",
  tumor_location        = "Tumor location",
  recurrent_glioma      = "Recurrent glioma"
)

shap_long <- shap_long %>%
  mutate(
    feature_label = ifelse(feature %in% names(feature_labels),
                           feature_labels[feature],
                           feature)
  )

#---------------------------------------------------------------
# Global SHAP importance
#---------------------------------------------------------------
shap_importance <- shap_long %>%
  group_by(feature, feature_label) %>%
  summarise(mean_abs_shap = mean(abs(shap), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abs_shap))

write.csv(shap_importance, "7shap_importance.csv", row.names = FALSE)
write.csv(shap_long, "7shap_long.csv", row.names = FALSE)

#---------------------------------------------------------------
# Plot 1: mean |SHAP|
#---------------------------------------------------------------
p_bar <- shap_importance %>%
  mutate(feature_label = factor(feature_label, levels = rev(feature_label))) %>%
  ggplot(aes(x = mean_abs_shap, y = feature_label)) +
  geom_col(width = 0.7) +
  labs(
    x = "Mean absolute SHAP value",
    y = NULL,
    title = paste0("Global feature importance: ", algorithm_name, " ", dataset_name)
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

pdf("7SHAP_bar.pdf", width = 7, height = 5)
print(p_bar)
dev.off()

#---------------------------------------------------------------
# Plot 2: beeswarm-like SHAP plot
#---------------------------------------------------------------
top_features <- head(shap_importance$feature, 15)

plot_df <- shap_long %>%
  filter(feature %in% top_features) %>%
  left_join(
    shap_importance %>% select(feature, mean_abs_shap),
    by = "feature"
  ) %>%
  group_by(feature) %>%
  mutate(
    value_plot = if (all(is.na(value))) {
      NA_real_
    } else if (dplyr::n_distinct(value[!is.na(value)]) <= 1) {
      rep(0.5, dplyr::n())
    } else {
      dplyr::percent_rank(value)
    }
  ) %>%
  ungroup() %>%
  arrange(mean_abs_shap) %>%
  mutate(
    feature_label = factor(
      feature_label,
      levels = unique(feature_label[order(mean_abs_shap)])
    )
  )

p_bee <- ggplot(plot_df, aes(x = shap, y = feature_label, color = value_plot)) +
  geom_point(
    alpha = 0.85,
    size = 2.1,
    position = position_jitter(height = 0.18, width = 0)
  ) +
  scale_color_gradientn(
    colours = c("#2C7BB6", "#ABD9E9", "#F7F7F7", "#FDAE61", "#D7191C"),
    values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = c("Low", "High")
  ) +
  labs(
    x = "SHAP value",
    y = NULL,
    color = "Relative feature value",
    title = paste0("SHAP summary: ", algorithm_name, " ", dataset_name)
  ) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    legend.key.height = grid::unit(1.2, "cm")
  )

pdf("7SHAP_beeswarm.pdf", width = 10, height = 6)
print(p_bee)
dev.off()
