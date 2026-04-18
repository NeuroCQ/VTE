#===============================================================
# 2forest.R
#===============================================================

rm(list = ls())

library(dplyr)
library(stringr)
library(ggplot2)

setwd("D:/XiangResearch/X2VTE/1data")
data <- read.csv("df_PND.csv", header = TRUE, row.names = 1)

data$PND <- data$postop_functional_complication
data$postop_functional_complication <- NULL
data$OS <- NULL
data$OSTime <- NULL
data$vte_class <- NULL
data$VTE_clinical <- NULL

data$PND <- factor(data$PND, levels = c("No", "Yes"))

data <- data %>%
  dplyr::select(PND, everything())



setwd("D:/XiangResearch/X2VTE/3pnd")

load('data_set_samples.Rdata')

data_set_samples <- subset(data_set_samples, data_set_samples$dataset == 'Train')

data <- data[which(rownames(data) %in% rownames(data_set_samples)),]


#===============================================================

#===============================================================



vars <- setdiff(names(data), "PND")

uni_res <- lapply(vars, function(v){
  fml <- as.formula(paste("PND ~", v))
  fit <- glm(fml, data = data, family = binomial)
  s <- summary(fit)$coef
  if (nrow(s) < 2) return(NULL)
  
  tibble(
    Variable = v,
    OR  = exp(s[2,1]),
    LCL = exp(s[2,1] - 1.96*s[2,2]),
    UCL = exp(s[2,1] + 1.96*s[2,2]),
    p   = s[2,4]
  )
})


uni_table <- bind_rows(uni_res) %>%
  arrange(p) %>%
  mutate(
    p_bonf = p.adjust(p, method = "bonferroni"),
    p_bh   = p.adjust(p, method = "BH")
  )

head(uni_table, 15)

uni_table <- subset(uni_table, uni_table$p < 0.05)


uni_table_fmt <- uni_table %>%
  mutate(
    OR_CI = sprintf("%.2f (%.2f-%.2f)", OR, LCL, UCL),
    p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  ) %>%
  select(Variable, OR_CI, p)

uni_table_fmt <- as.data.frame(uni_table_fmt)
uni_table_fmt

write.csv(uni_table_fmt, file = '2Uni.csv')

#===============================================================

#===============================================================


feature_labels <- c(
  age = "Age (per year)",
  sex = "Sex (Male vs Female)",
  bmi = "BMI (kg/m^2)",
  smoking_history = "Smoking history (Yes vs No)",
  alcohol_history = "Alcohol history (Yes vs No)",
  hypertension_history = "Hypertension history (Yes vs No)",
  diabetes_history = "Diabetes history (Yes vs No)",
  hyperlipidemia_history = "Hyperlipidemia history (Yes vs No)",
  coronary_heart_disease_history = "Coronary heart disease history (Yes vs No)",
  stroke_history = "Stroke history (Yes vs No)",
  prior_vte_history = "Prior VTE history (Yes vs No)",
  pre_anticoagulant_use = "Preop anticoagulant use (Yes vs No)",
  pre_caprini_score = "Preop Caprini score (per point)",
  pre_d_dimer_mg_l = "Preop D-dimer (mg/L)",
  pre_fibrinogen_g_l = "Preop fibrinogen (g/L)",
  pre_fdp_ug_ml = "Preop FDP (ug/mL)",
  pre_pt_s = "Preop PT (s)",
  pre_aptt_s = "Preop APTT (s)",
  pre_inr = "Preop INR",
  pre_plt_10e9_l = "Preop platelets (10^9/L)",
  pre_rbc_10e12_l = "Preop RBC (10^12/L)",
  pre_wbc_10e9_l = "Preop WBC (10^9/L)",
  pre_neutrophil_pct = "Preop neutrophil (%)",
  pre_neutrophil_abs_10e9_l = "Preop neutrophils abs (10^9/L)",
  pre_monocyte_pct = "Preop monocyte (%)",
  pre_lymphocyte_pct = "Preop lymphocyte (%)",
  pre_albumin_g_l = "Preop albumin (g/L)",
  pre_hb_g_l = "Preop hemoglobin (g/L)",
  pre_total_protein_g_l = "Preop total protein (g/L)",
  pre_creatinine_umol_l = "Preop creatinine (umol/L)",
  pre_urea_mmol_l = "Preop urea (mmol/L)",
  pre_uric_acid_umol_l = "Preop uric acid (umol/L)",
  pre_na_mmol_l = "Preop sodium (mmol/L)",
  pre_k_mmol_l = "Preop potassium (mmol/L)",
  pre_ca_mmol_l = "Preop calcium (mmol/L)",
  pre_cl_mmol_l = "Preop chloride (mmol/L)",
  pre_glucose_mmol_l = "Preop glucose (mmol/L)",
  pre_total_chol_mmol_l = "Preop total cholesterol (mmol/L)",
  pre_triglyceride_mmol_l = "Preop triglyceride (mmol/L)",
  pre_hdl_c_mmol_l = "Preop HDL-C (mmol/L)",
  pre_ldl_c_mmol_l = "Preop LDL-C (mmol/L)",
  pre_tbil_umol_l = "Preop total bilirubin (umol/L)",
  pre_dbil_umol_l = "Preop direct bilirubin (umol/L)",
  pre_alt_u_l = "Preop ALT (U/L)",
  pre_ast_u_l = "Preop AST (U/L)",
  tumor_location = "Tumor location",
  tumor_laterality = "Tumor laterality",
  tumor_spread = "Tumor spread",
  tumor_max_diameter_cm = "Tumor max diameter (cm)",
  recurrent_glioma = "Recurrent glioma (Yes vs No)"
)

df_forest <- uni_table_fmt %>%
  mutate(
    OR  = as.numeric(str_extract(OR_CI, "^[0-9.]+")),
    LCL = as.numeric(str_extract(OR_CI, "(?<=\\()[0-9.]+")),
    UCL = as.numeric(str_extract(OR_CI, "(?<=-)[0-9.]+(?=\\))")),
    Variable_label = ifelse(Variable %in% names(feature_labels), feature_labels[Variable], Variable)
  ) %>%
  arrange(OR) %>%
  mutate(Variable_label = factor(Variable_label, levels = Variable_label))

xmax_plot <- max(df_forest$UCL, na.rm = TRUE)
x_orci <- xmax_plot * 1.8
x_pval <- xmax_plot * 2.7

df_forest <- df_forest %>%
  mutate(
    text_orci = OR_CI,
    text_p = paste0("P = ", p)
  )

p_forest_text <- ggplot(df_forest, aes(x = OR, y = Variable_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.18, linewidth = 0.8) +
  geom_point(size = 2) +
  geom_text(aes(x = x_orci, label = text_orci), hjust = 0, size = 3.4) +
  geom_text(aes(x = x_pval, label = text_p), hjust = 0, size = 3.4) +
  annotate("text", x = x_orci, y = nrow(df_forest) + 1, label = "OR (95% CI)", hjust = 0, size = 3.6, fontface = "bold") +
  annotate("text", x = x_pval, y = nrow(df_forest) + 1, label = "P value", hjust = 0, size = 3.6, fontface = "bold") +
  scale_x_log10() +
  labs(x = "Odds ratio (log scale)", y = NULL) +
  coord_cartesian(
    xlim = c(min(df_forest$LCL, na.rm = TRUE) * 0.8, x_pval * 1.25),
    clip = "off"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    plot.margin = ggplot2::margin(t = 10, r = 140, b = 10, l = 10, unit = "pt")
  )

pdf(file = "2forest.pdf", height = 6.5, width = 15)
print(p_forest_text)
dev.off()
