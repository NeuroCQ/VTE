#===============================================================
# 1table.R
#===============================================================

rm(list = ls())

library(dplyr)
library(caret)
library(gtsummary)
library(officer)
library(flextable)

setwd("D:/XiangResearch/X2VTE/1data")
data <- read.csv("df_PND.csv", header = TRUE, row.names = 1)

data$PND <- data$postop_functional_complication
data$postop_functional_complication <- NULL
data$OS <- NULL
data$OSTime <- NULL
data$vte_class <- NULL
data$VTE_clinical <- NULL

data$PND <- factor(data$PND, levels = c("No", "Yes"))
data_model <- data %>% dplyr::relocate(PND)

data_model <- data_model %>%
  mutate(across(where(is.character), as.factor))

#===============================================================

#===============================================================


is_binary_vec <- function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) == 2
}

convert_binary_cols <- function(df,
                                exclude = character(0),
                                label_01 = TRUE,
                                yes_labels = c("1", "Yes", "YES", "yes", "Y", "True", "TRUE", "true"),
                                no_labels  = c("0", "No", "NO", "no", "N", "False", "FALSE", "false")) {
  cols <- setdiff(names(df), exclude)
  
  for (v in cols) {
    x <- df[[v]]
    if (is.list(x)) next
    if (!is_binary_vec(x)) next
    
    ux <- unique(x[!is.na(x)])
    ux_chr <- as.character(ux)
    
    if (label_01 && all(ux_chr %in% c("0", "1"))) {
      df[[v]] <- factor(as.character(x), levels = c("0", "1"), labels = c("No", "Yes"))
      next
    }
    
    if (all(ux_chr %in% c(yes_labels, no_labels))) {
      mapped <- ifelse(as.character(x) %in% yes_labels, "Yes",
                       ifelse(as.character(x) %in% no_labels, "No", NA))
      df[[v]] <- factor(mapped, levels = c("No", "Yes"))
      next
    }
    
    levs <- sort(unique(as.character(x[!is.na(x)])))
    df[[v]] <- factor(as.character(x), levels = levs)
  }
  
  df
}

#===============================================================

#===============================================================


set.seed(8)
train_index <- createDataPartition(y = data_model$PND, p = 0.8, list = FALSE)
train_data <- data_model[train_index, , drop = FALSE]
test_data  <- data_model[-train_index, , drop = FALSE]

data_model$dataset <- ifelse(
  rownames(data_model) %in% rownames(train_data),
  "Train",
  "Test"
) %>% factor(levels = c("Train", "Test"))

data_model <- convert_binary_cols(
  data_model,
  exclude = c("PND", "dataset"),
  label_01 = TRUE
)

data_set_samples <- data_model %>%
  dplyr::select(PND, dataset)

setwd("D:/XiangResearch/X2VTE/3pnd")
save(data_set_samples, file = "data_set_samples.Rdata")

#===============================================================

#===============================================================


feature_labels <- c(
  age = "Age (per year)",
  sex = "Sex",
  bmi = "BMI (kg/m^2)",
  smoking_history = "Smoking history",
  alcohol_history = "Alcohol history",
  hypertension_history = "Hypertension history",
  diabetes_history = "Diabetes history",
  hyperlipidemia_history = "Hyperlipidemia history",
  coronary_heart_disease_history = "Coronary heart disease history",
  stroke_history = "Stroke history",
  prior_vte_history = "Prior VTE history",
  pre_anticoagulant_use = "Preop anticoagulant use",
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
  recurrent_glioma = "Recurrent glioma"
)

tbl_train_test <- data_model %>%
  tbl_summary(
    by = dataset,
    label = as.list(
      feature_labels[names(feature_labels) %in% names(data_model)]
    ),
    type = list(all_categorical() ~ "categorical",
                pre_caprini_score ~ "continuous"),
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    percent = "column",
    missing = "no"
  ) %>%
  add_p(test = list(pre_caprini_score ~ "wilcox.test")) %>%
  bold_labels()

#===============================================================

#===============================================================


ft <- as_flex_table(tbl_train_test)
doc <- read_docx() %>%
  body_add_par("Table 1. Baseline characteristics (Train vs Test)",
               style = "heading 1") %>%
  body_add_flextable(ft)

print(doc, target = "1Table.docx")
