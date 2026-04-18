#===============================================================
# 4model.R
# Train final models using repeated CV in training set
# Evaluate once on the full held-out test set
# Export summary metrics, per-sample predictions, and final models
#===============================================================

rm(list = ls())

library(dplyr)
library(caret)
library(pROC)
library(randomForest)
library(gbm)
library(kernlab)

#===============================================================

#===============================================================


setwd("D:/XiangResearch/X2VTE/2vte/save")

df_top5_train    <- read.csv("df_top5_train.csv",    row.names = 1)
df_top5_test     <- read.csv("df_top5_test.csv",     row.names = 1)
df_top10_train   <- read.csv("df_top10_train.csv",   row.names = 1)
df_top10_test    <- read.csv("df_top10_test.csv",    row.names = 1)
df_top13_train   <- read.csv("df_top13_train.csv",   row.names = 1)
df_top13_test    <- read.csv("df_top13_test.csv",    row.names = 1)
df_caprini_train <- read.csv("df_caprini_train.csv", row.names = 1)
df_caprini_test  <- read.csv("df_caprini_test.csv",  row.names = 1)

outcome <- "VTE"
pos <- "Yes"
neg <- "No"

datasets <- list(
  TOP5    = list(train = df_top5_train,    test = df_top5_test),
  TOP10   = list(train = df_top10_train,   test = df_top10_test),
  TOP13   = list(train = df_top13_train,   test = df_top13_test),
  CAPRINI = list(train = df_caprini_train, test = df_caprini_test)
)

algorithms <- list(
  GLM = list(method = "glm",       tuneLength = NULL),
  RDF = list(method = "rf",        tuneLength = 5),
  GBA = list(method = "gbm",       tuneLength = 5),
  SVM = list(method = "svmRadial", tuneLength = 5)
)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

#===============================================================

#===============================================================


model_list <- list()
pred_list <- list()
metric_list <- list()

for (ds in names(datasets)) {

  cat("\nDataset:", ds, "\n")

  train_df <- datasets[[ds]]$train
  test_df  <- datasets[[ds]]$test

  train_df[[outcome]] <- factor(train_df[[outcome]], levels = c(neg, pos))
  test_df[[outcome]]  <- factor(test_df[[outcome]], levels = c(neg, pos))

  for (alg in names(algorithms)) {

    cat("  Algorithm:", alg, "\n")

    method_name <- algorithms[[alg]]$method
    tl <- algorithms[[alg]]$tuneLength

    set.seed(5)

    if (is.null(tl)) {
      fit <- train(
        as.formula(paste0(outcome, " ~ .")),
        data = train_df,
        method = method_name,
        trControl = ctrl,
        metric = "ROC"
      )
    } else {
      fit <- train(
        as.formula(paste0(outcome, " ~ .")),
        data = train_df,
        method = method_name,
        trControl = ctrl,
        metric = "ROC",
        tuneLength = tl
      )
    }
    
    model_name <- paste0(alg, "_", ds, "_model")
    model_list[[model_name]] <- fit
    
    prob_yes <- predict(fit, newdata = test_df, type = "prob")[, pos]
    pred_cls <- predict(fit, newdata = test_df, type = "raw")
    
    cm <- confusionMatrix(pred_cls, test_df[[outcome]], positive = pos)
    
    roc_obj <- roc(
      response = test_df[[outcome]],
      predictor = prob_yes,
      levels = c(neg, pos),
      direction = "<",
      quiet = TRUE
    )
    
    auc_val <- as.numeric(auc(roc_obj))
    auc_ci  <- as.numeric(ci.auc(roc_obj))
    
    best_cut <- coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    )
    
    if (is.data.frame(best_cut)) {
      cutoff    <- as.numeric(best_cut$threshold[1])
      sens_cut  <- as.numeric(best_cut$sensitivity[1])
      spec_cut  <- as.numeric(best_cut$specificity[1])
    } else if (is.list(best_cut) && !is.null(best_cut$threshold)) {
      cutoff    <- as.numeric(best_cut$threshold[[1]])
      sens_cut  <- as.numeric(best_cut$sensitivity[[1]])
      spec_cut  <- as.numeric(best_cut$specificity[[1]])
    } else {
      cutoff    <- as.numeric(best_cut["threshold"])
      sens_cut  <- as.numeric(best_cut["sensitivity"])
      spec_cut  <- as.numeric(best_cut["specificity"])
    }
    
    cat(
      sprintf(
        "    Youden cutoff = %.3f | Sensitivity = %.3f | Specificity = %.3f\n",
        cutoff, sens_cut, spec_cut
      )
    )
    
    pred_youden <- ifelse(prob_yes >= cutoff, pos, neg)
    pred_youden <- factor(pred_youden, levels = c(neg, pos))
    cm_youden <- confusionMatrix(pred_youden, test_df[[outcome]], positive = pos)
    
    
    y_num <- ifelse(test_df[[outcome]] == pos, 1, 0)
    brier <- mean((prob_yes - y_num)^2)
    
    eps <- 1e-6
    p_clip <- pmin(pmax(prob_yes, eps), 1 - eps)
    logit_p <- qlogis(p_clip)
    
    cal_intercept_fit <- glm(y_num ~ offset(logit_p), family = binomial)
    cal_slope_fit     <- glm(y_num ~ logit_p, family = binomial)
    
    cal_intercept <- coef(cal_intercept_fit)[1]
    cal_slope     <- coef(cal_slope_fit)[2]
    
    pred_list[[length(pred_list) + 1]] <- data.frame(
      Dataset = ds,
      Algorithm = alg,
      row_id = rownames(test_df),
      truth = as.character(test_df[[outcome]]),
      prob_yes = prob_yes,
      pred_class_default = as.character(pred_cls),
      pred_class_youden = as.character(pred_youden),
      threshold_youden = cutoff,
      stringsAsFactors = FALSE
    )
    
    metric_list[[length(metric_list) + 1]] <- data.frame(
      Dataset = ds,
      Algorithm = alg,
      AUC = auc_val,
      AUC_LCL = auc_ci[1],
      AUC_UCL = auc_ci[3],
      Accuracy_default = unname(cm$overall["Accuracy"]),
      Sensitivity_default = unname(cm$byClass["Sensitivity"]),
      Specificity_default = unname(cm$byClass["Specificity"]),
      Precision_default = unname(cm$byClass["Pos Pred Value"]),
      F1_default = unname(cm$byClass["F1"]),
      Threshold_Youden = cutoff,
      Accuracy_Youden = unname(cm_youden$overall["Accuracy"]),
      Sensitivity_Youden = sens_cut,
      Specificity_Youden = spec_cut,
      Precision_Youden = unname(cm_youden$byClass["Pos Pred Value"]),
      F1_Youden = unname(cm_youden$byClass["F1"]),
      Brier = brier,
      Calibration_Intercept = unname(cal_intercept),
      Calibration_Slope = unname(cal_slope),
      stringsAsFactors = FALSE
    )
  }
}

#===============================================================

#===============================================================


pred_long <- bind_rows(pred_list)
metrics_table <- bind_rows(metric_list)

setwd("D:/XiangResearch/X2VTE/2vte")

write.csv(metrics_table, "4result_full_test.csv", row.names = FALSE)
write.csv(pred_long, "4pred_long.csv", row.names = FALSE)
write.csv(metrics_table, "4metrics.csv", row.names = FALSE)

threshold_table <- metrics_table %>%
  select(
    Dataset, Algorithm,
    Threshold_Youden,
    Sensitivity_Youden,
    Specificity_Youden,
    Accuracy_Youden
  ) %>%
  arrange(Dataset, Algorithm)

write.csv(threshold_table, "4thresholds.csv", row.names = FALSE)


for (nm in names(model_list)) {
  saveRDS(model_list[[nm]], file.path("save", paste0(nm, ".rds")))
}
