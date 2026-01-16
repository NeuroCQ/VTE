#===============================================================

#===============================================================
rm(list = ls())
library(dplyr)
library(caret)
library(pROC)
library(randomForest)
library(gbm)
library(kernlab)

df_top5_train    <- read.csv("df_top5_train.csv",    row.names = 1)
df_top5_test     <- read.csv("df_top5_test.csv",     row.names = 1)
df_top10_train   <- read.csv("df_top10_train.csv",   row.names = 1)
df_top10_test    <- read.csv("df_top10_test.csv",    row.names = 1)
df_top15_train   <- read.csv("df_top15_train.csv",   row.names = 1)
df_top15_test    <- read.csv("df_top15_test.csv",    row.names = 1)
df_top20_train   <- read.csv("df_top20_train.csv",   row.names = 1)
df_top20_test    <- read.csv("df_top20_test.csv",    row.names = 1)
df_caprini_train <- read.csv("df_caprini_train.csv", row.names = 1)
df_caprini_test  <- read.csv("df_caprini_test.csv",  row.names = 1)

#===============================================================

#===============================================================

ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)


datasets <- list(
  TOP5    = list(train = df_top5_train,    test = df_top5_test),
  TOP10   = list(train = df_top10_train,   test = df_top10_test),
  TOP15   = list(train = df_top15_train,   test = df_top15_test),
  TOP20   = list(train = df_top20_train,   test = df_top20_test),
  CAPRINI = list(train = df_caprini_train, test = df_caprini_test)
)

outcome <- "VTE"
pos <- "Yes"
neg <- "No"

algorithms <- list(
  RDF  = list(method = "rf",        tuneLength = 5),
  GLM = list(method = "glm",       tuneLength = NULL),
  GBA = list(method = "gbm",       tuneLength = 5),
  SVM = list(method = "svmRadial", tuneLength = 5)
)

#===============================================================

#===============================================================

models <- list()
metrics_list <- list()

set.seed(5)

n_splits <- 20
test_frac <- 0.8

for (ds in names(datasets)) {
  
  cat("\nDataset:", ds, "\n")
  
  train_df <- datasets[[ds]]$train
  test_df  <- datasets[[ds]]$test
  
  train_df[[outcome]] <- factor(train_df[[outcome]], levels = c(neg, pos))
  test_df[[outcome]]  <- factor(test_df[[outcome]],  levels = c(neg, pos))
  
  for (alg in names(algorithms)) {
    
    cat("  Algorithm:", alg, "\n")
    
    method_name <- algorithms[[alg]]$method
    tl <- algorithms[[alg]]$tuneLength

    if (is.null(tl)) {
      fit <- train(as.formula(paste0(outcome, " ~ .")),
                   data = train_df,
                   method = method_name,
                   trControl = ctrl,
                   metric = "ROC")
    } else {
      fit <- train(as.formula(paste0(outcome, " ~ .")),
                   data = train_df,
                   method = method_name,
                   trControl = ctrl,
                   metric = "ROC",
                   tuneLength = tl)
    }
    
    models[[paste(ds, alg, sep = "_")]] <- fit
    
    n_sub <- floor(nrow(test_df) * test_frac)
    
    for (split_id in 1:n_splits) {
      
      idx <- sample(seq_len(nrow(test_df)), size = n_sub, replace = FALSE)
      test_sub <- test_df[idx, , drop = FALSE]
      test_sub[[outcome]] <- factor(test_sub[[outcome]], levels = c(neg, pos))
      

      prob_yes <- predict(fit, newdata = test_sub, type = "prob")[, pos]
      pred_cls <- predict(fit, newdata = test_sub, type = "raw")
      
      cm <- confusionMatrix(pred_cls, test_sub[[outcome]], positive = pos)
      
  
      auc_val <- tryCatch({
        roc_obj <- pROC::roc(test_sub[[outcome]], prob_yes,
                             levels = c(neg, pos), direction = "<", quiet = TRUE)
        as.numeric(pROC::auc(roc_obj))
      }, error = function(e) NA_real_)
      
      metrics_list[[length(metrics_list) + 1]] <- data.frame(
        Dataset = ds,
        Algorithm = alg,
        split_id = split_id,
        AUC = auc_val,
        Accuracy = unname(cm$overall["Accuracy"]),
        Sensitivity = unname(cm$byClass["Sensitivity"]),
        Specificity = unname(cm$byClass["Specificity"]),
        row.names = NULL
      )
    }
  }
}

metrics_table <- do.call(rbind, metrics_list)
print(metrics_table)

getwd()
                          
write.csv(metrics_table, file = '3result.csv')

