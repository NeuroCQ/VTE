#===============================================================
# 3prepare.R
#===============================================================

rm(list = ls())

library(dplyr)
library(stringr)
library(caret)
library(randomForest)
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
data <- data %>% dplyr::relocate(PND)


#===============================================================
# 2. Variable preprocessing
#===============================================================

convert_binary_column <- function(x, pos = "Yes", neg = "No") {
  x <- as.character(x)
  x[x == pos] <- 1
  x[x == neg] <- 0
  x[!(x %in% c("0", "1"))] <- NA
  as.integer(x)
}

binary_cols <- c(
  "smoking_history",
  "alcohol_history",
  "hypertension_history",
  "diabetes_history",
  "hyperlipidemia_history",
  "coronary_heart_disease_history",
  "stroke_history",
  "prior_vte_history",
  "pre_anticoagulant_use",
  "recurrent_glioma"
)

binary_cols <- intersect(binary_cols, names(data))
cat("Columns to convert:\n")
print(binary_cols)

data[binary_cols] <- lapply(data[binary_cols], convert_binary_column)

cat("\nSummary after conversion:\n")
for (col in binary_cols) {
  cat(sprintf("\n%s: ", col))
  print(table(data[[col]], useNA = "always"))
}

#---------------------------------------------------------------
# sex
#---------------------------------------------------------------
table(data$sex)
data$sex[data$sex == "Female"] <- 1
data$sex[data$sex == "Male"] <- 0
data$sex <- as.numeric(as.character(data$sex))
table(data$sex)

#---------------------------------------------------------------
# tumor_spread
#---------------------------------------------------------------
table(data$tumor_spread)
data$tumor_spread[data$tumor_spread == "localized"] <- 0
data$tumor_spread[data$tumor_spread == "regional"] <- 1
data$tumor_spread <- as.numeric(as.character(data$tumor_spread))
table(data$tumor_spread)

#---------------------------------------------------------------
# tumor_laterality
#---------------------------------------------------------------
table(data$tumor_laterality)
data$tumor_laterality[data$tumor_laterality == "Left"] <- 1
data$tumor_laterality[data$tumor_laterality == "Right"] <- 2
data$tumor_laterality[data$tumor_laterality == "Bilateral"] <- 3
data$tumor_laterality <- as.numeric(as.character(data$tumor_laterality))
table(data$tumor_laterality)

#---------------------------------------------------------------
# tumor_location
#---------------------------------------------------------------
table(data$tumor_location)
data$tumor_location[data$tumor_location == "Frontal"] <- 1
data$tumor_location[data$tumor_location == "Deep/Other"] <- 2
data$tumor_location[data$tumor_location == "Non-frontal lobar"] <- 3
data$tumor_location[data$tumor_location == "Ventricle"] <- 4
data$tumor_location <- as.numeric(as.character(data$tumor_location))
table(data$tumor_location)

str(data)

#===============================================================
# 3. Keep selected variables for modeling
#===============================================================

setwd("D:/XiangResearch/X2VTE/3pnd")
dir.create("save", showWarnings = FALSE, recursive = TRUE)

sele_variables <- read.csv("2Uni.csv", header = TRUE, row.names = 1)
sele_variables <- c("PND", sele_variables$Variable)

data_model <- data[, which(colnames(data) %in% sele_variables)]

#===============================================================
# 4. Convert character to factor if any
#===============================================================

data_model <- data_model %>%
  mutate(across(where(is.character), as.factor))

#===============================================================
# 5. Standardization on full dataset
#===============================================================

outcome <- "PND"

binary_vars_no_scale <- c(
  "sex",
  "smoking_history",
  "alcohol_history",
  "hypertension_history",
  "diabetes_history",
  "hyperlipidemia_history",
  "coronary_heart_disease_history",
  "stroke_history",
  "prior_vte_history",
  "pre_anticoagulant_use",
  "tumor_spread",
  "recurrent_glioma"
)

binary_vars_no_scale <- intersect(binary_vars_no_scale, names(data_model))

numeric_vars <- names(data_model)[sapply(data_model, is.numeric)]
continuous_vars <- setdiff(numeric_vars, c(outcome, binary_vars_no_scale))

cat("\nContinuous variables to standardize:\n")
print(continuous_vars)

if (length(continuous_vars) > 0) {
  pp_all <- caret::preProcess(
    data_model[, continuous_vars, drop = FALSE],
    method = c("center", "scale")
  )
  
  data_model[, continuous_vars] <- predict(
    pp_all,
    data_model[, continuous_vars, drop = FALSE]
  )
  
  saveRDS(pp_all, file = file.path("save", "scaler_all_data.rds"))
}

cat("\nData structure after standardization:\n")
str(data_model)

#===============================================================
# 6. Train / test split
#===============================================================

load("data_set_samples.Rdata")

trainSamples <- subset(data_set_samples, dataset == "Train")
trainSamples <- rownames(trainSamples)

testSamples <- subset(data_set_samples, dataset == "Test")
testSamples <- rownames(testSamples)

train_data <- data_model[which(rownames(data_model) %in% trainSamples), ]
test_data  <- data_model[which(rownames(data_model) %in% testSamples), ]

cat("\nTrain/Test dimensions:\n")
print(dim(train_data))
print(dim(test_data))

#===============================================================
# 7. Random forest for feature importance
#===============================================================

ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(8)

rf_model <- train(
  PND ~ .,
  data = train_data,
  method = "rf",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  importance = TRUE
)

#===============================================================
# 8. Variable importance plot
#===============================================================


varImp_rf <- varImp(rf_model, scale = TRUE)
imp <- varImp_rf$importance

if (all(c("No", "Yes") %in% colnames(imp))) {
  imp$Overall <- rowMeans(imp[, c("No", "Yes")], na.rm = TRUE)
}

imp <- imp[order(imp$Overall, decreasing = TRUE), , drop = FALSE]
imp$Variable <- rownames(imp)

p <- ggplot(imp, aes(x = reorder(Variable, Overall), y = Overall)) +
  geom_bar(stat = "identity", width = 0.7, fill = "grey30") +
  coord_flip() +
  labs(x = NULL, y = "Variable importance") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8)
  )

pdf(file = "3importance.pdf", height = 5, width = 5)
print(p)
dev.off()


#===============================================================
# 9. Select top variables
#===============================================================


feat_raw <- rownames(imp)

feat_base <- feat_raw %>%
  stringr::str_replace("^tumor_location.*", "tumor_location") %>%
  stringr::str_replace("^tumor_laterality.*", "tumor_laterality") %>%
  stringr::str_replace("^tumor_spread.*", "tumor_spread")

feat_base <- feat_base[!duplicated(feat_base)]
feat_base <- intersect(feat_base, names(data_model))

top5_vars  <- feat_base[1:5]
top10_vars <- feat_base[1:10]

cat("\nTop 5 variables:\n")
print(top5_vars)
cat("\nTop 10 variables:\n")
print(top10_vars)

df_top5_train  <- train_data %>% dplyr::select(PND, all_of(top5_vars))
df_top5_test   <- test_data  %>% dplyr::select(PND, all_of(top5_vars))
df_top10_train <- train_data %>% dplyr::select(PND, all_of(top10_vars))
df_top10_test  <- test_data  %>% dplyr::select(PND, all_of(top10_vars))

#===============================================================
# 10. Build train/test datasets
#===============================================================

setwd("D:/XiangResearch/X2VTE/3pnd/save")

write.csv(df_top5_train,  file = "df_top5_train.csv", row.names = TRUE)
write.csv(df_top5_test,   file = "df_top5_test.csv", row.names = TRUE)
write.csv(df_top10_train, file = "df_top10_train.csv", row.names = TRUE)
write.csv(df_top10_test,  file = "df_top10_test.csv", row.names = TRUE)

