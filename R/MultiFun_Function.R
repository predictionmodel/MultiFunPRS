
# Part 1: Construct prediction model with Functional annotations and Omics summary statistics####

#' This function trains a machine learning model using provided labels (`Y_Data`) and features (`X_Data),
#' and predicts the probability of each SNP being causal.
#'
#' @param Y_Data A data frame containing summary statistics with at least columns: POS, Index (Y).
#' @param X_Data A data frame containing annotation features with at least column: POS.
#' @param Prediction_Model A string indicating the prediction algorithm to use (e.g., "randomforest").
#' @return A list containing:
#' \item{Probability}{Predicted probability for each SNP.}
#' \item{ROC}{ROC curve object from pROC::roc().}
#' @importFrom dplyr select mutate filter bind_cols
#' @importFrom stats na.pass
#' @export
PriorCausalProbability <- function(Y_Data, X_Data, Prediction_Model) {
  library(pROC)
  
  # Add label Y based on Index
  Y_Data$Y <- Y_Data$Index
  
  # Merge features and labels by SNP position
  data <- merge(X_Data, as.data.frame(Y_Data[, c("POS", "Y")]), by = "POS", all.x = TRUE)
  
  # Remove unnecessary columns
  cols_to_remove <- c("POS", "CHR", "SNP", "CM")
  testdata <- data[, !(colnames(data) %in% cols_to_remove), drop = FALSE]
  traindata <- data[!is.na(data$Y), ]
  traindata <- traindata[, !(colnames(traindata) %in% cols_to_remove), drop = FALSE]
  
  # Balance the dataset
  traindata <- balance_samples(traindata)
  
  # Convert Y to factor
  traindata$Y <- as.factor(traindata$Y)
  testdata$Y <- as.factor(testdata$Y)
  
  # Train prediction model and obtain predictive probability based on specified algorithm
  PredictModel <- switch(tolower(Prediction_Model),
                         logisticregression = Logistic,
                         randomforest = RandomForest,
                         svm = SVM,
                         gradientboosing = GBM,
                         adaboost = AdaBoost)
  
  predict_result <- PredictModel(Y, traindata, testdata)
  names(predict_result) <- c("Probability", "ROC")
  return(predict_result)
}

#' Train a Random Forest classifier and predict causal SNP probabilities.
#'
#' @param Y Response variable name (ignored but kept for interface consistency).
#' @param traindata Training data including the response variable "Y".
#' @param testdata Test data without the response variable.
#' @return A list containing:
#' \item{Probability}{Predicted probability of being causal (class 1).}
#' \item{ROC}{ROC curve object from pROC::roc().}
#' @import randomForest
#' @importFrom pROC roc
#' @keywords internal
#' @export
RandomForest <- function(Y, traindata, testdata) {
  library(randomForest)
  library(caret)
  
  set.seed(20230204)
  RandomForest<-randomForest(Y~.,data=traindata,importance=TRUE,na.action = na.pass)
  
  testdata2<-testdata[,-which(colnames(testdata)=="Y")]
  RandomForest_Probability<-predict(RandomForest,testdata2,type="prob")[,2]
  RandomForest_roc<-roc(testdata$Y,RandomForest_Probability)
  list<-list(RandomForest_Probability,RandomForest_roc)
  return(list)
}

#' Helper function to balance positive and negative samples in the training dataset.
#'
#' @param df A data frame containing a binary response variable "Y" (0 or 1).
#' @param pos_limit Maximum number of positive samples allowed.
#' @param neg_ratio Ratio of negative to positive samples after downsampling.
#' @return A balanced data frame with reduced sample size.
#' @keywords internal
#' @export

balance_samples <- function(df, pos_limit = 1000, neg_ratio = 4) {
  library(dplyr)
  
  df_pos <- df %>% filter(Y == 1)
  df_neg <- df %>% filter(Y == 0)
  
  if (nrow(df_pos) > pos_limit) {
    set.seed(20230204)
    df_pos <- slice_sample(df_pos, n = pos_limit)
  }
  
  if (nrow(df_neg) > nrow(df_pos) * neg_ratio) {
    set.seed(20230204)
    df_neg <- slice_sample(df_neg, n = nrow(df_pos) * neg_ratio)
  }
  
  bind_rows(df_pos, df_neg)
}


###########################################################################################
# Part 2: Compute weights for multiple omics layers using Multiple Kernel Learning (MKL)####
#'
#' This function integrates multiple omics summary datasets and computes weights
#'
#' @param Multi_Y_Data A list of data frames, each representing one omics layer's summary statistics.
#' @param X_Data A data frame containing annotation features with at least column: POS.
#' @return A numeric vector of weights corresponding to each omics layer.
#' @import RMKL
#' @import mixKernel
#' @import dplyr
#' @import data.table
#' @export
MultiOmics_Weight_MKL <- function(Multi_Y_Data, X_Data) {
  library(RMKL)
  library(mixKernel) 
  library(dplyr)
  library(data.table)
  
  Traindata_MultiOmics <- list()
  ImportanceVar_MultiOmics <- list()
  
  for (i in seq_along(Multi_Y_Data)) {
    summarystatistic <- Multi_Y_Data[[i]]
    #Y:Define whether a causal SNP is a Y,which is by summary statistics
    summarystatistic$Y <- summarystatistic$Index
    
    data <- left_join(X_Data, summarystatistic[, c("POS", "Y")], by = "POS")
    
    traindata <- data %>%
      filter(!is.na(Y)) %>%
      balance_samples()
    
    Traindata_MultiOmics[[i]] <- traindata
    traindata <- traindata[,-which(colnames(traindata) %in% c("POS","CHR","SNP","CM"))]
    ImportanceVar_MultiOmics[[i]] <- RandomForest_ImportanceVar(Y, traindata)
  }
  
  # Further sample balancing for the combined dataset
  traindata_all <- NULL
  for (i in seq_along(Multi_Y_Data)) {
    temp <- as.data.frame(Traindata_MultiOmics[[i]]) %>%
      filter(Y %in% c(0, 1)) %>%
      balance_samples(pos_limit = round(1000 / (2 * length(Multi_Y_Data))), neg_ratio = 1)
    
    traindata_all <- bind_rows(traindata_all, temp)
  }
  
  traindata_all <- distinct(traindata_all, POS, .keep_all = TRUE)
  
  # Convert Y labels from 0 to -1
  traindata_all <- mutate(traindata_all, Y = ifelse(Y == 0, -1, Y))
  
  # Compute kernel matrices and calculate weights
  SNPMatrix <- list()
  for (i in seq_along(Multi_Y_Data)) {
    temp <- traindata_all %>%
      dplyr::select(all_of(ImportanceVar_MultiOmics[[i]])) %>%
      dplyr::select_if(~ !is.character(.))
    
    temp <- temp[, colSums(is.na(temp) | temp == 0 | temp == "") < nrow(temp)]
    
    kernel <- compute.kernel(temp, kernel.func = "linear")
    SNPMatrix[[i]] <- kernel$kernel
  }
  
  C <- 1
  Weight <- SEMKL.classification(SNPMatrix, traindata_all$Y, C)$gamma
  return(Weight)
}

#' Extract important variables from Random Forest model.
#'
#' @param Y Response variable name (ignored but kept for interface consistency).
#' @param traindata Training data including the response variable "Y".
#' @return Character vector of selected feature names.
#' @import randomForest
#' @importFrom ggplot2 ggplot geom_line labs theme element_blank element_rect
#' @keywords internal
#' @export
RandomForest_ImportanceVar <- function(Y, traindata) {
  set.seed(20230204)
  traindata$Y <- as.factor(traindata$Y)
  rf <- randomForest(Y ~ ., data = traindata, importance = TRUE, na.action = na.pass)
  
  importance_df <- data.frame(importance(rf), check.names = FALSE)
  importance_df <- importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = TRUE), ]
  
  Select_Number <- 100
  selected_vars <- rownames(importance_df[1:Select_Number, ])
  
  return(selected_vars)
}

#############################################################################
# Part 3: MultiFun: Integrate multi-omics predictions into a unified causal SNP probability score using weighted fusion####
#'
#' @param MultiOmicsSummary A list of data frames, each representing one omics layer's summary statistics.
#' @param PredictionData Annotation data containing SNP positions and features.
#' @param PredictionModel Machine learning method used for prior probability estimation.
#' @return A list containing:
#' \item{MultiOmics_Weight}{Numeric vector of weights for each omics layer.}
#' \item{MultiOmics_PriorProbability}{Data frame of predicted probabilities per omics and combined scores.}
#' @importFrom dplyr select
#' @export
MultiFunR <- function(MultiOmicsSummary, PredictionData, PredictionModel) {
  
  # Step 1: Estimate prior probabilities for each omics layer
  PriorProbability <- list()
  for (i in seq_along(MultiOmicsSummary)) {
    PriorProbability[[i]] <- PriorCausalProbability(MultiOmicsSummary[[i]], PredictionData, PredictionModel)
  }
  
  # Combine prior probabilities into a single data frame
  Omics_PriorProbability <- PredictionData[, c(1:4), drop = FALSE]
  for (i in seq_along(MultiOmicsSummary)) {
    omics_name <- names(MultiOmicsSummary)[i]
    Omics_PriorProbability[[omics_name]] <- PriorProbability[[i]][["Probability"]]
    cat(paste0(omics_name, "'s AUC: ", PriorProbability[[i]][["ROC"]]$auc, "\n"))
  }
  
  # Step 2: Estimate weights using MKL
  MultiOmics_Weight <- MultiOmics_Weight_MKL(MultiOmicsSummary, PredictionData)
  
  # Step 3: Apply weights to calculate weighted average (MKL-based)
  Omics_PriorProbability$PriorCausalProbability_MKL <- 0
  for (i in seq_along(MultiOmicsSummary)) {
    omics_name <- names(MultiOmicsSummary)[i]
    Omics_PriorProbability$PriorCausalProbability_MKL <-
      Omics_PriorProbability$PriorCausalProbability_MKL +
      MultiOmics_Weight[i] * Omics_PriorProbability[[omics_name]]
  }
  
  # Step 4: Use MAX across omics as an alternative integration method
  Omics_PriorProbability$PriorCausalProbability_MAX <-
    apply(Omics_PriorProbability[names(MultiOmicsSummary)], 1, max)
  
  # Return results
  list(
    MultiOmics_Weight = MultiOmics_Weight,
    MultiOmics_PriorProbability = Omics_PriorProbability
  )
}