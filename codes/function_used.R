# GSVA function  --------------------------------------------------------

do_gsva_Poisson <- function(tpm_matrix,gsva_feature_list,annot_matrix,disease_name,control_name){
  
  list_name <- list()
  
  list_name$gsva_score <- gsva(expr = tpm_matrix %>% as.matrix() ,
                               gset.idx.list = gsva_feature_list, method = "gsva",
                               kcdf = "Poisson", parallel.sz = 8L)
  
  list_name$feature_score <- list_name$gsva_score %>% 
    t %>% as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    left_join(annot_matrix %>% dplyr::select(ID, class)) %>% 
    relocate(class, .after = ID) %>% 
    mutate(class = factor(class,levels = c(disease_name,control_name)))
}

do_gsva_Gaussian <- function(tpm_matrix,gsva_feature_list,annot_matrix,disease_name,control_name){
  
  list_name <- list()
  
  list_name$gsva_score <- gsva(expr = tpm_matrix %>% as.matrix() ,
                               gset.idx.list = gsva_feature_list, method = "gsva",
                               kcdf = "Gaussian", parallel.sz = 8L)
  
  list_name$feature_score <- list_name$gsva_score %>% 
    t %>% as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    left_join(annot_matrix %>% dplyr::select(ID, class)) %>% 
    relocate(class, .after = ID) %>% 
    mutate(class = factor(class,levels = c(disease_name,control_name)))
}

# random forest --------------------------------------------------------
do_ml_rf <- function(training_data,testing_data_1,testing_data_2,preprocess_step,disease_name,seed = NULL){
  require(tidyverse)
  require(tidymodels)
  require(doParallel)
  
  if (!is.null(seed)) {
    set.seed(seed)  
  }
  
  cl<-makeCluster(8)
  registerDoParallel(cl)
  
  list_name <- list()
  
  # Choose preprocessing based on input parameter
  if(preprocess_step == "PCA"){
    list_name$recipe <- 
      recipe(class ~ ., data = training_data[,-1]) %>% 
      step_normalize(all_numeric()) %>%
      step_pca(all_numeric(), num_comp = 10) %>% 
      step_nzv(all_predictors()) %>%
      step_zv(all_predictors())
  }else if(preprocess_step == "NOR"){ 
    list_name$recipe <- 
      recipe(class ~ ., data = training_data[,-1]) %>% 
      step_normalize(all_numeric()) %>%
      step_nzv(all_predictors()) %>%
      step_zv(all_predictors())
  } else if (preprocess_step == "NON"){
    list_name$recipe <- 
      recipe(class ~ ., data = training_data[,-1]) %>% 
      step_nzv(all_predictors()) %>%
      step_zv(all_predictors())
  }
  
  # Cross-validation setup
  list_name$cv_folds <- 
    vfold_cv(training_data[,-1],v = 10, strata = class, repeats = 3)
  
  # Random Forest Model setup
  list_name$rf_model <- rand_forest(trees = 1000,
                                    mtry = tune(),
                                    min_n = tune()) %>% 
    set_engine("ranger",importance = "impurity") %>% 
    set_mode("classification")
  
  # Workflow setup
  list_name$rf_workflow <- 
    workflow() %>% 
    add_recipe(list_name$recipe) %>% 
    add_model(list_name$rf_model)
  
  # Grid setup and tuning
  list_name$rf_grid <- 
    grid_max_entropy(extract_parameter_set_dials(list_name$rf_model) %>% 
                       finalize(x = training_data[,-1] %>% dplyr::select(-class)), size = 20)
  
  list_name$rf_tuned <- 
    tune_grid(object = list_name$rf_workflow, 
              resamples =  list_name$cv_folds, 
              grid = list_name$rf_grid, 
              metrics = metric_set(accuracy, bal_accuracy),
              control=tune::control_grid(verbose=TRUE))
  
  # Select the best hyperparameters
  list_name$rf_best_hyper <- list_name$rf_tuned %>%
    select_best(metric = "accuracy")
  
  # Finalize and fit the workflow
  list_name$rf_final_workflow <- 
    list_name$rf_workflow %>% 
    finalize_workflow(list_name$rf_best_hyper)
  
  list_name$rf_fit <- 
    list_name$rf_final_workflow %>% 
    fit(data = training_data[,-1])
  
  # results matrix
  
  list_name$result <- matrix(0,nrow = 1,ncol = 8,dimnames = list(c("rf"),
                                                                 c("cv_accuracy","cv_bal_accuracy",
                                                                   "pre_accuracy_1","pre_bal_accuracy_1","auc_1",
                                                                   "pre_accuracy_2","pre_bal_accuracy_2","auc_2")))
  
  list_name$result[1,1] <- (list_name$rf_tuned %>% collect_metrics() %>% filter(.metric == "accuracy") %>% arrange(-mean) %>% dplyr::select(mean))[1,] %>% pull
  list_name$result[1,2] <- (list_name$rf_tuned %>% collect_metrics() %>% filter(.metric == "bal_accuracy") %>% arrange(-mean) %>% dplyr::select(mean))[1,] %>% pull
  
  # test_data_1
  
  list_name$rf_pred_1 <-  
    predict(list_name$rf_fit, testing_data_1[,-1]) %>% 
    bind_cols(predict(list_name$rf_fit, testing_data_1[,-1], type = "prob")) %>% 
    bind_cols(testing_data_1 %>% dplyr::select(class))
  
  list_name$rf_accuracy_1 <- 
    list_name$rf_pred_1 %>% 
    accuracy(truth = class, .pred_class) %>% 
    dplyr::select(.estimate) %>% pull
  
  list_name$rf_bal_accuracy_1  <- 
    list_name$rf_pred_1 %>% 
    bal_accuracy(truth = class, .pred_class) %>% 
    dplyr::select(.estimate) %>% pull
  
  list_name$rf_auc_1 <- 
    list_name$rf_pred_1 %>% 
    roc_auc(truth = class, paste0(".pred_",disease_name)) %>% 
    dplyr::select(.estimate) %>% pull
  
  list_name$result[1,3] <- list_name$rf_accuracy_1
  list_name$result[1,4] <- list_name$rf_bal_accuracy_1
  list_name$result[1,5] <- list_name$rf_auc_1
  
  # test_data_2
  if(!is.null(testing_data_2)){
    list_name$rf_pred_2 <-  
      predict(list_name$rf_fit, testing_data_2[,-1]) %>% 
      bind_cols(predict(list_name$rf_fit, testing_data_2[,-1], type = "prob")) %>% 
      bind_cols(testing_data_2 %>% dplyr::select(class))
    
    list_name$rf_accuracy_2 <- 
      list_name$rf_pred_2 %>% 
      accuracy(truth = class, .pred_class) %>% 
      dplyr::select(.estimate) %>% pull
    
    list_name$rf_bal_accuracy_2  <- 
      list_name$rf_pred_2 %>% 
      bal_accuracy(truth = class, .pred_class) %>% 
      dplyr::select(.estimate) %>% pull
    
    list_name$rf_auc_2 <- 
      list_name$rf_pred_2 %>% 
      roc_auc(truth = class, paste0(".pred_",disease_name)) %>% 
      dplyr::select(.estimate) %>% pull
    
    list_name$result[1,6] <- list_name$rf_accuracy_2
    list_name$result[1,7] <- list_name$rf_bal_accuracy_2
    list_name$result[1,8] <- list_name$rf_auc_2
  }
  
  return(list_name)
}

# calib_curve  --------------------------------------------------------
do_calib_curve <- function(prediction_res,disease_name){
  output_list <- list()
  calib_data <- prediction_res %>%
    dplyr::select(class, paste0(".pred_",disease_name)) %>%
    dplyr::rename(prob = paste0(".pred_",disease_name)) %>% 
    mutate(observed = as.numeric(class == disease_name))
  
  output_list$brier_score <- mean((calib_data$observed - calib_data$prob)^2)
  
  output_list$curve_data <- calib_data %>%
    mutate(prob_bin = cut(prob, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
    group_by(prob_bin) %>%
    summarise(mean_prob = mean(prob), mean_observed = mean(observed), .groups = 'drop')
  
  return(output_list)
}