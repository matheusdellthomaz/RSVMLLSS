n_cores <- detectCores() - 2
cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)
calculate_metrics <- function(pred_df) {
  if (nrow(pred_df) == 0) return(NULL)
  
  pred_df %>%
    mutate(
      biais_brut = .pred - AUCt,
      bias_rel = (.pred - AUCt) / AUCt,
      bias_rel_square = bias_rel^2,
      biais_brut_sqr = biais_brut^2
    ) %>%
    summarise(
      biais_brut = mean(biais_brut),
      rmse = sqrt(mean(biais_brut_sqr)),
      biais_rel = mean(bias_rel),
      relative_rmse = sqrt(mean(bias_rel_square)),
      biais_out_20percent = mean((bias_rel) > 0.2),
      nb_out_20percent = sum((bias_rel) > 0.2),
      n = n()
    ) %>%
    mutate(
      R2 = summary(lm(AUCt ~ .pred, data = pred_df))$r.squared
    )
}
test_combinations <- function(data, time_predictors, fixed_predictors = NULL, response,
                              external_dfs = list(), combo_sizes = c(2, 3),
                              include_fixed_predictors = TRUE,
                              include_ratio_predictors = FALSE, 
                              include_diff_predictors = FALSE,
                              algorithms = c("rf", "xgb", "svm", "glmnet"),
                              return_fits = FALSE,
                              return_workflows = FALSE) {
  
  # Lists
  all_results <- list()
  all_best_params <- list()
  all_tabs <- list()
  all_fits <- list()
  all_workflows <- list()
  
  # Control
  ctrl <- control_resamples(save_pred = TRUE, save_workflow = TRUE)
  race_ctrl <- control_race(save_pred = TRUE)
  
  # Algorithms
  valid_algs <- c("rf", "xgb", "svm", "glmnet")
  algorithms <- intersect(algorithms, valid_algs)
  
  # Derived Predictors
  create_derived_predictors <- function(df, predictors, 
                                        include_ratio, include_diff) {
    new_df <- df
    ratio_names <- c()
    diff_names <- c()
    
    if (include_ratio && length(predictors) >= 2) {
      ratios <- combn(predictors, 2, simplify = FALSE) %>%
        map(function(pair) {
          name <- paste0(pair[2], "_div_", pair[1])
          new_df <<- new_df %>%
            mutate(!!name := ifelse(!!sym(pair[1]) == 0, NA, !!sym(pair[2]) / !!sym(pair[1])))
          name
        })
      ratio_names <- unlist(ratios)
    }
    
    if (include_diff && length(predictors) >= 2) {
      diffs <- combn(predictors, 2, simplify = FALSE) %>%
        map(function(pair) {
          name <- paste0(pair[2], "_sub_", pair[1])
          new_df <<- new_df %>%
            mutate(!!name := !!sym(pair[2]) - !!sym(pair[1]))
          name
        })
      diff_names <- unlist(diffs)
    }
    
    list(data = new_df, ratio_names = ratio_names, diff_names = diff_names)
  }
  
  # Algs specs
  model_specs <- list(
    rf = rand_forest(
      mtry = tune(),
      min_n = tune(),
      trees = tune()
    ) %>% 
      set_engine("ranger", nthread = 10) %>% 
      set_mode("regression"),
    
    xgb = boost_tree(
      mode = "regression",
      mtry = tune(),
      trees = tune(),
      min_n = tune(),
      sample_size = tune(),
      tree_depth = tune(),
      learn_rate = tune()
    ) %>% 
      set_engine("xgboost", nthread = 10),
    
    svm = svm_linear(
      mode = "regression",
      cost = tune(),
      margin = tune()
    ) %>% 
      set_engine("kernlab", nthread = 10),
    
    glmnet = linear_reg(
      penalty = tune(),
      mixture = tune()
    ) %>% 
      set_engine("glmnet", nthread = 10)
  )
  
  # Combination process
  for (n in combo_sizes) {
    size_results <- list()
    size_best_params <- list()
    size_tabs <- list()
    size_fits <- list()
    size_workflows <- list()
    
    # Time combinations
    time_combos <- combn(time_predictors, n, simplify = FALSE)
    
    for (combo in time_combos) {
      # Comb names
      combo_name <- paste(combo, collapse = "_")
      suffix <- ""
      
      # Suffix
      if (include_fixed_predictors && length(fixed_predictors) > 0) {
        suffix <- paste0(suffix, "_fixed(", paste(fixed_predictors, collapse = "+"), ")")
      }
      
      derived <- create_derived_predictors(
        data, combo, include_ratio_predictors, include_diff_predictors
      )
      main_data <- derived$data
      all_predictors <- c(combo, 
                          if(include_fixed_predictors) fixed_predictors,
                          derived$ratio_names, 
                          derived$diff_names)
      
      # Derived predictor suffix
      if (include_ratio_predictors && length(derived$ratio_names) > 0) {
        suffix <- paste0(suffix, "_ratio")
      }
      if (include_diff_predictors && length(derived$diff_names) > 0) {
        suffix <- paste0(suffix, "_diff")
      }
      
      combo_base_name <- paste0(combo_name, suffix)
      cat("\n--- Processando combinação:", combo_base_name, "---\n")
      
      # Data
      df_subset <- main_data %>% 
        select(all_of(c(all_predictors, response))) %>%
        na.omit()
      
      # Train/Test split
      set.seed(123)
      split <- initial_split(df_subset, strata = !!sym(response), prop = 0.75)
      train <- training(split)
      test <- testing(split)
      
      # Recipe
      rec <- recipe(as.formula(paste(response, "~ .")), data = train) %>% 
        step_normalize(all_numeric_predictors()) %>%
        step_zv(all_numeric_predictors())
      
      # Alg selection
      for (alg in algorithms) {
        alg_name <- paste0(combo_base_name,alg)
        cat("  Algoritmo:", alg, "\n")
        
        tryCatch({
          # WF
          wf <- workflow() %>% 
            add_recipe(rec) %>% 
            add_model(model_specs[[alg]])
          
          # CV10F
          set.seed(1234)
          folds <- vfold_cv(train, strata = !!sym(response))
          
          # Tuning
          tune_res <- tune_race_anova(
            wf, 
            resamples = folds, 
            grid = 60,  
            metrics = metric_set(rmse),
            control = race_ctrl
          )
          
          # Best hyperparameters
          best_params <- select_best(tune_res, metric = "rmse")
          final_wf <- finalize_workflow(wf, best_params)
          
          # Fit
          final_fit <- fit(final_wf, train)
          
          # Save WF 
          if (return_fits) {
            size_fits[[alg_name]] <- final_fit
          }
          if (return_workflows) {
            size_workflows[[alg_name]] <- final_wf
          }
          
          # CV Metrics
          cv_res <- fit_resamples(final_wf, folds, control = ctrl)
          cv_preds <- collect_predictions(cv_res)
          cv_metrics <- calculate_metrics(cv_preds) %>% 
            mutate(val = "CV10F", Combinacao = alg_name, Algoritmo = alg)
          
          # Test Metrics
          test_preds <- predict(final_fit, test) %>% bind_cols(test)
          test_metrics <- calculate_metrics(test_preds) %>% 
            mutate(val = "Test", Combinacao = alg_name, Algoritmo = alg)
          
          # External Metrics
          external_metrics <- list()
          for (df_name in names(external_dfs)) {
            df <- external_dfs[[df_name]]
            
            # Derived predictors for external data
            derived_ext <- create_derived_predictors(
              df, combo, include_ratio_predictors, include_diff_predictors
            )
            df_ext <- derived_ext$data
            
            # Verify
            missing_vars <- setdiff(all_predictors, colnames(df_ext))
            if (length(missing_vars) > 0) {
              warning(paste("Variáveis faltantes em", df_name, ":", paste(missing_vars, collapse = ", ")))
              next
            }
            
            preds <- predict(final_fit, df_ext) %>% 
              bind_cols(df_ext %>% select(all_of(response)))
            metrics <- calculate_metrics(preds) %>% 
              mutate(val = df_name, Combinacao = alg_name, Algoritmo = alg)
            
            external_metrics[[df_name]] <- metrics
          }
          
          # Save results
          size_results[[alg_name]] <- cv_metrics$rmse
          size_best_params[[alg_name]] <- best_params
          all_metrics <- c(list(cv_metrics, test_metrics), external_metrics)
          size_tabs[[alg_name]] <- bind_rows(all_metrics)
          
        }, error = function(e) {
          message("Erro no algoritmo ", alg, " para combinação ", combo_base_name, ": ", e$message)
        })
      }
    }
    
    all_results[[paste0("size_", n)]] <- size_results
    all_best_params[[paste0("size_", n)]] <- size_best_params
    all_tabs[[paste0("size_", n)]] <- bind_rows(size_tabs)
    if (return_fits) {
      all_fits[[paste0("size_", n)]] <- size_fits
    }
    if (return_workflows) {
      all_workflows[[paste0("size_", n)]] <- size_workflows
    }
  }
  
  if (return_fits && return_workflows) {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs,
      fits = all_fits,
      workflows = all_workflows
    ))
  } else if (return_fits) {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs,
      fits = all_fits
    ))
  } else if (return_workflows) {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs,
      workflows = all_workflows
    ))
  } else {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs
    ))
  }
}
