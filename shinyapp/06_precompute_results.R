# =============================================================================
#  06_precompute_results.R - cross-validation and test results of the models of
#  the app, for the "All results" tab
#
#  Run once (it takes a while; Block Forest and DIABLO are the slowest):
#    source("shinyapp/06_precompute_results.R")      # from the repo root
#  Output: shinyapp/results/all_results.rds (one row per model). Each model is
#  saved when it finishes: if the script stops, run it again and it continues
#  with the models that are missing (delete the file to start from scratch).
#
#  Models
#    - Original models (random forest, XGBoost, Lasso, elastic net; DIABLO
#      components + classifier; Block Forest) on each omics layer alone (RNA,
#      protein, DNA methylation) and on the three layers together.
#    - DIABLO analysis: sparse DIABLO with its own prediction (RNA + protein +
#      methylation, the tuned ncomp and keepX of that analysis).
#
#  Settings of XGBoost, Lasso and elastic net (direct join)
#    - the RNA and protein models of the original XGBoost / elastic-net report
#      use the settings selected in that report ("tuned (report)");
#    - the other ones are tuned here with the same procedure and grids as the
#      report ("tuned (same grid)"): 5-fold training CV (same folds), mean
#      validation balanced accuracy; XGBoost depth 1 or 3 and 50, 100 or 200
#      rounds; glmnet alpha 0.25, 0.5 or 0.75 (Lasso: 1) and 12 lambdas; for a
#      single high-dimensional layer (RNA or methylation) also the number of
#      features (216, 500, 1000 or 2000). Ties: fewer features, then the
#      simpler setting, as in the report.
#    - random forest, Block Forest and the classifiers on DIABLO components
#      keep the settings of the workbench ("app default").
#  Because the settings were chosen with these same folds, the CV value of a
#  tuned model is slightly optimistic; the test value is not.
#
#  For every model
#    - CV: 5-fold cross-validation on the training patients (folds stratified
#      by stage, seed 20261110 as in the original report). The whole model is
#      fitted inside each fold (feature ranking, DIABLO, Block Forest weights).
#    - Test: the model trained on all training patients, evaluated on the test
#      patients with conformal prediction (40% of the holdout patients
#      calibrate, the rest are the test patients), as in the workbench.
# =============================================================================

if (file.exists("shinyapp/modeling.R")) setwd("shinyapp")
suppressPackageStartupMessages(library(ggplot2))
source("modeling.R", encoding = "UTF-8")

out_file <- "results/all_results.rds"
dir.create("results", showWarnings = FALSE)

# ---- Models -----------------------------------------------------------------
join_datasets  <- names(dataset_layers)            # each layer alone and the three together
multi_omics    <- multi_datasets
classifiers    <- c("rf", "xgb", "lasso", "enet")

new_config <- function(dataset, integration, model, origin = "Original models",
                       rna = 216, protein = 216, methylation = 216, params = NULL,
                       settings = "app default") {
  if (length(dataset_layers[[dataset]]) == 1) integration <- "none"
  list(dataset = dataset, integration = integration, model = model, origin = origin,
       rna = rna, protein = protein, methylation = methylation, params = params,
       settings = settings)
}
configs <- c(
  unlist(lapply(join_datasets, function(d) lapply(classifiers, function(m) new_config(d, "join", m))),
         recursive = FALSE),
  unlist(lapply(multi_omics, function(d) lapply(classifiers, function(m) new_config(d, "diablo", m))),
         recursive = FALSE),
  lapply(multi_omics, function(d) new_config(d, "block", "rf")),
  list(new_config("rna_protein_methylation", "diablo", "diablo", origin = "DIABLO analysis",
                  rna = 2000, protein = 200, methylation = 5000, settings = "tuned (analysis)"))
)

# Settings selected in the original XGBoost / elastic-net report (the lambdas
# are values of the same grid, printed rounded in the report)
closest_lambda <- function(value) enet_lambdas[which.min(abs(enet_lambdas - value))]
report_settings <- list(
  "rna|xgb"         = list(rna = 2000, params = list(depth = 1, rounds = 200)),
  "protein|xgb"     = list(params = list(depth = 1, rounds = 200)),
  "rna|enet"        = list(rna = 216, params = list(alpha = 0.25, lambda = closest_lambda(0.1110))),
  "protein|enet"    = list(params = list(alpha = 0.50, lambda = closest_lambda(0.0390)))
)

config_id <- function(cfg) paste(cfg$dataset, cfg$integration, cfg$model, sep = "|")

# ---- Helpers ----------------------------------------------------------------
settings_of <- function(cfg, conformal) {
  layers <- dataset_layers[[cfg$dataset]]
  list(dataset = cfg$dataset, integration = cfg$integration, model = cfg$model,
       rna = if ("RNA" %in% layers) cfg$rna else NULL,
       protein = if ("Protein" %in% layers) cfg$protein else NULL,
       methylation = if ("Methylation" %in% layers) cfg$methylation else NULL,
       conformal = conformal, params = cfg$params)
}

# Data of one fold: the training patients split into fit and validation
# patients; the feature ranking is learned again on the fit patients
fold_data <- function(data, fit, valid) {
  d <- data
  d$train <- lapply(data$train, function(x) x[fit, , drop = FALSE])
  d$holdout <- lapply(data$train, function(x) x[valid, , drop = FALSE])
  d$y_train <- data$y_train[fit]
  d$y_holdout <- data$y_train[valid]
  d$holdout_ids <- rownames(data$train[[1]])[valid]
  d$calibration <- integer()
  d$ranked <- lapply(d$train, function(x) {
    variances <- sort(matrixStats::colVars(x), decreasing = TRUE)
    names(variances)[variances > 0]
  })
  d
}

# Feature counts no larger than the features available in the fold
fit_counts <- function(settings, data) {
  keys <- c(RNA = "rna", Protein = "protein", Methylation = "methylation")
  for (layer in names(keys)) {
    key <- keys[[layer]]
    if (!is.null(settings[[key]])) settings[[key]] <- min(settings[[key]], length(data$ranked[[layer]]))
  }
  settings
}

balanced_accuracy <- function(observed, scores) {
  predicted <- factor(stages[max.col(scores[, stages, drop = FALSE], ties.method = "first")], levels = stages)
  counts <- table(observed, predicted)
  mean(diag(counts) / rowSums(counts))
}

# ---- Tuning (same procedure and grids as the original report) ---------------
tune_config <- function(data, cfg, folds) {
  layers <- dataset_layers[[cfg$dataset]]
  count_key <- if (identical(layers, "RNA")) "rna" else if (identical(layers, "Methylation")) "methylation" else NULL
  count_options <- if (is.null(count_key)) NA else c(216, 500, 1000, 2000)
  scores <- do.call(rbind, lapply(sort(unique(folds)), function(k) {
    d <- fold_data(data, which(folds != k), which(folds == k))
    do.call(rbind, lapply(count_options, function(count) {
      option <- cfg
      if (!is.null(count_key)) option[[count_key]] <- count
      blocks <- prepare_blocks(d, fit_counts(settings_of(option, FALSE), d))
      x <- do.call(cbind, blocks$train)
      x_valid <- do.call(cbind, blocks$holdout)
      weights <- stage_weights(d$y_train)
      if (cfg$model == "xgb") {
        training <- xgboost::xgb.DMatrix(x, label = as.integer(d$y_train) - 1L, weight = weights, nthread = 2)
        valid <- xgboost::xgb.DMatrix(x_valid, nthread = 2)
        do.call(rbind, lapply(c(1L, 3L), function(depth) {
          fit <- xgboost::xgb.train(
            params = list(objective = "multi:softprob", num_class = length(stages), max_depth = depth,
                          eta = 0.05, tree_method = "hist", nthread = 2),
            data = training, nrounds = 200, verbose = 0)
          do.call(rbind, lapply(c(50L, 100L, 200L), function(rounds) {
            p <- predict(fit, valid, iterationrange = c(1L, rounds), strict_shape = TRUE)
            colnames(p) <- stages
            data.frame(count = count, depth = depth, rounds = rounds,
                       balanced = balanced_accuracy(d$y_holdout, p))
          }))
        }))
      } else {
        alphas <- if (cfg$model == "lasso") 1 else c(0.25, 0.50, 0.75)
        do.call(rbind, lapply(alphas, function(alpha) {
          fit <- glmnet::glmnet(x, d$y_train, family = "multinomial", type.multinomial = "grouped",
                                weights = weights, alpha = alpha, lambda = enet_lambdas)
          p <- predict(fit, newx = x_valid, s = enet_lambdas, type = "response")
          do.call(rbind, lapply(seq_along(enet_lambdas), function(i) {
            data.frame(count = count, alpha = alpha, lambda = enet_lambdas[i],
                       balanced = balanced_accuracy(d$y_holdout, p[, , i]))
          }))
        }))
      }
    }))
  }))
  scores$count[is.na(scores$count)] <- 0
  if (cfg$model == "xgb") {
    means <- aggregate(balanced ~ count + depth + rounds, scores, mean)
    best <- means[order(-means$balanced, means$count, means$depth, means$rounds)[1], ]
    cfg$params <- list(depth = best$depth, rounds = best$rounds)
  } else {
    means <- aggregate(balanced ~ count + alpha + lambda, scores, mean)
    best <- means[order(-means$balanced, means$count, -means$lambda, abs(means$alpha - 0.5))[1], ]
    cfg$params <- list(alpha = best$alpha, lambda = best$lambda)
  }
  if (!is.null(count_key)) cfg[[count_key]] <- best$count
  cfg$settings <- "tuned (same grid)"
  cfg
}

describe_params <- function(cfg) {
  p <- cfg$params
  if (is.null(p)) return(cfg$settings)
  text <- if (!is.null(p$depth)) {
    sprintf("depth %d, %d rounds", as.integer(p$depth), as.integer(p$rounds))
  } else {
    sprintf("alpha %.2f, lambda %.4g", if (is.null(p$alpha)) 1 else p$alpha, p$lambda)
  }
  paste0(cfg$settings, ": ", text)   # the app shows only the values (settings_text())
}

# ---- Evaluation -------------------------------------------------------------
evaluate_model <- function(data, cfg, folds) {
  started <- Sys.time()
  cv <- do.call(rbind, lapply(sort(unique(folds)), function(k) {
    d <- fold_data(data, which(folds != k), which(folds == k))
    run <- run_analysis(d, fit_counts(settings_of(cfg, FALSE), d))
    data.frame(fold = k, balanced = run$balanced, accuracy = run$accuracy, macro_f1 = run$macro_f1)
  }))
  final <- run_analysis(data, settings_of(cfg, TRUE))
  auc <- vapply(final$roc, function(curve) as.numeric(pROC::auc(curve)), numeric(1))
  data.frame(
    id = config_id(cfg), Origin = cfg$origin, Data = dataset_labels[[cfg$dataset]],
    Integration = if (cfg$integration == "none") "-" else integration_labels[[cfg$integration]],
    Model = if (cfg$integration == "block") "Block Forest" else model_labels[[cfg$model]],
    Features = paste(sprintf("%s %d", names(final$features), final$features), collapse = ", "),
    Settings = describe_params(cfg),
    cv_balanced = mean(cv$balanced), cv_balanced_sd = sd(cv$balanced),
    cv_accuracy = mean(cv$accuracy), cv_macro_f1 = mean(cv$macro_f1),
    cv_folds = paste(round(cv$balanced, 3), collapse = " / "),
    test_balanced = final$balanced, test_accuracy = final$accuracy, test_macro_f1 = final$macro_f1,
    test_auc = mean(auc), coverage = mean(final$conformal$covered),
    set_size = mean(final$conformal$size),
    test_patients = final$n_test, minutes = as.numeric(difftime(Sys.time(), started, units = "mins")),
    error = NA_character_, stringsAsFactors = FALSE)
}

# Stack two result tables with possibly different columns
bind_results <- function(a, b) {
  if (is.null(a)) return(b)
  columns <- union(names(a), names(b))
  for (column in setdiff(columns, names(a))) a[[column]] <- NA
  for (column in setdiff(columns, names(b))) b[[column]] <- NA
  rbind(a[columns], b[columns])
}

# ---- Run --------------------------------------------------------------------
data <- load_stage_data()
set.seed(20261110)                                  # folds as in the original report
folds <- integer(length(data$y_train))
for (rows in split(seq_along(data$y_train), data$y_train)) {
  folds[sample(rows)] <- rep(1:5, length.out = length(rows))
}

results <- if (file.exists(out_file)) readRDS(out_file) else NULL
# Keep only the models of the current list
if (!is.null(results)) {
  results <- results[results$id %in% vapply(configs, config_id, character(1)), , drop = FALSE]
  saveRDS(results, out_file)
}
done <- if (is.null(results)) character() else results$id[is.na(results$error)]
for (i in seq_along(configs)) {
  cfg <- configs[[i]]
  id <- config_id(cfg)
  if (id %in% done) next
  message(sprintf("[%s] %d/%d %s", format(Sys.time(), "%H:%M"), i, length(configs), id))
  row <- tryCatch({
    key <- paste(cfg$dataset, cfg$model, sep = "|")
    if (cfg$integration %in% c("none", "join") && key %in% names(report_settings)) {
      cfg <- modifyList(cfg, report_settings[[key]])
      cfg$settings <- "tuned (report)"
    } else if (cfg$integration %in% c("none", "join") && cfg$model %in% c("xgb", "lasso", "enet")) {
      message("  tuning")
      cfg <- tune_config(data, cfg, folds)
    }
    evaluate_model(data, cfg, folds)
  }, error = function(e) {
    message("  error: ", conditionMessage(e))
    data.frame(id = id, Origin = cfg$origin, Data = dataset_labels[[cfg$dataset]],
               Integration = cfg$integration, Model = cfg$model, error = conditionMessage(e),
               stringsAsFactors = FALSE)
  })
  if (!is.null(results)) results <- results[results$id != id, , drop = FALSE]
  results <- bind_results(results, row)
  saveRDS(results, out_file)
}
message("Done: ", out_file)
