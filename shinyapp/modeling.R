stages <- c("Stage I", "Stage II", "Stage III/IV")
dataset_labels <- c(rna = "RNA", protein = "Protein", integrated = "RNA + protein")
integration_labels <- c(join = "Direct join", block = "Block Forest", diablo = "DIABLO")
model_labels <- c(rf = "Random forest", xgb = "XGBoost",
                  lasso = "Lasso", enet = "Elastic net (glmnet)")
model_labels <- c(model_labels, diablo = "DIABLO (prediction)")
dataset_labels <- c(dataset_labels, methylation = "DNA methylation",
                    rna_protein_methylation = "RNA + protein + methylation")
layer_labels <- c(RNA = "RNA", Protein = "Protein", Methylation = "DNA methylation")
dataset_layers <- list(rna = "RNA", protein = "Protein", methylation = "Methylation",
                       rna_protein_methylation = c("RNA", "Protein", "Methylation"))
dataset_choices <- setNames(names(dataset_layers), dataset_labels[names(dataset_layers)])
multi_datasets <- "rna_protein_methylation"
max_methylation <- 10000
# keepX per component tuned in the DIABLO report.
diablo_keepx <- list(RNA = c(50, 100), Protein = c(5, 40), Methylation = c(25, 100))
# Lambda grid of the XGBoost / elastic-net report (tuned settings).
enet_lambdas <- 10^seq(-0.5, -3, length.out = 12)

# Fill missing values with training medians.
fill_medians <- function(x, medians) {
  missing <- which(is.na(x), arr.ind = TRUE)
  if (nrow(missing)) x[missing] <- medians[missing[, "col"]]
  x
}

# Load the project data once per app process. These files are read-only inputs.
load_stage_data <- function() {
  # App data (05_make_app_data.R) if present, otherwise the project data.
  data_dir <- if (dir.exists("data/processed")) "data" else if (dir.exists("../../data")) "../../data" else "../data"
  read_input <- function(path) readRDS(file.path(data_dir, path))
  read_assay <- function(path) {
    data <- read_input(path)
    stopifnot(identical(data$patientID, rownames(data)), !anyNA(data$patientID),
              !anyDuplicated(data$patientID), !anyDuplicated(names(data)))
    values <- as.matrix(data[, names(data) != "patientID", drop = FALSE])
    stopifnot(is.numeric(values), !any(is.infinite(values)))
    values
  }

  train <- list(RNA = read_assay("processed/rna.train.filtered.rds"),
                Protein = read_assay("processed/prot.train.filtered.rds"))
  holdout <- list(RNA = read_assay("processed/rna.test.filtered.rds"),
                  Protein = read_assay("processed/prot.test.filtered.rds"))
  protein_original <- read_assay("data_split/prot.train.rds")
  qc <- read_input("processed/prot_qc_params.rds")
  screen <- colnames(protein_original)[colMeans(is.na(protein_original)) <= 0.20]
  screen <- intersect(screen, colnames(holdout$Protein))
  eligible <- rowMeans(is.na(holdout$Protein[, screen, drop = FALSE])) <= qc$max_na_sample

  train_ids <- intersect(rownames(train$RNA), rownames(train$Protein))
  holdout_ids <- intersect(rownames(holdout$RNA), rownames(holdout$Protein)[eligible])
  stopifnot(!any(train_ids %in% holdout_ids),
            identical(colnames(train$Protein), qc$kept_proteins))
  for (assay in names(train)) {
    train[[assay]] <- train[[assay]][train_ids, , drop = FALSE]
    holdout[[assay]] <- holdout[[assay]][holdout_ids, colnames(train[[assay]]), drop = FALSE]
  }
  stopifnot(!anyNA(train$RNA), !anyNA(holdout$RNA))

  # Learn protein medians and variance rankings from the filtered training set.
  medians <- matrixStats::colMedians(train$Protein, na.rm = TRUE)
  stopifnot(all(is.finite(medians)))
  fill_protein <- function(x) {
    missing <- which(is.na(x), arr.ind = TRUE)
    x[missing] <- medians[missing[, "col"]]
    x
  }
  train$Protein <- fill_protein(train$Protein)
  holdout$Protein <- fill_protein(holdout$Protein)
  # Methylation: complete training regions, the most variable ones kept.
  meth_train <- read_assay("processed/meth.train.region_M.rds")
  meth_test <- read_assay("processed/meth.test.region_M.rds")
  stopifnot(all(train_ids %in% rownames(meth_train)), all(holdout_ids %in% rownames(meth_test)))
  meth_train <- meth_train[train_ids, colSums(is.na(meth_train)) == 0, drop = FALSE]
  meth_var <- matrixStats::colVars(meth_train)
  keep <- colnames(meth_train)[order(meth_var, decreasing = TRUE)]
  keep <- head(keep, min(max_methylation, sum(meth_var > 0)))
  train$Methylation <- meth_train[, keep, drop = FALSE]
  holdout$Methylation <- fill_medians(meth_test[holdout_ids, keep, drop = FALSE],
                                      matrixStats::colMedians(train$Methylation))
  rm(meth_train, meth_test)
  ranked <- lapply(train, function(x) {
    variances <- sort(matrixStats::colVars(x), decreasing = TRUE)
    names(variances)[variances > 0]
  })

  read_stage <- function(path, ids) {
    clinical <- read_input(path)
    stopifnot(identical(clinical$patientID, rownames(clinical)),
              all(ids %in% rownames(clinical)))
    groups <- c("stage i", "stage ii", "stage iii")
    stage <- match(clinical[ids, "outcome_p.stage"], groups)
    stopifnot(!anyNA(stage))
    factor(stages[stage], levels = stages)
  }
  y_train <- read_stage("data_split/clin.train.rds", train_ids)
  y_holdout <- read_stage("data_split/clin.test.rds", holdout_ids)
  stopifnot(all(table(y_train) >= 5), all(table(y_holdout) > 0))

  # Reproduce the calibration and test split used in the reports.
  set.seed(20261109)
  calibration <- sample(seq_along(holdout_ids), size = round(0.4 * length(holdout_ids)))
  list(train = train, holdout = holdout, ranked = ranked, y_train = y_train,
       y_holdout = y_holdout, calibration = calibration, holdout_ids = holdout_ids,
       layers = names(train))
}

require_model_packages <- function(settings) {
  model_package <- c(rf = "randomForest", xgb = "xgboost", lasso = "glmnet", enet = "glmnet",
                     diablo = "mixOmics")
  packages <- if (settings$integration == "block") "blockForest" else model_package[[settings$model]]
  packages <- c(packages, "pROC")
  if (settings$integration == "diablo") packages <- c(packages, "mixOmics")
  if (settings$conformal) packages <- c(packages, "predictset")
  available <- vapply(packages, function(package) {
    suppressWarnings(suppressPackageStartupMessages(requireNamespace(package, quietly = TRUE)))
  }, logical(1))
  if (!all(available)) stop("Install these R packages before running: ",
                            paste(packages[!available], collapse = ", "), call. = FALSE)
}

# Select the requested number of predictors in each assay by training variance.
prepare_blocks <- function(data, settings) {
  counts <- vapply(dataset_layers[[settings$dataset]], function(layer) {
    count <- switch(layer, RNA = settings$rna, Protein = settings$protein,
                    Methylation = settings$methylation)
    if (is.null(count)) NA_real_ else as.numeric(count)
  }, numeric(1))
  train <- holdout <- list()
  for (assay in names(counts)) {
    count <- counts[[assay]]
    available <- length(data$ranked[[assay]])
    if (length(count) != 1L || !is.finite(count) || count != floor(count) ||
        count < 2L || count > available) {
      stop(sprintf("%s feature count must be a whole number from 2 to %d.", assay, available),
           call. = FALSE)
    }
    features <- head(data$ranked[[assay]], count)
    train[[assay]] <- data$train[[assay]][, features, drop = FALSE]
    holdout[[assay]] <- data$holdout[[assay]][, features, drop = FALSE]
    colnames(train[[assay]]) <- colnames(holdout[[assay]]) <- paste0(assay, "__", features)
  }
  list(train = train, holdout = holdout, counts = counts)
}

# Fit supervised components on training patients and project new patients.
diablo_components <- function(train, holdout, y) {
  design <- matrix(0.1, nrow = length(train), ncol = length(train),
                   dimnames = list(names(train), names(train)))
  diag(design) <- 0
  model <- mixOmics::block.splsda(
    train, y, ncomp = 2, design = design,
    keepX = lapply(train, function(block) rep(ncol(block), 2)),
    max.iter = 500
  )
  join_components <- function(variates) {
    do.call(cbind, lapply(names(train), function(assay) {
      x <- variates[[assay]]
      colnames(x) <- paste0(assay, "_component_", seq_len(ncol(x)))
      x
    }))
  }
  list(train = join_components(model$variates),
       holdout = join_components(predict(model, newdata = holdout)$variates))
}

stage_weights <- function(y) {
  weights <- 1 / as.numeric(table(y)[as.character(y)])
  weights / mean(weights)
}

# The package documents binary classification. Fit one forest per stage and
# normalize the three one-vs-rest scores so they share a common scale.
predict_block_forest <- function(x, y, new_x, counts) {
  ends <- cumsum(counts)
  blocks <- lapply(seq_along(counts), function(i) (ends[i] - counts[i] + 1):ends[i])
  # blockfor() constructs a data frame, so use the same safe names for prediction.
  colnames(x) <- colnames(new_x) <- make.names(colnames(x), unique = TRUE)
  test_frame <- as.data.frame(new_x)
  scores <- vapply(stages, function(stage) {
    binary <- factor(ifelse(y == stage, "Target", "Other"), levels = c("Other", "Target"))
    fit <- blockForest::blockfor(
      x, binary, blocks = blocks, block.method = "BlockForest", num.trees = 500,
      nsets = 20, num.trees.pre = 100, probability = TRUE, case.weights = stage_weights(y),
      num.threads = 2, verbose = FALSE
    )
    predict(fit$forest, data = test_frame, num.threads = 2, verbose = FALSE)$predictions[, "Target"]
  }, numeric(nrow(new_x)))
  scores[rowSums(scores) == 0, ] <- 1
  scores / rowSums(scores)
}

predict_classifier <- function(x, y, new_x, model, integration, params = NULL) {
  weights <- stage_weights(y)
  if (model == "rf") {
    fit <- randomForest::randomForest(x, y, weights = weights, ntree = 1000, nodesize = 10)
    return(list(scores = predict(fit, newdata = new_x, type = "prob"),
                details = "1,000 trees · node size 10"))
  }
  if (model == "xgb") {
    depth <- if (is.null(params$depth)) 2 else params$depth
    rounds <- if (is.null(params$rounds)) 100 else params$rounds
    training <- xgboost::xgb.DMatrix(x, label = as.integer(y) - 1L,
                                     weight = weights, nthread = 2)
    fit <- xgboost::xgb.train(
      params = list(objective = "multi:softprob", num_class = length(stages),
                    max_depth = depth, eta = 0.05, tree_method = "hist", nthread = 2),
      data = training, nrounds = rounds, verbose = 0
    )
    scores <- matrix(predict(fit, xgboost::xgb.DMatrix(new_x, nthread = 2), reshape = TRUE), nrow = nrow(new_x))
    colnames(scores) <- stages
    return(list(scores = scores, details = sprintf("%d rounds · depth %d · learning rate 0.05",
                                                   as.integer(rounds), as.integer(depth))))
  }

  alpha <- if (model == "lasso") 1 else 0.5
  if (!is.null(params$lambda) && integration != "diablo") {
    # Settings chosen beforehand by 5-fold training CV.
    if (!is.null(params$alpha)) alpha <- params$alpha
    fit <- glmnet::glmnet(x, y, family = "multinomial", type.multinomial = "grouped",
                          alpha = alpha, weights = weights, lambda = enet_lambdas)
    scores <- predict(fit, newx = new_x, s = params$lambda, type = "response")[, , 1]
    return(list(scores = scores, details = sprintf("Alpha %.2f · lambda %.4g · tuned by 5-fold training CV",
                                                   alpha, params$lambda)))
  }
  if (integration == "diablo") {
    # The supervised DIABLO projection uses all training labels. Use a fixed
    # penalty here so inner validation labels do not enter the projection.
    lambda <- 0.03
    fit <- glmnet::glmnet(x, y, family = "multinomial", type.multinomial = "grouped",
                          alpha = alpha, weights = weights, lambda = c(0.3, 0.1, lambda))
    tuning <- "fixed penalty on DIABLO components"
  } else {
    fold_id <- integer(length(y))
    for (rows in split(seq_along(y), y)) {
      fold_id[sample(rows)] <- rep(1:5, length.out = length(rows))
    }
    fit <- glmnet::cv.glmnet(x, y, family = "multinomial", type.multinomial = "grouped",
                             alpha = alpha, weights = weights, foldid = fold_id)
    lambda <- fit$lambda.min
    tuning <- "5-fold training CV"
  }
  scores <- predict(fit, newx = new_x, s = lambda, type = "response")[, , 1]
  list(scores = scores, details = sprintf("Alpha %.2f · lambda %.4g · %s", alpha, lambda, tuning))
}

# Sparse DIABLO of the DIABLO report. The weighted prediction of the layers is
# rescaled to scores that sum to 1, as needed for ROC curves and LAC sets.
predict_diablo <- function(train, y, holdout) {
  design <- matrix(0.1, nrow = length(train), ncol = length(train),
                   dimnames = list(names(train), names(train)))
  diag(design) <- 0
  keepx <- lapply(names(train), function(layer) {
    k <- if (layer %in% names(diablo_keepx)) diablo_keepx[[layer]] else rep(ncol(train[[layer]]), 2)
    pmin(k, ncol(train[[layer]]))
  })
  names(keepx) <- names(train)
  model <- mixOmics::block.splsda(train, y, ncomp = 2, design = design, keepX = keepx, max.iter = 500)
  prediction <- predict(model, newdata = holdout, dist = "centroids.dist")
  scores <- matrix(prediction$WeightedPredict[, , 2], ncol = length(stages))
  scores <- pmax(scores, 0)
  scores[rowSums(scores) == 0, ] <- 1
  scores <- scores / rowSums(scores)
  colnames(scores) <- stages
  stopifnot(nrow(scores) == nrow(holdout[[1]]), all(is.finite(scores)))
  keep_text <- paste(sprintf("%s %s", names(keepx), vapply(keepx, paste, character(1), collapse = "/")),
                     collapse = ", ")
  list(scores = scores, details = sprintf("Sparse DIABLO · 2 components · design 0.1 · keepX %s", keep_text))
}

point_results <- function(scores, truth, y_train) {
  predicted <- factor(stages[max.col(scores, ties.method = "first")], levels = stages)
  cm <- table(Observed = truth, Predicted = predicted)
  correct <- diag(cm)
  recall <- correct / rowSums(cm)
  precision <- correct / colSums(cm)
  f1 <- 2 * correct / (rowSums(cm) + colSums(cm))
  majority <- stages[which.max(table(y_train))]

  list(predicted = predicted, confusion = cm, accuracy = mean(predicted == truth),
       balanced = mean(recall), macro_f1 = mean(f1), baseline = mean(truth == majority),
       majority = majority,
       stage_metrics = data.frame(Stage = stages, Patients = as.integer(rowSums(cm)),
                                  Recall = unname(recall), Precision = unname(precision),
                                  F1 = unname(f1)))
}

# Score each observed stage against all other stages using its own score column.
stage_roc <- function(scores, truth) {
  setNames(lapply(stages, function(stage) {
    pROC::roc(response = as.integer(truth == stage), predictor = scores[, stage],
              levels = c(0, 1), direction = "<", quiet = TRUE)
  }), stages)
}

# LAC ranks each candidate against the nonconformity scores of the true
# stages in the calibration subset. The test labels are used only for coverage.
conformal_results <- function(cal_scores, cal_truth, test_scores, test_truth) {
  label_index <- function(y) cbind(seq_along(y), as.integer(y))
  calibration_error <- 1 - cal_scores[label_index(cal_truth)]
  p_values <- predictset::conformal_pvalue(calibration_error, as.vector(1 - test_scores))
  included <- matrix(p_values > 0.10, nrow = nrow(test_scores),
                     dimnames = dimnames(test_scores))
  covered <- included[label_index(test_truth)]
  size <- rowSums(included)
  sets <- apply(included, 1, function(row) {
    if (any(row)) paste(stages[row], collapse = " / ") else "(empty)"
  })
  groups <- c(list(Overall = seq_along(test_truth)), split(seq_along(test_truth), test_truth))
  coverage <- do.call(rbind, lapply(names(groups), function(group) {
    rows <- groups[[group]]
    interval <- binom.test(sum(covered[rows]), length(rows))$conf.int
    data.frame(Group = group, Patients = length(rows), Coverage = mean(covered[rows]),
               Lower = interval[1], Upper = interval[2])
  }))
  list(sets = sets, covered = covered, size = size, coverage = coverage)
}

run_analysis <- function(data, settings, progress = function(detail) NULL) {
  require_model_packages(settings)
  progress("Selecting training features")
  blocks <- prepare_blocks(data, settings)

  set.seed(20261110)
  if (settings$integration == "diablo" && settings$model != "diablo") {
    progress("Learning DIABLO components")
    matrices <- diablo_components(blocks$train, blocks$holdout, data$y_train)
  } else {
    matrices <- list(train = do.call(cbind, blocks$train), holdout = do.call(cbind, blocks$holdout))
  }

  progress("Fitting the classifier and scoring the holdout")
  if (settings$integration == "block") {
    fit <- list(scores = predict_block_forest(matrices$train, data$y_train,
                                              matrices$holdout, blocks$counts),
                details = "Three one-vs-rest forests · 500 trees each · 20 tuning sets")
  } else if (settings$model != "diablo") {
    fit <- predict_classifier(matrices$train, data$y_train, matrices$holdout,
                              settings$model, settings$integration, settings$params)
  }
  if (settings$model == "diablo") {
    if (settings$integration != "diablo") {
      stop("'DIABLO (prediction)' needs the DIABLO integration.", call. = FALSE)
    }
    fit <- predict_diablo(blocks$train, data$y_train, blocks$holdout)
  }
  scores <- fit$scores[, stages, drop = FALSE]
  stopifnot(nrow(scores) == length(data$y_holdout), all(is.finite(scores)),
            all(abs(rowSums(scores) - 1) < 1e-6))

  calibration <- if (settings$conformal) data$calibration else integer()
  test <- setdiff(seq_along(data$y_holdout), calibration)
  truth <- data$y_holdout[test]
  stopifnot(all(table(truth) > 0))
  test_scores <- scores[test, , drop = FALSE]
  result <- point_results(test_scores, truth, data$y_train)
  result$roc <- stage_roc(test_scores, truth)
  result$predictions <- data.frame(
    Patient = data$holdout_ids[test], Observed = truth, Predicted = result$predicted,
    round(test_scores, 3), row.names = NULL, check.names = FALSE
  )
  result$settings <- settings
  result$features <- blocks$counts
  result$model_inputs <- ncol(matrices$train)
  result$details <- fit$details
  result$n_train <- length(data$y_train)
  result$n_cal <- length(calibration)
  result$n_test <- length(test)

  if (settings$conformal) {
    progress("Calibrating LAC prediction sets")
    result$conformal <- conformal_results(scores[calibration, , drop = FALSE],
                                         data$y_holdout[calibration], test_scores, truth)
    result$predictions[["Prediction set"]] <- result$conformal$sets
    result$predictions$Covered <- result$conformal$covered
  }
  result
}

# New patients: refit the model, score holdout and new patients together, and
# calibrate LAC sets on all holdout patients. New data are not imputed.
# Features (original names) that the settings use, by layer.
needed_features <- function(data, settings) {
  blocks <- prepare_blocks(data, settings)
  lapply(blocks$train, function(x) sub("^[^_]+__", "", colnames(x)))
}

# Check that the new data are complete; return the patients in every layer.
check_new_data <- function(data, settings, new_layers) {
  features <- needed_features(data, settings)
  problems <- character()
  for (layer in names(features)) {
    x <- new_layers[[layer]]
    if (is.null(x)) {
      problems <- c(problems, paste0(layer_labels[[layer]], ": file missing."))
      next
    }
    absent <- setdiff(features[[layer]], colnames(x))
    if (length(absent)) {
      problems <- c(problems, sprintf("%s: %d of the %d features used by the model are not in the file (e.g. %s).",
                                      layer_labels[[layer]], length(absent), length(features[[layer]]),
                                      paste(head(absent, 5), collapse = ", ")))
    } else if (anyNA(x[, features[[layer]], drop = FALSE])) {
      n <- sum(rowSums(is.na(x[, features[[layer]], drop = FALSE])) > 0)
      problems <- c(problems, sprintf("%s: missing or unknown values in %d patient(s).",
                                      layer_labels[[layer]], n))
    }
  }
  if (length(problems)) {
    stop(paste(c("The data are not complete for this model (the app does not impute):", problems),
               collapse = "\n"), call. = FALSE)
  }
  ids <- Reduce(intersect, lapply(new_layers[names(features)], rownames))
  if (!length(ids)) stop("No patient is present in every file.", call. = FALSE)
  ids
}

score_new_patients <- function(data, settings, new_layers, level = 0.90) {
  settings$conformal <- TRUE
  if (settings$integration == "block") settings$model <- "rf"
  require_model_packages(settings)
  if (settings$model == "diablo" && settings$integration != "diablo") {
    stop("'DIABLO (prediction)' needs the DIABLO integration.", call. = FALSE)
  }
  ids <- check_new_data(data, settings, new_layers)

  # Add the new patients after the holdout; unused columns get training medians.
  extended <- data
  for (layer in dataset_layers[[settings$dataset]]) {
    old <- data$holdout[[layer]]
    new <- matrix(NA_real_, length(ids), ncol(old), dimnames = list(ids, colnames(old)))
    common <- intersect(colnames(old), colnames(new_layers[[layer]]))
    new[, common] <- new_layers[[layer]][ids, common, drop = FALSE]
    new <- fill_medians(new, matrixStats::colMedians(data$train[[layer]]))
    extended$holdout[[layer]] <- rbind(old, new)
  }

  blocks <- prepare_blocks(extended, settings)
  set.seed(20261110)
  if (settings$integration == "diablo" && settings$model != "diablo") {
    matrices <- diablo_components(blocks$train, blocks$holdout, data$y_train)
  } else {
    matrices <- list(train = do.call(cbind, blocks$train), holdout = do.call(cbind, blocks$holdout))
  }
  scores <- if (settings$integration == "block") {
    predict_block_forest(matrices$train, data$y_train, matrices$holdout, blocks$counts)
  } else if (settings$model == "diablo") {
    predict_diablo(blocks$train, data$y_train, blocks$holdout)$scores
  } else {
    predict_classifier(matrices$train, data$y_train, matrices$holdout,
                       settings$model, settings$integration, settings$params)$scores
  }
  scores <- matrix(scores[, stages], ncol = length(stages), dimnames = list(NULL, stages))
  stopifnot(nrow(scores) == length(data$y_holdout) + length(ids), all(is.finite(scores)))

  n <- length(data$y_holdout)
  calibration_error <- 1 - scores[cbind(seq_len(n), as.integer(data$y_holdout))]
  calibration_error <- calibration_error[!data$holdout_ids %in% ids]   # a test patient does not calibrate itself
  new_scores <- scores[n + seq_along(ids), , drop = FALSE]
  rownames(new_scores) <- ids
  uncertainty_table(new_scores, calibration_error, level)
}

# LAC sets with the same rule as conformal_results().
conformal_sets <- function(calibration_error, scores, level = 0.90) {
  p_values <- predictset::conformal_pvalue(calibration_error, as.vector(1 - scores))
  included <- matrix(p_values > 1 - level, nrow = nrow(scores), dimnames = dimnames(scores))
  sets <- apply(included, 1, function(row) if (any(row)) paste(stages[row], collapse = " / ") else "(empty)")
  list(sets = unname(sets), size = rowSums(included))
}

# Score a stage needs to enter the set (dashed line in the patient plot).
conformal_threshold <- function(calibration_error, level = 0.90) {
  n <- length(calibration_error)
  k <- ceiling((n + 1) * level)
  1 - (if (k > n) 1 else sort(calibration_error)[k])
}

uncertainty_table <- function(scores, calibration_error, level = 0.90) {
  sets <- conformal_sets(calibration_error, scores, level)
  data.frame(Patient = rownames(scores),
             Predicted = stages[max.col(scores, ties.method = "first")],
             Probability = round(apply(scores, 1, max), 3),
             round(scores, 3),
             `Prediction set` = sets$sets,
             Uncertainty = ifelse(sets$size == 1, "Low: one stage",
                           ifelse(sets$size == 2, "Medium: two stages possible",
                           ifelse(sets$size == 3, "High: any stage possible",
                                  "Unusual patient: no stage fits"))),
             `Coverage target` = sprintf("%.0f%%", 100 * level),
             `Set threshold` = round(conformal_threshold(calibration_error, level), 3),
             row.names = NULL, check.names = FALSE)
}
