# =============================================================================
#  app_helpers.R - small helpers of the interface: (i) help icons, dataset
#  choices and the JavaScript conditions that show the inputs of each layer
# =============================================================================

repo_url <- "https://github.com/ssayolsb/Group2-tcga-brca-multiomics"

# (i) icon with a help text (hover, focus or tap)
info <- function(text) {
  bslib::tooltip(tags$span(class = "info-icon", tabindex = "0", icon("circle-info")),
                 text, placement = "right")
}
with_info <- function(label, text) tags$span(label, " ", info(text))

help_text <- list(
  dataset = "Which omics data the model uses: one layer (RNA, protein or DNA methylation) or the three together. The three layers together need an integration method.",
  integration = "Join directly: all features in one matrix. Block Forest: a random forest with a weight per layer, tuned on the out-of-bag error. DIABLO: supervised components per layer (2 per layer) used as classifier inputs, or DIABLO's own prediction.",
  model = "Random forest and XGBoost are tree ensembles; Lasso and elastic net are penalised multinomial regressions; DIABLO (prediction) lets the sparse DIABLO model predict the stage (only with the DIABLO integration).",
  counts = "The most variable features of this layer in the training patients are used.",
  conformal = "Uses 40% of the holdout patients to calibrate prediction sets (LAC): each test patient gets the set of stages compatible with 90% coverage. The other holdout patients are the test patients."
)

# JavaScript conditions on a dataset input (id "dataset" by default)
# Written out because shiny sources R/ before modeling.R (see multi_datasets).
multi_dataset_js <- "(input.dataset === 'rna_protein_methylation')"
layer_js <- function(layer, id = "dataset") {
  datasets <- names(dataset_layers)[vapply(dataset_layers, function(l) layer %in% l, logical(1))]
  sprintf("[%s].indexOf(input.%s) > -1", paste0("'", datasets, "'", collapse = ", "), id)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# One row of the "Compare runs" table
run_summary <- function(run, number) {
  s <- run$settings
  auc <- vapply(run$roc, function(curve) as.numeric(pROC::auc(curve)), numeric(1))
  data.frame(Run = number, Data = dataset_labels[[s$dataset]],
             Integration = if (s$integration == "none") "-" else integration_labels[[s$integration]],
             Model = if (s$integration == "block") "Block Forest" else model_labels[[s$model]],
             Features = paste(sprintf("%s %d", names(run$features), run$features), collapse = ", "),
             `Balanced accuracy` = round(run$balanced, 3),
             Accuracy = round(run$accuracy, 3),
             `Macro F1` = round(run$macro_f1, 3), `Mean AUC` = round(mean(auc), 3),
             `Test patients` = run$n_test,
             `Conformal coverage` = if (is.null(run$conformal)) NA else round(mean(run$conformal$covered), 3),
             check.names = FALSE)
}
