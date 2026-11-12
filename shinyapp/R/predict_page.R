# =============================================================================
#  predict_page.R - "Predict new patients" page
#  The model is set up on this page (same choices as the workbench) and trained
#  on the training patients; a chosen test patient is then predicted, and the
#  other holdout patients calibrate its conformal prediction set
#  (score_new_patients() in modeling.R). Uploading a data set: coming soon.
# =============================================================================

predict_page <- function(data) {
  feature_limits <- lengths(data$ranked)
  layer_files <- lapply(data$layers, function(layer) {
    fileInput(paste0("p_file_", layer),
              with_info(paste(layer_labels[[layer]], "file"),
                        "Rows = patients (patientID column or first column), columns = features with the names of the processed files. Same scale as the training data, no missing values."),
              accept = c(".csv", ".tsv", ".txt", ".rds"))
  })
  count_input <- function(layer, id, label, value) {
    if (!layer %in% data$layers) return(NULL)
    conditionalPanel(layer_js(layer, "p_dataset"),
                     numericInput(id, label, min = 2, max = feature_limits[[layer]],
                                  value = min(value, feature_limits[[layer]]), step = 1))
  }
  multi_js <- gsub("input.dataset", "input.p_dataset", multi_dataset_js, fixed = TRUE)

  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      title = "1. Set up the model",
      width = 340,
      selectInput("p_dataset", with_info("Training data", help_text$dataset), choices = dataset_choices),
      conditionalPanel(multi_js,
                       selectInput("p_integration", with_info("Integration method", help_text$integration),
                                   choices = c("Join directly" = "join", "Block Forest" = "block", "DIABLO" = "diablo"))),
      conditionalPanel(sprintf("!(%s && input.p_integration === 'block')", multi_js),
                       selectInput("p_model", with_info("Classifier", help_text$model),
                                   choices = setNames(names(model_labels), model_labels))),
      count_input("RNA", "p_rna", "Number of RNA genes", 216),
      count_input("Protein", "p_protein", "Number of protein markers", feature_limits[["Protein"]]),
      count_input("Methylation", "p_methylation", "Number of methylation regions", 216),
      sliderInput("p_conf_level",
                  with_info("Target coverage of the prediction sets", "Proportion of patients whose true stage should be inside their prediction set. Higher coverage gives larger sets."),
                  min = 0.70, max = 0.95, value = 0.90, step = 0.05),
      helpText("The model is trained on the training patients when you predict; the holdout patients calibrate the prediction sets.")
    ),
    bslib::card(
      class = "run-card",
      bslib::card_body(
        tags$p(class = "stage-kicker mb-1", "TCGA-BRCA  /  NEW PATIENTS"),
        tags$h3(class = "mb-2", "Predict the stage of new patients"),
        tags$p(class = "stage-note mb-0",
               "1. Set up the model on the left. 2. Choose a test patient and predict its stage. Every prediction comes with its probability and its conformal prediction set.")
      )
    ),
    bslib::navset_card_tab(
      selected = "2. Single patient",
      bslib::nav_panel(
        "2. Upload a data set", icon = icon("upload"),
        tags$div(class = "alert alert-danger",
                 icon("circle-info"), tags$b(" Coming soon. "),
                 "Uploading a data set will be available in a future version of the app. Use 'Single patient' to try the predictions."),
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(
            width = 330,
            tags$fieldset(disabled = NA,                   # shown but not active yet
              uiOutput("p_needed"),
              layer_files,
              actionButton("p_predict_upload", "Predict uploaded patients", class = "btn-primary w-100",
                           icon = icon("wand-magic-sparkles")))
          ),
          tags$p(class = "stage-note", "Upload one file per data layer (patients in rows, features in columns).")
        )
      ),
      bslib::nav_panel(
        "2. Single patient", icon = icon("user"),
        bslib::layout_sidebar(
          sidebar = bslib::sidebar(
            width = 330,
            selectizeInput("p_sp_patient",
                           with_info("Test patient", "One of the test patients (its true stage is known, so the prediction can be checked). It is not used to calibrate its own prediction set."),
                           choices = data$holdout_ids),
            actionButton("p_sp_predict", "Predict this patient", class = "btn-primary w-100",
                         icon = icon("wand-magic-sparkles"))
          ),
          uiOutput("p_sp_result")
        )
      )
    )
  )
}

# Bar plot of one patient with the conformal inclusion threshold
plot_patient <- function(pred) {
  probs <- data.frame(Stage = factor(stages, levels = stages), Probability = as.numeric(unlist(pred[1, stages])))
  set <- strsplit(pred[["Prediction set"]][1], " / ", fixed = TRUE)[[1]]
  probs$In_set <- as.character(probs$Stage) %in% set
  threshold <- pred[["Set threshold"]][1]
  level <- pred[["Coverage target"]][1]
  subtitle <- if (sum(probs$In_set) == 1) sprintf("Prediction set (%s coverage): %s", level, set)
    else if (sum(probs$In_set) == 0) "No stage reaches the threshold: the patient looks unlike the calibration patients"
    else sprintf("The model cannot choose between %s at %s coverage", paste(set, collapse = " and "), level)
  ggplot(probs, aes(Stage, Probability, fill = In_set)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = scales::label_percent(accuracy = 0.1)(Probability)), vjust = -0.4, size = 4.5) +
    geom_hline(yintercept = threshold, linetype = "dashed", colour = "#A65028") +
    annotate("text", x = 3.45, y = threshold, vjust = -0.5, hjust = 1, size = 3.8, colour = "#A65028",
             label = sprintf("inclusion threshold (%s sets)", level)) +
    scale_fill_manual(values = c("TRUE" = "#176084", "FALSE" = "#c5ced6"), guide = "none") +
    scale_y_continuous(limits = c(0, 1.1), labels = scales::label_percent()) +
    labs(title = "Predicted stage probabilities", subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 13)
}

predict_server <- function(input, output, session, data) {
  # Settings of this page
  p_settings <- reactive({
    dataset <- req(input$p_dataset)
    layers <- dataset_layers[[dataset]]
    integration <- if (dataset %in% multi_datasets) req(input$p_integration) else "none"
    model <- if (integration == "block") "rf" else req(input$p_model)
    list(dataset = dataset, integration = integration, model = model,
         rna = if ("RNA" %in% layers) input$p_rna else NULL,
         protein = if ("Protein" %in% layers) input$p_protein else NULL,
         methylation = if ("Methylation" %in% layers) input$p_methylation else NULL,
         conformal = TRUE)
  })

  # "DIABLO (prediction)" only with the DIABLO integration
  observe({
    diablo_ok <- isTRUE(input$p_dataset %in% multi_datasets) && identical(input$p_integration, "diablo")
    choices <- if (diablo_ok) names(model_labels) else setdiff(names(model_labels), "diablo")
    selected <- if (isolate(input$p_model) %in% choices) isolate(input$p_model) else "rf"
    updateSelectInput(session, "p_model", choices = setNames(choices, model_labels[choices]),
                      selected = selected)
  })
  p_features <- reactive(tryCatch(needed_features(data, p_settings()), error = function(e) NULL))

  output$p_needed <- renderUI({
    layers <- dataset_layers[[req(input$p_dataset)]]
    helpText(tags$b("Files needed: "), paste(layer_labels[layers], collapse = ", "),
             ". A patient must be present in every file. The data must be complete (the app does not impute).")
  })
  show_error <- function(msg) {
    lines <- strsplit(msg, "\n", fixed = TRUE)[[1]]
    tags$div(class = "alert alert-danger", tags$b(lines[1]),
             if (length(lines) > 1) tags$ul(lapply(lines[-1], tags$li)))
  }
  run_prediction <- function(new_layers) {
    tryCatch(
      withProgress(message = "Training the model and predicting", value = 0.1, {
        score_new_patients(data, p_settings(), new_layers, level = input$p_conf_level)
      }),
      error = function(e) data.frame(error = conditionMessage(e)))
  }

  # ---- Single patient -----------------------------------------------------------
  sp_pred <- reactiveVal(NULL)
  observeEvent(input$p_sp_predict, {
    patient <- req(input$p_sp_patient)
    f <- req(p_features())
    layers <- lapply(setNames(names(f), names(f)), function(l) {
      data$holdout[[l]][patient, f[[l]], drop = FALSE]
    })
    pred <- run_prediction(layers)
    if (is.null(pred$error)) pred$Observed <- as.character(data$y_holdout[match(patient, data$holdout_ids)])
    sp_pred(pred)
  })
  output$p_sp_result <- renderUI({
    pred <- sp_pred()
    if (is.null(pred)) return(NULL)
    if (!is.null(pred$error)) return(show_error(pred$error[1]))
    tagList(
      bslib::layout_columns(
        bslib::value_box("Predicted stage", sprintf("%s (%s)", pred$Predicted[1],
                                                    scales::label_percent(accuracy = 0.1)(pred$Probability[1])), fill = FALSE),
        bslib::value_box("True stage", pred$Observed[1], fill = FALSE),
        bslib::value_box(with_info("Uncertainty", "From the size of the conformal prediction set: one stage = low, two = medium, three = high; an empty set means the patient looks unlike the calibration patients."),
                         pred$Uncertainty[1], fill = FALSE),
        bslib::value_box(with_info("Prediction set", "Stages compatible with this patient's data at the target coverage."),
                         pred[["Prediction set"]][1], fill = FALSE),
        fill = FALSE),
      plotOutput("p_sp_plot", height = "300px"),
      tags$hr())
  })
  output$p_sp_plot <- renderPlot({
    pred <- req(sp_pred()); req(is.null(pred$error))
    plot_patient(pred)
  })
}
