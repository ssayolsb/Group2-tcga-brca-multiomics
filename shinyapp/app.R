suppressPackageStartupMessages(suppressWarnings({
  library(shiny)
  library(ggplot2)
}))
source("modeling.R", local = TRUE, encoding = "UTF-8")   # UTF-8 symbols also on Windows
# Helpers ((i) icons, dataset choices), About page and Predict new patients page
for (f in c("R/app_helpers.R", "R/about_page.R", "R/predict_page.R", "R/results_page.R")) {
  source(f, local = TRUE, encoding = "UTF-8")
}
stage_data <- load_stage_data()
format_rate <- scales::label_percent(accuracy = 0.1)
feature_limits <- lengths(stage_data$ranked)
holdout_n <- length(stage_data$y_holdout)
calibration_n <- length(stage_data$calibration)

workbench_ui <- bslib::layout_sidebar(
  sidebar = bslib::sidebar(
    title = "Model setup",
    width = 330,
    selectInput("dataset", with_info("Training data", help_text$dataset), choices = dataset_choices),
    conditionalPanel(
      multi_dataset_js,
      selectInput("integration", with_info("Integration method", help_text$integration),
                  choices = c("Join directly" = "join", "Block Forest" = "block",
                              "DIABLO" = "diablo"))
    ),
    conditionalPanel(
      sprintf("!(%s && input.integration === 'block')", multi_dataset_js),
      selectInput("model", with_info("Classifier", help_text$model),
                  choices = setNames(names(model_labels), model_labels))
    ),
    conditionalPanel(
      sprintf("%s && input.integration === 'block'", multi_dataset_js),
      tags$div(class = "alert alert-info py-2",
               tags$strong("Classifier"), tags$br(), "Random forest (Block Forest engine)")
    ),
    conditionalPanel(
      sprintf("%s && input.integration === 'diablo'", multi_dataset_js),
      helpText("DIABLO learns two components per layer from the training patients. 'DIABLO (prediction)' lets the sparse DIABLO model predict the stage itself.")
    ),
    conditionalPanel(
      layer_js("RNA"),
      numericInput("rna", with_info("Number of RNA genes", help_text$counts), min = 2, max = feature_limits[["RNA"]],
                   value = min(216, feature_limits[["RNA"]]), step = 1)
    ),
    conditionalPanel(
      layer_js("Protein"),
      numericInput("protein", with_info("Number of protein markers", help_text$counts), min = 2, max = feature_limits[["Protein"]],
                   value = feature_limits[["Protein"]], step = 1)
    ),
    conditionalPanel(
      layer_js("Methylation"),
      numericInput("methylation", with_info("Number of methylation regions", help_text$counts),
                   min = 2, max = feature_limits[["Methylation"]],
                   value = min(216, feature_limits[["Methylation"]]), step = 1)
    ),
    helpText("Features are ranked by variance in the training data."),
    tags$hr(),
    checkboxInput("conformal", with_info("Add conformal prediction", help_text$conformal), value = FALSE),
    uiOutput("split_note"),
    actionButton("run", "Train and test", class = "btn-primary w-100")
  ),
  tags$style(HTML("
    .run-card { margin-bottom: 1rem; }
    .stage-kicker { font-size: .75rem; font-weight: 700; letter-spacing: .09em; color: #176084; }
    .stage-note { color: #4a6174; }
    .result-tab { padding: .6rem .15rem 0; }
    .metric-row { margin-bottom: 1rem; }
    .summary-table { margin-top: .8rem; }
  ")),
  tags$style(HTML("
    .info-icon { color: #176084; cursor: pointer; font-size: .9em; }
    .about-card h5 { color: #176084; }
    .download-row .btn { margin: 0 .4rem .4rem 0; }
  ")),
  bslib::card(
    class = "run-card",
    bslib::card_body(
      tags$p(class = "stage-kicker mb-1", "TCGA-BRCA  /  MODEL WORKBENCH"),
      tags$h3(class = "mb-2", "Predict pathological stage from tumour assays"),
      tags$p(class = "stage-note mb-0",
             "Train on the saved RNA, protein and DNA methylation data and evaluate Stage I, Stage II, and Stage III/IV.")
    )
  ),
  uiOutput("run_status"),
  uiOutput("result_panels"),
  uiOutput("compare_runs")
)

ui <- bslib::page_navbar(
  title = "TCGA-BRCA stage prediction",
  theme = bslib::bs_theme(version = 5, bg = "#f4f7fa", fg = "#19344a", primary = "#176084"),
  fillable = FALSE,
  bslib::nav_panel("About", icon = icon("house"), about_page(stage_data)),
  bslib::nav_panel("Model workbench", icon = icon("flask"), workbench_ui),
  bslib::nav_panel("Predict new patients", icon = icon("user-plus"), predict_page(stage_data)),
  bslib::nav_panel("All results", icon = icon("ranking-star"), results_page()),
  bslib::nav_spacer(),
  bslib::nav_item(tags$a(href = repo_url, target = "_blank", icon("github"), " Code"))
)

workbench_server <- function(input, output, session) {
  result <- reactiveVal(NULL)
  run_error <- reactiveVal(NULL)
  runs <- reactiveVal(NULL)   # one row per run of this session ("Compare runs")

  settings <- reactive({
    dataset <- req(input$dataset)
    integration <- if (dataset %in% multi_datasets) req(input$integration) else "none"
    layers <- dataset_layers[[dataset]]
    model <- if (integration == "block") "rf" else req(input$model)
    list(dataset = dataset, integration = integration, model = model,
         rna = if ("RNA" %in% layers) input$rna else NULL,
         protein = if ("Protein" %in% layers) input$protein else NULL,
         methylation = if ("Methylation" %in% layers) input$methylation else NULL,
         conformal = isTRUE(input$conformal))
  })

  output$split_note <- renderUI({
    if (isTRUE(input$conformal)) {
      helpText(sprintf("%d patients for calibration and %d for final testing.",
                       calibration_n, holdout_n - calibration_n))
    } else {
      helpText(sprintf("All %d holdout patients are used for testing.", holdout_n))
    }
  })

  # Capture a configuration only when the button is clicked.
  observeEvent(input$run, {
    config <- settings()
    result(NULL)
    run_error(NULL)
    tryCatch(
      withProgress(message = "Training and evaluating", value = 0, {
        answer <- run_analysis(stage_data, config, progress = function(detail) {
          incProgress(0.18, detail = detail)
        })
        result(answer)
        runs(rbind(runs(), run_summary(answer, NROW(runs()) + 1)))
      }),
      error = function(error) run_error(conditionMessage(error))
    )
  }, ignoreInit = TRUE)

  output$run_status <- renderUI({
    error <- run_error()
    if (!is.null(error)) return(tags$div(class = "alert alert-danger", role = "alert", error))
    run <- result()
    if (is.null(run)) {
      return(bslib::card(
        bslib::card_body(
          tags$h5("Ready to train"),
          tags$p(class = "stage-note mb-0",
                 sprintf("%d training patients and %d holdout patients are available.",
                         length(stage_data$y_train), holdout_n))
        )
      ))
    }

    config <- run$settings
    title <- if (config$integration == "block") "Block Forest" else model_labels[[config$model]]
    if (config$integration %in% c("join", "diablo")) {
      title <- paste(title, "with", integration_labels[[config$integration]])
    }
    counts <- paste(sprintf("%s %d", names(run$features), run$features), collapse = " · ")
    tagList(
      if (!identical(settings(), config)) {
        tags$div(class = "alert alert-info",
                 "Settings have changed. Train and test again to evaluate the new configuration.")
      },
      bslib::card(
        class = "run-card",
        bslib::card_body(
          tags$p(class = "stage-kicker mb-1", "LAST RUN"),
          tags$h5(class = "mb-1", paste(title, "·", dataset_labels[[config$dataset]])),
          tags$p(class = "mb-1", sprintf("%d trained · %d calibrated · %d tested",
                                        run$n_train, run$n_cal, run$n_test)),
          tags$p(class = "stage-note mb-0", paste(counts, "·", run$model_inputs,
                                                 "classifier inputs ·", run$details))
        )
      )
    )
  })

  output$result_panels <- renderUI({
    run <- req(result())
    tabs <- list(
      bslib::nav_panel(
        "Performance",
        tags$div(class = "result-tab",
          bslib::layout_columns(
            bslib::value_box("Accuracy", format_rate(run$accuracy), fill = FALSE),
            bslib::value_box("Balanced accuracy", format_rate(run$balanced), fill = FALSE),
            bslib::value_box("Macro F1", format_rate(run$macro_f1), fill = FALSE),
            col_widths = c(4, 4, 4), fill = FALSE, class = "metric-row"
          ),
          tags$p(class = "stage-note",
                 paste0("The majority-stage baseline predicts ", run$majority,
                        " for every patient and reaches ", format_rate(run$baseline), " accuracy.")),
          plotOutput("confusion", height = "355px"),
          tags$div(class = "summary-table", tableOutput("class_metrics")),
          tags$hr(),
          tags$p(class = "stage-note",
                 "Each ROC curve compares one stage with all other stages as the score threshold varies. The legend gives the area under each curve (AUC)."),
          plotOutput("roc", height = "390px")
        )
      ),
      bslib::nav_panel(
        "Patient predictions",
        tags$div(class = "result-tab",
                 tags$p(class = "stage-note", "One row per test patient. Stage columns show model scores."),
                 DT::DTOutput("predictions"))
      )
    )
    if (!is.null(run$conformal)) {
      tabs <- c(tabs, list(bslib::nav_panel(
        "Conformal prediction",
        tags$div(class = "result-tab",
          tags$p("Least Ambiguous Classifier (LAC) compares each candidate stage score with the calibration scores. A stage enters the set when its conformal p-value exceeds 0.10."),
          tags$p(class = "stage-note",
                 "The 90% target applies across patients as a group when calibration and test patients are exchangeable."),
          bslib::layout_columns(
            bslib::value_box("Observed coverage", format_rate(mean(run$conformal$covered)),
                             fill = FALSE),
            bslib::value_box("Mean stages per set", sprintf("%.2f", mean(run$conformal$size)),
                             fill = FALSE),
            col_widths = c(6, 6), fill = FALSE, class = "metric-row"
          ),
          bslib::layout_columns(
            plotOutput("coverage", height = "295px"),
            plotOutput("set_size", height = "295px"),
            col_widths = c(6, 6), fill = FALSE
          )
        )
      )))
    }
    do.call(bslib::navset_card_tab, tabs)
  })

  output$confusion <- renderPlot({
    cm <- req(result())$confusion
    cells <- as.data.frame(cm, responseName = "Count")
    cells$Fraction <- as.vector(prop.table(cm, margin = 1))
    ggplot(cells, aes(Predicted, Observed, fill = Fraction)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = Count, colour = Fraction > 0.60), size = 5) +
      scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "#19344a"), guide = "none") +
      scale_fill_gradient(low = "#eaf3f8", high = "#176084", limits = c(0, 1),
                          labels = scales::label_percent(), name = "Within stage") +
      scale_y_discrete(limits = rev(stages)) +
      labs(title = "Confusion matrix", x = "Predicted stage", y = "Observed stage") +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank())
  })

  output$class_metrics <- renderTable({
    metrics <- req(result())$stage_metrics
    for (metric in c("Recall", "Precision", "F1")) {
      values <- metrics[[metric]]
      metrics[[metric]] <- ifelse(is.finite(values), format_rate(values), "—")
    }
    metrics
  }, width = "100%")

  output$roc <- renderPlot({
    run <- req(result())
    curves <- run$roc
    auc <- vapply(curves, function(curve) as.numeric(pROC::auc(curve)), numeric(1))
    names(curves) <- sprintf("%s (AUC %.3f)", names(curves), auc)

    pROC::ggroc(curves, legacy.axes = TRUE, linewidth = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
      scale_colour_manual(values = c("#176084", "#A65028", "#22765B")) +
      coord_equal() +
      labs(title = "Held-out ROC curves",
           subtitle = sprintf("One-vs-rest on %d test patients", run$n_test),
           x = "False positive rate", y = "True positive rate", colour = NULL) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
  })

  output$predictions <- DT::renderDT({
    DT::datatable(req(result())$predictions, rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  }, server = FALSE)

  output$coverage <- renderPlot({
    coverage <- req(req(result())$conformal)$coverage
    coverage$Group <- factor(coverage$Group, levels = rev(coverage$Group))
    ggplot(coverage, aes(Coverage, Group)) +
      geom_vline(xintercept = 0.90, colour = "grey50", linetype = "dashed") +
      geom_linerange(aes(xmin = Lower, xmax = Upper), colour = "#176084") +
      geom_point(colour = "#176084", size = 2.8) +
      scale_x_continuous(limits = c(0, 1), labels = scales::label_percent()) +
      labs(title = "Coverage by observed stage", subtitle = "Exact 95% binomial intervals",
           x = "Observed coverage", y = NULL) +
      theme_minimal(base_size = 12)
  })

  output$set_size <- renderPlot({
    size <- req(req(result())$conformal)$size
    sizes <- data.frame(Size = factor(size, levels = 0:3))
    ggplot(sizes, aes(Size)) +
      geom_bar(fill = "#176084", width = 0.65) +
      scale_x_discrete(drop = FALSE) +
      labs(title = "Prediction set sizes", x = "Stages in set", y = "Test patients") +
      theme_minimal(base_size = 12)
  })

  # "DIABLO (prediction)" only with the DIABLO integration
  observe({
    diablo_ok <- isTRUE(input$dataset %in% multi_datasets) && identical(input$integration, "diablo")
    choices <- if (diablo_ok) names(model_labels) else setdiff(names(model_labels), "diablo")
    selected <- if (isolate(input$model) %in% choices) isolate(input$model) else "rf"
    updateSelectInput(session, "model", choices = setNames(choices, model_labels[choices]),
                      selected = selected)
  })

  output$compare_runs <- renderUI({
    req(runs())
    bslib::card(
      class = "mt-3",
      bslib::card_header("Compare runs (this session)"),
      bslib::card_body(
        tags$p(class = "stage-note",
               "Every run is added here. All runs are tested on holdout patients; with conformal prediction on, 40% of them are used for calibration and the rest for testing."),
        plotOutput("runs_plot", height = "300px"),
        DT::DTOutput("runs_table"),
        tags$div(class = "download-row mt-2",
                 downloadButton("dl_runs", "Download table (CSV)", class = "btn-sm"),
                 actionButton("clear_runs", "Clear table", class = "btn-sm btn-outline-secondary"))
      )
    )
  })
  output$runs_plot <- renderPlot({
    h <- req(runs())
    labels <- paste0(h$Run, ". ", h$Model, ifelse(h$Integration == "-", "", paste0(" (", h$Integration, ")")))
    h$Label <- factor(labels, levels = rev(labels))
    ggplot(h, aes(`Balanced accuracy`, Label, fill = Data)) +
      geom_col(width = 0.6) +
      geom_vline(xintercept = 1 / 3, linetype = "dashed", colour = "grey50") +
      scale_fill_manual(values = c("RNA" = "#176084", "Protein" = "#A65028", "DNA methylation" = "#22765B",
                                   "RNA + protein + methylation" = "#7A4E9C")) +
      scale_x_continuous(limits = c(0, 1), labels = scales::label_percent()) +
      labs(title = "Balanced accuracy on the test patients", subtitle = "Dashed line: chance level",
           x = NULL, y = NULL, fill = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
  })
  output$runs_table <- DT::renderDT({
    DT::datatable(req(runs()), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  }, server = FALSE)
  output$dl_runs <- downloadHandler(
    filename = function() "compare_runs.csv",
    content = function(file) utils::write.csv(req(runs()), file, row.names = FALSE))
  observeEvent(input$clear_runs, runs(NULL))
}

server <- function(input, output, session) {
  workbench_server(input, output, session)
  predict_server(input, output, session, stage_data)
  results_server(input, output, session)
}

shinyApp(ui, server)
