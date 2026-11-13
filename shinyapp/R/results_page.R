# =============================================================================
#  results_page.R - "All results" tab: cross-validation and test results of
#  the models, computed beforehand with 06_precompute_results.R
# =============================================================================

results_file <- "results/all_results.rds"

results_page <- function() {
  tagList(
    bslib::card(
      class = "run-card",
      bslib::card_body(
        tags$p(class = "stage-kicker mb-1", "TCGA-BRCA  /  ALL RESULTS"),
        tags$h3(class = "mb-2", "Which model works best?"),
        tags$p(class = "stage-note mb-0",
               "Every model was evaluated twice: 5-fold cross-validation on the training patients, and the test patients (never used to fit or choose the model).")
      )
    ),
    uiOutput("results_summary"),
    bslib::card(bslib::card_body(
      tags$div(class = "d-flex flex-wrap gap-5 align-items-end",
               sliderInput("results_n", "Models shown", min = 1, max = 30, value = 5, step = 1,
                           ticks = FALSE, width = "320px"),
               radioButtons("results_order", "Best models by", inline = TRUE,
                            choices = c("Training (cross-validation)" = "cv", "Test patients" = "test"),
                            selected = "cv")),
      plotOutput("results_plot", height = "auto"),
      tags$p(class = "stage-note mb-0",
             "Models are compared by cross-validation; test results are shown for reference only."))),
    bslib::card(bslib::card_body(
      fillable = FALSE,
      tags$div(class = "download-row mb-2",
               downloadButton("dl_results", "Download table (CSV)", class = "btn-sm")),
      tags$p(class = "stage-note",
             "The Settings column shows the settings used by each model."),
      DT::DTOutput("results_table")
    ))
  )
}

# Settings of a model as values only (tuned values, or the fixed workbench settings)
settings_text <- function(r) {
  tuned <- grepl(": ", r$Settings)
  fixed <- ifelse(r$Model == "Random forest", "1000 trees, node size 10",
           ifelse(r$Model == "XGBoost", "depth 2, 100 rounds",
           ifelse(r$Model == "Lasso", ifelse(r$Integration == "DIABLO", "alpha 1, lambda 0.03", "alpha 1, lambda by 5-fold CV"),
           ifelse(r$Model == "Elastic net (glmnet)", ifelse(r$Integration == "DIABLO", "alpha 0.5, lambda 0.03", "alpha 0.5, lambda by 5-fold CV"),
           ifelse(r$Model == "Block Forest", "500 trees per forest, 20 tuning sets",
           ifelse(r$Model == "DIABLO (prediction)", "2 components, keepX RNA 50/100, protein 5/40, methylation 25/100",
                  r$Settings))))))
  ifelse(tuned, sub(".*: ", "", r$Settings), fixed)
}

results_server <- function(input, output, session) {
  results <- reactive({
    if (!file.exists(results_file)) return(NULL)
    r <- readRDS(results_file)
    r[is.na(r$error) & sub("\\|.*", "", r$id) %in% names(dataset_layers), , drop = FALSE]
  })
  label <- function(r) paste(r$Model, "·", r$Data, if (r$Integration != "-") paste("·", r$Integration) else "")

  output$results_summary <- renderUI({
    r <- results()
    if (is.null(r) || !nrow(r)) {
      return(tags$div(class = "alert alert-info",
                      "No results yet. Run 06_precompute_results.R (see the README) to fill this tab."))
    }
    best_cv <- r[which.max(r$cv_balanced), ]
    best_test <- r[which.max(r$test_balanced), ]
    bslib::layout_columns(
      bslib::value_box(with_info("Chosen by cross-validation", "The model with the highest mean balanced accuracy in the 5-fold cross-validation of the training patients. Its test result shows how well it works on patients that were not used to choose it."),
                       label(best_cv), tags$p(sprintf("CV %s ± %s · test %s",
                                                      format_rate(best_cv$cv_balanced), format_rate(best_cv$cv_balanced_sd),
                                                      format_rate(best_cv$test_balanced))), fill = FALSE),
      bslib::value_box(with_info("Best on the test patients", "Only for reference: picking the model with the best test result uses the test patients to choose, so this value is optimistic."),
                       label(best_test), tags$p(sprintf("test %s · CV %s", format_rate(best_test$test_balanced),
                                                        format_rate(best_test$cv_balanced))), fill = FALSE),
      bslib::value_box("Models compared", nrow(r), tags$p("Balanced accuracy by chance: 33.3%"), fill = FALSE),
      col_widths = c(5, 5, 2), fill = FALSE)
  })

  observe(updateSliderInput(session, "results_n", max = nrow(req(results()))))
  shown_n <- reactive(max(1, min(nrow(req(results())), as.integer(input$results_n %||% 5), na.rm = TRUE)))
  output$results_plot <- renderPlot({
    r <- req(results()); req(nrow(r) > 0)
    by_cv <- identical(input$results_order, "cv")
    r <- head(r[order(-(if (by_cv) r$cv_balanced else r$test_balanced)), , drop = FALSE], shown_n())
    r$Label <- vapply(seq_len(nrow(r)), function(i) label(r[i, ]), character(1))
    r$Label <- factor(r$Label, levels = rev(r$Label))                  # best model on top
    long <- rbind(
      data.frame(Label = r$Label, Data = r$Data, Set = "Training (5-fold CV)",
                 Value = r$cv_balanced, SD = r$cv_balanced_sd),
      data.frame(Label = r$Label, Data = r$Data, Set = "Test patients",
                 Value = r$test_balanced, SD = NA))
    long$Set <- factor(long$Set, levels = c("Training (5-fold CV)", "Test patients"))
    ggplot(long, aes(Value, Label, fill = Data)) +
      geom_col(width = 0.65) +
      geom_errorbar(aes(xmin = Value - SD, xmax = Value + SD), width = 0.3, colour = "grey30",
                    orientation = "y", na.rm = TRUE) +
      geom_vline(xintercept = 1 / 3, linetype = "dashed", colour = "grey50") +
      facet_wrap(~ Set) +
      scale_fill_manual(values = c("RNA" = "#176084", "Protein" = "#A65028", "DNA methylation" = "#22765B",
                                   "RNA + protein + methylation" = "#7A4E9C")) +
      scale_x_continuous(labels = scales::label_percent()) +
      coord_cartesian(xlim = c(0, 0.75)) +
      labs(title = "Balanced accuracy of every model",
           subtitle = sprintf("Ordered by the %s · error bars: ± 1 SD over the 5 folds · dashed line: chance (33%%)",
                              if (by_cv) "cross-validation result" else "test result"),
           x = NULL, y = NULL, fill = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
            strip.text = element_text(size = 15, face = "bold"))
  }, height = function() 160 + 45 * shown_n())

  results_display <- reactive({
    r <- req(results())
    data.frame(Origin = r$Origin, Model = r$Model, Data = r$Data, Integration = r$Integration,
               Features = r$Features, Settings = settings_text(r),
               `CV balanced accuracy` = sprintf("%s ± %s", format_rate(r$cv_balanced), format_rate(r$cv_balanced_sd)),
               `CV folds` = r$cv_folds,
               `Test balanced accuracy` = format_rate(r$test_balanced),
               `Test accuracy` = format_rate(r$test_accuracy), `Test macro F1` = format_rate(r$test_macro_f1),
               `Test mean AUC` = round(r$test_auc, 3), `Conformal coverage` = format_rate(r$coverage),
               `Mean set size` = round(r$set_size, 2), check.names = FALSE)[order(-(if (identical(input$results_order, "cv")) r$cv_balanced else r$test_balanced)), ]
  })
  output$results_table <- DT::renderDT({
    DT::datatable(results_display(), rownames = FALSE,
                  options = list(pageLength = 5, lengthMenu = c(3, 5, 10, 25), scrollX = TRUE))
  }, server = FALSE)
  output$dl_results <- downloadHandler(
    filename = function() "all_results.csv",
    content = function(file) utils::write.csv(req(results()), file, row.names = FALSE))
}
