# =============================================================================
#  about_page.R - About page: what the app does, data, methods, team, links
# =============================================================================

about_page <- function(data) {
  card <- function(title, ...) bslib::card(class = "about-card", bslib::card_body(tags$h5(title), ...))
  n_feat <- lengths(data$ranked)
  tagList(
    tags$div(class = "alert alert-warning",
             icon("triangle-exclamation"), tags$b(" Not for clinical use. "),
             "This is a course project for teaching and research."),
    bslib::card(
      class = "run-card",
      bslib::card_body(
        tags$p(class = "stage-kicker mb-1", "TCGA-BRCA  /  ABOUT THIS APP"),
        tags$h3(class = "mb-2", "Predict the pathological stage of breast tumours from omics data"),
        tags$p(class = "stage-note mb-0",
               "Train machine-learning models on RNA, protein and DNA methylation data of TCGA-BRCA tumours, evaluate them on held-out patients, and predict new patients with the uncertainty of every prediction.")
      )
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      card("What the app does",
           tags$ul(
             tags$li(tags$b("Model workbench: "), "choose the data, the integration method and the classifier. The model is trained on ",
                     length(data$y_train), " training patients and tested on ", length(data$y_holdout), " test patients."),
             tags$li(tags$b("Uncertainty: "), "optional conformal prediction gives, for each patient, the set of stages compatible with 90% coverage."),
             tags$li(tags$b("Predict new patients: "), "choose a test patient and get its predicted stage, with its uncertainty, next to its true stage. Uploading a data set of new patients will be available in a future version.")
           )),
      card("How to use it",
           tags$ol(
             tags$li("Open ", tags$b("Model workbench"), ", choose the settings (each one has an (i) icon) and click ", tags$b("Train and test"), "."),
             tags$li("Open ", tags$b("Predict new patients"), " to predict the stage of one patient (uploading a data set: coming soon)."),
             tags$li("Open ", tags$b("All results"), " to compare all the models: cross-validation on the training patients and results on the test patients.")
           )),
      card("Data",
           tags$p("TCGA-BRCA tumours split once into training (70%) and test (30%) patients. All preprocessing decisions were made on the training patients and then applied to the test patients."),
           tags$ul(
             tags$li(tags$b("RNA: "), "RNA-seq counts normalised as log2(CPM + 1), filtered in the quality control (genes with zero variance, mostly zeros or very low expression removed). The app ranks the genes by their variance in the training patients and uses the most variable ones (",
                     n_feat[["RNA"]], " available)."),
             tags$li(tags$b("Protein: "), "RPPA values (log2, median-centred) of the ", n_feat[["Protein"]],
                     " proteins kept by the quality control."),
             tags$li(tags$b("DNA methylation: "), "450k beta values normalised with BMIQ, averaged by region (promoter, gene body, CpG island) and converted to M-values. The app uses the most variable regions in the training patients (",
                     n_feat[["Methylation"]], " available).")
           ),
           tags$p(tags$b("Missing values: "), "the RNA data have none. In the protein and methylation data, features with many missing values were removed during quality control, and any value still missing is replaced by the median of the training patients. New patients are not imputed, so their data must be complete."),
           tags$p(class = "stage-note mb-0", "Outcome: AJCC pathological stage I, II and III/IV.")),
      card("Methods",
           tags$ul(
             tags$li("Features ranked by training variance in each layer."),
             tags$li("Integration of the three layers: direct join, Block Forest or DIABLO (components as classifier inputs, or DIABLO's own prediction)."),
             tags$li("Classifiers: random forest, XGBoost, Lasso, elastic net and DIABLO prediction. Class weights for the imbalanced stages (except DIABLO prediction)."),
             tags$li("Evaluation: 5-fold cross-validation on the training patients (All results) and, on the test patients, accuracy, balanced accuracy, macro F1, per-stage metrics and ROC curves; conformal prediction sets (LAC).")
           )),
      card("Team",
           tags$p("Group 2: Multiomics - Parniyan, Sergi, Baizhen, Paula, Yueling."),
           tags$p(class = "mb-0", "Machine Learning for Life Sciences, group project.")),
      card("Code and documentation",
           tags$p("All the code is open on GitHub, with the steps to reproduce the analysis."),
           tags$p(icon("github"), " ", tags$a(href = repo_url, target = "_blank", repo_url)),
           tags$p(class = "stage-note mb-0", "Licence: MIT. Developed with AI-assisted coding and writing tools."))
    )
  )
}
