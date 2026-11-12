# Shiny app: TCGA-BRCA stage prediction

Interactive app to train and compare stage prediction models on RNA, protein
and DNA methylation data, and to predict the pathological stage of a test
patient with the uncertainty of the prediction.

> Course project for teaching and research. Not for clinical use.

## Files

| File | What it contains |
|---|---|
| `app.R` | Model workbench (original app) and the navigation bar with all the tabs |
| `modeling.R` | Data loading and models: original code plus DNA methylation, the three layers together, DIABLO prediction and prediction of single patients with uncertainty |
| `R/app_helpers.R` | (i) help texts, dataset conditions, "Compare runs" rows |
| `R/about_page.R` | About tab |
| `R/predict_page.R` | Predict new patients tab |
| `R/results_page.R` | All results tab |
| `05_make_app_data.R` | Makes `data/`, the reduced copy of the project data that the app reads |
| `data/` | Output of `05_make_app_data.R` (RNA: 5000 most variable genes; protein; DNA methylation: 10000 most variable regions; stage) |
| `06_precompute_results.R` | Cross-validation and test results of all models (for the All results tab) |
| `results/all_results.rds` | Output of `06_precompute_results.R` |

## Tabs

- **About** – what the app does, data, methods, team and links.
- **Model workbench** – choose the data (RNA, protein or DNA methylation alone,
  or the three together), the number of features, the integration method
  (direct join, Block Forest, DIABLO) and the classifier (random forest,
  XGBoost, Lasso, elastic net, DIABLO prediction). The model is trained on the
  training patients and evaluated on the test patients. Option: conformal
  prediction. Every run is added to "Compare runs" (table and bar plot).
- **Predict new patients** – set up a model, choose a test patient and get its
  predicted stage next to its true stage, with the probability of each stage,
  the conformal prediction set and an uncertainty level. Uploading a data set
  of new patients: coming soon.
- **All results** – 5-fold cross-validation (training patients) and test
  results of all models, computed beforehand with `06_precompute_results.R`.

## Install and run

Package versions are recorded with `renv`. The first time, in the `shinyapp/`
project:

```r
renv::restore()
```

Run the app from the repository root:

```r
shiny::runApp("shinyapp")
```

The app reads its own reduced data in `data/` (made by `05_make_app_data.R`
from the project data). If `data/` does not exist, it reads the project data in
`../data/` (not in the repository; see the main README). To recompute the All
results tab: `source("shinyapp/06_precompute_results.R")`.
