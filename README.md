# Group2-tcga-brca-multiomics

**Group 2: Multiomics** · Machine Learning for Life Sciences · Parniyan, Sergi, Baizhen, Paula, Yueling

Prediction of the pathological stage (I, II, III/IV) of breast tumours from
multi-omics data of TCGA-BRCA (RNA-seq, RPPA protein and DNA methylation), and
a Shiny app to train and compare models and to predict new patients with the
uncertainty of every prediction.

- **Code:** <https://github.com/ssayolsb/Group2-tcga-brca-multiomics>
- **Licence:** MIT (see [LICENSE](LICENSE))

> Course project for teaching and research. The models are not validated for
> clinical use.

---

## Analysis steps

The files are numbered in the order of the analysis; each report in
`reports/` has the number of the code that produces it.

| Step | Code | Report | Output |
|---|---|---|---|
| 01 Data preparation | `01_data_preprocess/01_data_preparation.R` / `.Rmd` | `reports/01_data_preparation.html` | `data/data_split/*.train.rds`, `*.test.rds` (70/30 split) |
| 02 Quality control (training set) | `01_data_preprocess/02_qc_rnaseq.R` / `.Rmd`<br>`02_qc_prot.R` / `.Rmd`<br>`02_qc_meth.R` / `.Rmd` | `reports/02_qc_rnaseq.html`<br>`reports/02_qc_prot.html`<br>`reports/02_qc_meth.html` | `data/processed/*.train.*`, `*_qc_params.rds` |
| 03 Same QC applied to the test set | `01_data_preprocess/03_apply_preprocessing.R` | – | `data/processed/*.test.*`, `test_qc_log.csv` |
| 04 Models: XGBoost and elastic net | `02_ML_methods/xgboost_elastic_net/` | `reports/04_xgboost_elastic_net.html` | Cross-validated settings, test results, conformal sets |
| 04 Models: DIABLO | `02_ML_methods/diablo/` | `reports/04_diablo.html` | Tuned sparse DIABLO, test results |
| 05 Reduced data for the app | `shinyapp/05_make_app_data.R` | – | `shinyapp/data/` (the files the app reads; in the repository) |
| 06 Results of all the app models | `shinyapp/06_precompute_results.R` | "All results" tab of the app | `shinyapp/results/all_results.rds` |
| App | `shinyapp/` | – | Shiny app |

**General rule:** every data-driven step (QC filters, feature ranking, model
fitting and tuning) is done on the training patients only and then applied
unchanged to the test patients and to new patients.

## Folder structure

```
Group2-tcga-brca-multiomics/
├── 01_data_preprocess/     data preparation, split and QC (scripts + Rmd reports)
├── 02_ML_methods/          one folder per method (source of the model reports)
│   ├── xgboost_elastic_net/
│   └── diablo/
├── reports/                rendered HTML reports, numbered like the code
├── shinyapp/               Shiny app, results of all models, own renv environment
└── data/                   NOT in the repository (see below)
```

## How to reproduce

Open `Group2-tcga-brca-multiomics.Rproj` in RStudio (paths are relative to the
project root) and run in this order:

1. `01_data_preprocess/01_data_preparation.R` – download, inclusion criteria and split.
2. `01_data_preprocess/02_qc_rnaseq.R`, `02_qc_prot.R`, `02_qc_meth.R` – QC on the
   training set (the BMIQ normalisation of the methylation data takes several hours).
   The `.Rmd` versions produce the reports and reuse the saved BMIQ result:
   `rmarkdown::render("01_data_preprocess/02_qc_prot.Rmd", output_dir = "reports")`.
3. `01_data_preprocess/03_apply_preprocessing.R` – the same QC applied to the test set.
4. The model reports in `02_ML_methods/` (see the README of each folder).
5. `source("shinyapp/05_make_app_data.R")` – reduced copy of the data used by the app (`shinyapp/data/`).
6. `source("shinyapp/06_precompute_results.R")` – cross-validation and test results
   of all the app models.
7. `shiny::runApp("shinyapp")` – the app.

## Data

The data are public (TCGA-BRCA, Bioconductor package `curatedTCGAData`
2.1.1) but are **not stored in this repository** because of their size
(`data/` and `*.rds` are in `.gitignore`). The app has its own reduced copy
in `shinyapp/data/` (step 5), which is in the repository.
Step 1 regenerates them.

## Software

- **Data preparation and QC:** R (4.3 or later) with the packages loaded at the
  top of each script (e.g. `curatedTCGAData`, `IlluminaHumanMethylation450kanno.ilmn12.hg19`,
  `wateRmelon`, `matrixStats`, `ggplot2`).
- **App and 06_precompute_results.R:** exact package versions recorded with
  `renv` in `shinyapp/renv.lock`; `renv::restore()` in the `shinyapp/` project
  installs them (see [`shinyapp/README.md`](shinyapp/README.md)).
