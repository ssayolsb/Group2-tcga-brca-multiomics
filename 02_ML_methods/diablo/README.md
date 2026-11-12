# 04 Models: DIABLO

**Report:** [`reports/04_diablo.html`](../../reports/04_diablo.html)
**Code:** `04_diablo.Rmd` (R Markdown) – *to be added*

## What it does

Supervised multi-omics integration with DIABLO (mixOmics `block.splsda`) of
RNA, protein and DNA methylation. The number of components and the number of
features per layer and component (keepX) are tuned by 5-fold cross-validation
(5 repeats) on the training set; the final sparse model predicts the stage of
the test patients (weighted vote, centroids distance).

## Input

- `data/processed/rna.train.filtered.rds`, `prot.train.filtered.rds`, `meth.train.region_M.rds`
- `data/processed/*.test.*` (test set processed with the training QC)
- `data/data_split/clin.train.rds`, `clin.test.rds`

## Output

The report: tuned number of components and keepX, cross-validated and test
balanced error rate, confusion matrix, sample, loading and circos plots.

## How to run

From the project root: `rmarkdown::render("02_ML_methods/diablo/04_diablo.Rmd", output_dir = "reports")`
(the parameters of the report set the data folder and the number of features).

## Settings selected (used again in `shinyapp/06_precompute_results.R`)

- Input: the 2000 most variable genes, 200 proteins and 5000 methylation regions.
- Design 0.1 between layers, 2 components.
- keepX (component 1 / 2): RNA 50 / 100, protein 5 / 40, methylation 25 / 100.
