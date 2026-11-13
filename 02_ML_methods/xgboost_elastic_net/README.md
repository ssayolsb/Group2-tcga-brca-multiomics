# 04 Models: XGBoost and elastic net

**Report:** [`reports/04_xgboost_elastic_net.html`](../../reports/04_xgboost_elastic_net.html)
**Code:** `04_xgboost_elastic_net.qmd` (Quarto)

## What it does

XGBoost and multinomial elastic-net logistic regression on RNA, protein and
RNA + protein. Settings and feature counts are chosen by 5-fold training
cross-validation (folds stratified by stage), and the six final models are
evaluated on the same test patients, with bootstrap intervals, ROC curves and
split conformal prediction sets (LAC, 90% target).

**Note:** this report uses two data layers (RNA and protein). In the app, these
models were extended to DNA methylation as a third layer, alone or together
with RNA and protein.

## Input

- `data/processed/rna.train.filtered.rds`, `prot.train.filtered.rds`, `prot_qc_params.rds`
- `data/data_split/rna.test.rds`, `prot.test.rds`, `prot.train.rds`, `clin.train.rds`, `clin.test.rds`

## Output

The report: selected settings, test metrics (accuracy, balanced accuracy,
macro F1, AUC), confusion matrices and conformal coverage.

## How to run

From the project root:

```
quarto render 02_ML_methods/xgboost_elastic_net/04_xgboost_elastic_net.qmd --execute-dir . --output-dir ../../reports
```

## Settings selected (used again in `shinyapp/06_precompute_results.R`)

| Model | Data | Features | Settings |
|---|---|---|---|
| XGBoost | RNA | 2000 | depth 1, 200 rounds |
| XGBoost | Protein | 216 | depth 1, 200 rounds |
| XGBoost | RNA + protein | 432 | depth 1, 200 rounds |
| Elastic net | RNA | 216 | alpha 0.25, lambda 0.111 |
| Elastic net | Protein | 216 | alpha 0.50, lambda 0.039 |
| Elastic net | RNA + protein | 432 | alpha 0.75, lambda 0.023 |
