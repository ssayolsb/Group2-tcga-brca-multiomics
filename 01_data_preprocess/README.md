# 01–03 Data preparation and quality control

Scripts to download the TCGA-BRCA data, split them into training and test
sets, and run the quality control (QC). Every QC decision is learned on the
**training set** and saved, and the test set is processed with the saved
decisions only.

| Script | Input | Output | What it does |
|---|---|---|---|
| `01_data_preparation.R` / `.Rmd` | `curatedTCGAData` (2.1.1) | `data/data_split/*.train.rds`, `*.test.rds` | Download, inclusion criteria, outcome, RNA-seq counts to log2(CPM + 1), 70/30 random split |
| `02_qc_rnaseq.R` / `.Rmd` | `rna.train.rds` | `processed/rna.train.filtered.rds`, `rna_qc_params.rds` | Gene filters and sample checks (incl. PCA) |
| `02_qc_prot.R` / `.Rmd` | `prot.train.rds` | `processed/prot.train.filtered.rds`, `prot_qc_params.rds` | Proteins and patients with too many missing values (≤ 5% per patient) |
| `02_qc_meth.R` / `.Rmd` | `meth.train.rds` | `processed/meth.train.region_M.rds`, `meth_qc_params.rds` | Probe filters, patients with > 5% missing CpGs, BMIQ normalisation, mean beta per region (promoter, gene body, CpG island), M-values |
| `03_apply_preprocessing.R` | `*.test.rds` + `*_qc_params.rds` | `processed/*.test.*` and `test_qc_log.csv` | The training QC applied to the test set (nothing is learned from the test data) |

The `.Rmd` versions are the QC reports (HTML in `reports/`, with the same
number as the code); they read the saved results and do not rerun the slow
steps.

## Inclusion criteria and outcome

- Primary solid tumours, female patients, data in all three omics, known
  pathological stage (stage X excluded).
- Outcome: AJCC pathological stage in three classes: Stage I (I, IA, IB),
  Stage II (II, IIA, IIB), Stage III/IV (IIIA, IIIB, IIIC, IV).
- Split: 70% training, 30% test, at the patient level, before any
  data-driven preprocessing.

## Data format

All files: patients in rows, features in columns, plus a `patientID` column
(12-character TCGA barcode, e.g. `TCGA-A1-A0SB`).

## Runtime

BMIQ normalisation of the methylation data is done sample by sample and takes
several hours for the training set (about 40% of that for the test set).
