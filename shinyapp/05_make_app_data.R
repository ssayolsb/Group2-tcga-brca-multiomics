# =============================================================================
#  05_make_app_data.R - reduced copy of the project data for the app
#
#  When the app is published only the shinyapp/ folder is uploaded, so the app
#  needs its own data. This script copies to shinyapp/data/ the files that
#  load_stage_data() (modeling.R) reads, with the same names and folders:
#    - RNA: the 5000 most variable genes of the training set (train and test)
#    - Protein: all QC proteins (train and test), the original training file
#      (sample screen) and the QC parameters
#    - DNA methylation: the 10000 most variable regions without missing
#      training values (train and test)
#    - Clinical: only the patient ID and the stage
#  The rankings are made on the training patients only, as in the app, so
#  the app gives the same results with these files for up to 5000 genes and
#  10000 regions.
#
#  Run once, after 01_data_preprocess/03_apply_preprocessing.R:
#    source("shinyapp/05_make_app_data.R")      # from the repo root
# =============================================================================

if (file.exists("shinyapp/app.R")) setwd("shinyapp")
src <- "../data"       # project data (not in the repository)
out <- "data"          # app data (in the repository)
n_genes   <- 5000
n_regions <- 10000

dir.create(file.path(out, "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "data_split"), recursive = TRUE, showWarnings = FALSE)

read_layer <- function(path) readRDS(file.path(src, path))
save_layer <- function(x, path) saveRDS(x, file.path(out, path), compress = "xz")

# Keep the patientID column and the given features, in this order
keep_columns <- function(df, features) df[, c("patientID", features), drop = FALSE]

# Most variable features of the training data (features with missing values
# in training are excluded, as in the app)
top_features <- function(train, n) {
  x <- as.matrix(train[, names(train) != "patientID", drop = FALSE])
  x <- x[, colSums(is.na(x)) == 0, drop = FALSE]
  v <- matrixStats::colVars(x)
  names(v) <- colnames(x)
  names(sort(v[v > 0], decreasing = TRUE))[seq_len(min(n, sum(v > 0)))]
}

# ---- RNA ----------------------------------------------------------------------
rna_train <- read_layer("processed/rna.train.filtered.rds")
rna_test  <- read_layer("processed/rna.test.filtered.rds")
genes <- top_features(rna_train, n_genes)
save_layer(keep_columns(rna_train, genes), "processed/rna.train.filtered.rds")
save_layer(keep_columns(rna_test, genes), "processed/rna.test.filtered.rds")
rm(rna_train, rna_test); invisible(gc())

# ---- Protein (small: copied unchanged) ----------------------------------------
for (path in c("processed/prot.train.filtered.rds", "processed/prot.test.filtered.rds",
               "processed/prot_qc_params.rds", "data_split/prot.train.rds")) {
  save_layer(read_layer(path), path)
}

# ---- DNA methylation ------------------------------------------------------------
meth_train <- read_layer("processed/meth.train.region_M.rds")
regions <- top_features(meth_train, n_regions)
save_layer(keep_columns(meth_train, regions), "processed/meth.train.region_M.rds")
rm(meth_train); invisible(gc())
meth_test <- read_layer("processed/meth.test.region_M.rds")
save_layer(keep_columns(meth_test, regions), "processed/meth.test.region_M.rds")
rm(meth_test); invisible(gc())

# ---- Clinical: patient ID and stage only -----------------------------------------
for (path in c("data_split/clin.train.rds", "data_split/clin.test.rds")) {
  clin <- read_layer(path)
  save_layer(clin[, c("patientID", "pathologic_stage", "outcome_p.stage"), drop = FALSE], path)
}

# ---- Summary ----------------------------------------------------------------------
files <- list.files(out, recursive = TRUE, full.names = TRUE)
print(data.frame(file = sub(paste0("^", out, "/"), "", files),
                 MB = round(file.size(files) / 1e6, 2)))
message("App data written to shinyapp/", out, " (", round(sum(file.size(files)) / 1e6, 1), " MB)")
