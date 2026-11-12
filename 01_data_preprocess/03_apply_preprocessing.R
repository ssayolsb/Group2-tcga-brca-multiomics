# =============================================================================
#  APPLY THE TRAINING QC TO THE TEST SET
#
#  The QC scripts (02_qc_*.R) learned their decisions on the training set and
#  saved them in data/processed/*_qc_params.rds. This script applies exactly
#  the same decisions to the test set. Nothing is learned from the test data.
#
#  Input : data/data_split/rna.test.rds, meth.test.rds, prot.test.rds
#          data/processed/rna_qc_params.rds   (kept_genes)
#          data/processed/meth_qc_params.rds  (kept_probes, region_lists, min_obs)
#          data/processed/prot_qc_params.rds  (kept_proteins, max_na_sample)
#  Output: data/processed/rna.test.filtered.rds
#          data/processed/meth.test.region_beta.rds   region-level beta values
#          data/processed/meth.test.region_M.rds      region-level M-values
#          data/processed/prot.test.filtered.rds
#          data/processed/test_qc_log.csv             patients removed and why
#  Same format as the training files: patients in rows, patientID column first,
#  features in the same order as in training.
#
#  Steps per omics (same as in the training QC)
#    RNA (log2(CPM + 1), computed per sample in 01_data_preparation.R)
#        keep the training genes
#    Protein (RPPA)
#        keep the training proteins; remove patients with > max_na_sample
#        missing proteins (same rule as in training)
#    Methylation (450k beta values)
#        keep the training CpGs (probe filters) -> remove patients with > 5%
#        missing CpGs (same rule as in training) -> BMIQ per sample ->
#        mean beta per region (same region definitions) -> M-values
#
#  Runtime: BMIQ is applied sample by sample, as in training. With ~190 test
#  patients it takes roughly 40% of the training time (several hours on a
#  laptop). Set do_bmiq <- FALSE only to test the script.
#
#  Sections
#    0. Setup
#    1. RNA-seq
#    2. Proteomics (RPPA)
#    3. Methylation: CpG filtering and sample QC
#    4. Methylation: BMIQ normalisation
#    5. Methylation: region-level features and M-values
#    6. Save the log of removed patients
# =============================================================================

rm(list = ls())


# =============================================================================
# 0. SETUP
# =============================================================================
library(matrixStats)

split_dir <- "data/data_split"
proc_dir  <- "data/processed"

do_bmiq  <- TRUE                                   # FALSE = skip BMIQ (testing only)
n_cores  <- max(1, parallel::detectCores() - 1)    # BMIQ in parallel
max_na_meth_sample <- 0.05                         # same threshold as in 02_qc_meth.R

# Read a split file as a numeric matrix (patients x features)
read_split <- function(file) {
  df <- readRDS(file)
  x  <- as.matrix(df[, colnames(df) != "patientID", drop = FALSE])
  rownames(x) <- df$patientID
  x
}

# Save in the same format as the training files (patientID as first column)
save_processed <- function(x, file) {
  df <- data.frame(patientID = rownames(x), x, check.names = FALSE)
  rownames(df) <- df$patientID
  saveRDS(df, file)
}

removed_log <- list()                             # patients removed, per omics


# =============================================================================
# 1. RNA-SEQ: keep the training genes
# =============================================================================
rna_params <- readRDS(file.path(proc_dir, "rna_qc_params.rds"))
stopifnot("rna_qc_params.rds must contain kept_genes" = !is.null(rna_params$kept_genes))
rna_test   <- read_split(file.path(split_dir, "rna.test.rds"))

missing_genes <- setdiff(rna_params$kept_genes, colnames(rna_test))
length(missing_genes)                              # should be 0 (same platform)
rna_test <- rna_test[, intersect(rna_params$kept_genes, colnames(rna_test)), drop = FALSE]

dim(rna_test)
save_processed(rna_test, file.path(proc_dir, "rna.test.filtered.rds"))
rm(rna_test); gc()


# =============================================================================
# 2. PROTEOMICS (RPPA): keep the training proteins; same sample rule
# =============================================================================
prot_params <- readRDS(file.path(proc_dir, "prot_qc_params.rds"))
stopifnot("prot_qc_params.rds must contain kept_proteins" = !is.null(prot_params$kept_proteins))
if (is.null(prot_params$max_na_sample)) prot_params$max_na_sample <- 0.05   # default of 02_qc_prot.R
prot_test   <- read_split(file.path(split_dir, "prot.test.rds"))

prot_test <- prot_test[, intersect(prot_params$kept_proteins, colnames(prot_test)), drop = FALSE]

na_prot_sample <- rowMeans(is.na(prot_test))
bad_prot <- names(na_prot_sample)[na_prot_sample > prot_params$max_na_sample]
bad_prot
if (length(bad_prot)) {
  removed_log$prot <- data.frame(patientID = bad_prot, omics = "prot",
                                 reason = sprintf("%.1f%% missing proteins",
                                                  100 * na_prot_sample[bad_prot]))
  prot_test <- prot_test[!rownames(prot_test) %in% bad_prot, , drop = FALSE]
}

dim(prot_test)
save_processed(prot_test, file.path(proc_dir, "prot.test.filtered.rds"))
rm(prot_test); gc()


# =============================================================================
# 3. METHYLATION: CpG FILTERING AND SAMPLE QC
# The probe filters (non-CpG, sex chromosomes, SNPs, cross-reactive, > 20%
# missing in training) are summarised in kept_probes: the test set keeps
# exactly the same CpGs.
# =============================================================================
meth_params <- readRDS(file.path(proc_dir, "meth_qc_params.rds"))
stopifnot("meth_qc_params.rds must contain kept_probes, region_lists and min_obs" =
            all(c("kept_probes", "region_lists", "min_obs") %in% names(meth_params)))

# Read and keep only the training CpGs before converting to a matrix (memory)
df <- readRDS(file.path(split_dir, "meth.test.rds"))
keep_cols <- intersect(meth_params$kept_probes, colnames(df))
length(setdiff(meth_params$kept_probes, keep_cols))   # training CpGs absent in test (should be 0)
meth <- as.matrix(df[, keep_cols, drop = FALSE])
rownames(meth) <- df$patientID
rm(df); gc()
dim(meth)

# Same sample rule as in training: > 5% missing CpGs = low quality
na_meth_sample <- rowMeans(is.na(meth))
bad_meth <- names(na_meth_sample)[na_meth_sample > max_na_meth_sample]
bad_meth
if (length(bad_meth)) {
  removed_log$meth <- data.frame(patientID = bad_meth, omics = "meth",
                                 reason = sprintf("%.1f%% missing CpGs",
                                                  100 * na_meth_sample[bad_meth]))
  meth <- meth[!rownames(meth) %in% bad_meth, , drop = FALSE]
}


# =============================================================================
# 4. METHYLATION: BMIQ NORMALISATION (per sample, as in training)
# =============================================================================
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(RPMM)
library(parallel)

anno     <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
type     <- anno[colnames(meth), "Type"]
design_v <- ifelse(type == "I", 1, 2)

if (do_bmiq) {
  # Normalise one block of samples (rows); failed samples are left unchanged
  bmiq_block <- function(m, design_v) {
    for (i in seq_len(nrow(m))) {
      ok  <- !is.na(m[i, ])
      res <- tryCatch(
        wateRmelon::BMIQ(m[i, ok], design_v[ok], plots = FALSE, pri = FALSE)$nbeta,
        error = function(e) NULL
      )
      if (!is.null(res)) m[i, ok] <- res
    }
    m
  }

  blocks <- split(seq_len(nrow(meth)), cut(seq_len(nrow(meth)), n_cores, labels = FALSE))
  chunks <- lapply(blocks, function(i) meth[i, , drop = FALSE])

  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, { library(wateRmelon); library(RPMM) })
  t0  <- Sys.time()
  res <- parLapply(cl, chunks, bmiq_block, design_v = design_v)
  stopCluster(cl)
  print(Sys.time() - t0)

  meth <- do.call(rbind, res)[rownames(meth), ]
  rm(chunks, res); gc()
}


# =============================================================================
# 5. METHYLATION: REGION-LEVEL FEATURES AND M-VALUES
# Same region definitions (region_lists) and the same rule (min_obs) as in
# training; the columns are built in the same order as the training files:
# promoters, gene bodies, CpG islands.
# =============================================================================
library(Matrix)

# Mean beta of the observed CpGs of each region (same function as 02_qc_meth.R)
region_mean <- function(m, rl, min_obs) {
  rl <- lapply(rl, function(p) intersect(p, colnames(m)))  # CpGs present in the test data
  A <- sparseMatrix(i    = match(unlist(rl), colnames(m)),
                    j    = rep(seq_along(rl), lengths(rl)),
                    x    = 1,
                    dims = c(ncol(m), length(rl)),
                    dimnames = list(colnames(m), names(rl)))
  obs <- !is.na(m)
  m0  <- m
  m0[!obs] <- 0
  sums   <- as.matrix(m0  %*% A)
  counts <- as.matrix(obs %*% A)
  out <- sums / counts
  out[counts < min_obs] <- NA
  out
}

rl <- meth_params$region_lists
min_obs <- meth_params$min_obs

meth_promoter <- region_mean(meth, rl$promoter, min_obs)
meth_body     <- region_mean(meth, rl$body,     min_obs)
meth_island   <- region_mean(meth, rl$island,   min_obs)
colnames(meth_promoter) <- paste0(colnames(meth_promoter), ".promoter")
colnames(meth_body)     <- paste0(colnames(meth_body),     ".genebody")
colnames(meth_island)   <- paste0(colnames(meth_island),   ".island")

methB_test <- cbind(meth_promoter, meth_body, meth_island)
rm(meth_promoter, meth_body, meth_island, meth); gc()

# Same column order as the training file
train_cols <- setdiff(colnames(readRDS(file.path(proc_dir, "meth.train.region_M.rds"))), "patientID")
c(train = length(train_cols), test = ncol(methB_test),
  missing_in_test = length(setdiff(train_cols, colnames(methB_test))))
methB_test <- methB_test[, intersect(train_cols, colnames(methB_test)), drop = FALSE]

# M-values with the same offset as in training
to_M <- function(b, offset = 1e-6) log2((b + offset) / (1 - b + offset))
methM_test <- to_M(methB_test)

dim(methM_test)
mean(is.na(methM_test))                             # proportion of missing region values

save_processed(methB_test, file.path(proc_dir, "meth.test.region_beta.rds"))
save_processed(methM_test, file.path(proc_dir, "meth.test.region_M.rds"))


# =============================================================================
# 6. SAVE THE LOG OF REMOVED PATIENTS
# =============================================================================
removed <- do.call(rbind, removed_log)
if (is.null(removed)) removed <- data.frame(patientID = character(0), omics = character(0),
                                            reason = character(0))
removed
write.csv(removed, file.path(proc_dir, "test_qc_log.csv"), row.names = FALSE)
