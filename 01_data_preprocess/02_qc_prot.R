# =============================================================================
#  PROTEOMICS (RPPA) - QC AND PREPROCESSING
#  TCGA-BRCA, training set
#
#  Input : prot.train.rds  (samples in rows, proteins in columns, plus patientID)
#          clin.train.rds  (clinical data, one row per patient)
#          rna.train.rds   (optional, for the RNA-protein concordance check)
#  Output: prot.train.filtered.rds  (filtered RPPA data)
#          prot_qc_params.rds       (parameters learned on train, to apply to test)
#
#  RPPA values are relative protein levels: log2 scale, normalised and
#  median-centred per antibody by TCGA. Negative values are expected
#  (below the cohort median). Values are comparable across patients for the
#  same protein, but not between different proteins.
#
#  Convention: numeric MATRIX with SAMPLES IN ROWS and PROTEINS IN COLUMNS.
#
#  Sections
#    0. Setup: packages, plot theme and helper functions
#    1. Read data
#    2. Sanity checks
#    3. Missing values (per protein and per sample)
#    4. Data-driven filtering (learned on train)
#    5. Distributions per sample and outlier samples
#    6. Protein variability
#    7. Duplicated or mixed-up samples
#    8. PCA (structure, outcome and batch)
#    9. RNA-protein concordance (sample identity check)
#   10. Save
# =============================================================================

rm(list = ls())


# =============================================================================
# 0. SETUP
# =============================================================================
library(matrixStats)   # fast row/column statistics on matrices
library(ggplot2)

# ---- Common theme for all plots (same as the other QC reports) ----------------
theme_qc <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40", size = 9),
          legend.position = "top")
}

# ---- Helper: barplot of missing-value proportions grouped in bins -------------
#   step     : bin width (0.05 = 5% bins)
#   max_frac : last bin edge; values above it are grouped in a "> x%" bin
# Colours: green = no missing, blue = some missing, red = high missingness
plot_na_bins <- function(na_frac, title, x_lab, step = 0.05, max_frac = 1) {
  upper  <- seq(step, max_frac, by = step)
  breaks <- c(-Inf, 0, upper, if (max_frac < 1) Inf)
  labels <- c("0%", paste0(round(100 * (upper - step)), "-", round(100 * upper), "%"),
              if (max_frac < 1) paste0("> ", round(100 * max_frac), "%"))
  df <- as.data.frame(table(bin = cut(na_frac, breaks = breaks, labels = labels)))
  df$pct <- df$Freq / sum(df$Freq)
  
  n_bins   <- nrow(df)
  high_lab <- if (max_frac < 1) paste0("> ", round(100 * max_frac), "% missing") else "> 90% missing"
  is_high  <- if (max_frac < 1) seq_len(n_bins) == n_bins else c(0, upper) > 0.9
  df$group <- factor(ifelse(seq_len(n_bins) == 1, "No missing",
                            ifelse(is_high, high_lab, "Some missing")),
                     levels = c("No missing", "Some missing", high_lab))
  
  ggplot(df, aes(bin, Freq, fill = group)) +
    geom_col(width = 0.8) +
    geom_text(data = subset(df, Freq > 0),
              aes(label = sprintf("%d\n(%.0f%%)", Freq, 100 * pct)),
              vjust = -0.2, size = 2.6, lineheight = 0.85) +
    scale_fill_manual(values = c("#1A9850", "#2C7FB8", "#D73027"), name = NULL, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = title, subtitle = sprintf("n = %d", length(na_frac)),
         x = x_lab, y = "Count") +
    theme_qc() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# =============================================================================
# 1. READ DATA
# =============================================================================
path.x <- "./data/data_split"

df <- readRDS(file.path(path.x, "prot.train.rds"))
prot <- as.matrix(df[, colnames(df) != "patientID"])   # drop the ID column (not numeric)
rownames(prot) <- df$patientID                         # patient IDs as row names
rm(df)

clin <- readRDS(file.path(path.x, "clin.train.rds"))
clin <- clin[rownames(prot), ]                         # same patients, same order

class(prot); dim(prot)
head(colnames(prot), 20)                               # antibody names (total and phospho forms)


# =============================================================================
# 2. SANITY CHECKS
# Expected: log2 relative values centred around 0 (negatives are normal),
# no duplicated patients or antibodies.
# =============================================================================
data.frame(
  check = c("Patients", "Proteins (antibodies)", "Phospho-specific antibodies",
            "Missing values", "Duplicated patients", "Duplicated antibodies",
            "Minimum value", "Maximum value", "Median of protein medians"),
  value = c(nrow(prot), ncol(prot),
            sum(grepl("_p[STY]", colnames(prot))),          # e.g. AKT_pS473
            sum(is.na(prot)),
            sum(duplicated(rownames(prot))), sum(duplicated(colnames(prot))),
            round(min(prot, na.rm = TRUE), 2), round(max(prot, na.rm = TRUE), 2),
            round(median(colMedians(prot, na.rm = TRUE)), 3))
)


# =============================================================================
# 3. MISSING VALUES
# =============================================================================

# ---- 3a. Missing values per protein --------------------------------------------
# Some antibodies were added to or removed from the TCGA RPPA panel between
# batches, so they can be missing in a subset of samples.
na_prot <- colMeans(is.na(prot))
plot_na_bins(na_prot, "Missing values per protein",
             "Percentage of samples with missing value")

# ---- 3b. Missing values per sample (proteins measured in most samples) -------
qc_prot   <- names(na_prot)[na_prot <= 0.2]
na_sample <- rowMeans(is.na(prot[, qc_prot]))
plot_na_bins(na_sample, "Missing values per sample",
             "Percentage of proteins with missing value", step = 0.01, max_frac = 0.10)

# Samples with more than 5% missing proteins are flagged as low quality
bad_samples <- names(na_sample)[na_sample > 0.05]
bad_samples


# =============================================================================
# 4. DATA-DRIVEN FILTERING (learned on the training set)
# =============================================================================

# ---- 4a. Remove low-quality samples -------------------------------------------
prot <- prot[!rownames(prot) %in% bad_samples, , drop = FALSE]
clin <- clin[rownames(prot), ]
# If samples are removed, also remove them from rna and meth to keep all datasets aligned

# ---- 4b. Remove proteins missing in more than 20% of samples ------------------
na_prot <- colMeans(is.na(prot))
sum(na_prot > 0.2)                       # number of proteins removed
prot    <- prot[, na_prot <= 0.2, drop = FALSE]

# =============================================================================
# 5. DISTRIBUTIONS PER SAMPLE AND OUTLIER SAMPLES
# After TCGA normalisation, all samples should have similar distributions
# centred around 0. A shifted or much wider distribution may indicate a
# technical problem (e.g. poor protein loading).
# =============================================================================

# ---- 5a. Boxplot of protein values per sample (ordered by median) -------------
box_df <- data.frame(patientID = rep(rownames(prot), ncol(prot)),
                     value     = as.vector(prot))
box_df$patientID <- factor(box_df$patientID,
                           levels = rownames(prot)[order(rowMedians(prot))])

ggplot(box_df, aes(x = patientID, y = value)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.2, fill = "grey85") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#D73027") +
  labs(title = "Distribution of protein levels per sample",
       subtitle = "Samples ordered by median; dashed line = cohort median (0)",
       x = "Samples", y = "log2 relative protein level") +
  theme_qc() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ---- 5b. Sample-level summary statistics and Tukey outliers -------------------
sample_stats <- data.frame(patientID = rownames(prot),
                           median    = rowMedians(prot),
                           iqr       = rowIQRs(prot))

# Tukey fences (Q1 - 3 IQR, Q3 + 3 IQR) on the sample medians and IQRs
tukey_out <- function(x, k = 3) {
  q <- quantile(x, c(0.25, 0.75)); r <- diff(q)
  x < q[1] - k * r | x > q[2] + k * r
}
sample_stats$outlier <- tukey_out(sample_stats$median) | tukey_out(sample_stats$iqr)

ggplot(sample_stats, aes(x = median, y = iqr, colour = outlier)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "#D73027"),
                      labels = c("Other samples", "Tukey outlier (3 x IQR)"), name = NULL) +
  labs(title = "Median vs. spread of protein levels per sample",
       subtitle = "Samples far from the cloud may have technical problems",
       x = "Median log2 protein level", y = "Interquartile range") +
  theme_qc()

outlier_samples <- sample_stats$patientID[sample_stats$outlier]
outlier_samples   # flagged for inspection, not removed automatically


# =============================================================================
# 6. PROTEIN VARIABILITY
# Proteins with very low variability carry little information to distinguish
# patients. With ~200 proteins, all are usually kept; this plot is for QC.
# =============================================================================
# ---- 6a. Variability of RPPA proteins (SD across samples) ---------------------
# Standard deviation of each protein across samples, for total and phospho-specific
# antibodies. Proteins with near-zero SD would carry no information to distinguish
# patients; with ~200 proteins, all are usually kept.
prot_stats <- data.frame(protein = colnames(prot),
                         sd      = colSds(prot),
                         phospho = ifelse(grepl("_p[STY]", colnames(prot)),
                                          "Phospho-specific", "Total protein"))

ggplot(prot_stats, aes(x = sd, fill = phospho)) +
  geom_histogram(bins = 40, colour = "white") +
  scale_fill_manual(values = c(`Total protein` = "#2C7FB8", `Phospho-specific` = "#D73027"),
                    name = NULL) +
  labs(title = "Variability of RPPA proteins",
       subtitle = "Standard deviation of log2 relative levels across samples",
       x = "Standard deviation across samples", y = "Proteins") +
  theme_qc()

# Most and least variable proteins
head(prot_stats[order(prot_stats$sd, decreasing = TRUE), ], 10)
head(prot_stats[order(prot_stats$sd), ], 10)

# ---- 6b. Distribution across patients of random proteins ----------------------
# RPPA values are median-centred (all medians ~0), so proteins are stratified
# by variability (tertiles of SD) instead of by mean level
set.seed(1)
strata <- cut(prot_stats$sd, quantile(prot_stats$sd, c(0, 1/3, 2/3, 1)),
              labels = c("Low", "Medium", "High"), include.lowest = TRUE)
names(strata) <- prot_stats$protein
pick <- unlist(lapply(split(prot_stats$protein, strata), sample, size = 4))

prot_df_long <- data.frame(protein = rep(pick, each = nrow(prot)),
                           stratum = rep(strata[pick], each = nrow(prot)),
                           value   = as.vector(prot[, pick]))
prot_df_long$protein <- factor(prot_df_long$protein, levels = pick)

ggplot(prot_df_long, aes(value, fill = stratum)) +
  geom_histogram(bins = 25, colour = "white") +
  facet_wrap(~ protein, scales = "free", ncol = 4) +
  scale_fill_manual(values = c(Low = "#7FCDBB", Medium = "#2C7FB8", High = "#253494"),
                    name = "Variability (SD)") +
  labs(title = "Distribution across patients of randomly selected proteins",
       subtitle = "4 random proteins per tertile of standard deviation",
       x = "log2 relative protein level", y = "Patients") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40"),
        strip.text    = element_text(face = "bold"),
        legend.position = "top")

# =============================================================================
# 7. DUPLICATED OR MIXED-UP SAMPLES
# Two different patients should never have near-identical profiles (r > 0.99).
# =============================================================================
cor_mat <- cor(t(prot))                 # sample x sample correlation
diag(cor_mat) <- NA                     # ignore self-correlation

dup_pairs <- which(cor_mat > 0.99, arr.ind = TRUE)
dup_pairs <- dup_pairs[dup_pairs[, 1] < dup_pairs[, 2], , drop = FALSE]   # each pair once
data.frame(sample1 = rownames(cor_mat)[dup_pairs[, 1]],
           sample2 = colnames(cor_mat)[dup_pairs[, 2]],
           r       = cor_mat[dup_pairs])

max_cor <- data.frame(patientID = rownames(cor_mat),
                      max_r     = apply(cor_mat, 1, max, na.rm = TRUE))

ggplot(max_cor, aes(x = max_r)) +
  geom_histogram(bins = 50, fill = "grey60", colour = "white") +
  labs(title = "Highest correlation with another sample",
       subtitle = "All RPPA proteins; values near 1 suggest duplicated samples",
       x = "Pearson r", y = "Samples") +
  theme_qc()


# =============================================================================
# 8. PCA
# Centred, not scaled: all proteins are on the same (log2 relative) scale.
# =============================================================================
pc      <- prcomp(prot, center = TRUE, scale. = FALSE)
var_exp <- round(100 * summary(pc)$importance[2, 1:2], 1)

pca_df <- data.frame(PC1     = pc$x[, 1],
                     PC2     = pc$x[, 2],
                     stage   = clin[rownames(prot), "outcome_p.stage"],
                     outlier = rownames(prot) %in% outlier_samples)

# ---- 8a. PCA coloured by pathological stage -----------------------------------
ggplot(pca_df, aes(x = PC1, y = PC2, colour = stage)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA of RPPA protein levels",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "Stage") +
  theme_qc() +
  theme(legend.position = "right")

# ---- 8b. PCA highlighting the outlier samples from step 5 ----------------------
ggplot(pca_df, aes(x = PC1, y = PC2, colour = outlier)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "grey50", `TRUE` = "#D73027"),
                      labels = c("Other samples", "Tukey outlier"), name = NULL) +
  labs(title = "PCA of RPPA protein levels: outlier samples",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_qc() +
  theme(legend.position = "right")

# ---- 8c. Proteins that drive PC1 and PC2 (top absolute loadings) --------------
head(sort(abs(pc$rotation[, 1]), decreasing = TRUE), 10)
head(sort(abs(pc$rotation[, 2]), decreasing = TRUE), 10)

# ---- 8d. Optional: PCA coloured by RPPA batch (plate) --------------------------
# The RPPA plate is part of the full TCGA barcode (characters 22-25), which was
# trimmed when the data were prepared. It is recovered from curatedTCGAData.
do_batch <- FALSE
if (do_batch) {
  library(curatedTCGAData)
  library(MultiAssayExperiment)
  
  mae <- curatedTCGAData(diseaseCode = "BRCA", assays = "RPPAArray",
                         version = "2.1.1", dry.run = FALSE)
  barcodes <- colnames(mae[["BRCA_RPPAArray-20160128"]])
  
  batch_info <- data.frame(patientID   = substr(barcodes, 1, 12),
                           sample_type = substr(barcodes, 14, 15),   # 01 = primary tumour
                           plate       = substr(barcodes, 22, 25))
  batch_info <- batch_info[batch_info$sample_type == "01", ]
  batch_info <- batch_info[!duplicated(batch_info$patientID), ]
  rownames(batch_info) <- batch_info$patientID
  rm(mae)
  
  pca_df$plate <- batch_info[rownames(pca_df), "plate"]
  table(pca_df$plate, useNA = "ifany")
  
  ggplot(pca_df, aes(x = PC1, y = PC2, colour = plate)) +
    geom_point(size = 2, alpha = 0.7) +
    labs(title = "PCA of RPPA protein levels coloured by plate",
         x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
         colour = "Plate") +
    theme_qc() +
    theme(legend.position = "right")
}


# =============================================================================
# 9. RNA-PROTEIN CONCORDANCE (sample identity check)
# For well-known breast cancer markers, gene expression and protein level should
# be positively correlated. A sample with discordant values in several markers
# may be a sample swap. Antibody names differ between panel versions: check the
# exact names with grep() and adjust the mapping below if needed.
# =============================================================================
do_rna_check <- TRUE
if (do_rna_check) {
  df  <- readRDS(file.path(path.x, "rna.train.rds"))
  rna <- as.matrix(df[, colnames(df) != "patientID"])
  rownames(rna) <- df$patientID
  rm(df)
  
  # Find the antibody names of the markers in this dataset
  grep("ER|PR|HER2|AR|GATA3|CDH1|E-Cad", colnames(prot), value = TRUE, ignore.case = TRUE)
  
  # Antibody -> gene mapping (adjust the antibody names to those found above)
  markers <- c(`ER-alpha`   = "ESR1",
               `PR`         = "PGR",
               `HER2`       = "ERBB2",
               `AR`         = "AR",
               `GATA3`      = "GATA3",
               `E-Cadherin` = "CDH1")
  markers <- markers[names(markers) %in% colnames(prot) & markers %in% colnames(rna)]
  
  common <- intersect(rownames(prot), rownames(rna))
  concord <- data.frame(antibody = names(markers),
                        gene     = markers,
                        r        = sapply(names(markers), function(a) {
                          cor(prot[common, a], rna[common, markers[a]], method = "spearman")
                        }))
  concord   # expected: clearly positive correlations (e.g. ER, PR, HER2 > 0.6)
  
  # Scatter plots of protein vs. RNA for each marker
  conc_df <- do.call(rbind, lapply(names(markers), function(a) {
    data.frame(marker  = sprintf("%s (protein) vs %s (RNA)", a, markers[a]),
               protein = prot[common, a],
               rna     = rna[common, markers[a]])
  }))
  
  ggplot(conc_df, aes(x = rna, y = protein)) +
    geom_point(size = 1, alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, colour = "#D73027", linewidth = 0.6) +
    facet_wrap(~ marker, scales = "free") +
    labs(title = "RNA-protein concordance of breast cancer markers",
         x = "log2(RSEM + 1)", y = "log2 relative protein level") +
    theme_qc()
}


# =============================================================================
# 10. SAVE
# Samples in rows, with a patientID column like the other datasets
# =============================================================================
path.x <- "./data/processed"
dir.create(path.x, showWarnings = FALSE, recursive = TRUE)

prot_df <- data.frame(prot, check.names = FALSE)
prot_df$patientID <- rownames(prot_df)
prot_df <- prot_df[, c(ncol(prot_df), 1:(ncol(prot_df) - 1))]   # patientID as first column

# QC parameters learned on train, to apply the same QC to test / Shiny data
prot_qc_params <- list(
  kept_proteins = colnames(prot),   # Proteins kept after filtering (step 4b)
  max_na_sample = 0.05              # Max fraction of missing proteins per sample (step 3b)
)

saveRDS(prot_df, file.path(path.x, "prot.train.filtered.rds"))
saveRDS(prot_qc_params, file.path(path.x, "prot_qc_params.rds"))   # QC parameters for test / Shiny

