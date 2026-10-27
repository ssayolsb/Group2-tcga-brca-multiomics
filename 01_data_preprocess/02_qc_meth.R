# =============================================================================
#  DNA METHYLATION (Illumina 450k) - QC, NORMALISATION AND REGION-LEVEL FEATURES
#  TCGA-BRCA, training set
#
#  Input : meth.train.rds  (samples in rows, probes in columns, plus patientID)
#          clin.train.rds  (clinical data, one row per patient)
#  Output: meth.train.region_beta.rds  region-level beta values
#          meth.train.region_M.rds     region-level M-values
#          meth.train.cpg_beta.rds     CpG-level beta values (filtered, BMIQ)
#          meth_qc_params.rds          QC parameters (kept CpGs, region definitions)
#                                      to apply the same QC to test / Shiny data
# 
#  Convention used throughout: numeric MATRIX with SAMPLES IN ROWS and
#  PROBES / FEATURES IN COLUMNS.
#
#  Sections
#    0. Setup: packages, plot theme and helper functions
#    1. Read data
#    2. Probe annotation
#    3. Sample-level QC (missing values, sex check)
#    4. Probe filtering based on annotation
#    5. Data-driven probe filtering (missing values)
#    6. QC plots on the filtered data (variability, duplicates, densities, PCA)
#    7. Batch information (array plate and centre)
#    8. Probe type bias and BMIQ normalisation
#    9. Region-level methylation (promoters, gene bodies, CpG islands)
#   10. QC of region-level features
#   11. M-values
#   12. Save
# =============================================================================

rm(list = ls())


# =============================================================================
# 0. SETUP
# =============================================================================

# ---- Packages (run once) -----------------------------------------------------
# BiocManager::install(c("IlluminaHumanMethylation450kanno.ilmn12.hg19",
#                        "matrixStats", "wateRmelon"))
# install.packages("remotes")
# remotes::install_version("RPMM", repos = "https://cran.r-project.org")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # probe annotation (chr, type, SNPs, genes)
library(matrixStats)                                   # fast row/column statistics on matrices
library(ggplot2)

# ---- Common theme for all plots (same as the transcriptomics QC report) -------
theme_qc <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40", size = 9),
          legend.position = "top")
}

# ---- Helper: barplot of missing-value proportions grouped in bins -------------
# Same style as the zero barplots in the transcriptomics QC report.
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
# matrixStats functions (colVars, rowMedians...) need a matrix, not a data
# frame, so the data are converted right after reading.
# =============================================================================
path.x <- "./data/data_split"

df <- readRDS(file.path(path.x, "meth.train.rds"))
meth <- as.matrix(df[, colnames(df) != "patientID"])   # drop the ID column (not numeric)
rownames(meth) <- df$patientID                         # patient IDs as row names
rm(df); gc()                                           # free memory (the data are large)

clin <- readRDS(file.path(path.x, "clin.train.rds"))
clin <- clin[rownames(meth), ]                         # same patients, same order

class(meth); dim(meth)


# =============================================================================
# 2. PROBE ANNOTATION
# Illumina annotation for each probe: chromosome, probe type (I/II), SNPs and
# genes. Here probes are in ROWS (one row per probe).
# =============================================================================
anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
anno <- anno[intersect(colnames(meth), rownames(anno)), ]   # only probes present in our data

x_probes <- rownames(anno)[which(anno$chr == "chrX")]
y_probes <- rownames(anno)[which(anno$chr == "chrY")]


# =============================================================================
# 3. SAMPLE-LEVEL QC
# Per-sample checks: each sample is assessed with its own data only.
# =============================================================================

# ---- 3a. Missing values per probe (all probes) -------------------------------
# Proportion of samples with NA for each probe. TCGA level 3 data contain a block
# of probes that are NA in (almost) all samples (masked during TCGA processing).
# These probes would inflate the missing rate of every sample by the same amount
# (~18%), so they are excluded before assessing samples (step 3b).
na_probe <- colMeans(is.na(meth))
plot_na_bins(na_probe, "Missing values per probe",
             "Percentage of samples with missing value")

# ---- 3b. Missing values per sample (probes measured in most samples) ---------
# Only probes measured in at least 80% of samples are used, so the per-sample
# missing rate reflects sample quality (e.g. poor DNA, failed hybridisation).
qc_probes <- names(na_probe)[na_probe <= 0.2]
na_sample <- rowMeans(is.na(meth[, qc_probes]))
plot_na_bins(na_sample, "Missing values per sample",
             "Percentage of probes with missing value", step = 0.01, max_frac = 0.10)

# Samples with more than 5% missing probes are flagged as low quality
bad_samples <- names(na_sample)[na_sample > 0.05]
bad_samples

# ---- 3c. Sex check (before removing X/Y probes) ------------------------------
# All patients are female. Females should form a single cloud: Y probes give only
# background signal and X probes show intermediate methylation (X inactivation).
# A sample far from the cloud could be a sample swap. Low chrX methylation can
# also be biological (loss of the inactive X, frequent in basal-like tumours).
sex_check <- data.frame(patientID = rownames(meth),
                        x_median  = rowMedians(meth[, x_probes], na.rm = TRUE),
                        y_mean    = rowMeans(meth[, y_probes], na.rm = TRUE))

ggplot(sex_check, aes(x = x_median, y = y_mean)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Sex check",
       subtitle = "Median chrX vs. mean chrY beta per sample",
       x = "Median chrX beta", y = "Mean chrY beta") +
  theme_qc()

# ---- 3d. Remove low-quality samples ------------------------------------------
meth <- meth[!rownames(meth) %in% bad_samples, , drop = FALSE]
clin <- clin[rownames(meth), ]
# If samples are removed, also remove them from rna and prot to keep all datasets aligned


# =============================================================================
# 4. PROBE FILTERING BASED ON ANNOTATION
# These filters only use the probe annotation (not the sample values), so they
# are identical for training and test.
# =============================================================================

# ---- 4a. Non-CpG probes ------------------------------------------------------
# "ch." probes measure CpH (non-CpG) methylation; "rs." probes are SNP controls
non_cpg <- colnames(meth)[!grepl("^cg", colnames(meth))]

# ---- 4b. Sex-chromosome probes -----------------------------------------------
# X/Y methylation differs strongly by sex and X inactivation; Y probes are only
# noise in females
sex_chr <- c(x_probes, y_probes)

# ---- 4c. SNP-affected probes (MAF >= 1%) -------------------------------------
# A SNP at the CpG site or at the single-base extension (SBE) position changes
# the measured signal, so the genotype would look like a methylation change
snp_probes <- rownames(anno)[which(anno$CpG_maf >= 0.01 | anno$SBE_maf >= 0.01)]

# ---- 4d. Cross-reactive probes (Chen et al. 2013) ----------------------------
# Probes that also hybridise to other genomic locations (non-specific probes)
url  <- "https://raw.githubusercontent.com/sirselim/illumina450k_filtering/master/48639-non-specific-probes-Illumina450k.csv"
chen <- read.csv(url)
xreactive <- if ("TargetID" %in% names(chen)) chen$TargetID else chen[, 1]
length(xreactive)

# ---- 4e. Remove annotation-based probes --------------------------------------
remove_probes <- unique(c(non_cpg, sex_chr, snp_probes, xreactive))

# Summary table: probes flagged by each reason (a probe can have several reasons)
data.frame(reason = c("Non-CpG (ch./rs.)", "Sex chromosomes", "SNP at CpG/SBE (MAF >= 1%)",
                      "Cross-reactive", "Total removed (unique)"),
           probes = c(length(non_cpg), length(sex_chr), length(snp_probes),
                      sum(xreactive %in% colnames(meth)),
                      sum(colnames(meth) %in% remove_probes)))

meth <- meth[, !colnames(meth) %in% remove_probes]          # probes are COLUMNS
anno <- anno[intersect(colnames(meth), rownames(anno)), ]   # keep annotation in sync
gc()
dim(meth)


# =============================================================================
# 5. DATA-DRIVEN PROBE FILTERING (learned on the training set)
# =============================================================================

# ---- 5a. Probes with many missing values -------------------------------------
# Remove probes with missing values in more than 20% of samples (includes the
# block of probes masked by TCGA seen in step 3a)
na_probe <- colMeans(is.na(meth))
sum(na_probe > 0.2)                      # number of probes removed
meth     <- meth[, na_probe <= 0.2]
gc()


# =============================================================================
# 6. QC PLOTS ON THE FILTERED DATA
# =============================================================================

# ---- 6a. Most variable probes ------------------------------------------------
# The 5,000 most variable probes capture the main structure of the data and make
# the correlation between samples much faster than using all probes
probe_var <- colVars(meth, na.rm = TRUE)
top_var   <- colnames(meth)[order(probe_var, decreasing = TRUE)[1:5000]]

# Distribution of probe variability (SD is on the same scale as beta values)
probe_sd <- sqrt(probe_var)
cut_5k   <- sort(probe_sd, decreasing = TRUE)[5000]    # SD of the 5,000th probe
cut_10k  <- sort(probe_sd, decreasing = TRUE)[10000]   # SD of the 10,000th probe

ggplot(data.frame(sd = probe_sd), aes(x = sd)) +
  geom_histogram(bins = 100, fill = "grey60", colour = "white") +
  geom_vline(xintercept = cut_5k,  linetype = "dashed", colour = "#D73027") +
  geom_vline(xintercept = cut_10k, linetype = "dotted", colour = "#2C7FB8") +
  labs(title = "Variability of methylation probes",
       subtitle = sprintf("Red dashed: top 5,000 (SD > %.3f, QC plots); blue dotted: top 10,000 (SD > %.3f)",
                          cut_5k, cut_10k),
       x = "Standard deviation of beta values across samples", y = "Probes") +
  theme_qc()

# ---- 6b. Duplicated or mixed-up samples --------------------------------------
# Correlation between SAMPLES (cor() correlates columns, hence the transpose).
# Two different patients should never have near-identical profiles (r > 0.99).
cor_mat <- cor(t(meth[, top_var]), use = "pairwise.complete.obs")   # sample x sample
diag(cor_mat) <- NA                                                 # ignore self-correlation

dup_pairs <- which(cor_mat > 0.99, arr.ind = TRUE)
dup_pairs <- dup_pairs[dup_pairs[, 1] < dup_pairs[, 2], , drop = FALSE]   # each pair once
data.frame(sample1 = rownames(cor_mat)[dup_pairs[, 1]],
           sample2 = colnames(cor_mat)[dup_pairs[, 2]],
           r       = cor_mat[dup_pairs])

# Highest correlation of each sample with any other sample:
# near 1 = possible duplicate; far left = sample unlike any other (check it)
max_cor <- data.frame(patientID = rownames(cor_mat),
                      max_r     = apply(cor_mat, 1, max, na.rm = TRUE))

ggplot(max_cor, aes(x = max_r)) +
  geom_histogram(bins = 50, fill = "grey60", colour = "white") +
  labs(title = "Highest correlation with another sample",
       subtitle = "5,000 most variable probes; values near 1 suggest duplicated samples",
       x = "Pearson r", y = "Samples") +
  theme_qc()

# ---- 6c. Beta value distributions (all probes) -------------------------------
# Beta values should be bimodal (peaks near 0 and 1). A sample with a flattened
# or shifted curve may have technical problems.
grid <- seq(0, 1, length.out = 256)   # common x positions for all densities

# Density of each sample on the common grid (one column per sample)
dens_mat <- apply(meth, 1, function(x) {
  density(x, na.rm = TRUE, from = 0, to = 1, n = 256)$y
})

# Long format for ggplot: one row per (sample, beta value)
dens_df <- data.frame(beta      = rep(grid, ncol(dens_mat)),
                      density   = as.vector(dens_mat),
                      patientID = rep(colnames(dens_mat), each = length(grid)))

# Height of the unmethylated peak per sample: maximum density for beta <= 0.1
low_zone  <- grid <= 0.1
peak_low  <- apply(dens_mat[low_zone, ], 2, max)

# Flag samples with a low (shifted) unmethylated peak
peak_cut <- 2.5
flagged  <- names(peak_low)[peak_low < peak_cut]
length(flagged)
flagged

# Highlight the flagged samples in the density plot
dens_df$flagged <- dens_df$patientID %in% flagged

ggplot(dens_df, aes(x = beta, y = density, group = patientID, colour = flagged)) +
  geom_line(aes(alpha = flagged), linewidth = 0.4) +
  scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "#D73027"),
                      labels = c("Other samples", "Flagged samples"), name = NULL) +
  scale_alpha_manual(values = c(`FALSE` = 0.15, `TRUE` = 0.8), guide = "none") +
  labs(title = "Beta value densities per sample",
       subtitle = sprintf("%d samples with unmethylated peak density < %.1f highlighted",
                          length(flagged), peak_cut),
       x = "Beta", y = "Density") +
  theme_qc()

# ---- 6d. PCA -----------------------------------------------------------------
# Main structure of the data and outliers. prcomp() does not accept NAs, so only
# probes without any missing value are used.
pc_dat <- meth[, colSums(is.na(meth)) == 0]       # all probes without any NA
dim(pc_dat)
pc <- prcomp(pc_dat)
rm(pc_dat); gc()

var_exp <- round(100 * summary(pc)$importance[2, 1:2], 1)  # % variance of PC1 and PC2

# PCA coloured by pathological stage (outcome)
pca_df <- data.frame(PC1   = pc$x[, 1],
                     PC2   = pc$x[, 2],
                     stage = clin[rownames(meth), "outcome_p.stage"])

ggplot(pca_df, aes(x = PC1, y = PC2, colour = stage)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA of methylation",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "Stage") +
  theme_qc()

# PCA coloured by race
pca_df <- data.frame(PC1   = pc$x[, 1],
                     PC2   = pc$x[, 2],
                     stage = clin[rownames(meth), "race"])

ggplot(pca_df, aes(x = PC1, y = PC2, colour = stage)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA of methylation",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "Race") +
  theme_qc()

# PCA coloured by radiation therapy
pca_df <- data.frame(PC1   = pc$x[, 1],
                     PC2   = pc$x[, 2],
                     stage = clin[rownames(meth), "radiation_therapy"])

ggplot(pca_df, aes(x = PC1, y = PC2, colour = stage)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA of methylation",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "radiation_therapy") +
  theme_qc()

# PCA coloured by histological type
pca_df <- data.frame(PC1   = pc$x[, 1],
                     PC2   = pc$x[, 2],
                     stage = clin[rownames(meth), "histological_type"])

ggplot(pca_df, aes(x = PC1, y = PC2, colour = stage)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA of methylation",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "histological_type") +
  theme_qc()

# PCA coloured by age (years_to_birth)
pca_df <- data.frame(PC1   = pc$x[, 1],
                     PC2   = pc$x[, 2],
                     stage = clin[rownames(meth), "years_to_birth"])

ggplot(pca_df, aes(x = PC1, y = PC2, colour = stage)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA of methylation",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "years_to_birth") +
  theme_qc()


# =============================================================================
# 7. BATCH INFORMATION (array plate and tissue source site)
# The plate is part of the full TCGA sample barcode, which was trimmed when the
# data were prepared. It is recovered from curatedTCGAData (cached download).
# =============================================================================
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)

# Load only the methylation assay (from cache) to recover the full sample barcodes
mae <- curatedTCGAData(diseaseCode = "BRCA", assays = "Methylation_methyl450",
                       version = "2.1.1", dry.run = FALSE)
barcodes <- colnames(mae[["BRCA_Methylation_methyl450-20160128"]])

batch_info <- data.frame(
  barcode     = barcodes,
  patientID   = substr(barcodes, 1, 12),
  sample_type = substr(barcodes, 14, 15),   # 01 = primary solid tumour
  plate       = substr(barcodes, 22, 25),   # methylation array plate
  centre      = substr(barcodes, 6, 7)      # tissue source site
)
batch_info <- batch_info[batch_info$sample_type == "01", ]           # primary tumours only
batch_info <- batch_info[!duplicated(batch_info$patientID), ]         # one per patient
rownames(batch_info) <- batch_info$patientID

batch_info <- batch_info[rownames(meth), ]   # same patients and order as meth
table(batch_info$plate)

rm(mae); gc()

# Are the samples with a shifted unmethylated peak concentrated in some plates?
tab <- table(plate = batch_info$plate, flagged = rownames(meth) %in% flagged)
tab[tab[, "TRUE"] > 0, , drop = FALSE]   # plates containing flagged samples

# PCA coloured by plate
pca_df$plate <- batch_info[rownames(pca_df), "plate"]

ggplot(pca_df, aes(x = PC1, y = PC2, colour = plate)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA of methylation coloured by array plate",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "Plate") +
  theme_qc() +
  theme(legend.position = "right")

# PCA coloured by tissue source site (centre)
pca_df$centre <- batch_info[rownames(pca_df), "centre"]

ggplot(pca_df, aes(x = PC1, y = PC2, colour = centre)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA of methylation coloured by tissue source site",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "Centre") +
  theme_qc() +
  theme(legend.position = "right")

rm(cor_mat, dens_mat, dens_df); gc()


# =============================================================================
# 8. PROBE TYPE BIAS AND BMIQ NORMALISATION
#
# !!!! BMIQ is slow: 9.331187 hours on the test computer !!!!
#   CPU : 11th Gen Intel Core i7-1185G7 @ 3.00GHz (4 cores / 8 threads, 7 used in parallel)
#   RAM : 16 GB
#   OS  : Windows, R 4.3.2
#   Set do_bmiq <- FALSE to skip it (for testing the rest of the script)
# 
# The 450k has two probe designs. Type II probes have a smaller dynamic range
# (peaks closer to the middle). BMIQ adjusts type II probes to the type I
# distribution. It is applied per sample, so it does not leak information
# between samples (the test set must also be normalised with BMIQ).
# =============================================================================
type <- anno[colnames(meth), "Type"]

# Median beta per probe BEFORE normalisation
type_before <- data.frame(median_beta = colMedians(meth, na.rm = TRUE),
                          type        = paste("Type", type),
                          step        = "Before BMIQ")

library(RPMM)
library(parallel)

# BMIQ normalisation in parallel (one sample per task; samples are independent)
do_bmiq <- TRUE
if (do_bmiq) {
  design_v <- ifelse(type == "I", 1, 2)
  
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
  
  # Split the samples into one block per core and run the blocks in parallel
  n_cores <- max(1, detectCores() - 1)           # leave one core free
  blocks  <- split(seq_len(nrow(meth)), cut(seq_len(nrow(meth)), n_cores, labels = FALSE))
  chunks  <- lapply(blocks, function(i) meth[i, , drop = FALSE])
  
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, { library(wateRmelon); library(RPMM) })
  t0 <- Sys.time()
  res <- parLapply(cl, chunks, bmiq_block, design_v = design_v)
  stopCluster(cl)
  print(Sys.time() - t0)
  
  meth <- do.call(rbind, res)[rownames(meth), ]   # same sample order as before
  rm(chunks, res); gc()
}

# Median beta per probe AFTER normalisation
type_after <- data.frame(median_beta = colMedians(meth, na.rm = TRUE),
                         type        = paste("Type", type),
                         step        = "After BMIQ")

# Combine before and after, and plot side by side
type_df <- rbind(type_before, type_after)
type_df$step <- factor(type_df$step, levels = c("Before BMIQ", "After BMIQ"))

ggplot(type_df, aes(x = median_beta, colour = type)) +
  geom_density(linewidth = 0.8) +
  facet_wrap(~ step) +
  labs(title = "Probe type bias before and after BMIQ",
       subtitle = "Median beta per probe across samples",
       x = "Median beta", y = "Density", colour = "Probe type") +
  theme_qc()


# =============================================================================
# 9. REGION-LEVEL METHYLATION
# Mean beta of the CpGs in each gene promoter, gene body and CpG island
# (missing values ignored, i.e. na.rm = TRUE).
#   Promoter  : TSS1500, TSS200, 5'UTR, 1stExon (Illumina annotation)
#   Gene body : Body, 3'UTR
#   CpG island: probes with Relation_to_Island == "Island", grouped by island
# =============================================================================
library(Matrix)

min_probes <- 1   # minimum number of CpGs in a region (from the annotation)
min_obs    <- 1   # minimum number of non-missing CpGs in a sample for that region

anno_r <- anno[colnames(meth), c("UCSC_RefGene_Name", "UCSC_RefGene_Group",
                                 "Relation_to_Island", "Islands_Name")]

# ---- 9a. Gene promoters and gene bodies (one region per gene) ----------------
# A probe can map to several genes/transcripts: fields are ";"-separated and paired
anno_g  <- anno_r[anno_r$UCSC_RefGene_Name != "", ]
genes   <- strsplit(anno_g$UCSC_RefGene_Name,  ";")
regions <- strsplit(anno_g$UCSC_RefGene_Group, ";")
probe_gene <- data.frame(probe  = rep(rownames(anno_g), lengths(genes)),
                         gene   = unlist(genes),
                         region = unlist(regions))

probe_gene$region_type <- ifelse(probe_gene$region %in% c("TSS1500", "TSS200", "5'UTR", "1stExon"), "promoter",
                                 ifelse(probe_gene$region %in% c("Body", "3'UTR"), "body", NA))
probe_gene <- probe_gene[!is.na(probe_gene$region_type), ]

# If a probe is in the promoter of one transcript and the body of another
# transcript of the SAME gene, it is counted only in the promoter
probe_gene <- probe_gene[order(probe_gene$region_type != "promoter"), ]
probe_gene <- probe_gene[!duplicated(probe_gene[, c("probe", "gene")]), ]

# ---- 9b. CpG islands (one region per island) ---------------------------------
probe_island <- data.frame(probe  = rownames(anno_r),
                           island = anno_r$Islands_Name)
probe_island <- probe_island[anno_r$Relation_to_Island == "Island" &
                               probe_island$island != "", ]

# ---- 9c. Region lists: region name -> its CpGs --------------------------------
make_list <- function(probes, region_ids) {
  rl <- split(probes, region_ids)
  rl[lengths(rl) >= min_probes]
}

region_lists <- list(
  promoter = with(probe_gene[probe_gene$region_type == "promoter", ], make_list(probe, gene)),
  body     = with(probe_gene[probe_gene$region_type == "body", ],     make_list(probe, gene)),
  island   = with(probe_island,                                       make_list(probe, island))
)
sapply(region_lists, length)   # number of promoters, gene bodies and CpG islands

# ---- 9d. Mean beta per region (na.rm = TRUE), via sparse matrix products ------
# Equivalent to rowMeans(m[, probes], na.rm = TRUE) for every region, but all
# regions are computed at once:
#   A      : probe x region indicator matrix (1 if the probe belongs to the region)
#   sums   : sum of the observed betas of each region, per sample
#   counts : number of observed CpGs of each region, per sample
#   mean   = sums / counts
region_mean <- function(m, rl) {
  A <- sparseMatrix(i    = match(unlist(rl), colnames(m)),
                    j    = rep(seq_along(rl), lengths(rl)),
                    x    = 1,
                    dims = c(ncol(m), length(rl)),
                    dimnames = list(colnames(m), names(rl)))
  
  obs  <- !is.na(m)          # TRUE where the CpG was measured
  m0   <- m
  m0[!obs] <- 0              # missing values contribute 0 to the sum
  
  sums   <- as.matrix(m0  %*% A)   # sum of observed betas per sample and region
  counts <- as.matrix(obs %*% A)   # number of observed CpGs per sample and region
  
  out <- sums / counts             # mean of observed CpGs (= rowMeans with na.rm = TRUE)
  out[counts < min_obs] <- NA      # too few observed CpGs: not reliable
  out
}

# check.names = FALSE keeps the original region names (e.g. "HLA-A", "chr1:123-456")
meth_promoter <- data.frame(region_mean(meth, region_lists$promoter), check.names = FALSE)
meth_body     <- data.frame(region_mean(meth, region_lists$body),     check.names = FALSE)
meth_island   <- data.frame(region_mean(meth, region_lists$island),   check.names = FALSE)
gc()

# ---- 9e. Summary --------------------------------------------------------------
data.frame(region    = c("Gene promoters", "Gene bodies", "CpG islands"),
           features  = c(ncol(meth_promoter), ncol(meth_body), ncol(meth_island)),
           samples   = c(nrow(meth_promoter), nrow(meth_body), nrow(meth_island)),
           pct_na    = round(100 * c(mean(is.na(meth_promoter)),
                                     mean(is.na(meth_body)),
                                     mean(is.na(meth_island))), 3))

# ---- 9f. Join the three region types in one dataset ---------------------------
meth_promoter$patientID <- rownames(meth_promoter)
meth_body$patientID     <- rownames(meth_body)
meth_island$patientID   <- rownames(meth_island)

# Add a suffix to every feature (all columns except patientID), so that the
# promoter and the gene body of the same gene can be distinguished
add_suffix <- function(df, suffix) {
  feat <- colnames(df) != "patientID"
  colnames(df)[feat] <- paste0(colnames(df)[feat], suffix)
  df
}

meth2 <- merge(add_suffix(meth_promoter, ".promoter"),
               add_suffix(meth_body,     ".genebody"), by = "patientID")
meth2 <- merge(meth2, add_suffix(meth_island, ".island"), by = "patientID")

meth2 <- meth2[match(rownames(meth), meth2$patientID), ]
rownames(meth2) <- meth2$patientID

# =============================================================================
# 10. QC OF REGION-LEVEL FEATURES
# =============================================================================

# ---- 10a. Matrix of region-level features without NAs (for the PCA) -----------
X <- as.matrix(meth2[, colnames(meth2) != "patientID"])
rownames(X) <- meth2$patientID
X <- X[, colSums(is.na(X)) == 0]          # PCA does not accept NAs

feature_type <- sub(".*\\.", "", colnames(X))   # promoter / genebody / island
table(feature_type)

# ---- 10b. PCA loadings of region-level features -------------------------------
# Centred, not scaled: all features are on the same scale (mean beta, 0-1)
pc      <- prcomp(X, center = TRUE, scale. = FALSE)
var_exp <- round(100 * summary(pc)$importance[2, 1:2], 1)

# Loadings plot: each point is a feature, coloured by region type
load_df <- data.frame(PC1  = pc$rotation[, 1],
                      PC2  = pc$rotation[, 2],
                      type = feature_type)

ggplot(load_df, aes(x = PC1, y = PC2, colour = type)) +
  geom_point(size = 0.6, alpha = 0.4) +
  scale_colour_manual(values = c(promoter = "#1A9850", genebody = "#2C7FB8", island = "#D73027")) +
  labs(title = "PCA of region-level methylation features",
       subtitle = "Each point is a feature (loadings); features far from the origin drive the component",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "Region type") +
  theme_qc() +
  theme(legend.position = "right")

# ---- 10c. Mean methylation and variability of region-level features -----------
feat_stats <- data.frame(mean = colMeans(X, na.rm = TRUE),
                         sd   = colSds(X, na.rm = TRUE),
                         type = sub(".*\\.", "", colnames(X)))

# Mean beta of each feature across samples
ggplot(feat_stats, aes(x = mean, fill = type)) +
  geom_histogram(bins = 100, colour = "white") +
  facet_wrap(~ type, ncol = 1) +
  scale_fill_manual(values = c(promoter = "#1A9850", genebody = "#2C7FB8", island = "#D73027"),
                    guide = "none") +
  labs(title = "Mean methylation of region-level features",
       subtitle = "Mean beta value of each feature across samples",
       x = "Mean beta value across samples", y = "Features") +
  theme_qc()

# SD of each feature across samples (same scale as beta values)
ggplot(feat_stats, aes(x = sd, fill = type)) +
  geom_histogram(bins = 100, colour = "white") +
  facet_wrap(~ type, ncol = 1) +
  scale_fill_manual(values = c(promoter = "#1A9850", genebody = "#2C7FB8", island = "#D73027"),
                    guide = "none") +
  labs(title = "Variability of region-level methylation features",
       subtitle = "Standard deviation of mean beta values across samples",
       x = "Standard deviation of beta values across samples", y = "Features") +
  theme_qc()


# =============================================================================
# 11. M-VALUES
# M = log2(beta / (1 - beta)): better statistical properties than beta values
# (more symmetric distribution, more constant variance). The offset avoids
# infinite values at 0 and 1.
# =============================================================================
to_M <- function(b, offset = 1e-6) log2((b + offset) / (1 - b + offset))

# All region-level features (NAs kept)
X <- as.matrix(meth2[, names(meth2) != "patientID"])
rownames(X) <- meth2$patientID

methM <- to_M(X)

# ---- 11a. Distribution of random features: beta vs M-values -------------------
# Stratify features by mean beta (tertiles) and pick 4 random features per stratum
set.seed(1)
feat_mean <- colMeans(X, na.rm = TRUE)
strata    <- cut(feat_mean, quantile(feat_mean, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                 labels = c("Low", "Medium", "High"), include.lowest = TRUE)
pick <- unlist(lapply(split(colnames(methM), strata), sample, size = 4))

# Plot the distribution across patients of the selected features
plot_feature_dist <- function(mat, x_lab, title) {
  feat_df <- data.frame(feature = rep(pick, each = nrow(mat)),
                        stratum = rep(strata[match(pick, colnames(mat))], each = nrow(mat)),
                        value   = as.vector(mat[, pick]))
  feat_df$feature <- factor(feat_df$feature, levels = pick)
  
  ggplot(feat_df, aes(value, fill = stratum)) +
    geom_histogram(bins = 25, colour = "white", na.rm = TRUE) +
    facet_wrap(~ feature, scales = "free", ncol = 4) +
    scale_fill_manual(values = c(Low = "#7FCDBB", Medium = "#2C7FB8", High = "#253494"),
                      name = "Mean methylation") +
    labs(title = title,
         subtitle = "4 random features per tertile of mean beta value",
         x = x_lab, y = "Patients") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40"),
          strip.text    = element_text(face = "bold", size = 7),
          legend.position = "top")
}

plot_feature_dist(X,             "Beta value",        "Distribution across patients of randomly selected features (beta values)")
plot_feature_dist(methM,         "M-value",           "Distribution across patients of randomly selected features (M-values)")
plot_feature_dist(scale(methM),  "z-transf(M-value)", "Distribution across patients of randomly selected features z-transf(M-values)")

# =============================================================================
# 12. SAVE
# Samples in rows, with a patientID column like the other datasets
# =============================================================================
path.x <- "./data/processed"
dir.create(path.x, showWarnings = FALSE, recursive = TRUE)

# QC parameters learned on train, to apply the same QC to test / Shiny data
# (computed before adding patientID, so only CpG names are stored)
meth_qc_params <- list(
  kept_probes  = colnames(meth),   # CpGs kept after probe filtering (before aggregation)
  region_lists = region_lists,     # CpG-to-region mapping (promoter, body, island)
  min_obs      = min_obs           # Minimum non-NA CpGs per region and sample
)

methM <- data.frame(methM, check.names = FALSE)
methM$patientID <- rownames(methM)
methM <- methM[, c(ncol(methM), 1:(ncol(methM) - 1))]   # patientID as first column

meth_df <- data.frame(meth, check.names = FALSE)
meth_df$patientID <- rownames(meth_df)
meth_df <- meth_df[, c(ncol(meth_df), 1:(ncol(meth_df) - 1))]

saveRDS(meth2,          file.path(path.x, "meth.train.region_beta.rds"))  # region-level beta values
saveRDS(methM,          file.path(path.x, "meth.train.region_M.rds"))     # region-level M-values
saveRDS(meth_df,        file.path(path.x, "meth.train.cpg_beta.rds"))     # CpG-level beta values (filtered, BMIQ)
saveRDS(meth_qc_params, file.path(path.x, "meth_qc_params.rds"))          # QC parameters for test / Shiny

gc()



