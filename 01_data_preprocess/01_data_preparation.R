# =============================================================================
#  DATA DOWNLOAD, ALIGNMENT AND TRAIN/TEST SPLIT
#  TCGA-BRCA multi-omics (RNA-seq, DNA methylation 450k, RPPA) + clinical data
#
#  Input : downloaded from curatedTCGAData (version 2.1.1)
#  Output: data_split/rna.train.rds,  data_split/rna.test.rds
#          data_split/meth.train.rds, data_split/meth.test.rds
#          data_split/prot.train.rds, data_split/prot.test.rds
#          data_split/clin.train.rds, data_split/clin.test.rds
#
#  General principle: the split is done BEFORE any data-driven preprocessing
#  (filtering, imputation, scaling), so that every QC decision is learned on the
#  training set only and then applied to the test set. The only transformation
#  applied here is log2(CPM + 1), which is computed sample by sample and
#  therefore does not leak information between train and test.
#
#  Inclusion criteria: primary solid tumours, female patients, data available in
#  all three omics, and a known pathological stage (stage x excluded).
#
#  Outcome: pathological stage in three categories
#    stage i   = stage i, ia, ib
#    stage ii  = stage ii, iia, iib
#    stage iii = stage iiia, iiib, iiic, iv
#
#  Convention (output): data.frames with PATIENTS IN ROWS and FEATURES IN
#  COLUMNS, plus a patientID column. Row names are 12-character patient
#  barcodes (e.g. TCGA-A1-A0SB).
#
#  Sections
#    0. Install packages (run once)
#    1. Load libraries
#    2. Check available assays
#    3. Download the selected assays
#    4. Keep primary tumours and patients with all three omics
#    5. Extract matrices and clinical data
#   5b. RNA-seq: raw counts to log2(CPM + 1)
#    6. Trim sample barcodes to patient IDs
#    7. Filter clinical data, define outcome and align all datasets
#    8. Train/test split (70/30)
#    9. Save
# =============================================================================

rm(list = ls())

# =============================================================================
# 0. INSTALL PACKAGES (run once)
# =============================================================================
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install(c("curatedTCGAData", "MultiAssayExperiment", "TCGAutils"))


# =============================================================================
# 1. LOAD LIBRARIES
# =============================================================================
library(curatedTCGAData)        # curated TCGA data as MultiAssayExperiment objects
library(MultiAssayExperiment)   # container linking several omics to the same patients
library(TCGAutils)              # helpers for TCGA barcodes and sample types


# =============================================================================
# 2. CHECK AVAILABLE ASSAYS (no download)
# =============================================================================
# dry.run = TRUE only lists the assays available for BRCA in this version
curatedTCGAData(diseaseCode = "BRCA", assays = "*",
                version = "2.1.1", dry.run = TRUE)


# =============================================================================
# 3. DOWNLOAD THE SELECTED ASSAYS
# =============================================================================
#   RNASeq2Gene           : RSEM gene-level counts (raw, non-negative)
#   Methylation_methyl450 : Illumina 450k beta values (0-1)
#   RPPAArray             : reverse-phase protein array (log2, median-centred)
# The first download can take a while; files are cached locally afterwards.
brca <- curatedTCGAData(
  diseaseCode = "BRCA",
  assays      = c("RNASeq2Gene", "Methylation_methyl450", "RPPAArray"),
  version     = "2.1.1",
  dry.run     = FALSE
)

brca
names(experiments(brca))   # exact assay names used below


# =============================================================================
# 4. KEEP PRIMARY TUMOURS AND PATIENTS WITH ALL THREE OMICS
# =============================================================================
# Sample type is encoded in barcode positions 14-15
brca <- TCGAprimaryTumors(brca) # keep: 01 = primary solid tumour; remove: 06 = metastatic, 11 = solid tissue normal (adjacent normal), etc.
brca <- intersectColumns(brca)  # keep: participants with data in the three datasets

lengths(colnames(brca))         # samples per assay (should be equal)


# =============================================================================
# 5. EXTRACT MATRICES AND CLINICAL DATA
# =============================================================================
# Check names with names(experiments(brca)) and adjust if they differ
# All omics matrices are FEATURES IN ROWS x SAMPLES IN COLUMNS at this point
counts <- assay(brca[["BRCA_RNASeq2Gene-20160128"]])
meth   <- as.matrix(assay(brca[["BRCA_Methylation_methyl450-20160128"]]))  # HDF5 -> in-memory matrix
prot   <- assay(brca[["BRCA_RPPAArray-20160128"]])
clin   <- as.data.frame(colData(brca))                                     # one row per patient

dim(counts); dim(meth); dim(prot); dim(clin)


# -----------------------------------------------------------------------------
# 5b. RNA-SEQ: CHECK RAW COUNTS AND NORMALISE TO log2(CPM + 1)
# -----------------------------------------------------------------------------
# Raw counts should be non-negative, without NAs, and (almost) integers.
# RSEM expected counts can contain some non-integer values, which is expected.
data.frame(check = c("Minimum value", "Maximum value", "Negative values",
                     "Missing values", "Proportion of integer values"),
           value = c(min(counts, na.rm = TRUE), max(counts, na.rm = TRUE),
                     sum(counts < 0, na.rm = TRUE), sum(is.na(counts)),
                     round(mean(counts == round(counts), na.rm = TRUE), 3)))

# CPM corrects for sequencing depth; it uses only each sample's own total,
# so it is safe to compute before the split (no train/test leakage).
lib_size <- colSums(counts)                          # total counts per sample
rna      <- log2(t(t(counts) / lib_size * 1e6) + 1)  # log2(CPM + 1), genes x samples
range(rna)                                           # expected: 0 to ~20, no negatives
rm(counts); gc()                                     # free memory

# Optional: inspect the sample-to-patient mapping
smap <- as.data.frame(sampleMap(brca))
head(smap)


# =============================================================================
# 6. TRIM SAMPLE BARCODES TO PATIENT IDs
# =============================================================================
# Omics columns use full sample barcodes (e.g. TCGA-A1-A0SB-01A-11R-A144-07),
# while clinical rows use patient IDs (first 12 characters: TCGA-A1-A0SB).
colnames(rna)  <- TCGAbarcode(colnames(rna))
colnames(meth) <- TCGAbarcode(colnames(meth))
colnames(prot) <- TCGAbarcode(colnames(prot))

# Check duplicated patients (multiple aliquots of the same tumour)
sum(duplicated(colnames(rna)))
sum(duplicated(colnames(meth)))
sum(duplicated(colnames(prot)))

# Keep the first aliquot per patient
rna  <- rna[,  !duplicated(colnames(rna))]
meth <- meth[, !duplicated(colnames(meth))]
prot <- prot[, !duplicated(colnames(prot))]


# =============================================================================
# 7. FILTER CLINICAL DATA, DEFINE OUTCOME AND ALIGN ALL DATASETS
# =============================================================================
# Remove male participants (male breast cancer is rare and biologically distinct)
clin <- clin[which(clin$gender %in% "female"), ]

# Remove participants without pathological stage (NA or "stage x" = not assessable)
clin <- clin[!is.na(clin$pathologic_stage) & !clin$pathologic_stage %in% "stage x", ]

# Define outcome: pathological stage in three categories
#   i   --> i + ia + ib
#   ii  --> ii + iia + iib
#   iii --> iiia + iiib + iiic + iv   (all remaining stages)
clin$outcome_p.stage <- ifelse(clin$pathologic_stage %in% c("stage i", "stage ia", "stage ib"), "stage i",
                               ifelse(clin$pathologic_stage %in% c("stage ii", "stage iia", "stage iib"), "stage ii",
                                      "stage iii"))

# Check the recoding: each original stage must map to exactly one category
table(clin$pathologic_stage, clin$outcome_p.stage, useNA = "ifany")

# Patients present in all four datasets (after the clinical filters)
common <- Reduce(intersect, list(colnames(rna), colnames(meth),
                                 colnames(prot), rownames(clin)))
length(common)

# Same patients, same order, in every dataset
rna  <- rna[,  common]
meth <- meth[, common]
prot <- prot[, common]
clin <- clin[common, ]

# Sanity checks: all should be TRUE
all(colnames(rna)  == rownames(clin))
all(colnames(meth) == rownames(clin))
all(colnames(prot) == rownames(clin))

# Transpose to PATIENTS IN ROWS x FEATURES IN COLUMNS (standard ML layout)
# check.names = FALSE keeps original feature names (e.g. "HER2_pY1248", "TP53")
rna  <- data.frame(t(rna),  check.names = FALSE, stringsAsFactors = FALSE)
meth <- data.frame(t(meth), check.names = FALSE, stringsAsFactors = FALSE)
prot <- data.frame(t(prot), check.names = FALSE, stringsAsFactors = FALSE)

# Add patient ID as an explicit column (useful for merging later)
rna$patientID  <- rownames(rna)
meth$patientID <- rownames(meth)
prot$patientID <- rownames(prot)


# =============================================================================
# 8. TRAIN/TEST SPLIT (70/30)
# =============================================================================
# Seed for a reproducible train/test split
set.seed(20261014)
# Random split at patient level; the same patients go to train or test in all
# four datasets. Stage proportions are checked below.
train_ids <- sample(rownames(clin), size = round(0.7 * nrow(clin)))
test_ids  <- setdiff(rownames(clin), train_ids)

# Apply the same split to every dataset
rna.train  <- rna[train_ids, ]
meth.train <- meth[train_ids, ]
prot.train <- prot[train_ids, ]
clin.train <- clin[train_ids, ]

rna.test   <- rna[test_ids, ]
meth.test  <- meth[test_ids, ]
prot.test  <- prot[test_ids, ]
clin.test  <- clin[test_ids, ]

# Sanity checks: no overlap and identical row order across blocks (all TRUE)
length(intersect(train_ids, test_ids)) == 0

all(rownames(rna.train) == rownames(clin.train))
all(rownames(rna.train) == rownames(meth.train))
all(rownames(rna.train) == rownames(prot.train))

all(rownames(rna.test) == rownames(clin.test))
all(rownames(rna.test) == rownames(meth.test))
all(rownames(rna.test) == rownames(prot.test))

# Stage proportions in each set (should be similar)
prop.table(table(clin.train$outcome_p.stage))
prop.table(table(clin.test$outcome_p.stage))


# =============================================================================
# 9. SAVE
# =============================================================================
path.x <- "./data/data_split"
dir.create(path.x, showWarnings = FALSE, recursive = TRUE)

# Save each object to its own .rds file (e.g. rna.train.rds)
objs <- list(rna.train  = rna.train,  meth.train = meth.train,
             prot.train = prot.train, clin.train = clin.train,
             rna.test   = rna.test,   meth.test  = meth.test,
             prot.test  = prot.test,  clin.test  = clin.test)

for (nm in names(objs)) {
  saveRDS(objs[[nm]], file = file.path(path.x, paste0(nm, ".rds")))
}
