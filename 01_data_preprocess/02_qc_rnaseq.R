# =============================================================================
#  TRANSCRIPTOMICS (RNA-seq) - EXPLORATORY DATA ANALYSIS AND QUALITY CONTROL
#  TCGA-BRCA, training set
#
#  Input : rna.train.rds   (samples in rows, genes in columns, plus patientID)
#                          values: RSEM counts normalised as log2(CPM + 1)
#          clin.train.rds  (clinical data, one row per patient)
#  Output: rna.train.filtered.rds  (filtered genes, log2(CPM + 1))
#          rna_qc_params.rds       (parameters learned on train, to apply to test)
#
#  General principle: the goal of QC is not to remove what is *different*, but to
#  distinguish what is technical (to be removed or adjusted) from what is
#  biological (to be kept). 
#
#  Convention: numeric MATRIX with SAMPLES IN ROWS and GENES IN COLUMNS.
#
#  Sections
#    0. Setup: packages, plot theme and helper functions
#    1. Read data
#    2. Sanity checks
#    3. Zeros and uninformative genes (filtering)
#    4. Per-patient distributions
#    5. Patient-patient correlation (connectivity)
#    6. Gene distributions (random genes, skewness)
#    7. Gene-gene correlation (subset of most variable genes)
#    8. Extreme values per patient
#    9. Principal component analysis
#   10. QC summary per patient
#   11. Save
# =============================================================================

rm(list = ls())


# =============================================================================
# 0. SETUP
# =============================================================================
library(matrixStats)   # fast row/column statistics on matrices
library(ggplot2)

# ---- Common theme for all plots ------------------------------------------------
theme_qc <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40", size = 9),
          legend.position = "top")
}

# ---- Helper: barplot of zero proportions grouped in bins -----------------------
#   step     : bin width (0.05 = 5% bins)
#   max_frac : last bin edge; values above it are grouped in a "> x%" bin
# Colours: green = no zeros, blue = some zeros, red = mostly zeros
plot_zero_bins <- function(zero_frac, title, x_lab, step = 0.05, max_frac = 1) {
  upper  <- seq(step, max_frac, by = step)
  breaks <- c(-Inf, 0, upper, if (max_frac < 1) Inf)
  labels <- c("0%", paste0(round(100 * (upper - step)), "-", round(100 * upper), "%"),
              if (max_frac < 1) paste0("> ", round(100 * max_frac), "%"))
  df <- as.data.frame(table(bin = cut(zero_frac, breaks = breaks, labels = labels)))
  df$pct <- df$Freq / sum(df$Freq)
  
  n_bins   <- nrow(df)
  high_lab <- if (max_frac < 1) paste0("> ", round(100 * max_frac), "% zeros") else "> 90% zeros"
  is_high  <- if (max_frac < 1) seq_len(n_bins) == n_bins else c(0, upper) > 0.9
  df$group <- factor(ifelse(seq_len(n_bins) == 1, "No zeros",
                            ifelse(is_high, high_lab, "Some zeros")),
                     levels = c("No zeros", "Some zeros", high_lab))
  
  ggplot(df, aes(bin, Freq, fill = group)) +
    geom_col(width = 0.8) +
    geom_text(data = subset(df, Freq > 0),
              aes(label = sprintf("%d\n(%.0f%%)", Freq, 100 * pct)),
              vjust = -0.2, size = 2.6, lineheight = 0.85) +
    scale_fill_manual(values = c("#1A9850", "#2C7FB8", "#D73027"), name = NULL, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = title, subtitle = sprintf("n = %d", length(zero_frac)),
         x = x_lab, y = "Count") +
    theme_qc() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ---- Helper: robust z-score (median and MAD instead of mean and SD) ------------
robust_z <- function(x) (x - median(x)) / mad(x)


# =============================================================================
# 1. READ DATA
# =============================================================================
path.x <- "./data/data_split"

df <- readRDS(file.path(path.x, "rna.train.rds"))
X  <- as.matrix(df[, colnames(df) != "patientID"])   # samples x genes
rownames(X) <- df$patientID
rm(df); gc()

clin <- readRDS(file.path(path.x, "clin.train.rds"))
clin <- clin[rownames(X), ]                         # same patients, same order

class(X); dim(X)


# =============================================================================
# 2. SANITY CHECKS
# Expected for log2(CPM + 1): no missing values, no negatives, range 0 to ~20,
# non-integer values (except zeros, since log2(0 + 1) = 0).
# =============================================================================
data.frame(
  check = c("Patients", "Genes", "Missing values", "Duplicated patients",
            "Duplicated genes", "Minimum value", "Maximum value", "Negative values",
            "Proportion of non-integer values (among non-zero values)",
            "Overall proportion of zeros"),
  value = c(nrow(X), ncol(X), sum(is.na(X)), sum(duplicated(rownames(X))),
            sum(duplicated(colnames(X))), round(min(X), 2), round(max(X), 2), sum(X < 0),
            round(mean(X[X != 0] != round(X[X != 0])), 3), round(mean(X == 0), 3))
)


# =============================================================================
# 3. ZEROS AND UNINFORMATIVE GENES
# Since there are no missing values, zeros play the role of "not detected".
# =============================================================================
zero_gene   <- colMeans(X == 0)   # proportion of patients with zero, per gene
zero_sample <- rowMeans(X == 0)   # proportion of genes with zero, per patient

# ---- 3a. Zeros per gene ---------------------------------------------------------
# Expected: asymmetric U shape (most genes without zeros, a group of genes not
# expressed in breast tissue, and "on/off" genes in between)
plot_zero_bins(zero_gene, "Zeros per gene", "Percentage of patients with value = 0")

# ---- 3b. Zeros per patient ------------------------------------------------------
# All patients should have a similar proportion of zeros
plot_zero_bins(zero_sample, "Zeros per patient", "Percentage of genes with value = 0",
               step = 0.01, max_frac = 0.5)

# ---- 3c. Expression level of detected values -----------------------------------
# Mean of the non-zero values of each gene: how strongly a gene is expressed
# when it is detected. A bimodal distribution separates low-expression noise from
# clearly expressed genes; the valley between the two modes is used as threshold.
nz_mean <- apply(X, 2, function(x) if (any(x > 0)) mean(x[x > 0]) else 0)

# Data-driven valley: minimum of the density between the two modes.
# Check it in the plot and adjust the search range or the value if needed.
d      <- density(nz_mean)
search <- d$x > 1 & d$x < 7
nz_cut <- round(d$x[search][which.min(d$y[search])], 2)
nz_cut

ggplot(data.frame(nz_mean), aes(nz_mean)) +
  geom_histogram(bins = 60, fill = "#2C7FB8", colour = "white") +
  geom_vline(xintercept = nz_cut, linetype = 2, colour = "#D73027") +
  labs(title = "Mean expression of non-zero values per gene",
       subtitle = sprintf("Dashed red line = valley of the distribution (%.2f)", nz_cut),
       x = "Mean of non-zero log2(CPM + 1)", y = "Genes") +
  theme_qc()

# ---- 3d. Filtering --------------------------------------------------------------
# Genes are removed if they carry no information or look like detection noise:
#   Zero variance                       : SD = 0
#   Mostly zeros                        : zeros in > 90% of patients
#   Low expression AND frequent zeros   : non-zero mean < nz_cut and zeros in >= 10%
# Genes with low but consistent expression (detected in > 90% of patients) are
# kept. No filtering on variance is applied at this step.
gene_sd      <- colSds(X)
rm_zero_var  <- gene_sd == 0
rm_high_zero <- zero_gene > 0.9
rm_low_expr  <- nz_mean < nz_cut & zero_gene >= 0.1
remove       <- rm_zero_var | rm_high_zero | rm_low_expr

# Record of removed genes and the reason
filter_log <- data.frame(gene = colnames(X), nz_mean = nz_mean, zero_prop = zero_gene,
                         zero_var = rm_zero_var, high_zero = rm_high_zero,
                         low_expr = rm_low_expr, removed = remove)

X_filt <- X[, !remove]

data.frame(rule  = c("Zero variance", "Zeros in > 90% of patients",
                     "Low expression and frequent zeros", "Total removed (union)",
                     "Genes kept"),
           genes = c(sum(rm_zero_var), sum(rm_high_zero), sum(rm_low_expr),
                     sum(remove), ncol(X_filt)))


# =============================================================================
# 4. PER-PATIENT DISTRIBUTIONS
# With normalised data, all patients should have similar distributions. A
# patient whose median or spread is far from the rest (|robust z| > 3) is flagged.
# =============================================================================
sample_stats <- data.frame(sample        = rownames(X_filt),
                           sample_median = rowMedians(X_filt),
                           iqr           = rowIQRs(X_filt))
sample_stats$z_median <- robust_z(sample_stats$sample_median)
sample_stats$z_iqr    <- robust_z(sample_stats$iqr)
sample_stats$flag     <- abs(sample_stats$z_median) > 3 | abs(sample_stats$z_iqr) > 3

long <- data.frame(sample = rep(rownames(X_filt), times = ncol(X_filt)),
                   value  = as.vector(X_filt))
long$sample <- factor(long$sample,
                      levels = sample_stats$sample[order(sample_stats$sample_median)])
long$flag   <- sample_stats$flag[match(long$sample, sample_stats$sample)]

# ---- 4a. Boxplot per patient ---------------------------------------------------
ggplot(long, aes(sample, value, fill = flag)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = c(`FALSE` = "grey80", `TRUE` = "#D73027"), guide = "none") +
  labs(title = "Expression distribution per patient",
       subtitle = "Red = unusual median or spread (|robust z| > 3)",
       x = "Patient (ordered by median)", y = "log2(CPM + 1)") +
  theme_qc() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ---- 4b. Density per patient ---------------------------------------------------
ggplot(long, aes(value, group = sample, colour = flag)) +
  geom_density(linewidth = 0.3) +
  scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#D73027"), guide = "none") +
  labs(title = "Density of expression values per patient",
       x = "log2(CPM + 1)", y = "Density") +
  theme_qc()

subset(sample_stats, flag)
rm(long); gc()


# =============================================================================
# 5. PATIENT-PATIENT CORRELATION
# Standardized connectivity (Z.k, Oldham et al. 2012): correlations are turned
# into connection strengths ((1 + r) / 2)^2, summed per patient and standardised.
# Patients with Z.k < -2 are weakly connected to the rest (candidate outliers).
# =============================================================================
C  <- cor(t(X_filt))                 # patient x patient
A  <- ((1 + C) / 2)^2
k  <- rowSums(A) - 1
zk <- (k - mean(k)) / sd(k)
med_cor <- apply(C, 1, function(x) median(x[x < 1]))

cor_stats <- data.frame(sample = rownames(X_filt), med_cor = med_cor, zk = zk)
cor_stats$group  <- ifelse(cor_stats$zk < -2, "Z.k < -2", "OK")
cor_stats$sample <- factor(cor_stats$sample, levels = cor_stats$sample[order(cor_stats$zk)])

ggplot(cor_stats, aes(sample, zk, colour = group)) +
  geom_point(size = 2) +
  geom_hline(yintercept = -2, linetype = 2, colour = "#D73027") +
  scale_colour_manual(values = c("OK" = "#2C7FB8", "Z.k < -2" = "#D73027"), name = NULL) +
  labs(title = "Standardized connectivity per patient",
       subtitle = "Patients below -2 are weakly correlated with the rest",
       x = "Patient (ordered by Z.k)", y = "Z.k") +
  theme_qc() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

summary(med_cor)
subset(cor_stats, zk < -2, select = c(sample, med_cor, zk))


# ---- 5. Correlation heatmap (hierarchical clustering order) --------------------
hc  <- hclust(as.dist(1 - C), method = "average")
ord <- rownames(C)[hc$order]
heat_df <- data.frame(s1 = factor(rep(rownames(C), times = ncol(C)), levels = ord),
                      s2 = factor(rep(colnames(C), each = nrow(C)), levels = ord),
                      r  = as.vector(C))

ggplot(heat_df, aes(s1, s2, fill = r)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Pearson r") +
  labs(title = "Patient-patient correlation (hierarchical clustering order)", x = NULL, y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text = element_blank(), panel.grid = element_blank())

rm(heat_df, A); gc()


# =============================================================================
# 6. GENE DISTRIBUTIONS
# =============================================================================

# ---- 6a. Randomly selected genes ------------------------------------------------
# Four random genes from each tertile of mean expression
set.seed(20261013)
gene_mean <- colMeans(X_filt)
strata    <- cut(gene_mean, quantile(gene_mean, c(0, 1/3, 2/3, 1)),
                 labels = c("Low", "Medium", "High"), include.lowest = TRUE)
pick <- unlist(lapply(split(colnames(X_filt), strata), sample, size = 4))

gene_df <- data.frame(gene    = factor(rep(pick, each = nrow(X_filt)), levels = pick),
                      stratum = rep(strata[match(pick, colnames(X_filt))], each = nrow(X_filt)),
                      value   = as.vector(X_filt[, pick]))

ggplot(gene_df, aes(value, fill = stratum)) +
  geom_histogram(bins = 25, colour = "white") +
  facet_wrap(~ gene, scales = "free", ncol = 4) +
  scale_x_continuous(n.breaks = 3) +
  scale_fill_manual(values = c(Low = "#7FCDBB", Medium = "#2C7FB8", High = "#253494"),
                    name = "Mean expression") +
  labs(title = "Distribution across patients of randomly selected genes",
       subtitle = "4 random genes per expression tertile",
       x = "log2(CPM + 1)", y = "Patients") +
  theme_qc(11) +
  theme(strip.text = element_text(face = "bold"))

# ---- 6b. Skewness -----------------------------------------------------------------
# Reported descriptively. Red lines = Tukey "far out" fences (Q1 - 3 IQR,
# Q3 + 3 IQR) on the skewness of all genes.
gene_skew <- apply(X_filt, 2, function(x) mean((x - mean(x))^3) / sd(x)^3)
q    <- quantile(gene_skew, c(0.25, 0.75))
lims <- c(q[1] - 3 * diff(q), q[2] + 3 * diff(q))

ggplot(data.frame(skew = gene_skew), aes(skew)) +
  geom_histogram(bins = 60, fill = "#2C7FB8", colour = "white") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_vline(xintercept = lims, linetype = 2, colour = "#D73027") +
  labs(title = "Skewness of gene distributions",
       subtitle = sprintf("Grey = symmetric; red = Tukey far-out fences (%.2f, %.2f)",
                          lims[1], lims[2]),
       x = "Skewness", y = "Genes") +
  theme_qc()

head(round(sort(gene_skew[gene_skew < lims[1] | gene_skew > lims[2]]), 2), 30)


# =============================================================================
# 7. GENE-GENE CORRELATION
# A gene x gene matrix with all genes (~15,000-20,000) would need several GB of
# memory, so the 2,000 most variable genes are used. A permutation null (each
# gene shuffled independently) gives the chance level of correlation.
# =============================================================================
n_top     <- 2000
top_genes <- colnames(X_filt)[order(colVars(X_filt), decreasing = TRUE)[1:n_top]]
Xg        <- X_filt[, top_genes]

Cg <- cor(Xg)
r_pairs <- abs(Cg[upper.tri(Cg)])
diag(Cg) <- NA
max_r <- apply(abs(Cg), 1, max, na.rm = TRUE)

set.seed(1)
X_perm <- apply(Xg, 2, sample)
Cp <- cor(X_perm)
r_pairs_null <- abs(Cp[upper.tri(Cp)])
diag(Cp) <- NA
max_r_null <- apply(abs(Cp), 1, max, na.rm = TRUE)

cut_pairs <- quantile(r_pairs_null, 0.95)
cut_max   <- quantile(max_r_null, 0.95)

# Barplot of |r| in bins of 0.1, with the chance level as a dashed line
plot_cor_bins <- function(r, null_cut, title, y_lab) {
  breaks <- seq(0, 1, by = 0.1)
  df <- data.frame(mid  = head(breaks, -1) + 0.05,
                   Freq = as.vector(table(cut(r, breaks = breaks, include.lowest = TRUE))))
  df$pct   <- df$Freq / sum(df$Freq)
  df$level <- factor(ifelse(df$mid < 0.5, "|r| < 0.5", ifelse(df$mid < 0.8, "0.5 - 0.8", "> 0.8")),
                     levels = c("|r| < 0.5", "0.5 - 0.8", "> 0.8"))
  
  ggplot(df, aes(mid, Freq, fill = level)) +
    geom_col(width = 0.08) +
    geom_text(data = subset(df, Freq > 0),
              aes(label = sprintf("%d\n(%.1f%%)", Freq, 100 * pct)),
              vjust = -0.2, size = 2.6, lineheight = 0.85) +
    geom_vline(xintercept = null_cut, linetype = 2, colour = "#D73027") +
    scale_fill_manual(values = c("#1A9850", "#2C7FB8", "#D73027"), name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = breaks, limits = c(0, 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = title,
         subtitle = sprintf("n = %d | dashed red line = chance level (95th percentile of permuted data = %.2f)",
                            length(r), null_cut),
         x = "|Pearson r|", y = y_lab) +
    theme_qc()
}

plot_cor_bins(r_pairs, cut_pairs,
              sprintf("Correlation between gene pairs (%d most variable genes)", n_top), "Gene pairs")
plot_cor_bins(max_r, cut_max,
              "Strongest correlation of each gene with any other gene", "Genes")

# Summary: proportion above chance and theoretical chance level 1.96 / sqrt(n - 3)
data.frame(measure = c("% gene pairs above chance level",
                       "% genes with a partner above chance level",
                       "Chance level (permutation)", "Chance level (theoretical)"),
           value   = round(c(100 * mean(r_pairs > cut_pairs), 100 * mean(max_r > cut_max),
                             cut_pairs, 1.96 / sqrt(nrow(Xg) - 3)), 3))

# Most correlated gene pairs
top_pairs <- which(abs(Cg) > 0.9 & upper.tri(Cg), arr.ind = TRUE)
head(data.frame(gene1 = rownames(Cg)[top_pairs[, 1]],
                gene2 = colnames(Cg)[top_pairs[, 2]],
                r     = round(Cg[top_pairs], 3)), 20)

rm(Cg, Cp, X_perm, Xg); gc()


# =============================================================================
# 8. EXTREME VALUES PER PATIENT
# For each gene a robust z-score is computed; values with |robust z| > 5 are
# extreme. The number of extreme values per patient is compared with a Tukey
# far-out fence (Q3 + 3 IQR). Extreme expression can be biological
# (amplifications, deletions, silencing), so values are flagged, not removed.
# =============================================================================
med <- colMedians(X_filt)
mdv <- colMads(X_filt)
ok  <- mdv > 0                      # genes with MAD = 0 cannot be standardised
Z   <- sweep(sweep(X_filt[, ok], 2, med[ok]), 2, mdv[ok], "/")

extreme   <- abs(Z) > 5
n_extreme <- rowSums(extreme)
q         <- quantile(n_extreme, c(0.25, 0.75))
fence     <- q[2] + 3 * diff(q)

ext_df <- data.frame(sample = names(n_extreme), n = n_extreme)
ext_df$sample <- factor(ext_df$sample, levels = ext_df$sample[order(ext_df$n)])
ext_df$flag   <- ext_df$n > fence

ggplot(ext_df, aes(sample, n, fill = flag)) +
  geom_col() +
  geom_hline(yintercept = fence, linetype = 2, colour = "#D73027") +
  scale_fill_manual(values = c(`FALSE` = "#2C7FB8", `TRUE` = "#D73027"), guide = "none") +
  labs(title = "Extreme values per patient",
       subtitle = sprintf("Genes with |robust z| > 5; dashed red line = Q3 + 3 IQR (%.1f)", fence),
       x = "Patient (ordered)", y = "Number of genes") +
  theme_qc() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

flagged_ext <- as.character(ext_df$sample[ext_df$flag])
data.frame(patient   = flagged_ext,
           n_extreme = n_extreme[flagged_ext],
           row.names = NULL)

# Overall proportion of extreme values
c(n_extreme_values = sum(extreme), pct = round(100 * mean(extreme), 3))

rm(Z, extreme); gc()


# =============================================================================
# 9. PRINCIPAL COMPONENT ANALYSIS
# Centred but not scaled, so the most variable genes (where most biological
# signal lies) carry more weight.
# =============================================================================
pca     <- prcomp(X_filt, center = TRUE, scale. = FALSE)
var_exp <- pca$sdev^2 / sum(pca$sdev^2)

# ---- 9a. Variance explained -----------------------------------------------------
ggplot(data.frame(pc = 1:15, var = 100 * var_exp[1:15]), aes(pc, var)) +
  geom_col(fill = "#2C7FB8") +
  geom_text(aes(label = sprintf("%.1f", var)), vjust = -0.4, size = 3) +
  scale_x_continuous(breaks = 1:15) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Variance explained per principal component",
       x = "Principal component", y = "Variance explained (%)") +
  theme_qc()

# ---- 9b. PCA data frame with QC metrics and biological scores --------------------

pca_df <- data.frame(sample        = rownames(X_filt),
                     PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
                     sample_median = sample_stats$sample_median[match(rownames(X_filt), sample_stats$sample)],
                     zeros         = rowMeans(X_filt == 0),
                     zk            = zk[rownames(X_filt)],
                     n_extreme     = n_extreme[rownames(X_filt)],
                     centre        = substr(rownames(X_filt), 6, 7),   # tissue source site
                     stage         = clin[rownames(X_filt), "outcome_p.stage"])
pca_df$flag <- pca_df$zk < -2 | pca_df$n_extreme > fence

# Helper: PCA scatter coloured by any variable (red circles = flagged patients)
plot_pca_by <- function(df, var, x = "PC1", y = "PC2") {
  p <- ggplot(df, aes(.data[[x]], .data[[y]], colour = .data[[var]])) +
    geom_point(size = 2.3, alpha = 0.85) +
    geom_point(data = subset(df, flag), shape = 21, size = 4.3,
               colour = "#D73027", stroke = 1) +
    labs(title = sprintf("PCA coloured by %s", var),
         subtitle = "Red circles = patients flagged in previous steps",
         x = sprintf("%s (%.1f%%)", x, 100 * var_exp[as.integer(sub("PC", "", x))]),
         y = sprintf("%s (%.1f%%)", y, 100 * var_exp[as.integer(sub("PC", "", y))]),
         colour = var) +
    theme_qc()
  if (is.numeric(df[[var]])) p + scale_colour_viridis_c() else p
}

# ---- 9c. What do the main components represent? ---------------------------------
plot_pca_by(pca_df, "stage")                              # outcome

# Top 10 genes (absolute loadings) of the first three components
top_load <- function(pc, n = 10) {
  l <- pca$rotation[, pc]
  names(l[order(-abs(l))][1:n])
}
data.frame(PC1 = top_load(1), PC2 = top_load(2), PC3 = top_load(3))

# ---- 9d. Association with QC metrics and centre ---------------------------------
n_pc <- 5
round(cor(pca$x[, 1:n_pc],
          pca_df[, c("sample_median", "zeros", "zk", "n_extreme")],
          method = "spearman"), 2)

# Kruskal-Wallis test of each PC by centre (centres with at least 3 patients)
big_centres <- names(which(table(pca_df$centre) >= 3))
in_big      <- pca_df$centre %in% big_centres
data.frame(PC = paste0("PC", 1:n_pc),
           p_centre = sapply(1:n_pc, function(i)
             signif(kruskal.test(pca$x[in_big, i] ~ factor(pca_df$centre[in_big]))$p.value, 2)))

ggplot(subset(pca_df, in_big), aes(reorder(centre, PC3, FUN = stats::median), PC3)) +
  geom_boxplot(fill = "grey85", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, colour = "#2C7FB8") +
  labs(title = "PC3 by tissue source site",
       subtitle = "Centres with at least 3 patients",
       x = "Centre", y = "PC3") +
  theme_qc() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- 9e. Sensitivity to the most atypical patient --------------------------------
# PCA with and without the patient with the lowest connectivity: if the
# loadings are almost identical, this patient does not drive the PCA
most_atypical <- names(which.min(zk))
most_atypical

pca_wo <- prcomp(X_filt[rownames(X_filt) != most_atypical, ], center = TRUE, scale. = FALSE)
data.frame(PC = paste0("PC", 1:n_pc),
           abs_cor_loadings = round(abs(diag(cor(pca$rotation[, 1:n_pc],
                                                 pca_wo$rotation[, 1:n_pc]))), 3))
rm(pca_wo); gc()


# =============================================================================
# 10. QC SUMMARY PER PATIENT
# Patients flagged in at least one step. They are flagged, not removed:
# downstream analyses can be repeated with and without them.
# =============================================================================
qc_table <- data.frame(
  patient           = rownames(X_filt),
  flag_distribution = sample_stats$flag[match(rownames(X_filt), sample_stats$sample)],
  flag_connectivity = zk[rownames(X_filt)] < -2,
  flag_extreme      = n_extreme[rownames(X_filt)] > fence
)
qc_table$n_flags <- rowSums(qc_table[, -1])
qc_flagged <- qc_table[qc_table$n_flags > 0, ]
qc_flagged[order(-qc_flagged$n_flags), ]


# =============================================================================
# 11. SAVE
# Samples in rows, with a patientID column like the other datasets
# =============================================================================
path.x <- "./data/processed"
dir.create(path.x, showWarnings = FALSE, recursive = TRUE)

rna_df <- data.frame(X_filt, check.names = FALSE)
rna_df$patientID <- rownames(rna_df)
rna_df <- rna_df[, c(ncol(rna_df), 1:(ncol(rna_df) - 1))]   # patientID as first column

saveRDS(rna_df, file.path(path.x, "rna.train.filtered.rds"))          # filtered genes, log2(CPM + 1)

# Parameters learned on the training set, to apply the same filtering to the test set
saveRDS(list(kept_genes = colnames(X_filt),   # genes kept after filtering
             nz_cut     = nz_cut,             # threshold of the non-zero mean
             filter_log = filter_log,         # reason for removal of each gene
             qc_flagged = qc_flagged),        # patients flagged in the QC
        file.path(path.x, "rna_qc_params.rds"))

gc()
