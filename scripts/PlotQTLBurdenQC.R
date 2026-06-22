library(data.table)
library(tidyverse)
library(optparse)

safe_ratio <- function(numerator, denominator) {
  numerator <- as.numeric(numerator)
  denominator <- as.numeric(denominator)
  out <- rep(NA_real_, length(numerator))
  valid <- is.finite(numerator) & is.finite(denominator) & denominator != 0
  out[valid] <- numerator[valid] / denominator[valid]
  out
}

option_list <- list(
  make_option(c("--QCFile"), type = "character", default = NULL, help = "Gene QC table (QTLGeneBurdenQC.tsv.gz)"),
  make_option(c("--MedianGenesPerBinByGeneSet"), type = "character", default = NULL,
              help = "Per-bin median gene counts for different removed gene-set scenarios"),
  make_option(c("--MissingnessWarnThreshold"), type = "double", default = 0.05,
              help = "Warn threshold for missingness fraction"),
  make_option(c("--MissingnessFailThreshold"), type = "double", default = 0.10,
              help = "Fail threshold for missingness fraction"),
  make_option(c("--ContextMissingnessFailThreshold"), type = "double", default = 0.05,
              help = "Fail threshold for context-missingness fraction"),
  make_option(c("--VarianceRatioWarnLower"), type = "double", default = 0.35,
              help = "Warn lower bound for empirical/population variance ratio"),
  make_option(c("--VarianceRatioWarnUpper"), type = "double", default = 3,
              help = "Warn upper bound for empirical/population variance ratio"),
  make_option(c("--VarianceRatioLower"), type = "double", default = 0.20,
              help = "Fail lower bound for empirical/population variance ratio"),
  make_option(c("--VarianceRatioUpper"), type = "double", default = 5,
              help = "Fail upper bound for empirical/population variance ratio"),
  make_option(c("--TailZWarnThreshold"), type = "double", default = 15,
              help = "Warn threshold for max centered burden z"),
  make_option(c("--TailZFailThreshold"), type = "double", default = 25,
              help = "Fail threshold for max centered burden z"),
  make_option(c("--DominantWarnThreshold"), type = "double", default = 0.98,
              help = "Warn threshold for dominant variant fraction"),
  make_option(c("--Prefix"), type = "character", default = "QTLGeneBurdenQC", help = "Output filename prefix")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (is.null(opt$QCFile)) {
  stop("--QCFile is required")
}

qc <- fread(opt$QCFile) %>%
  as_tibble()

required_cols <- c(
  "gene_type", "QC_Status", "fraction_samples_with_missing_genotype",
  "fraction_samples_without_context", "median_variance_ratio",
  "max_abs_centered_z_population", "dominant_variant_fraction_effect",
  "QC_FailReasons", "QC_WarnReasons"
)
missing_cols <- setdiff(required_cols, names(qc))
if (length(missing_cols) > 0) {
  stop("QC table is missing required columns: ", paste(missing_cols, collapse = ", "))
}

qc <- qc %>%
  mutate(
    QC_Status = factor(QC_Status, levels = c("PASS", "WARN", "FAIL"))
  )

plot_dir <- getwd()

theme_burden <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

out <- function(name) {
  file.path(plot_dir, paste0(opt$Prefix, "_", name))
}

status_counts <- qc %>%
  count(gene_type, QC_Status, name = "n") %>%
  mutate(gene_type = replace_na(as.character(gene_type), "unknown"))

g_status <- ggplot(status_counts, aes(x = gene_type, y = n, fill = QC_Status)) +
  geom_col(position = "fill", color = "white") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Gene type",
    y = "Proportion",
    fill = "QC status",
    title = "Gene-level QC status by gene type"
  ) +
  theme_burden()
ggsave(out("StatusByGeneType.pdf"), plot = g_status, width = 9, height = 4)

g_status_overall <- qc %>%
  count(QC_Status, name = "n") %>%
  mutate(QC_Status = forcats::fct_relevel(QC_Status, "PASS", "WARN", "FAIL"))

g_status2 <- ggplot(g_status_overall, aes(x = QC_Status, y = n, fill = QC_Status)) +
  geom_col(color = "white") +
  labs(
    x = "QC status",
    y = "Number of genes",
    title = "Gene-level QC status"
  ) +
  theme_burden()
ggsave(out("StatusOverall.pdf"), plot = g_status2, width = 6, height = 4)

missing_long <- qc %>%
  mutate(
    gene_type = replace_na(as.character(gene_type), "unknown"),
    missingness = fraction_samples_with_missing_genotype,
    context_missingness = fraction_samples_without_context
  ) %>%
  select(gene_type, missingness, context_missingness) %>%
  pivot_longer(cols = c(missingness, context_missingness), names_to = "metric", values_to = "value")

g_missing <- ggplot(missing_long, aes(x = value, fill = metric)) +
  geom_histogram(alpha = 0.65, bins = 60, position = "identity") +
  facet_wrap(~metric, scales = "free_y") +
  geom_vline(xintercept = opt$MissingnessWarnThreshold, linetype = "dashed", color = "orange") +
  geom_vline(xintercept = opt$MissingnessFailThreshold, linetype = "dashed", color = "red") +
  geom_vline(xintercept = opt$ContextMissingnessFailThreshold, linetype = "dotted", color = "red") +
  labs(
    x = "Fraction",
    y = "Number of genes",
    fill = "Metric",
    title = "Missingness and context coverage QC"
  ) +
  theme_burden()
ggsave(out("Missingness.pdf"), plot = g_missing, width = 10, height = 5)

g_ratio <- ggplot(qc, aes(x = median_variance_ratio)) +
  geom_histogram(bins = 60, fill = "#4c72b0", color = "white") +
  geom_vline(xintercept = c(opt$VarianceRatioWarnLower, opt$VarianceRatioLower), linetype = "dashed") +
  geom_vline(xintercept = c(opt$VarianceRatioWarnUpper, opt$VarianceRatioUpper), linetype = "dashed") +
  labs(
    x = "Median empirical/population variance ratio",
    y = "Number of genes",
    title = "Variance calibration by gene"
  ) +
  theme_burden()
ggsave(out("VarianceRatio.pdf"), plot = g_ratio, width = 8, height = 4)

g_tail <- ggplot(qc, aes(x = max_abs_centered_z_population, fill = QC_Status)) +
  geom_histogram(bins = 80, color = "white") +
  geom_vline(xintercept = c(opt$TailZWarnThreshold, opt$TailZFailThreshold), linetype = "dashed") +
  labs(
    x = "max |CenteredEffectZPopulation|",
    y = "Number of genes",
    title = "Burden tail diagnostics"
  ) +
  theme_burden()
ggsave(out("TailZ.pdf"), plot = g_tail, width = 8, height = 4)

g_dom <- ggplot(qc, aes(x = dominant_variant_fraction_effect)) +
  geom_histogram(bins = 80, fill = "#55a868", color = "white") +
  geom_vline(xintercept = opt$DominantWarnThreshold, linetype = "dashed", color = "orange") +
  labs(
    x = "Dominant variant fraction (abs burden)",
    y = "Number of genes",
    title = "Dominant variant fraction"
  ) +
  theme_burden()
ggsave(out("DominantVariantFraction.pdf"), plot = g_dom, width = 8, height = 4)

fail_counts <- qc %>%
  filter(!is.na(QC_FailReasons), QC_FailReasons != "") %>%
  mutate(QC_FailReasons = str_split(QC_FailReasons, ";")) %>%
  unnest(QC_FailReasons) %>%
  filter(QC_FailReasons != "") %>%
  count(QC_FailReasons, name = "count") %>%
  mutate(reason_type = "Fail", reason = QC_FailReasons) %>%
  select(reason_type, reason, count)

warn_counts <- qc %>%
  filter(!is.na(QC_WarnReasons), QC_WarnReasons != "") %>%
  mutate(QC_WarnReasons = str_split(QC_WarnReasons, ";")) %>%
  unnest(QC_WarnReasons) %>%
  filter(QC_WarnReasons != "") %>%
  count(QC_WarnReasons, name = "count") %>%
  mutate(reason_type = "Warn", reason = QC_WarnReasons) %>%
  select(reason_type, reason, count)

reason_counts <- bind_rows(fail_counts, warn_counts) %>%
  arrange(reason_type, desc(count))

if (nrow(reason_counts) > 0) {
  g_reason <- ggplot(reason_counts, aes(x = reorder(reason, count), y = count, fill = reason_type)) +
    geom_col(position = "dodge") +
    coord_flip() +
    labs(
      x = "Reason",
      y = "Number of genes",
      fill = "Type",
      title = "Top QC fail/warn reasons"
    ) +
    theme_burden()
  ggsave(out("ReasonCounts.pdf"), plot = g_reason, width = 9, height = 5)
} else {
  pdf(out("ReasonCounts.pdf"), width = 8, height = 4)
  plot.new()
  text(0.5, 0.5, "No fail/warn reasons present", cex = 1.4)
  dev.off()
}

if (!is.null(opt$MedianGenesPerBinByGeneSet)) {
  message("Generating gene-set removal impact plot: ", opt$MedianGenesPerBinByGeneSet)

  median_genes_by_set <- fread(opt$MedianGenesPerBinByGeneSet) %>%
    as_tibble()

  required_median_cols <- c(
    "PercentChangeBin", "gene_type", "median_genes", "q25_genes", "q75_genes",
    "median_abs_z", "GeneCategory", "n_removed_genes"
  )
  missing_median_cols <- setdiff(required_median_cols, names(median_genes_by_set))
  if (length(missing_median_cols) > 0) {
    stop("Median genes-by-gene-set file is missing required columns: ", paste(missing_median_cols, collapse = ", "))
  }

  median_genes_by_set <- median_genes_by_set %>%
    mutate(
      PercentChangeBin = factor(PercentChangeBin, levels = unique(PercentChangeBin)),
      gene_type = replace_na(as.character(gene_type), "unknown"),
      GeneCategory = replace_na(as.character(GeneCategory), "Other"),
      GeneCategory = factor(GeneCategory, levels = c("All", "CausalCodingVariantGenesRemoved", "DominantVariantGenesRemoved"))
    )

  g_median_genes_by_set <- ggplot(median_genes_by_set, aes(x = PercentChangeBin, y = median_genes, color = GeneCategory, group = GeneCategory)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = q25_genes, ymax = q75_genes, fill = GeneCategory), alpha = 0.15, color = NA, inherit.aes = FALSE) +
    facet_wrap(~gene_type, scales = "free_y") +
    labs(
      x = "Percent-change bin",
      y = "Median # genes per individual (within bin)",
      color = "Removed gene set",
      fill = "Removed gene set",
      title = "Effect of removing gene sets on median gene burden by bin",
      subtitle = "Gene categories = how many genes removed and median burden within each bin"
    ) +
    theme_burden() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(out("GeneSetMedianImpact.pdf"), plot = g_median_genes_by_set, width = 12, height = 7)
} else {
  pdf(out("GeneSetMedianImpact.pdf"), width = 8, height = 4)
  plot.new()
  text(0.5, 0.5, "Median genes by gene set file not supplied", cex = 1.4)
  dev.off()
}

message("Saved QC diagnostic plots with prefix: ", out(""))
