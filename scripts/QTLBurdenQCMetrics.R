library(data.table)
library(tidyverse)
library(optparse)

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }
  if (!is.null(sys.frame(1)$ofile)) {
    return(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")))
  }
  return(getwd())
}

SCRIPT_DIR <- get_script_dir()
source(file.path(SCRIPT_DIR, "qtl_burden_functions", "math_utils.R"))

extract_cli_numeric_arg <- function(arg_name, args, default_value = NULL) {
  long_form <- paste0("--", arg_name)
  eq_form <- paste0("^", long_form, "=")
  idx <- which(args == long_form)
  if (length(idx) == 1) {
    if (idx < length(args) && !startsWith(args[idx + 1], "--")) {
      return(list(
        value = suppressWarnings(as.numeric(args[idx + 1])),
        args = args[-c(idx, idx + 1)]
      ))
    }
    stop(sprintf("Missing value for %s", long_form))
  }

  idx_eq <- which(grepl(eq_form, args))
  if (length(idx_eq) == 1) {
    split_val <- sub(sprintf("^%s=", fixed = TRUE), "", args[idx_eq])
    return(list(
      value = suppressWarnings(as.numeric(split_val)),
      args = args[-idx_eq]
    ))
  }

  list(value = default_value, args = args)
}

raw_cli_args <- commandArgs(trailingOnly = TRUE)
tail_z_warn_sample_proportion_arg <- extract_cli_numeric_arg(
  "TailZWarnSampleProportionThreshold",
  raw_cli_args,
  default_value = NULL
)
clean_cli_args <- tail_z_warn_sample_proportion_arg$args

option_list <- list(
  optparse::make_option(c("--CleanedQTLBurden"), type = "character", default = NULL,
                        help = "QC-ready cleaned burden file"),
  optparse::make_option(c("--aFC"), type = "character", default = NULL,
                        help = "aFC input used to derive variant-level QC descriptors"),
  optparse::make_option(c("--MissingnessFailThreshold"), type = "double", default = 0.10,
                        help = "Fail genes if > this fraction of samples per gene have missing genotypes"),
  optparse::make_option(c("--MissingnessWarnThreshold"), type = "double", default = 0.05,
                        help = "Warn when fraction of samples with missing genotypes is above this value"),
  optparse::make_option(c("--ContextMissingnessFailThreshold"), type = "double", default = 0.05,
                        help = "Fail when allele-frequency/ancestry context is missing above this fraction"),
  optparse::make_option(c("--VarianceRatioLower"), type = "double", default = 0.20,
                        help = "Fail if empirical/population variance ratio is below this value"),
  optparse::make_option(c("--VarianceRatioUpper"), type = "double", default = 5,
                        help = "Fail if empirical/population variance ratio is above this value"),
  optparse::make_option(c("--VarianceRatioWarnLower"), type = "double", default = 0.35,
                        help = "Warn if empirical/population variance ratio is below this value"),
  optparse::make_option(c("--VarianceRatioWarnUpper"), type = "double", default = 3,
                        help = "Warn if empirical/population variance ratio is above this value"),
  optparse::make_option(c("--TailZFailThreshold"), type = "double", default = 25,
                        help = "Fail genes with extreme max |CenteredEffectZPopulation| above this value"),
  optparse::make_option(c("--TailZWarnThreshold"), type = "double", default = 15,
                        help = "Warn threshold for |CenteredEffectZPopulation| samples included in proportion check"),
  optparse::make_option(c("--TailZWarnSampleProportionThreshold"), type = "double", default = 0.05,
                        help = "Warn when fraction of samples with |CenteredEffectZPopulation| above TailZWarnThreshold exceeds this value"),
  optparse::make_option(c("--DominantVariantWarnThreshold"), type = "double", default = 0.98,
                        help = "Warn when a single variant explains more than this fraction of abs burden")
)

parser <- optparse::OptionParser(option_list = option_list)
opt <- tryCatch(
  {
    optparse::parse_args(parser, args = clean_cli_args)
  },
  error = function(e) {
    message("Argument parsing failed after stripping compatibility options: ", conditionMessage(e))
    stop(e)
  }
)

if (is.null(opt$CleanedQTLBurden) || is.null(opt$aFC)) {
  stop("--CleanedQTLBurden and --aFC are required")
}

MissingnessFailThreshold <- opt$MissingnessFailThreshold
MissingnessWarnThreshold <- opt$MissingnessWarnThreshold
ContextMissingnessFailThreshold <- opt$ContextMissingnessFailThreshold
VarianceRatioLower <- opt$VarianceRatioLower
VarianceRatioUpper <- opt$VarianceRatioUpper
VarianceRatioWarnLower <- opt$VarianceRatioWarnLower
VarianceRatioWarnUpper <- opt$VarianceRatioWarnUpper
TailZFailThreshold <- opt$TailZFailThreshold
TailZWarnThreshold <- opt$TailZWarnThreshold
TailZWarnSampleProportionThreshold <- if (is.null(tail_z_warn_sample_proportion_arg$value)) {
  if (is.null(opt$TailZWarnSampleProportionThreshold)) {
    0.05
  } else {
    opt$TailZWarnSampleProportionThreshold
  }
} else {
  tail_z_warn_sample_proportion_arg$value
}
DominantVariantWarnThreshold <- opt$DominantVariantWarnThreshold

message('Loading cleaned burden')
QTLBurdenZscores <- fread(opt$CleanedQTLBurden)
if (!all(c("ObservedZ", "UpOutlier", "DownOutlier") %in% names(QTLBurdenZscores))) {
  message("Expression outlier annotations are missing; proceeding with placeholder NA columns for QC summaries.")
  QTLBurdenZscores <- QTLBurdenZscores %>%
    mutate(
      ObservedZ = if ("ObservedZ" %in% names(.)) ObservedZ else NA_real_,
      UpOutlier = if ("UpOutlier" %in% names(.)) UpOutlier else as.logical(NA),
      DownOutlier = if ("DownOutlier" %in% names(.)) DownOutlier else as.logical(NA)
    )
}
if (!all(c("ObservedProteomicsZ", "UpProteomicsOutlier", "DownProteomicsOutlier") %in% names(QTLBurdenZscores))) {
  message("Proteomics outlier annotations are missing; proceeding with placeholder NA columns for QC summaries.")
  QTLBurdenZscores <- QTLBurdenZscores %>%
    mutate(
      ObservedProteomicsZ = if ("ObservedProteomicsZ" %in% names(.)) ObservedProteomicsZ else NA_real_,
      UpProteomicsOutlier = if ("UpProteomicsOutlier" %in% names(.)) UpProteomicsOutlier else as.logical(NA),
      DownProteomicsOutlier = if ("DownProteomicsOutlier" %in% names(.)) DownProteomicsOutlier else as.logical(NA)
    )
}

message('Loading aFC data')
aFC <- fread(opt$aFC)

message('Summarizing variant-level QC metrics from aFC')
aFCGeneQC <- aFC %>%
  group_by(pid) %>%
  summarise(
    n_afc_variants = n_distinct(sid),
    n_positive_effect_variants = sum(log2_aFC > 0, na.rm = TRUE),
    n_negative_effect_variants = sum(log2_aFC < 0, na.rm = TRUE),
    n_zero_effect_variants = sum(log2_aFC == 0, na.rm = TRUE),
    max_abs_log2_aFC = safe_max_abs(log2_aFC),
    mean_abs_log2_aFC = mean(abs(log2_aFC), na.rm = TRUE),
    total_abs_log2_aFC = sum(abs(log2_aFC), na.rm = TRUE),
    dominant_variant_fraction_effect = if_else(
      total_abs_log2_aFC > 0,
      max_abs_log2_aFC / total_abs_log2_aFC,
      NA_real_
    ),
    .groups = "drop"
  )

message('Creating gene QC table')
QTLGeneBurdenQC <- QTLBurdenZscores %>%
  as_tibble() %>%
  group_by(pid, gene_id, gene_name, gene_type, CausalCodingVariantPresent) %>%
  summarise(
    n_samples = n_distinct(sample),
    n_samples_nonzero_burden = n_distinct(sample[!is.na(predicted_effect) & predicted_effect != 0]),
    fraction_samples_nonzero_burden = n_samples_nonzero_burden / n_samples,
    n_samples_with_noisy_burden_calls = n_distinct(sample[is_noisy_extreme_call %in% TRUE]),
    fraction_samples_with_noisy_burden_calls = n_samples_with_noisy_burden_calls / n_samples,
    n_samples_with_missing_genotype = n_distinct(sample[has_missing_genotype %in% TRUE]),
    fraction_samples_with_missing_genotype = n_samples_with_missing_genotype / n_samples,
    max_missing_genotypes = safe_max(n_missing_genotypes),
    mean_missing_genotypes = mean(n_missing_genotypes, na.rm = TRUE),
    max_abs_predicted_effect = safe_max_abs(predicted_effect),
    mean_abs_predicted_effect = mean(abs(predicted_effect), na.rm = TRUE),
    median_abs_predicted_effect = median(abs(predicted_effect), na.rm = TRUE),
    max_abs_centered_effect = safe_max_abs(CenteredEffectPopulation),
    mean_abs_centered_effect = mean(abs(CenteredEffectPopulation), na.rm = TRUE),
    median_abs_centered_effect = median(abs(CenteredEffectPopulation), na.rm = TRUE),
    max_abs_centered_z_population = safe_max_abs(CenteredEffectZPopulation),
    median_abs_centered_z_population = median(abs(CenteredEffectZPopulation), na.rm = TRUE),
    fraction_samples_with_large_centered_z_population = if_else(
      n_samples > 0,
      mean(abs(CenteredEffectZPopulation) > TailZWarnThreshold, na.rm = TRUE),
      NA_real_
    ),
    max_abs_centered_z_empirical_population = safe_max_abs(CenteredEffectZEmpiricalPopulation),
    median_abs_centered_z_empirical_population = median(abs(CenteredEffectZEmpiricalPopulation), na.rm = TRUE),
    fraction_samples_without_context = mean(
      is.na(CenteredEffectZPopulation) |
        is.na(CenteredEffectZEmpiricalPopulation) |
        is.na(GeneVariance_Population) |
        is.na(EmpiricalVariance_Population),
      na.rm = TRUE
    ),
    median_variance_ratio = median(safe_ratio(EmpiricalVariance_Population, GeneVariance_Population), na.rm = TRUE),
    max_variance_ratio = safe_max(safe_ratio(EmpiricalVariance_Population, GeneVariance_Population)),
    max_abs_observed_z = safe_max_abs(ObservedZ),
    median_abs_observed_z = median(abs(ObservedZ), na.rm = TRUE),
    n_up_expression_outliers = sum(UpOutlier, na.rm = TRUE),
    n_down_expression_outliers = sum(DownOutlier, na.rm = TRUE),
    max_abs_observed_proteomics_z = safe_max_abs(ObservedProteomicsZ),
    median_abs_observed_proteomics_z = median(abs(ObservedProteomicsZ), na.rm = TRUE),
    n_up_proteomics_outliers = sum(UpProteomicsOutlier, na.rm = TRUE),
    n_down_proteomics_outliers = sum(DownProteomicsOutlier, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(aFCGeneQC, by = "pid") %>%
  mutate(
    qc_flag_missingness_fail = fraction_samples_with_missing_genotype > MissingnessFailThreshold,
    qc_flag_missingness_warn = fraction_samples_with_missing_genotype > MissingnessWarnThreshold,
    qc_flag_context_fail = fraction_samples_without_context > ContextMissingnessFailThreshold,
    qc_flag_variance_fail = is.na(median_variance_ratio) |
      (median_variance_ratio < VarianceRatioLower) |
      (median_variance_ratio > VarianceRatioUpper),
    qc_flag_variance_warn = !is.na(median_variance_ratio) &
      ((median_variance_ratio >= VarianceRatioLower & median_variance_ratio < VarianceRatioWarnLower) |
       (median_variance_ratio > VarianceRatioWarnUpper & median_variance_ratio <= VarianceRatioUpper)),
    qc_flag_tail_fail = !is.na(max_abs_centered_z_population) & (max_abs_centered_z_population > TailZFailThreshold),
    qc_flag_tail_warn = !is.na(max_abs_centered_z_population) &
      (fraction_samples_with_large_centered_z_population > TailZWarnSampleProportionThreshold) &
      !qc_flag_tail_fail,
    qc_flag_dominant_variant_warn = dominant_variant_fraction_effect >= DominantVariantWarnThreshold &
      n_afc_variants < 3,
    qc_fail_count = qc_flag_missingness_fail + qc_flag_context_fail + qc_flag_variance_fail + qc_flag_tail_fail,
    qc_warn_count = qc_flag_missingness_warn + qc_flag_variance_warn + qc_flag_tail_warn + qc_flag_dominant_variant_warn,
    QC_Status = case_when(
      qc_fail_count > 0 ~ "FAIL",
      qc_warn_count > 0 ~ "WARN",
      TRUE ~ "PASS"
    )
  ) %>%
  rowwise() %>%
  mutate(
    QC_FailReasons = {
      reasons <- c(
        if (qc_flag_missingness_fail) "high_sample_missingness" else NULL,
        if (qc_flag_context_fail) "missing_allelefreq_or_ancestry_context" else NULL,
        if (qc_flag_variance_fail) "empirical_vs_population_variance_mismatch" else NULL,
        if (qc_flag_tail_fail) "extreme_centered_burden_tail" else NULL
      )
      if (length(reasons) == 0) NA_character_ else paste(reasons, collapse = ";")
    },
    QC_WarnReasons = {
      reasons <- c(
        if (qc_flag_missingness_warn) "elevated_sample_missingness" else NULL,
        if (qc_flag_variance_warn) "marginal_variance_mismatch" else NULL,
        if (qc_flag_tail_warn) "large_centered_burden_tail" else NULL,
        if (qc_flag_dominant_variant_warn) "single_variant_dominance" else NULL
      )
      if (length(reasons) == 0) NA_character_ else paste(reasons, collapse = ";")
    }
  ) %>%
  ungroup()

QTLGeneBurdenQC %>% write_tsv('QTLGeneBurdenQC.tsv.gz')

QTLGeneBurdenStatusList <- QTLGeneBurdenQC %>%
  select(
    pid,
    gene_id,
    gene_name,
    gene_type,
    CausalCodingVariantPresent,
    QC_Status,
    qc_fail_count,
    qc_warn_count,
    QC_FailReasons,
    QC_WarnReasons,
    fraction_samples_with_missing_genotype,
    fraction_samples_without_context,
    median_variance_ratio,
    max_abs_centered_z_population
  ) %>%
  mutate(
    is_pass = QC_Status == "PASS",
    is_warn = QC_Status == "WARN",
    is_fail = QC_Status == "FAIL"
  )

QTLGeneBurdenStatusList %>% write_tsv('QTLGeneBurdenStatusList.tsv.gz')
