library(tidyverse)
library(data.table)
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
source(file.path(SCRIPT_DIR, "qtl_burden_functions", "median_summary_functions.R"))

option_list <- list(
  optparse::make_option(c("--CleanedQTLBurden"), type = "character", default = NULL,
                        help = "Cleaned QTL burden file from CleanQTLBurden.R"),
  optparse::make_option(c("--aFC"), type = "character", default = NULL,
                        help = "aFC input used to identify dominant-variant-gene outliers"),
  optparse::make_option(c("--DominantVariantWarnThreshold"), type = "double", default = 0.98,
                        help = "Warn when a single variant explains more than this fraction of abs burden")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (is.null(opt$CleanedQTLBurden) || is.null(opt$aFC)) {
  stop("--CleanedQTLBurden and --aFC are required")
}

DominantVariantWarnThreshold <- opt$DominantVariantWarnThreshold

message('Loading cleaned burden')
QTLBurdenZscores <- fread(opt$CleanedQTLBurden)

message('Summarizing burden call models for median-gene metrics')
QTLBurdenFiltered_AllCalls <- QTLBurdenZscores %>%
  filter(!is.na(PercentChangeBin), !is.na(CenteredEffectZPopulation)) %>%
  mutate(
    individual_id = sample,
    gene_id = str_remove(pid, '\\..*')
  )

message("Filtering KO-removal cohort")
GeneBurdenCounts <- QTLBurdenFiltered_AllCalls %>%
  mutate(PercentChange = (2^CenteredEffectPopulation - 1) * 100) %>%
  mutate(PredictedKO = PercentChange < -50) %>%
  group_by(pid) %>%
  summarize(NumKO = sum(PredictedKO))

KOPassList <- GeneBurdenCounts %>%
  filter(NumKO < 100) %>%
  pull(pid)

QTLBurdenFiltered_HighKORemoved <- QTLBurdenFiltered_AllCalls %>%
  filter(pid %in% KOPassList)

QTLBurdenFiltered_HighConfidence <- QTLBurdenFiltered_AllCalls %>%
  filter(!(is_dosage_extreme_call & is_noisy_extreme_call))

message("Loading aFC and computing dominant-variant-gene list")
aFC <- fread(opt$aFC)

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

DominantVariantWarnGenes <- aFCGeneQC %>%
  filter(
    !is.na(dominant_variant_fraction_effect),
    dominant_variant_fraction_effect >= DominantVariantWarnThreshold
  ) %>%
  pull(pid)

CodingVariantGenes <- QTLBurdenFiltered_AllCalls %>%
  filter(CausalCodingVariantPresent) %>%
  distinct(pid) %>%
  pull(pid)

write_median_gene_summaries(
  df = QTLBurdenFiltered_AllCalls,
  model_label = "AllCalls",
  output_suffix = "AllCalls",
  coding_variant_genes = CodingVariantGenes,
  dominant_variant_warn_genes = DominantVariantWarnGenes
)
write_median_gene_summaries(
  df = QTLBurdenFiltered_HighConfidence,
  model_label = "HighConfidence",
  output_suffix = "HighConfidence",
  coding_variant_genes = CodingVariantGenes,
  dominant_variant_warn_genes = DominantVariantWarnGenes
)
write_median_gene_summaries(
  df = QTLBurdenFiltered_HighKORemoved,
  model_label = "HighKORemoved",
  output_suffix = "HighKORemoved",
  coding_variant_genes = CodingVariantGenes,
  dominant_variant_warn_genes = DominantVariantWarnGenes,
  write_legacy_outputs = TRUE
)
