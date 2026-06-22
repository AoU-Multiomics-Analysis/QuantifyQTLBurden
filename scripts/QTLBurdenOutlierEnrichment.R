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
source(file.path(SCRIPT_DIR, "qtl_burden_functions", "enrichment_functions.R"))

option_list <- list(
  optparse::make_option(c("--CleanedQTLBurden"), type = "character", default = NULL,
                        help = "Cleaned QTL burden file output from CleanQTLBurden.R"),
  optparse::make_option(c("--OutlierPermutationIterations"), type = "integer", default = 200,
                        help = "Number of permutations used for null enrichment benchmarking")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (is.null(opt$CleanedQTLBurden)) {
  stop("--CleanedQTLBurden is required")
}

OutlierPermutationIterations <- opt$OutlierPermutationIterations

message("Loading cleaned burden")
QTLBurdenFiltered_AllCalls <- fread(opt$CleanedQTLBurden) %>%
  filter(!is.na(PercentChangeBin), !is.na(CenteredEffectZPopulation))

message("Computing outlier enrichment")

thresholds <- QTLBurdenFiltered_AllCalls %>%
  filter(gene_type == "protein_coding") %>%
  group_by(PercentChangeBin) %>%
  dplyr::count(DownOutlier) %>%
  mutate(lower = extract_bin_lower_bound(PercentChangeBin)) %>%
  distinct(lower) %>%
  filter(!is.na(lower), lower != 0) %>%
  arrange(lower) %>%
  pull(lower)

OutlierEnrichmentsRefBin <- bind_rows(
  compute_enrichment_for_model(
    df = QTLBurdenFiltered_AllCalls,
    model_label = "AllCalls",
    thresholds = thresholds
  ),
  compute_enrichment_for_model(
    df = QTLBurdenFiltered_AllCalls %>% filter(!(is_dosage_extreme_call & is_noisy_extreme_call)),
    model_label = "HighConfidence",
    thresholds = thresholds
  ),
  compute_enrichment_for_model(
    df = QTLBurdenFiltered_AllCalls %>%
      group_by(pid) %>%
      filter(sum(PercentChangeCenteredEffectPopulation < -50, na.rm = TRUE) < 100) %>%
      ungroup(),
    model_label = "HighKORemoved",
    thresholds = thresholds
  )
) %>% 
  arrange(enrichment_model, type, focal_lower_bound)

OutlierEnrichmentsRefBin %>% write_tsv("QTLBurdenOutlierEnrichment.tsv")

if (OutlierPermutationIterations > 0) {
  message("Computing outlier enrichment permutation null benchmark for AllCalls only.")

  OutlierBenchmarkData_AllCalls <- QTLBurdenFiltered_AllCalls %>%
    filter(gene_type == "protein_coding") %>%
    select(PercentChangeBin, DownOutlier, UpOutlier)

  if (nrow(OutlierBenchmarkData_AllCalls) == 0) {
    stop("No protein_coding rows found for permutation benchmark.")
  }

  OutlierPermutationNull <- bind_rows(
    run_permutation_enrichment(
      benchmark_data = OutlierBenchmarkData_AllCalls,
      thresholds = thresholds,
      n_perm = OutlierPermutationIterations
    )
  ) %>%
    mutate(enrichment_model = "AllCalls")

  if (nrow(OutlierPermutationNull) == 0) {
    stop("Permutation engine returned zero rows.")
  }

  NullStats <- OutlierPermutationNull %>%
    group_by(type, focal_lower_bound, reference_lower_bound, enrichment_model) %>%
    summarize(
      permutation_count = dplyr::n(),
      median_log_odds_ratio = median(log_odds_ratio, na.rm = TRUE),
      q25_log_odds_ratio = quantile(log_odds_ratio, 0.25, na.rm = TRUE),
      q75_log_odds_ratio = quantile(log_odds_ratio, 0.75, na.rm = TRUE),
      se_null_log_odds_ratio = sd(log_odds_ratio, na.rm = TRUE),
      mean_log_odds_ratio = mean(log_odds_ratio, na.rm = TRUE),
      .groups = "drop"
    )

  EmpiricalP <- OutlierPermutationNull %>%
    group_by(type, focal_lower_bound, reference_lower_bound, enrichment_model) %>%
    summarize(
      empirical_p_greater = (sum(
        log_odds_ratio >= first(OutlierEnrichmentsRefBin$log_odds_ratio[
          OutlierEnrichmentsRefBin$enrichment_model == enrichment_model[1] &
            OutlierEnrichmentsRefBin$type == type[1] &
            OutlierEnrichmentsRefBin$focal_lower_bound == focal_lower_bound[1] &
            OutlierEnrichmentsRefBin$reference_lower_bound == reference_lower_bound[1]
        ]
      ), na.rm = TRUE) + 1) / (dplyr::n() + 1),
      empirical_p_two_sided = (sum(
        abs(log_odds_ratio) >= abs(first(OutlierEnrichmentsRefBin$log_odds_ratio[
          OutlierEnrichmentsRefBin$enrichment_model == enrichment_model[1] &
            OutlierEnrichmentsRefBin$type == type[1] &
            OutlierEnrichmentsRefBin$focal_lower_bound == focal_lower_bound[1] &
            OutlierEnrichmentsRefBin$reference_lower_bound == reference_lower_bound[1]
        ])),
        na.rm = TRUE
      ) + 1) / (dplyr::n() + 1),
      .groups = "drop"
    )

  OutlierEnrichmentsRefBin <- OutlierEnrichmentsRefBin %>%
    left_join(
      NullStats,
      by = c("type", "focal_lower_bound", "reference_lower_bound", "enrichment_model")
    ) %>%
    left_join(
      EmpiricalP,
      by = c("type", "focal_lower_bound", "reference_lower_bound", "enrichment_model")
    )

  OutlierEnrichmentPermutation <- OutlierEnrichmentsRefBin %>%
    select(
      enrichment_model,
      type,
      focal_lower_bound,
      reference_lower_bound,
      permutation_count = permutation_count,
      median_log_odds_ratio = median_log_odds_ratio,
      q25_log_odds_ratio = q25_log_odds_ratio,
      q75_log_odds_ratio = q75_log_odds_ratio,
      se_null_log_odds_ratio = se_null_log_odds_ratio,
      mean_log_odds_ratio = mean_log_odds_ratio,
      empirical_p_greater = empirical_p_greater,
      empirical_p_two_sided = empirical_p_two_sided
    )
} else {
  OutlierEnrichmentPermutation <- OutlierEnrichmentsRefBin %>%
    mutate(
      permutation_count = 0L,
      median_log_odds_ratio = NA_real_,
      q25_log_odds_ratio = NA_real_,
      q75_log_odds_ratio = NA_real_,
      se_null_log_odds_ratio = NA_real_,
      mean_log_odds_ratio = NA_real_,
      empirical_p_greater = NA_real_,
      empirical_p_two_sided = NA_real_
    ) %>%
    select(
      enrichment_model,
      type,
      focal_lower_bound,
      reference_lower_bound,
      permutation_count,
      median_log_odds_ratio,
      q25_log_odds_ratio,
      q75_log_odds_ratio,
      se_null_log_odds_ratio,
      mean_log_odds_ratio,
      empirical_p_greater,
      empirical_p_two_sided
    )
}

if (nrow(OutlierEnrichmentPermutation) != nrow(OutlierEnrichmentsRefBin)) {
  stop("Permutation output rows do not match enrichment rows; refusing to write incomplete output.")
}

outlier_enrichment_permutation_file <- "QTLBurdenOutlierEnrichmentPermutation.tsv"
OutlierEnrichmentPermutation %>%
  write_tsv(outlier_enrichment_permutation_file)

if (!file.exists(outlier_enrichment_permutation_file)) {
  stop("Failed to write permutation output file: ", outlier_enrichment_permutation_file)
}

if (file.info(outlier_enrichment_permutation_file)$size == 0) {
  stop("Permutation output file was written but is empty: ", outlier_enrichment_permutation_file)
}
