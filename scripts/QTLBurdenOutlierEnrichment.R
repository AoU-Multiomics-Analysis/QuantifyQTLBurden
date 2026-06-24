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

safe_first <- function(x) {
  if (length(x) == 0) {
    NA_real_
  } else {
    x[[1]]
  }
}

message("Loading cleaned burden")
QTLBurdenFiltered_AllCalls <- fread(opt$CleanedQTLBurden) %>%
  filter(!is.na(PercentChangeBin), !is.na(CenteredEffectZPopulation))

message("Checking available outlier annotations")
has_expression_outlier_columns <- all(c("UpOutlier", "DownOutlier") %in% names(QTLBurdenFiltered_AllCalls))
has_proteomics_outlier_columns <- all(c("UpProteomicsOutlier", "DownProteomicsOutlier") %in% names(QTLBurdenFiltered_AllCalls))

if (!has_expression_outlier_columns && !has_proteomics_outlier_columns) {
  message("No outlier annotation columns are present; skipping outlier enrichment and permutation outputs.")

  OutlierEnrichmentsRefBin <- tibble(
    OutlierSource = character(),
    enrichment_model = character(),
    type = character(),
    focal_lower_bound = numeric(),
    reference_lower_bound = numeric(),
    focal_outliers = integer(),
    focal_non_outliers = integer(),
    ref_outliers = integer(),
    ref_non_outliers = integer(),
    focal_prop = double(),
    ref_prop = double(),
    odds_ratio = double(),
    log_odds_ratio = double(),
    se_log_or = double(),
    ci_low = double(),
    ci_high = double(),
    p_value = double()
  )
  OutlierEnrichmentsRefBin %>%
    write_tsv("QTLBurdenOutlierEnrichment.tsv")

  OutlierEnrichmentNoisyFilterImpact <- tibble(
    OutlierSource = character(),
    type = character(),
    focal_lower_bound = numeric(),
    reference_lower_bound = numeric(),
    all_calls_log_odds_ratio = double(),
    high_confidence_log_odds_ratio = double(),
    all_calls_fisher_p = double(),
    high_confidence_fisher_p = double(),
    noisy_filter_abs_log_or_change = double(),
    noisy_filter_log_or_change = double(),
    noisy_filter_fisher_p_change = double(),
    noisy_filter_enrichment_improved = logical(),
    noisy_filter_fisher_improved = logical(),
    all_calls_focal_outliers = integer(),
    high_confidence_focal_outliers = integer(),
    all_calls_focal_outlier_rate = double(),
    high_confidence_focal_outlier_rate = double(),
    all_calls_ref_outlier_rate = double(),
    high_confidence_ref_outlier_rate = double(),
    all_calls_focal_non_outliers = integer(),
    high_confidence_focal_non_outliers = integer(),
    all_calls_ref_outliers = integer(),
    high_confidence_ref_outliers = integer(),
    all_calls_ref_non_outliers = integer(),
    high_confidence_ref_non_outliers = integer()
  )
  OutlierEnrichmentNoisyFilterImpact %>%
    write_tsv("QTLBurdenOutlierEnrichmentNoisyFilterImpact.tsv")

  OutlierEnrichmentPermutation <- tibble(
    OutlierSource = character(),
    enrichment_model = character(),
    type = character(),
    focal_lower_bound = numeric(),
    reference_lower_bound = numeric(),
    permutation_count = integer(),
    median_log_odds_ratio = double(),
    q25_log_odds_ratio = double(),
    q75_log_odds_ratio = double(),
    se_null_log_odds_ratio = double(),
    mean_log_odds_ratio = double(),
    empirical_p_greater = double(),
    empirical_p_two_sided = double()
  )
  OutlierEnrichmentPermutation %>%
    write_tsv("QTLBurdenOutlierEnrichmentPermutation.tsv")

  quit(save = "no", status = 0)
}

message("Computing outlier enrichment")

thresholds <- QTLBurdenFiltered_AllCalls %>%
  filter(gene_type == "protein_coding") %>%
  distinct(PercentChangeBin) %>%
  filter(!is.na(PercentChangeBin)) %>%
  mutate(lower = extract_bin_lower_bound(PercentChangeBin)) %>%
  distinct(lower) %>%
  filter(!is.na(lower), lower != 0) %>%
  arrange(lower) %>%
  pull(lower)

QTLBurdenFiltered_HighConfidence <- QTLBurdenFiltered_AllCalls %>%
  filter(!(is_dosage_extreme_call & is_noisy_extreme_call))

QTLBurdenFiltered_HighKORemoved <- QTLBurdenFiltered_AllCalls %>%
  group_by(pid) %>%
  filter(sum(PercentChangeCenteredEffectPopulation < -50, na.rm = TRUE) < 100) %>%
  ungroup()

derive_enrichment <- function(df, model_label, outlier_source, down_col, up_col) {
  df <- df %>%
    mutate(
      DownOutlier = .data[[down_col]],
      UpOutlier = .data[[up_col]]
    )

  compute_enrichment_for_model(
    df = df,
    model_label = model_label,
    thresholds = thresholds
  ) %>%
    mutate(OutlierSource = outlier_source)
}

expression_enrichments <- if (has_expression_outlier_columns) {
  bind_rows(
    derive_enrichment(
      df = QTLBurdenFiltered_AllCalls,
      model_label = "AllCalls",
      outlier_source = "Expression",
      down_col = "DownOutlier",
      up_col = "UpOutlier"
    ),
    derive_enrichment(
      df = QTLBurdenFiltered_HighConfidence,
      model_label = "HighConfidence",
      outlier_source = "Expression",
      down_col = "DownOutlier",
      up_col = "UpOutlier"
    ),
    derive_enrichment(
      df = QTLBurdenFiltered_HighKORemoved,
      model_label = "HighKORemoved",
      outlier_source = "Expression",
      down_col = "DownOutlier",
      up_col = "UpOutlier"
    )
  )
} else {
  tibble()
}

proteomics_enrichments <- if (has_proteomics_outlier_columns) {
  bind_rows(
    derive_enrichment(
      df = QTLBurdenFiltered_AllCalls,
      model_label = "AllCalls",
      outlier_source = "Proteomics",
      down_col = "DownProteomicsOutlier",
      up_col = "UpProteomicsOutlier"
    ),
    derive_enrichment(
      df = QTLBurdenFiltered_HighConfidence,
      model_label = "HighConfidence",
      outlier_source = "Proteomics",
      down_col = "DownProteomicsOutlier",
      up_col = "UpProteomicsOutlier"
    ),
    derive_enrichment(
      df = QTLBurdenFiltered_HighKORemoved,
      model_label = "HighKORemoved",
      outlier_source = "Proteomics",
      down_col = "DownProteomicsOutlier",
      up_col = "UpProteomicsOutlier"
    )
  )
} else {
  tibble()
}

OutlierEnrichmentsRefBin <- bind_rows(
  expression_enrichments,
  proteomics_enrichments
) %>%
  arrange(OutlierSource, enrichment_model, type, focal_lower_bound)

OutlierEnrichmentsRefBin %>% write_tsv("QTLBurdenOutlierEnrichment.tsv")

message("Comparing impact of removing noisy extreme calls from BurdenTailProbability-driven filtering")

OutlierEnrichmentNoisyFilterImpact <- OutlierEnrichmentsRefBin %>%
  filter(enrichment_model %in% c("AllCalls", "HighConfidence")) %>%
  select(
    OutlierSource,
    type,
    focal_lower_bound,
    reference_lower_bound,
    enrichment_model,
    log_odds_ratio,
    p_value,
    focal_prop,
    ref_prop,
    focal_outliers,
    focal_non_outliers,
    ref_outliers,
    ref_non_outliers
  ) %>%
  arrange(OutlierSource, type, focal_lower_bound, reference_lower_bound, enrichment_model) %>%
  pivot_wider(
    id_cols = c(OutlierSource, type, focal_lower_bound, reference_lower_bound),
    names_from = enrichment_model,
    values_from = c(log_odds_ratio, p_value, focal_prop, ref_prop, focal_outliers, focal_non_outliers, ref_outliers, ref_non_outliers),
    names_sep = "__"
  ) %>%
  mutate(
    noisy_filter_abs_log_or_change = abs(log_odds_ratio__HighConfidence) - abs(log_odds_ratio__AllCalls),
    noisy_filter_log_or_change = log_odds_ratio__HighConfidence - log_odds_ratio__AllCalls,
    noisy_filter_fisher_p_change = p_value__AllCalls - p_value__HighConfidence,
    noisy_filter_enrichment_improved = noisy_filter_abs_log_or_change > 0 & sign(log_odds_ratio__HighConfidence) == sign(log_odds_ratio__AllCalls),
    noisy_filter_stronger_only_if_same_direction = noisy_filter_abs_log_or_change > 0,
    noisy_filter_fisher_improved = p_value__HighConfidence < p_value__AllCalls
  ) %>%
  select(
    OutlierSource,
    type,
    focal_lower_bound,
    reference_lower_bound,
    all_calls_log_odds_ratio = log_odds_ratio__AllCalls,
    high_confidence_log_odds_ratio = log_odds_ratio__HighConfidence,
    all_calls_fisher_p = p_value__AllCalls,
    high_confidence_fisher_p = p_value__HighConfidence,
    noisy_filter_abs_log_or_change,
    noisy_filter_log_or_change,
    noisy_filter_fisher_p_change,
    noisy_filter_enrichment_improved,
    noisy_filter_fisher_improved,
    all_calls_focal_outliers = focal_outliers__AllCalls,
    high_confidence_focal_outliers = focal_outliers__HighConfidence,
    all_calls_focal_outlier_rate = focal_prop__AllCalls,
    high_confidence_focal_outlier_rate = focal_prop__HighConfidence,
    all_calls_ref_outlier_rate = ref_prop__AllCalls,
    high_confidence_ref_outlier_rate = ref_prop__HighConfidence,
    all_calls_focal_non_outliers = focal_non_outliers__AllCalls,
    high_confidence_focal_non_outliers = focal_non_outliers__HighConfidence,
    all_calls_ref_outliers = ref_outliers__AllCalls,
    high_confidence_ref_outliers = ref_outliers__HighConfidence,
    all_calls_ref_non_outliers = ref_non_outliers__AllCalls,
    high_confidence_ref_non_outliers = ref_non_outliers__HighConfidence
  ) %>%
  arrange(OutlierSource, type, focal_lower_bound, reference_lower_bound)

OutlierEnrichmentNoisyFilterImpact %>%
  write_tsv("QTLBurdenOutlierEnrichmentNoisyFilterImpact.tsv")

run_modality_permutation <- function(df, outlier_source, down_col, up_col) {
  OutlierBenchmarkData_AllCalls <- df %>%
    filter(gene_type == "protein_coding") %>%
    select(PercentChangeBin, all_of(down_col), all_of(up_col)) %>%
    mutate(
      OutlierSource = outlier_source,
      DownOutlier = .data[[down_col]],
      UpOutlier = .data[[up_col]]
    )

  if (nrow(OutlierBenchmarkData_AllCalls) == 0) {
    stop(sprintf("No protein_coding rows found for %s permutation benchmark.", outlier_source))
  }

  run_permutation_enrichment(
    benchmark_data = OutlierBenchmarkData_AllCalls,
    thresholds = thresholds,
    n_perm = OutlierPermutationIterations
  ) %>%
    mutate(
      OutlierSource = outlier_source,
      enrichment_model = "AllCalls"
    )
}

if (OutlierPermutationIterations > 0) {
  message("Computing outlier enrichment permutation null benchmark for AllCalls only.")

  OutlierPermutationNull <- bind_rows(
    if (has_expression_outlier_columns) {
      run_modality_permutation(QTLBurdenFiltered_AllCalls, "Expression", "DownOutlier", "UpOutlier")
    } else {
      tibble()
    },
    if (has_proteomics_outlier_columns) {
      run_modality_permutation(QTLBurdenFiltered_AllCalls, "Proteomics", "DownProteomicsOutlier", "UpProteomicsOutlier")
    } else {
      tibble()
    }
  )

  if (nrow(OutlierPermutationNull) == 0) {
    stop("Permutation engine returned zero rows.")
  }

  NullStats <- OutlierPermutationNull %>%
    group_by(type, focal_lower_bound, reference_lower_bound, enrichment_model, OutlierSource) %>%
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
    group_by(type, focal_lower_bound, reference_lower_bound, enrichment_model, OutlierSource) %>%
    summarize(
      observed_log_or = safe_first(OutlierEnrichmentsRefBin$log_odds_ratio[
        OutlierEnrichmentsRefBin$enrichment_model == enrichment_model[1] &
          OutlierEnrichmentsRefBin$type == type[1] &
          OutlierEnrichmentsRefBin$OutlierSource == OutlierSource[1] &
          OutlierEnrichmentsRefBin$focal_lower_bound == focal_lower_bound[1] &
          OutlierEnrichmentsRefBin$reference_lower_bound == reference_lower_bound[1]
      ]),
      empirical_p_greater = (sum(
        log_odds_ratio >= observed_log_or,
        na.rm = TRUE
      ) + 1) / (dplyr::n() + 1),
      empirical_p_two_sided = (sum(
        abs(log_odds_ratio) >= abs(observed_log_or),
        na.rm = TRUE
      ) + 1) / (dplyr::n() + 1),
      .groups = "drop"
    ) %>%
    select(-observed_log_or)

  OutlierEnrichmentsRefBin <- OutlierEnrichmentsRefBin %>%
    left_join(
      NullStats,
      by = c("type", "focal_lower_bound", "reference_lower_bound", "enrichment_model", "OutlierSource")
    ) %>%
    left_join(
      EmpiricalP,
      by = c("type", "focal_lower_bound", "reference_lower_bound", "enrichment_model", "OutlierSource")
    )

  OutlierEnrichmentPermutation <- OutlierEnrichmentsRefBin %>%
    select(
      OutlierSource,
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
      OutlierSource,
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
