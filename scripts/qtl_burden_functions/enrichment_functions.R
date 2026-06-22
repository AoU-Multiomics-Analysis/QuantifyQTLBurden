# Functions for outlier enrichment and permutation analyses.
extract_bin_lower_bound <- function(percent_change_bin) {
  as.numeric(stringr::str_match(
    as.character(percent_change_bin),
    "^[\\[\\(]?(-?\\d+\\.?\\d*),"
  )[, 2])
}

compute_bin_enrichment <- function(df, focal_lower_bound, ref_lower_bound = -10, outlier_col) {
  df2 <- df %>%
    dplyr::mutate(lower = extract_bin_lower_bound(PercentChangeBin)) %>%
    dplyr::filter(!is.na(lower)) %>%
    dplyr::filter(lower %in% c(focal_lower_bound, ref_lower_bound)) %>%
    dplyr::mutate(
      bin_group = dplyr::case_when(
        lower == focal_lower_bound ~ "focal",
        lower == ref_lower_bound ~ "reference"
      )
    ) %>%
    dplyr::group_by(bin_group, .data[[outlier_col]]) %>%
    dplyr::summarise(n = sum(n), .groups = "drop")

  a <- df2 %>% dplyr::filter(bin_group == "focal", .data[[outlier_col]] == TRUE) %>% dplyr::pull(n)
  b <- df2 %>% dplyr::filter(bin_group == "focal", .data[[outlier_col]] == FALSE) %>% dplyr::pull(n)
  c <- df2 %>% dplyr::filter(bin_group == "reference", .data[[outlier_col]] == TRUE) %>% dplyr::pull(n)
  d <- df2 %>% dplyr::filter(bin_group == "reference", .data[[outlier_col]] == FALSE) %>% dplyr::pull(n)

  if (length(a) == 0) a <- 0
  if (length(b) == 0) b <- 0
  if (length(c) == 0) c <- 0
  if (length(d) == 0) d <- 0

  a <- as.numeric(a)
  b <- as.numeric(b)
  c <- as.numeric(c)
  d <- as.numeric(d)

  a_cc <- a + 0.5
  b_cc <- b + 0.5
  c_cc <- c + 0.5
  d_cc <- d + 0.5

  log_odds_ratio <- log(a_cc) + log(d_cc) - log(b_cc) - log(c_cc)
  odds_ratio <- exp(log_odds_ratio)
  se_log_or <- sqrt(1 / a_cc + 1 / b_cc + 1 / c_cc + 1 / d_cc)

  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)

  tibble::tibble(
    focal_lower_bound = focal_lower_bound,
    reference_lower_bound = ref_lower_bound,
    focal_outliers = a,
    focal_non_outliers = b,
    ref_outliers = c,
    ref_non_outliers = d,
    focal_prop = ifelse(a + b > 0, a / (a + b), NA_real_),
    ref_prop = ifelse(c + d > 0, c / (c + d), NA_real_),
    odds_ratio = odds_ratio,
    log_odds_ratio = log_odds_ratio,
    se_log_or = se_log_or,
    ci_low = exp(log_odds_ratio - 1.96 * se_log_or),
    ci_high = exp(log_odds_ratio + 1.96 * se_log_or),
    p_value = fisher.test(tab, simulate.p.value = TRUE)$p.value
  )
}

run_permutation_enrichment <- function(benchmark_data, thresholds, n_perm, type_label, outlier_col) {
  permuted_log_or <- vector("list", n_perm)
  start_time <- Sys.time()
  report_every <- max(1L, floor(n_perm / 10))

  for (iter_idx in seq_len(n_perm)) {
    permuted_data <- benchmark_data %>%
      mutate(
        !!sym(outlier_col) := sample(.data[[outlier_col]])
      ) %>%
      group_by(PercentChangeBin) %>%
      dplyr::count(!!sym(outlier_col))

    permuted_log_or[[iter_idx]] <- purrr::map_dfr(
      thresholds,
      function(.x) {
        compute_bin_enrichment(
          permuted_data,
          focal_lower_bound = .x,
          ref_lower_bound = -10,
          outlier_col = outlier_col
        )
      }
    ) %>%
      mutate(
        permutation = iter_idx,
        type = type_label
      )

    if (iter_idx %% report_every == 0 || iter_idx == n_perm) {
      elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      remaining_seconds <- if (iter_idx == 0) {
        NA_real_
      } else {
        (elapsed_seconds / iter_idx) * (n_perm - iter_idx)
      }
      message(sprintf(
        "[%s] permutation %d/%d (%.1f%%) complete; elapsed=%0.1fs; eta=%0.1fs",
        type_label,
        iter_idx,
        n_perm,
        (iter_idx / n_perm) * 100,
        elapsed_seconds,
        ifelse(is.nan(remaining_seconds), NA_real_, remaining_seconds)
      ))
    }
  }

  bind_rows(permuted_log_or)
}

compute_enrichment_for_model <- function(df, model_label, thresholds) {
  down_outlier_count <- df %>%
    filter(gene_type == 'protein_coding') %>%
    group_by(PercentChangeBin) %>%
    dplyr::count(DownOutlier)

  up_outlier_count <- df %>%
    filter(gene_type == 'protein_coding') %>%
    group_by(PercentChangeBin) %>%
    dplyr::count(UpOutlier)

  results_down_ref_bin_comparison <- purrr::map_dfr(
    thresholds,
    ~ compute_bin_enrichment(
      down_outlier_count,
      focal_lower_bound = .x,
      ref_lower_bound = -10,
      outlier_col = "DownOutlier"
    )
  ) %>%
    mutate(type = "Down", enrichment_model = model_label)

  results_up_ref_bin_comparison <- purrr::map_dfr(
    thresholds,
    ~ compute_bin_enrichment(
      up_outlier_count,
      focal_lower_bound = .x,
      ref_lower_bound = -10,
      outlier_col = "UpOutlier"
    )
  ) %>%
    mutate(type = "Up", enrichment_model = model_label)

  bind_rows(results_down_ref_bin_comparison, results_up_ref_bin_comparison) %>%
    filter(focal_lower_bound != -10)
}
