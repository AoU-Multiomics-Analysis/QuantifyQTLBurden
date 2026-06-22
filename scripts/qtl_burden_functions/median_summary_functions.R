# Functions for median gene summaries by percent-change bin.
compute_median_genes_per_bin <- function(df, remove_pids = character(0), remove_noisy_dosage_extremes = TRUE, gene_category = "All") {
  filtered <- if (length(remove_pids) == 0) {
    df
  } else {
    df %>% filter(!pid %in% remove_pids)
  }

  if (remove_noisy_dosage_extremes) {
    filtered <- filtered %>%
      filter(!(is_dosage_extreme_call & is_noisy_extreme_call))
  }

  filtered %>%
    mutate(
      burden_tail_weight = case_when(
        PercentChangeCenteredEffectPopulation <= -50 ~ coalesce(burden_probability_loss50, 1),
        PercentChangeCenteredEffectPopulation >= 50 ~ coalesce(burden_probability_gain50, 1),
        TRUE ~ 1
      )
    ) %>%
    group_by(individual_id, PercentChangeBin, gene_type) %>%
    summarise(
      n_genes = n(),
      weighted_n_genes = sum(burden_tail_weight, na.rm = TRUE),
      median_abs_z_individual = median(abs(ObservedZ), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(PercentChangeBin, gene_type) %>%
    summarise(
      median_genes = median(n_genes, na.rm = TRUE),
      q25_genes = quantile(n_genes, 0.25, na.rm = TRUE),
      q75_genes = quantile(n_genes, 0.75, na.rm = TRUE),
      median_weighted_genes = median(weighted_n_genes, na.rm = TRUE),
      q25_weighted_genes = quantile(weighted_n_genes, 0.25, na.rm = TRUE),
      q75_weighted_genes = quantile(weighted_n_genes, 0.75, na.rm = TRUE),
      median_abs_z = median(median_abs_z_individual, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(PercentChangeBin = fct_inorder(PercentChangeBin)) %>%
    mutate(
      GeneCategory = gene_category,
      n_removed_genes = n_distinct(remove_pids)
    )
}

write_median_gene_summaries <- function(
  df,
  model_label,
  output_suffix = NULL,
  write_legacy_outputs = FALSE,
  coding_variant_genes = NULL,
  dominant_variant_warn_genes = NULL
) {
  coding_variant_genes <- if (is.null(coding_variant_genes)) {
    if (exists("CodingVariantGenes", inherits = TRUE)) {
      CodingVariantGenes
    } else {
      character(0)
    }
  } else {
    coding_variant_genes
  }

  dominant_variant_warn_genes <- if (is.null(dominant_variant_warn_genes)) {
    if (exists("DominantVariantWarnGenes", inherits = TRUE)) {
      DominantVariantWarnGenes
    } else {
      character(0)
    }
  } else {
    dominant_variant_warn_genes
  }

  has_suffix <- !is.null(output_suffix) && nchar(output_suffix) > 0
  base_prefix <- if (has_suffix) {
    paste0("QTLBurdenMedianGenesPerBin_", output_suffix)
  } else {
    "QTLBurdenMedianGenesPerBin"
  }
  by_gene_set_prefix <- if (has_suffix) {
    paste0("QTLBurdenMedianGenesPerBinByGeneSet_", output_suffix)
  } else {
    "QTLBurdenMedianGenesPerBinByGeneSet"
  }

  all_samples <- compute_median_genes_per_bin(
    df = df,
    remove_pids = character(0),
    remove_noisy_dosage_extremes = FALSE,
    gene_category = "All"
  )

  all_samples_noisy <- compute_median_genes_per_bin(
    df = df,
    remove_pids = character(0),
    remove_noisy_dosage_extremes = TRUE,
    gene_category = "All"
  )

  coding_removed <- compute_median_genes_per_bin(
    df = df,
    remove_pids = coding_variant_genes,
    remove_noisy_dosage_extremes = FALSE,
    gene_category = "CausalCodingVariantGenesRemoved"
  )

  coding_removed_noisy <- compute_median_genes_per_bin(
    df = df,
    remove_pids = coding_variant_genes,
    remove_noisy_dosage_extremes = TRUE,
    gene_category = "CausalCodingVariantGenesRemoved"
  )

  dominant_removed <- compute_median_genes_per_bin(
    df = df,
    remove_pids = dominant_variant_warn_genes,
    remove_noisy_dosage_extremes = FALSE,
    gene_category = "DominantVariantGenesRemoved"
  )

  dominant_removed_noisy <- compute_median_genes_per_bin(
    df = df,
    remove_pids = dominant_variant_warn_genes,
    remove_noisy_dosage_extremes = TRUE,
    gene_category = "DominantVariantGenesRemoved"
  )

  all_summary <- all_samples %>% mutate(enrichment_model = model_label)
  all_summary_noisy <- all_samples_noisy %>% mutate(enrichment_model = model_label)

  all_by_gene_set <- bind_rows(
    all_samples,
    coding_removed,
    dominant_removed
  ) %>%
    arrange(GeneCategory, gene_type, PercentChangeBin) %>%
    mutate(enrichment_model = model_label)

  all_by_gene_set_noisy <- bind_rows(
    all_samples_noisy,
    coding_removed_noisy,
    dominant_removed_noisy
  ) %>%
    arrange(GeneCategory, gene_type, PercentChangeBin) %>%
    mutate(enrichment_model = model_label)

  all_summary %>% write_tsv(paste0(base_prefix, ".tsv"))
  all_summary_noisy %>% write_tsv(paste0(base_prefix, "_DosageNoisyFiltered.tsv"))
  all_by_gene_set %>% write_tsv(paste0(by_gene_set_prefix, ".tsv"))
  all_by_gene_set_noisy %>% write_tsv(paste0(by_gene_set_prefix, "_DosageNoisyFiltered.tsv"))

  if (write_legacy_outputs) {
    all_summary %>%
      select(-enrichment_model) %>%
      write_tsv('QTLBurdenMedianGenesPerBin.tsv')
    all_summary_noisy %>%
      select(-enrichment_model) %>%
      write_tsv('QTLBurdenMedianGenesPerBin_DosageNoisyFiltered.tsv')
    all_by_gene_set %>%
      select(-enrichment_model) %>%
      write_tsv('QTLBurdenMedianGenesPerBinByGeneSet.tsv')
    all_by_gene_set_noisy %>%
      select(-enrichment_model) %>%
      write_tsv('QTLBurdenMedianGenesPerBinByGeneSet_DosageNoisyFiltered.tsv')
  }
}
