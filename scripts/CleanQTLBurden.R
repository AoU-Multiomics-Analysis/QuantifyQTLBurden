library(tidyverse)
library(data.table)
library(optparse)
library(rtracklayer)
library(arrow)

safe_max <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  max(x, na.rm = TRUE)
}

safe_max_abs <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  max(abs(x), na.rm = TRUE)
}

safe_ratio <- function(numerator, denominator) {
  numerator <- as.numeric(numerator)
  denominator <- as.numeric(denominator)
  out <- rep(NA_real_, length(numerator))
  valid <- is.finite(numerator) & is.finite(denominator) & denominator != 0
  out[valid] <- numerator[valid] / denominator[valid]
  out
}

####### PARSE ARGUMENTS #########
option_list <- list(
    optparse::make_option(c("--QTLBurden"), type = "character", default = NULL,
                          help = "QTL burden summary file"),
    optparse::make_option(c("--AlleleFrequencies"), type = "character", default = NULL,
                          help = "Allele frequency file"),
    optparse::make_option(c("--ExpressionZscores"), type = "character", default = NULL,
                          help = "Expression z-score matrix"),
    optparse::make_option(c("--aFC"), type = "character", default = NULL,
                          help = "Allelic fold change file"),
    optparse::make_option(c("--AncestryAssignments"), type = "character", default = NULL,
                          help = "Ancestry assignments file"),
    optparse::make_option(c("--GTF"), type = "character", default = NULL,
                          help = "GTF file used to get annotate types of genes"),
    optparse::make_option(c("--eQTLSusie"), type = "character", default = NULL,
                          help = "Susie file containing variant annotations"),
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
                          help = "Warn when max |CenteredEffectZPopulation| is above this value"),
    optparse::make_option(c("--DominantVariantWarnThreshold"), type = "double", default = 0.98,
                          help = "Warn when a single variant explains more than this fraction of abs burden"),
    optparse::make_option(c("--OutlierPermutationIterations"), type = "integer", default = 200,
                          help = "Number of permutations used for null enrichment benchmarking")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

PathaFC <- opt$aFC
PathAncestryAssignments <- opt$AncestryAssignments
PathExpressionZscores <- opt$ExpressionZscores
PathAlleleFrequencies <- opt$AlleleFrequencies
BurdenPath <- opt$QTLBurden
GTFPath <- opt$GTF
SusiePath <- opt$eQTLSusie
MissingnessFailThreshold <- opt$MissingnessFailThreshold
MissingnessWarnThreshold <- opt$MissingnessWarnThreshold
ContextMissingnessFailThreshold <- opt$ContextMissingnessFailThreshold
VarianceRatioLower <- opt$VarianceRatioLower
VarianceRatioUpper <- opt$VarianceRatioUpper
VarianceRatioWarnLower <- opt$VarianceRatioWarnLower
VarianceRatioWarnUpper <- opt$VarianceRatioWarnUpper
TailZFailThreshold <- opt$TailZFailThreshold
TailZWarnThreshold <- opt$TailZWarnThreshold
DominantVariantWarnThreshold <- opt$DominantVariantWarnThreshold
OutlierPermutationIterations <- opt$OutlierPermutationIterations

####### LOAD DATA ###############
message('Loading GTF')
gene_annotations <- rtracklayer::import(GTFPath)
GeneTypes <- gene_annotations %>% data.frame()%>% select(gene_id,gene_type,gene_name) %>% distinct() %>% mutate(gene_id = str_remove(gene_id,'\\..*'))

message('Loading eQTL Susie')
eQTLSusie <- fread(SusiePath)



message('Extracting genes with coding variants')
if ('consequence' %in% colnames(eQTLSusie)) {
CodingVariantGenes <- eQTLSusie %>% 
    filter(pip > 0.9) %>% 
    filter(str_detect(consequence,'frameshift') | str_detect(consequence,'stop'))   %>%  
    distinct(molecular_trait_id) %>% 
    pull(molecular_trait_id)
} else {
CodingVariantGenes <- eQTLSusie %>% 
    filter(pip > 0.9) %>% 
    filter(frameshift == TRUE | stop_gained == TRUE)   %>%  
    distinct(molecular_trait_id) %>% 
    pull(molecular_trait_id)
}


message("Loading ancestry assignments")
AncestryDf <- fread(PathAncestryAssignments) %>%
    select(research_id, ancestry_pred_other)

message("Loading allele frequencies")
AlleleFrequencyDf <- fread(PathAlleleFrequencies)

message("Loading aFC data")
aFC <- fread(PathaFC)

message("Loading expression z scores")
ExpressionZscores <- fread(PathExpressionZscores) %>%
    pivot_longer(
        cols = -sample_id,
        names_to = "pid",
        values_to = "ObservedZ"
    )

message("Loading burden data and merging")
QTLBurdenMerge <- fread(BurdenPath) %>%
    left_join(AncestryDf, by = c("sample" = "research_id")) %>%
    mutate(gene_id = str_remove(pid,'\\..*')) %>%
    left_join(GeneTypes,by = 'gene_id') %>%
    mutate(CausalCodingVariantPresent = pid %in% CodingVariantGenes) %>% 
    left_join(ExpressionZscores, by = c("pid", "sample" = "sample_id")) %>%
    mutate(
            UpOutlier = ObservedZ > 4,
            DownOutlier = ObservedZ < -4
        ) 



####### COMPUTE EXPECTED VALUES #########
message("Computing expected mean and variance per population")
PopulationMOCExpectedValues <- aFC %>%
    left_join(AlleleFrequencyDf, by = c("sid" = "ID")) %>%
    distinct() %>%
    pivot_longer(
        cols = -c(sid, pid, log2_aFC,sid_chr,sid_pos),
        names_to = "ancestry_pred_other",
        values_to = "af"
    ) %>%
    mutate(
        af = as.numeric(af),
        WeightedEffect = 2 * af * log2_aFC,
        variance = 2 * af * (1 - af) * log2_aFC^2
    ) %>%
    filter(!is.na(af), !is.na(log2_aFC)) %>%
    group_by(pid, ancestry_pred_other) %>%
    summarize(
        GeneVariance_Population = sum(variance),
        ExpectedShift_Population = sum(WeightedEffect),
        .groups = "drop"
    )

EmpiricalVariance <- QTLBurdenMerge %>%
    group_by(ancestry_pred_other, pid) %>%
    summarize(
        EmpiricalVariance_Population = var(predicted_effect, na.rm = TRUE),
        .groups = "drop"
    )

###### CALCULATE QTL BURDEN Z SCORES #########
message('Merging all data')
QTLBurdenZscores <- QTLBurdenMerge %>%
    left_join(PopulationMOCExpectedValues, by = c("pid", "ancestry_pred_other")) %>%
    left_join(EmpiricalVariance, by = c("pid", "ancestry_pred_other")) %>%
    mutate(
        CenteredEffectPopulation = predicted_effect - ExpectedShift_Population,
        CenteredEffectZPopulation = (predicted_effect - ExpectedShift_Population) / sqrt(GeneVariance_Population),
        CenteredEffectZEmpiricalPopulation = (predicted_effect - ExpectedShift_Population) / sqrt(EmpiricalVariance_Population)
    ) %>%
    mutate(PercentChangeCenteredEffectPopulation = (2^CenteredEffectPopulation -1) *100) %>% 
    mutate(PercentChangeCenteredEffectPopulation = pmin(PercentChangeCenteredEffectPopulation,200)) %>%
    mutate(
      PercentChangeBin = case_when(
        PercentChangeCenteredEffectPopulation >= -100 & PercentChangeCenteredEffectPopulation <= -75 ~ "[-100,-75]",
        PercentChangeCenteredEffectPopulation > -75  & PercentChangeCenteredEffectPopulation <= -50 ~ "(-75,-50]",
        PercentChangeCenteredEffectPopulation > -50  & PercentChangeCenteredEffectPopulation <= -25 ~ "(-50,-25]",
        PercentChangeCenteredEffectPopulation > -25  & PercentChangeCenteredEffectPopulation <= -10 ~ "(-25,-10]",
        PercentChangeCenteredEffectPopulation > -10  & PercentChangeCenteredEffectPopulation < 10 ~ "(-10,10)",
        PercentChangeCenteredEffectPopulation >= 10  & PercentChangeCenteredEffectPopulation < 25  ~ "[10,25)",
        PercentChangeCenteredEffectPopulation >= 25  & PercentChangeCenteredEffectPopulation < 50  ~ "[25,50)",
        PercentChangeCenteredEffectPopulation >= 50  & PercentChangeCenteredEffectPopulation < 75  ~ "[50,75)",
        PercentChangeCenteredEffectPopulation >= 75  & PercentChangeCenteredEffectPopulation < 100 ~ "[75,100)",
        PercentChangeCenteredEffectPopulation >= 100 & PercentChangeCenteredEffectPopulation < 125 ~ "[100,125)",
        PercentChangeCenteredEffectPopulation >= 125 & PercentChangeCenteredEffectPopulation < 150 ~ "[125,150)",
        PercentChangeCenteredEffectPopulation >= 150 & PercentChangeCenteredEffectPopulation < 175 ~ "[150,175)",
        PercentChangeCenteredEffectPopulation >= 175 & PercentChangeCenteredEffectPopulation <= 200 ~ "[175,200]",
        TRUE ~ NA_character_
      )
    )

QTLBurdenZscores %>% write_tsv('QTLBurdenSummary.cleaned.tsv.gz')

message('Creating gene QC table')
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

QTLGeneBurdenQC <- QTLBurdenZscores %>%
  group_by(pid, gene_id, gene_name, gene_type, CausalCodingVariantPresent) %>%
  summarise(
    n_samples = n_distinct(sample),
    n_samples_nonzero_burden = n_distinct(sample[!is.na(predicted_effect) & predicted_effect != 0]),
    fraction_samples_nonzero_burden = n_samples_nonzero_burden / n_samples,
    n_samples_with_missing_genotype = n_distinct(sample[has_missing_genotype %in% TRUE]),
    fraction_samples_with_missing_genotype = n_samples_with_missing_genotype / n_samples,
    n_samples_without_context = n_distinct(sample[
      is.na(CenteredEffectZPopulation) |
      is.na(CenteredEffectZEmpiricalPopulation) |
      is.na(GeneVariance_Population) |
      is.na(EmpiricalVariance_Population)
    ]),
    fraction_samples_without_context = n_samples_without_context / n_samples,
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
    max_abs_centered_z_empirical_population = safe_max_abs(CenteredEffectZEmpiricalPopulation),
    median_abs_centered_z_empirical_population = median(abs(CenteredEffectZEmpiricalPopulation), na.rm = TRUE),
    median_variance_ratio = median(safe_ratio(EmpiricalVariance_Population, GeneVariance_Population), na.rm = TRUE),
    max_variance_ratio = safe_max(safe_ratio(EmpiricalVariance_Population, GeneVariance_Population)),
    max_abs_observed_z = safe_max_abs(ObservedZ),
    median_abs_observed_z = median(abs(ObservedZ), na.rm = TRUE),
    n_up_expression_outliers = sum(UpOutlier, na.rm = TRUE),
    n_down_expression_outliers = sum(DownOutlier, na.rm = TRUE),
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
    qc_flag_tail_warn = !is.na(max_abs_centered_z_population) & (max_abs_centered_z_population > TailZWarnThreshold) & !qc_flag_tail_fail,
    qc_flag_dominant_variant_warn = dominant_variant_fraction_effect >= DominantVariantWarnThreshold,
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



######  GET  BURDEN COUNTS  #########
KnockoutCountPerGene <- QTLBurdenZscores %>% 
    mutate(PercentChange = (2^CenteredEffectPopulation -1) *100) %>% 
    mutate(PredictedKO = PercentChange < -50) %>% 
    group_by(pid) %>% 
    summarize(NumKO = sum(PredictedKO))

message('Creating gene count table')
GeneBurdenCounts <- QTLBurdenZscores %>% 
  group_by(pid) %>% 
  summarise(
    max_abs_burden = max(abs(PercentChangeCenteredEffectPopulation), na.rm = TRUE),
    mean_abs_burden = mean(abs(PercentChangeCenteredEffectPopulation), na.rm = TRUE),
    n_deletion_like = sum(PercentChangeCenteredEffectPopulation <= -50, na.rm = TRUE),
    n_duplication_like = sum(PercentChangeCenteredEffectPopulation >= 150, na.rm = TRUE),
    n_validated_deletion = sum(
      PercentChangeCenteredEffectPopulation <= -50 &
        ObservedZ <= -3,
      na.rm = TRUE
    ),
    n_validated_duplication = sum(
      PercentChangeCenteredEffectPopulation >= 50 &
        ObservedZ >= 3,
      na.rm = TRUE
    ),
    has_validated_deletion = n_validated_deletion > 0,
    has_validated_duplication = n_validated_duplication > 0,
    .groups = "drop"
  ) %>% 
  left_join(KnockoutCountPerGene, by = "pid") %>% 
  mutate(CausalCodingVariantPresent = pid %in% CodingVariantGenes)
GeneBurdenCounts %>% write_tsv('QTLGeneBurdenCounts.tsv.gz')

###### SUMMARIZE BURDEN GENES PER INDIVIDUAL ########
message('Summarizing number genes per percent bin ')
KOPassList <- GeneBurdenCounts %>%
    filter(NumKO < 100) %>%
    pull(pid)

QTLBurdenFiltered <- QTLBurdenZscores %>%
    filter(!is.na(PercentChangeBin), !is.na(CenteredEffectZPopulation)) %>%
    filter(pid %in% KOPassList) %>%
    mutate(
        individual_id = sample,
        gene_id = str_remove(pid, '\\..*')
    ) %>%
    select(-any_of(c("gene_type", "gene_name"))) %>%
    left_join(GeneTypes, by = 'gene_id')

QTLBurdenPerSampleGene <- QTLBurdenFiltered %>%
    select(
      individual_id,
      sample,
      pid,
      gene_id,
      gene_name,
      gene_type,
      CausalCodingVariantPresent,
      PercentChangeBin,
      PercentChangeCenteredEffectPopulation,
      ObservedZ,
      CenteredEffectPopulation,
      CenteredEffectZPopulation,
      CenteredEffectZEmpiricalPopulation,
      ExpectedShift_Population,
      GeneVariance_Population,
      EmpiricalVariance_Population,
      has_missing_genotype,
      n_missing_genotypes
    )

arrow::write_parquet(
  QTLBurdenPerSampleGene,
  "QTLBurdenPerSampleGene.parquet"
)

compute_median_genes_per_bin <- function(df, remove_pids = character(0), gene_category = "All") {
  filtered <- if (length(remove_pids) == 0) {
    df
  } else {
    df %>% filter(!pid %in% remove_pids)
  }

  filtered %>%
    group_by(individual_id, PercentChangeBin, gene_type) %>%
    summarise(
      n_genes = n(),
      median_abs_z_individual = median(abs(ObservedZ), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(PercentChangeBin, gene_type) %>%
    summarise(
      median_genes = median(n_genes, na.rm = TRUE),
      q25_genes = quantile(n_genes, 0.25, na.rm = TRUE),
      q75_genes = quantile(n_genes, 0.75, na.rm = TRUE),
      median_abs_z = median(median_abs_z_individual, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(PercentChangeBin = fct_inorder(PercentChangeBin)) %>%
    mutate(
      GeneCategory = gene_category,
      n_removed_genes = n_distinct(remove_pids)
    )
}

# and summarize the range
MedianGenesPerIndividual <- compute_median_genes_per_bin(
  df = QTLBurdenPerSampleGene,
  remove_pids = character(0),
  gene_category = "All"
)

MedianGenesPerIndividualCodingVariantGenesRemoved <- compute_median_genes_per_bin(
  df = QTLBurdenPerSampleGene,
  remove_pids = CodingVariantGenes,
  gene_category = "CausalCodingVariantGenesRemoved"
)

DominantVariantWarnGenes <- aFCGeneQC %>%
  filter(
    !is.na(dominant_variant_fraction_effect),
    dominant_variant_fraction_effect >= DominantVariantWarnThreshold
  ) %>%
  pull(pid)

MedianGenesPerIndividualDominantVariantGenesRemoved <- compute_median_genes_per_bin(
  df = QTLBurdenPerSampleGene,
  remove_pids = DominantVariantWarnGenes,
  gene_category = "DominantVariantGenesRemoved"
)

MedianGenesSummary <- bind_rows(MedianGenesPerIndividual,MedianGenesPerIndividualCodingVariantGenesRemoved)
MedianGenesSummary %>% write_tsv('QTLBurdenMedianGenesPerBin.tsv')

MedianGenesByGeneSetSummary <- bind_rows(
  MedianGenesPerIndividual,
  MedianGenesPerIndividualCodingVariantGenesRemoved,
  MedianGenesPerIndividualDominantVariantGenesRemoved
) %>%
  arrange(GeneCategory, gene_type, PercentChangeBin)

MedianGenesByGeneSetSummary %>% write_tsv('QTLBurdenMedianGenesPerBinByGeneSet.tsv')

####### COMPUTE OUTLIER ENRICHMENT PER PERCENT CHANGE BIN ##########
message('Computing outlier enrichment')

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
  }

  bind_rows(permuted_log_or)
}

DownOutlierCount <- QTLBurdenFiltered %>% 
    filter(gene_type == 'protein_coding') %>% 
    group_by(PercentChangeBin) %>% 
    dplyr::count(DownOutlier)
UpOutlierCount <- QTLBurdenFiltered %>% 
    filter(gene_type == 'protein_coding') %>% 
    group_by(PercentChangeBin) %>% 
    dplyr::count(UpOutlier)

thresholds <- DownOutlierCount %>%
  mutate(lower = extract_bin_lower_bound(PercentChangeBin)) %>%
  distinct(lower) %>%
  filter(!is.na(lower), lower != 0) %>%
  arrange(lower) %>%
  pull(lower)

results_down_ref_bin_comparison <- purrr::map_dfr(
  thresholds,
  ~ compute_bin_enrichment(
    DownOutlierCount,
    focal_lower_bound = .x,
    ref_lower_bound = -10,
    outlier_col = "DownOutlier"
  )
) %>%
  mutate(type = "Down")

results_up_ref_bin_comparison <- purrr::map_dfr(
  thresholds,
  ~ compute_bin_enrichment(
    UpOutlierCount,
    focal_lower_bound = .x,
    ref_lower_bound = -10,
    outlier_col = "UpOutlier"
  )
) %>%
  mutate(type = "Up")

OutlierEnrichmentsRefBin <- bind_rows(results_down_ref_bin_comparison, 
                                      results_up_ref_bin_comparison
                                     ) %>% 
                            filter(focal_lower_bound != -10)

OutlierEnrichmentsRefBin %>% write_tsv('QTLBurdenOutlierEnrichment.tsv')

if (OutlierPermutationIterations > 0) {
  message('Computing outlier enrichment permutation null benchmark')

  OutlierBenchmarkData <- QTLBurdenFiltered %>%
    filter(gene_type == 'protein_coding') %>%
    select(PercentChangeBin, DownOutlier, UpOutlier)

  OutlierPermutationNull <- bind_rows(
    run_permutation_enrichment(
      benchmark_data = OutlierBenchmarkData,
      thresholds = thresholds,
      n_perm = OutlierPermutationIterations,
      type_label = "Down",
      outlier_col = "DownOutlier"
    ),
    run_permutation_enrichment(
      benchmark_data = OutlierBenchmarkData,
      thresholds = thresholds,
      n_perm = OutlierPermutationIterations,
      type_label = "Up",
      outlier_col = "UpOutlier"
    )
  )

  NullStats <- OutlierPermutationNull %>%
    group_by(type, focal_lower_bound, reference_lower_bound) %>%
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
    group_by(type, focal_lower_bound, reference_lower_bound) %>%
    summarize(
      empirical_p_greater = (sum(log_odds_ratio >= first(OutlierEnrichmentsRefBin$log_odds_ratio[
        OutlierEnrichmentsRefBin$type == type[1] &
          OutlierEnrichmentsRefBin$focal_lower_bound == focal_lower_bound[1] &
          OutlierEnrichmentsRefBin$reference_lower_bound == reference_lower_bound[1]
      ]), na.rm = TRUE) + 1) / (dplyr::n() + 1),
      empirical_p_two_sided = (sum(abs(log_odds_ratio) >= abs(first(OutlierEnrichmentsRefBin$log_odds_ratio[
        OutlierEnrichmentsRefBin$type == type[1] &
          OutlierEnrichmentsRefBin$focal_lower_bound == focal_lower_bound[1] &
          OutlierEnrichmentsRefBin$reference_lower_bound == reference_lower_bound[1]
      ])), na.rm = TRUE) + 1) / (dplyr::n() + 1),
      .groups = "drop"
    )

  OutlierEnrichmentsRefBin %>%
    left_join(NullStats, by = c("type", "focal_lower_bound", "reference_lower_bound")) %>%
    left_join(EmpiricalP, by = c("type", "focal_lower_bound", "reference_lower_bound")) %>%
    write_tsv('QTLBurdenOutlierEnrichmentPermutation.tsv')
} else {
  OutlierEnrichmentsRefBin %>%
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
    write_tsv('QTLBurdenOutlierEnrichmentPermutation.tsv')
}
