library(tidyverse)
library(data.table)
library(optparse)
library(rtracklayer)

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
                          help = "Susie file containing variant annotations")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

PathaFC <- opt$aFC
PathAncestryAssignments <- opt$AncestryAssignments
PathExpressionZscores <- opt$ExpressionZscores
PathAlleleFrequencies <- opt$AlleleFrequencies
BurdenPath <- opt$QTLBurden
GTFPath <- opt$GTF
SusiePath <- opt$eQTLSusie

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
KOPassList <- GeneBurdenCounts %>% filter(NumKO < 100) %>% 
QTLBurdenFiltered <- QTLBurdenZscores %>% 
    filter(!is.na(PercentChangeBin), !is.na(z_score)) %>% 
    filter(pid %In% KOPassList) %>% 
    mutate(gene_id = str_remove(pid,'\\..*')) %>% 
    left_join(gene_types,by = 'gene_id')

# Count genes per individual within each percent-change bin
# and summarize the range
MedianGenesPerIndividual <- QTLBurdenFiltered %>%
    group_by(individual_id, PercentChangeBin,gene_type) %>%
    summarise(
        n_genes = n(),
        median_abs_z_individual = median(abs(ObservedZ), na.rm = TRUE),
        .groups = "drop"
    ) %>% 
    group_by(PercentChangeBin,gene_type) %>%
    summarise(
        median_genes = median(n_genes, na.rm = TRUE),
        q25_genes = quantile(n_genes, 0.25, na.rm = TRUE),
        q75_genes = quantile(n_genes, 0.75, na.rm = TRUE),
        median_abs_z = median(median_abs_z_individual, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(PercentChangeBin = fct_inorder(PercentChangeBin))  %>% 
    mutate(GeneCategory = 'All')

MedianGenesPerIndividualCodingVariantGenesRemoved <- QTLBurdenFiltered %>%
    filter(CausalCodingVariantPresent == FALSE) %>% 
    group_by(individual_id, PercentChangeBin,gene_type) %>%
    summarise(
        n_genes = n(),
        median_abs_z_individual = median(abs(ObservedZ), na.rm = TRUE),
        .groups = "drop"
    ) %>% 
    group_by(PercentChangeBin,gene_type) %>%
    summarise(
        median_genes = median(n_genes, na.rm = TRUE),
        q25_genes = quantile(n_genes, 0.25, na.rm = TRUE),
        q75_genes = quantile(n_genes, 0.75, na.rm = TRUE),
        median_abs_z = median(median_abs_z_individual, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(PercentChangeBin = fct_inorder(PercentChangeBin)) %>% 
    mutate(GeneCategory = 'CausalCodingVariantGenesRemoved')


MedianGenesSummary <- bind_rows(MedianGenesPerIndividual,MedianGenesPerIndividualCodingVariantGenesRemoved)
MedianGenesSummary %>% write_tsv('QTLBurdenMedianGenesPerBin.tsv')

####### COMPUTE OUTLIER ENRICHMENT PER PERCENT CHANGE BIN ##########
message('Computing outlier enrichment')
compute_bin_enrichment <- function(df, focal_lower_bound, ref_lower_bound = -10, outlier_col) {
  df2 <- df %>%
    dplyr::mutate(
      lower = as.numeric(
        stringr::str_match(
          as.character(PercentChangeBin),
          "^[\\[\\(]?(-?\\d+\\.?\\d*),"
        )[, 2]
      )
    ) %>%
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
  
  a <- df2 %>% dplyr::filter(bin_group == "focal",     .data[[outlier_col]] == TRUE)  %>% dplyr::pull(n)
  b <- df2 %>% dplyr::filter(bin_group == "focal",     .data[[outlier_col]] == FALSE) %>% dplyr::pull(n)
  c <- df2 %>% dplyr::filter(bin_group == "reference", .data[[outlier_col]] == TRUE)  %>% dplyr::pull(n)
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



DownOutlierCount <- QTLBurdenFiltered %>% 
    filter(gene_type == 'protein_coding') %>% 
    group_by(PercentChangeBin) %>% 
    dplyr::count(DownOutlier)
UpOutlierCount <- QTLBurdenFiltered %>% 
    filter(gene_type == 'protein_coding') %>% 
    group_by(PercentChangeBin) %>% 
    dplyr::count(UpOutlier)

thresholds <- DownOutlierCount %>%
  mutate(
    lower = as.numeric(
      stringr::str_extract(as.character(PercentChangeBin), "-?\\d+\\.?\\d*")
    )
  ) %>%
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

OutlierEnrichments %>% write_tsv('QTLBurdenOutlierEnrichment.tsv')
