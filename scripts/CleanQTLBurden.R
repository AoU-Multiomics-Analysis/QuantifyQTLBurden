library(tidyverse)
library(data.table)
library(optparse)
library(rtracklayer)
library(arrow)

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

####### PARSE ARGUMENTS #########
option_list <- list(
    optparse::make_option(c("--QTLBurden"), type = "character", default = NULL,
                          help = "QTL burden summary file"),
    optparse::make_option(c("--AlleleFrequencies"), type = "character", default = NULL,
                          help = "Allele frequency file"),
    optparse::make_option(c("--ExpressionZscores"), type = "character", default = NULL,
                          help = "Expression z-score matrix"),
    optparse::make_option(c("--ProteomicsZscores"), type = "character", default = NULL,
                          help = "Proteomics z-score matrix"),
    optparse::make_option(c("--aFC"), type = "character", default = NULL,
                          help = "Allelic fold change file"),
    optparse::make_option(c("--AncestryAssignments"), type = "character", default = NULL,
                          help = "Ancestry assignments file"),
    optparse::make_option(c("--GTF"), type = "character", default = NULL,
                          help = "GTF file used to get annotate types of genes"),
    optparse::make_option(c("--eQTLSusie"), type = "character", default = NULL,
                          help = "Susie file containing variant annotations"),
    optparse::make_option(c("--BurdenTailProbability"), type = "double", default = NA_real_,
                          help = "Tail-probability cutoff for keeping confident loss/gain calls when computing noisy-filtered median counts.")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

PathaFC <- opt$aFC
PathAncestryAssignments <- opt$AncestryAssignments
PathExpressionZscores <- opt$ExpressionZscores
PathProteomicsZscores <- opt$ProteomicsZscores
PathAlleleFrequencies <- opt$AlleleFrequencies
BurdenPath <- opt$QTLBurden
GTFPath <- opt$GTF
SusiePath <- opt$eQTLSusie
BurdenTailProbability <- if (is.finite(opt$BurdenTailProbability)) {
  opt$BurdenTailProbability
} else {
  0.9
}
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

HasExpressionZscores <- !is.null(PathExpressionZscores) && PathExpressionZscores != ""
if (HasExpressionZscores) {
  message("Loading expression z scores")
  ExpressionZscores <- fread(PathExpressionZscores) %>%
      pivot_longer(
          cols = -sample_id,
          names_to = "pid",
          values_to = "ObservedZ"
      ) %>%
      mutate(ObservedZ = as.numeric(ObservedZ))
} else {
  message("No expression z scores provided; skipping observed-expression outlier annotations.")
}

HasProteomicsZscores <- !is.null(PathProteomicsZscores) && PathProteomicsZscores != ""
if (HasProteomicsZscores) {
  message("Loading proteomics z scores")
  ProteomicsZscores <- fread(PathProteomicsZscores) %>%
      pivot_longer(
          cols = -sample_id,
          names_to = "pid",
          values_to = "ObservedProteomicsZ"
      ) %>%
      mutate(ObservedProteomicsZ = as.numeric(ObservedProteomicsZ))
} else {
  message("No proteomics z scores provided; skipping proteomics outlier annotations.")
}

message("Loading burden data and merging")
QTLBurdenMerge <- fread(BurdenPath) %>%
    left_join(AncestryDf, by = c("sample" = "research_id")) %>%
    mutate(gene_id = str_remove(pid,'\\..*')) %>%
    left_join(GeneTypes,by = 'gene_id') %>%
    mutate(CausalCodingVariantPresent = pid %in% CodingVariantGenes) %>%
    {
      merged_expr <- if (HasExpressionZscores) {
        left_join(., ExpressionZscores, by = c("pid", "sample" = "sample_id")) %>%
          mutate(
            UpOutlier = ObservedZ > 4,
            DownOutlier = ObservedZ < -4
          )
      } else {
        mutate(
          .,
          ObservedZ = as.numeric(NA_real_),
          UpOutlier = as.logical(NA),
          DownOutlier = as.logical(NA)
        )
      }

      if (HasProteomicsZscores) {
        merged_expr %>%
          left_join(ProteomicsZscores, by = c("pid", "sample" = "sample_id")) %>%
          mutate(
            UpProteomicsOutlier = ObservedProteomicsZ > 4,
            DownProteomicsOutlier = ObservedProteomicsZ < -4
          )
      } else {
        merged_expr %>%
          mutate(
            ObservedProteomicsZ = as.numeric(NA_real_),
            UpProteomicsOutlier = as.logical(NA),
            DownProteomicsOutlier = as.logical(NA)
          )
      }
    }

if (!"burden_probability_loss50" %in% names(QTLBurdenMerge)) {
  QTLBurdenMerge$burden_probability_loss50 <- NA_real_
}
if (!"burden_probability_gain50" %in% names(QTLBurdenMerge)) {
  QTLBurdenMerge$burden_probability_gain50 <- NA_real_
}



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

if (is.finite(BurdenTailProbability) && BurdenTailProbability > 0) {
  QTLBurdenZscores <- QTLBurdenZscores %>%
    mutate(
      deletion_call_confident = ifelse(
        PercentChangeCenteredEffectPopulation <= -50,
        is.na(burden_probability_loss50) | (burden_probability_loss50 >= BurdenTailProbability),
        TRUE
      ),
      duplication_call_confident = ifelse(
        PercentChangeCenteredEffectPopulation >= 50,
        is.na(burden_probability_gain50) | (burden_probability_gain50 >= BurdenTailProbability),
        TRUE
      ),
      is_noisy_extreme_call = (!deletion_call_confident & PercentChangeCenteredEffectPopulation <= -50) |
        (!duplication_call_confident & PercentChangeCenteredEffectPopulation >= 50),
      is_dosage_extreme_call = PercentChangeCenteredEffectPopulation <= -50 |
        PercentChangeCenteredEffectPopulation >= 50
    ) %>%
    select(-deletion_call_confident, -duplication_call_confident)
} else {
  QTLBurdenZscores <- QTLBurdenZscores %>% mutate(is_noisy_extreme_call = FALSE, is_dosage_extreme_call = PercentChangeCenteredEffectPopulation <= -50 | PercentChangeCenteredEffectPopulation >= 50)
}

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
    n_duplication_like = sum(PercentChangeCenteredEffectPopulation >= 50, na.rm = TRUE),
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

QTLBurdenFiltered_AllCalls <- QTLBurdenZscores %>%
    filter(!is.na(PercentChangeBin), !is.na(CenteredEffectZPopulation)) %>%
    mutate(
        individual_id = sample,
        gene_id = str_remove(pid, '\\..*')
    ) %>%
    select(-any_of(c("gene_type", "gene_name"))) %>%
    left_join(GeneTypes, by = 'gene_id')

QTLBurdenFiltered_HighKORemoved <- QTLBurdenFiltered_AllCalls %>%
    filter(pid %in% KOPassList)

QTLBurdenFiltered_HighConfidence <- QTLBurdenFiltered_AllCalls %>%
    filter(!(is_dosage_extreme_call & is_noisy_extreme_call))

QTLBurdenFiltered <- QTLBurdenFiltered_HighKORemoved

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
      n_missing_genotypes,
      burden_probability_loss50,
      burden_probability_gain50,
      is_noisy_extreme_call,
      is_dosage_extreme_call
    )

arrow::write_parquet(
  QTLBurdenPerSampleGene,
  "QTLBurdenPerSampleGene.parquet"
)
