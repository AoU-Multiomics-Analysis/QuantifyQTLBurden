library(tidyverse)
library(data.table)
library(bedr)
library(optparse)

# build trigger
########### FUNCTIONS #############

resolve_optional_numeric_column <- function(df, requested_col, fallback_cols, default = NA_real_) {
  if (!is.null(requested_col) && requested_col != "" && requested_col != "auto") {
    if (!requested_col %in% names(df)) {
      return(rep(default, nrow(df)))
    }
    return(as.numeric(df[[requested_col]]))
  }

  for (col in fallback_cols) {
    if (col %in% names(df)) {
      return(as.numeric(df[[col]]))
    }
  }
  rep(default, nrow(df))
}

compute_burden_tail_probability <- function(beta, beta_se, dosage, threshold = log2(0.5), direction = "deletion") {
  beta <- as.numeric(beta)
  beta_se <- as.numeric(beta_se)
  dosage <- as.matrix(dosage)

  if (length(beta) != nrow(dosage)) {
    stop("length(beta) must equal nrow(dosage)")
  }
  if (length(beta_se) != length(beta)) {
    stop("length(beta_se) must equal length(beta)")
  }

  beta_filled <- beta
  beta_filled[is.na(beta_filled)] <- 0
  beta_se[is.na(beta_se)] <- 0

  direction <- match.arg(direction, c("deletion", "duplication"))

  # L = sum(dosage_i * beta_i), beta_i ~ Normal(mu_i, sigma_i)
  burden_mean <- drop(crossprod(beta_filled, dosage))
  burden_var <- drop(crossprod(beta_se^2, dosage^2))
  burden_sd <- sqrt(pmax(burden_var, 0))

  p_ko <- ifelse(
    burden_sd < 1e-12,
    as.numeric(if (direction == "deletion") burden_mean <= threshold else burden_mean >= threshold),
    if (direction == "deletion") {
      pnorm((threshold - burden_mean) / burden_sd)
    } else {
      1 - pnorm((threshold - burden_mean) / burden_sd)
    }
  )

  list(mean = burden_mean, sd = burden_sd, p = p_ko)
}

compute_load_metrics <- function(beta, G, effect_threshold = NULL, beta_se = NULL,
                                loss_threshold = log2(0.5), gain_threshold = log2(1.5),
                                burden_direction = c("deletion", "duplication", "both")) {
  G <- as.matrix(G)

  if (length(beta) != nrow(G)) {
    stop("length(beta) must equal nrow(G)")
  }
  if (is.null(beta_se)) {
    beta_se <- rep(0, length(beta))
  }

  beta <- as.numeric(beta)
  beta_se <- as.numeric(beta_se)

  if (length(beta_se) != length(beta)) {
    stop("length(beta_se) must equal nrow(G)")
  }

  samples <- colnames(G)
  if (is.null(samples)) {
    samples <- paste0("sample_", seq_len(ncol(G)))
  }

  # sample-level missing genotype flag
  n_missing_genotypes <- colSums(is.na(G))
  has_missing_genotype <- n_missing_genotypes > 0

  # replace missing values with 0 for load calculation
  G_filled <- G
  G_filled[is.na(G_filled)] <- 0

  beta_filled <- beta
  beta_filled[is.na(beta_filled)] <- 0

  beta_up <- pmax(beta_filled, 0)
  beta_down <- abs(pmin(beta_filled, 0))

  # directional loads
  up_load <- drop(crossprod(beta_up, G_filled))
  down_load <- drop(crossprod(beta_down, G_filled))
  net_load <- up_load - down_load
  total_load <- up_load + down_load

  # variant presence
  G_present <- G_filled > 0

  up_idx <- beta_filled > 0
  down_idx <- beta_filled < 0

  n_up_variants <- if (any(up_idx)) {
    colSums(G_present[up_idx, , drop = FALSE])
  } else {
    rep(0L, ncol(G_filled))
  }

  n_down_variants <- if (any(down_idx)) {
    colSums(G_present[down_idx, , drop = FALSE])
  } else {
    rep(0L, ncol(G_filled))
  }

  n_up_variants_above_threshold <- rep(NA_integer_, ncol(G_filled))
  n_down_variants_above_threshold <- rep(NA_integer_, ncol(G_filled))

  if (!is.null(effect_threshold)) {
    up_strong <- beta_filled > effect_threshold
    down_strong <- beta_filled < -effect_threshold

    n_up_variants_above_threshold <- if (any(up_strong)) {
      colSums(G_present[up_strong, , drop = FALSE])
    } else {
      rep(0L, ncol(G_filled))
    }

    n_down_variants_above_threshold <- if (any(down_strong)) {
      colSums(G_present[down_strong, , drop = FALSE])
    } else {
      rep(0L, ncol(G_filled))
    }
  }

  burden_direction <- match.arg(burden_direction, c("deletion", "duplication", "both"))

  burden_deletion <- compute_burden_tail_probability(
    beta = beta_filled,
    beta_se = beta_se,
    dosage = G_filled,
    threshold = loss_threshold,
    direction = "deletion"
  )
  burden_overexpression <- compute_burden_tail_probability(
    beta = beta_filled,
    beta_se = beta_se,
    dosage = G_filled,
    threshold = gain_threshold,
    direction = "duplication"
  )

  burden_probability <- if (burden_direction == "duplication") {
    burden_overexpression$p
  } else {
    burden_deletion$p
  }
  burden_probability_deletion50 <- burden_deletion$p
  burden_probability_duplication50 <- burden_overexpression$p
  burden_mean <- if (burden_direction == "duplication") {
    burden_overexpression$mean
  } else {
    burden_deletion$mean
  }
  burden_sd <- if (burden_direction == "duplication") {
    burden_overexpression$sd
  } else {
    burden_deletion$sd
  }

  # signed contributions
  contrib <- beta_filled * G_filled
  abs_contrib <- abs(contrib)

  up_contrib <- ifelse(contrib > 0, contrib, 0)
  down_contrib <- ifelse(contrib < 0, abs(contrib), 0)

  predicted_effect <- colSums(contrib)
  total_abs_effect <- colSums(abs_contrib)
  n_contributing_variants <- colSums(abs_contrib > 0)

  # dominant variant fraction: absolute
  dominant_variant_fraction_abs <- apply(abs_contrib, 2, function(x) {
    s <- sum(x, na.rm = TRUE)
    if (s == 0) return(NA_real_)
    max(x, na.rm = TRUE) / s
  })

  # dominant variant fraction: up
  dominant_variant_fraction_up <- apply(up_contrib, 2, function(x) {
    s <- sum(x, na.rm = TRUE)
    if (s == 0) return(NA_real_)
    max(x, na.rm = TRUE) / s
  })

  # dominant variant fraction: down
  dominant_variant_fraction_down <- apply(down_contrib, 2, function(x) {
    s <- sum(x, na.rm = TRUE)
    if (s == 0) return(NA_real_)
    max(x, na.rm = TRUE) / s
  })

  # effective number of variants
  denom_abs <- colSums(abs_contrib^2, na.rm = TRUE)
  N_eff_abs <- ifelse(denom_abs == 0, NA_real_, (total_abs_effect^2) / denom_abs)

  up_sum <- colSums(up_contrib, na.rm = TRUE)
  denom_up <- colSums(up_contrib^2, na.rm = TRUE)
  N_eff_up <- ifelse(denom_up == 0, NA_real_, (up_sum^2) / denom_up)

  down_sum <- colSums(down_contrib, na.rm = TRUE)
  denom_down <- colSums(down_contrib^2, na.rm = TRUE)
  N_eff_down <- ifelse(denom_down == 0, NA_real_, (down_sum^2) / denom_down)

  tibble::tibble(
    sample = samples,

    up_load = up_load,
    down_load = down_load,
    total_load = total_load,
    net_load = net_load,

    predicted_effect = predicted_effect,
    total_abs_effect = total_abs_effect,

    n_up_variants = n_up_variants,
    n_down_variants = n_down_variants,
    n_contributing_variants = n_contributing_variants,

    n_up_variants_above_threshold = n_up_variants_above_threshold,
    n_down_variants_above_threshold = n_down_variants_above_threshold,

    dominant_variant_fraction_abs = dominant_variant_fraction_abs,
    dominant_variant_fraction_up = dominant_variant_fraction_up,
    dominant_variant_fraction_down = dominant_variant_fraction_down,

    N_eff_abs = N_eff_abs,
    N_eff_up = N_eff_up,
    N_eff_down = N_eff_down,
    burden_mean = burden_mean,
    burden_sd = burden_sd,
    burden_probability = burden_probability,
    burden_probability_loss50 = burden_probability_deletion50,
    burden_probability_gain50 = burden_probability_duplication50,

    n_missing_genotypes = n_missing_genotypes,
    has_missing_genotype = has_missing_genotype
  )
}

ProcessGenotypes <- function(VCFGenotypes) {
VCFCleaned <- VCFGenotypes %>% 
                    select(-ID,-QUAL,-FILTER,-INFO,-FORMAT) %>% 
                    mutate(variant = paste0(CHROM,':',POS,'_',REF,'_',ALT)) %>% 
                    dplyr::select(-CHROM,-POS,-REF,-ALT) %>% 
                    mutate(across(c(!variant),~str_remove(.,':.*'))) %>% 
                    mutate(across(c(!variant),~ case_when(. == '0/0' ~ 0 ,
                                                          . == '1/0' ~ 1,
                                                          . == '0/1' ~ 1,
                                                          . == '1/1' ~ 2,
                                                          . == '1|1' ~ 2,
                                                          . == '1|0' ~ 1,
                                                          . == '0|1' ~ 1,
                                                          . == '0|0' ~ 0
                                                         )
                                 )
                          )  %>% 
                    column_to_rownames('variant')
    
VCFCleaned
    
}
ComputeQTLBurden <- function(ChromList, PosList, BetaList, VariantList, VCF, SEList = NULL,
                             loss_threshold = log2(0.5), gain_threshold = log2(1.5),
                             burden_direction = "both") {
    require(bedr)
    
    VariantBeta <- data.frame(variant = VariantList,Beta=BetaList)
    Chrom <- ChromList %>% unique()
    MaxPos <- max(PosList) + 1
    MinPos <- min(PosList) - 1
    TabixInput <- paste0(Chrom,':',MinPos,'-',MaxPos)
    message('Extracting genotype data')
    GenotypeData <- tabix(TabixInput,VCF) %>% 
                mutate(POS = as.numeric(POS)) %>% 
                filter(POS %in% PosList) %>% 
                ProcessGenotypes
    GenotypeDataSorted <- GenotypeData[VariantBeta$variant,]
    QTLBurden <- compute_load_metrics(
      beta = VariantBeta$Beta,
      G = GenotypeDataSorted,
      beta_se = SEList,
      loss_threshold = loss_threshold,
      gain_threshold = gain_threshold,
      burden_direction = burden_direction
    )
    QTLBurden
}


####### PARSE ARGUMENTS #########
option_list <- list(
    optparse::make_option(c("--AllelicFoldChangeData"), type="character", default=NULL,
                        help="Summary file containing alleleic fold change per variant", metavar = "type"),
    optparse::make_option(c("--VariantEffectSEColumn"), type="character", default="auto",
                        help="Optional aFC standard error column (e.g. log2_aFC_se). Default auto-detects common names.",
                        metavar = "type"),
    optparse::make_option(c("--LossThreshold"), type="double", default=NA,
                        help="Posterior loss threshold on total log2 burden. Default: log2(0.5) (50% decrease).", metavar = "type"),
    optparse::make_option(c("--GainThreshold"), type="double", default=NA,
                        help="Posterior gain threshold on total log2 burden. Default: log2(1.5) (+50% increase).", metavar = "type"),
    optparse::make_option(c("--BurdenDirection"), type="character", default="both",
                        help="Direction for posterior scoring: loss, gain, or both. Default: both", metavar = "type"),
    optparse::make_option(c("--GenotypeData"), type="character", default=NULL,
                        help="VCF file to be queried (ideally this is a vcf that has just been subset to the aFC variants)", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="VCF file to be queried (ideally this is a vcf that has just been subset to the aFC variants)", metavar = "type")
    )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
aFCPath <- opt$AllelicFoldChangeData
PathVCF <- opt$GenotypeData
OutputName <- paste0(opt$OutputPrefix, ".QTLBurdenSummary.tsv.gz")
VariantEffectSEColumn <- opt$VariantEffectSEColumn
LossThreshold <- if (is.finite(opt$LossThreshold)) {
  opt$LossThreshold
} else {
  log2(0.5)
}
GainThreshold <- if (is.finite(opt$GainThreshold)) {
  opt$GainThreshold
} else {
  log2(1.5)
}
normalize_direction <- function(x) {
  if (x %in% c("loss", "deletion")) {
    "deletion"
  } else if (x %in% c("gain", "duplication")) {
    "duplication"
  } else {
    "both"
  }
}

BurdenDirection <- if (opt$BurdenDirection %in% c("deletion", "duplication", "loss", "gain", "both")) {
  normalize_direction(opt$BurdenDirection)
} else {
  "both"
}
######## RUN QTL BURDEN ANALYSIS ##########
aFC <- fread(aFCPath)

QTLBurden <- aFC %>% 
    group_by(pid) %>% 
    group_modify(~ {
      se_vec <- resolve_optional_numeric_column(
        .x,
        requested_col = VariantEffectSEColumn,
        fallback_cols = c("log2_aFC_se", "se_log2_aFC", "aFC_se", "beta_se", "se")
      )

      ComputeQTLBurden(
        .x$sid_chr,
        .x$sid_pos,
        .x$log2_aFC,
        .x$sid,
        PathVCF,
        SEList = se_vec,
        loss_threshold = LossThreshold,
        gain_threshold = GainThreshold,
        burden_direction = BurdenDirection
      )
    })
QTLBurden %>% write_tsv(OutputName)
