# QTL Burden Workflows

The QTL burden workflows start from aFC weights and genotype dosages. They compute predicted genetically driven expression shifts for each individual and gene.

## Scientific Model

The primary burden statistic is:

```text
predicted_effect = sum(genotype_dosage_variant * log2_aFC_variant)
```
An optional analytic extreme-expression tail score is also reported:

```text
L = sum(genotype_dosage_variant * beta_variant)

For loss-like calls:
P(loss-like | data) = P(L < threshold_loss)
  - dosage-like perturbation state is fixed by the sample’s per-variant dosages

For over-expression/gain calls:
P(gain-like | data) = P(L > threshold_gain)
  - dosage-like perturbation state is fixed by the sample’s per-variant dosages

See [`burden-probability.md`](burden-probability.md) for full formulas, implementation details, and interpretation.
```

The cleaning step can also compute ancestry/population-centered values:

```text
ExpectedShift_Population = sum(2 * AF * log2_aFC)
GeneVariance_Population = sum(2 * AF * (1 - AF) * log2_aFC^2)
CenteredEffectPopulation = predicted_effect - ExpectedShift_Population
CenteredEffectZPopulation = CenteredEffectPopulation / sqrt(GeneVariance_Population)
PercentChangeCenteredEffectPopulation = (2^CenteredEffectPopulation - 1) * 100
```

## Main Workflow

[`../workflows/QuantifyQTLBurden.wdl`](../workflows/QuantifyQTLBurden.wdl) is the end-to-end burden workflow.

Stages:

1. `shard_afc_by_gene`: splits aFC weights into gene-based shards.
2. `QuantifyQTLBurden`: runs [`../scripts/QTLBurden.R`](../scripts/QTLBurden.R) on each shard.
3. `AggregateQTLBurden`: concatenates shard-level burden summaries.
4. `CleanBurdenData`: optionally runs [`../scripts/CleanQTLBurden.R`](../scripts/CleanQTLBurden.R) for annotation and summary outputs.

The Quantify branch (`QTLBurden.R`) itself does not require PIP and does not apply fine-mapping shrinkage; it is purely based on effect-size distributions and genotypes.

Inputs:

| Input | Type | Description |
| --- | --- | --- |
| `aFCWeights` | File | TSV or TSV.GZ containing fine-mapped QTL variants and allelic fold-change effects. |
| `VCF` | File | Genotype VCF, ideally subset to variants present in the aFC file. |
| `IndexVCF` | File | Index for `VCF`, typically `.tbi`. |
| `LossThreshold` | Float | Log2 cutoff for loss calls. Default: `-0.5849625007` (50% decrease). |
| `GainThreshold` | Float | Log2 cutoff for gain calls. Default: `0.5849625007` (+50%). |
| `BurdenTailProbability` | Float | Tail-probability cutoff used to flag low-confidence extreme calls (`is_noisy_extreme_call`) in outputs. Default: `0.9`. |
| `OutlierPermutationIterations` | Int | Number of permutations for outlier enrichment benchmarking (applies to the cleaning/annotation step). Default: `200`. |
| `VariantEffectSEColumn` | String | Optional column name for log2(aFC) standard error (or `auto`). |
| `AlleleFrequencies` | File | Allele-frequency table for expected burden and variance. |
| `ExpressionZscores` | File? | Matrix of observed expression z scores. If omitted, observed-expression outlier annotation and outlier enrichment outputs are disabled/sparse. |
| `ProteomicsZscores` | File? | Matrix of observed proteomics z scores. If omitted, proteomics outlier annotation and proteomics enrichment outputs are disabled/sparse. |
| `AncestryAssignments` | File | Table mapping individuals to ancestry groups. |
| `GTF` | File | Gene annotation file. |
| `eQTLSusie` | File | Fine-mapping/SuSiE annotation table. |
| `GenesPerShard` | Int | Number of genes per aFC shard. Default: `500`. |
| `AnnotateBurden` | Boolean | Whether to run cleaning and annotation. Default: `true`. |

## Cleaning-Only Workflow

[`../workflows/CleanQTLBurden.wdl`](../workflows/CleanQTLBurden.wdl) annotates a previously generated merged burden file.

Use it when `QTLBurdenSummary.AllGenes.tsv.gz` already exists and only the cleaning or annotation step needs to be rerun.

## Expected Input Columns

The aFC weights file should contain at least:

| Column | Description |
| --- | --- |
| `pid` | Gene or molecular trait ID. |
| `sid` | Variant ID matching the VCF variant format: `CHROM:POS_REF_ALT`. |
| `sid_chr` | Variant chromosome. |
| `sid_pos` | Variant position. |
| `log2_aFC` | Allelic fold-change effect size on the log2 scale. |
| `log2_aFC_se` | Optional standard error for log2_aFC (if available). |
| `log2_aFC_lower` | Optional lower 95% CI bound for log2_aFC (auto-converts to SE). |
| `log2_aFC_upper` | Optional upper 95% CI bound for log2_aFC (auto-converts to SE). |

The VCF should contain genotype fields for all samples to be scored. [`../scripts/QTLBurden.R`](../scripts/QTLBurden.R) parses GT values before the first `:` and supports:

```text
0/0, 0/1, 1/0, 1/1, 0|0, 0|1, 1|0, 1|1
```

Cleaning and annotation inputs should contain:

| File | Required columns or structure |
| --- | --- |
| `AlleleFrequencies` | `ID` plus one column per ancestry group. `ID` should match `sid`. |
| `ExpressionZscores` | `sample_id` plus one column per `pid` (optional). |
| `ProteomicsZscores` | `sample_id` plus one column per `<protein_id>_<gene_id[.version]>` (optional; may be a subset of measured proteins/genes). The gene portion is normalized by stripping version suffixes before join. |
| `AncestryAssignments` | `research_id`, `ancestry_pred_other`. |
| `GTF` | Gene-level metadata convertible to `gene_id`, `gene_type`, and `gene_name`. |
| `eQTLSusie` | `molecular_trait_id`, `pip`, and either `consequence` or boolean `frameshift` / `stop_gained` columns. |

## Outputs

The main aggregated output is `QTLBurdenSummary.AllGenes.tsv.gz`.

Important raw burden columns include:

| Column | Description |
| --- | --- |
| `pid` | Gene or molecular trait ID. |
| `sample` | Individual/sample ID from the VCF. |
| `predicted_effect` | Signed aggregate predicted expression effect. |
| `up_load` | Sum of positive-effect burden. |
| `down_load` | Sum of absolute negative-effect burden. |
| `net_load` | `up_load - down_load`. |
| `total_load` | `up_load + down_load`. |
| `total_abs_effect` | Sum of absolute variant contributions. |
| `is_noisy_extreme_call` | TRUE when a strong extreme call (<= -50% or >= 50% centered change) has low tail confidence below `BurdenTailProbability`. |
| `is_dosage_extreme_call` | TRUE when the centered effect is in an extreme dosage-driven bin (`<= -50%` or `>= 50%`). |
| `burden_mean` | Analytic mean total log2 burden for the active direction. |
| `burden_sd` | Analytic SD of the active directional burden. |
| `burden_probability` | Probabilistic extreme-expression score (one-sided tail), based on `LossThreshold` (loss-tail, consistent with fixed both-mode behavior). |
| `burden_probability_loss50` | Loss tail probability `P(L < LossThreshold)`. |
| `burden_probability_gain50` | Gain tail probability `P(L > GainThreshold)`. |
| `n_contributing_variants` | Number of variants with nonzero contribution. |
| `dominant_variant_fraction_abs` | Fraction of absolute burden explained by the largest contribution. |
| `N_eff_abs` | Effective number of contributing variants based on absolute effects. |
| `has_missing_genotype` | Whether any genotypes were missing for that individual/gene. |

When `AnnotateBurden = true`, the workflow can also emit:

Legacy unsuffixed median-bin outputs (`QTLBurdenMedianGenesPerBin.tsv`, `QTLBurdenMedianGenesPerBinByGeneSet.tsv`, and dosage-filtered variants) are no longer produced; use the suffixed per-model outputs below.

| Output | Description |
| --- | --- |
| `QTLBurdenSummary.cleaned.tsv.gz` | Annotated burden table. |
| `QTLGeneBurdenCounts.tsv.gz` | Gene-level summary of large predicted burden events. |
| `QTLGeneBurdenQC.tsv.gz` | Gene-level QC table. |
| `QTLGeneBurdenStatusList.tsv.gz` | Gene-level QC status list (`PASS`/`WARN`/`FAIL`) and reasons. |
| `QTLGeneBurdenQC_StatusByGeneType.pdf` | QC status distribution by gene type. |
| `QTLGeneBurdenQC_StatusOverall.pdf` | Overall QC status distribution across all genes. |
| `QTLGeneBurdenQC_Missingness.pdf` | Missingness/context coverage diagnostics. |
| `QTLGeneBurdenQC_VarianceRatio.pdf` | Variance-ratio QC diagnostics. |
| `QTLGeneBurdenQC_TailZ.pdf` | Tail behavior diagnostics for centered burden z scores. |
| `QTLGeneBurdenQC_DominantVariantFraction.pdf` | Dominant-variant fraction diagnostics. |
| `QTLGeneBurdenQC_ReasonCounts.pdf` | Frequency of fail/warn reasons by type. |
| `QTLBurdenMedianGenesPerBin_AllCalls.tsv` | Median genes per individual in each bin using all calls (unfiltered and noisy-filtered variants; all/noise-filtered gene-set versions also generated). |
| `QTLBurdenMedianGenesPerBin_HighConfidence.tsv` | Same as above, using high-confidence calls only (low-confidence dosage-extreme calls removed). |
| `QTLBurdenMedianGenesPerBin_HighKORemoved.tsv` | Same as above, excluding high-number KO genes (`pid` with `NumKO >= 100`). |
| `QTLBurdenMedianGenesPerBinByGeneSet_AllCalls.tsv` | Median genes per individual per bin by gene set (`All`, `CausalCodingVariantGenesRemoved`, `DominantVariantGenesRemoved`) for all calls. |
| `QTLBurdenMedianGenesPerBinByGeneSet_HighConfidence.tsv` | Same as above, using high-confidence calls. |
| `QTLBurdenMedianGenesPerBinByGeneSet_HighKORemoved.tsv` | Same as above, after high KO-count genes removed. |
| `QTLBurdenMedianGenesPerBin_AllCalls_DosageNoisyFiltered.tsv` | All-calls median by bin after removing dosage-extreme rows flagged noisy (`is_dosage_extreme_call & is_noisy_extreme_call`). |
| `QTLBurdenMedianGenesPerBin_HighConfidence_DosageNoisyFiltered.tsv` | High-confidence median by bin after noisy dosage-extreme filtering (already included in high-confidence model definition). |
| `QTLBurdenMedianGenesPerBin_HighKORemoved_DosageNoisyFiltered.tsv` | High-KO-removed median by bin after noisy dosage-extreme filtering. |
| `QTLBurdenMedianGenesPerBinByGeneSet_AllCalls_DosageNoisyFiltered.tsv` | By-gene-set noisy-filtered median table for all calls. |
| `QTLBurdenMedianGenesPerBinByGeneSet_HighConfidence_DosageNoisyFiltered.tsv` | By-gene-set noisy-filtered median table for high-confidence calls. |
| `QTLBurdenMedianGenesPerBinByGeneSet_HighKORemoved_DosageNoisyFiltered.tsv` | By-gene-set noisy-filtered median table after high-KO removal. |
| `QTLBurdenMedianGenesPerBin*` | Also includes `median_genes`, `q25_genes`, `q75_genes` (unweighted counts), and `median_weighted_genes`, `q25_weighted_genes`, `q75_weighted_genes` (weighted by `burden_tail_weight`: `P(loss-like)` for `<= -50`, `P(gain-like)` for `>= 50`, `1` otherwise). |
| `QTLBurdenPerSampleGene.parquet` | Per-sample, per-gene burden table used for downstream medians and enrichment analyses. |
| `QTLBurdenOutlierEnrichment.tsv` | Observed expression and proteomics outlier enrichment by burden bin for three models: `AllCalls`, `HighConfidence`, and `HighKORemoved`, with modality included in output rows. Empty (header-only) when no outlier annotations are provided. |
| `QTLBurdenOutlierEnrichmentNoisyFilterImpact.tsv` | Per-bin direct comparison of `AllCalls` vs `HighConfidence` (removing `is_noisy_extreme_call` within extreme bins), including change in log-odds ratio, Fisher p-value change, and contingency counts. Empty (header-only) when no outlier annotations are provided. |
| `QTLBurdenOutlierEnrichmentPermutation.tsv` | Permutation-based p-values and null odds-ratio summary for `AllCalls` only; other model rows are included with `NA` null fields. Empty (header-only) when no outlier annotations are provided. |

## Scripts

[`../scripts/QTLBurden.R`](../scripts/QTLBurden.R) computes raw QTL burden from aFC weights and genotype dosages.

[`../scripts/CleanQTLBurden.R`](../scripts/CleanQTLBurden.R) joins annotations, ancestry assignments, expression z scores, allele frequencies, and coding-variant summaries.
