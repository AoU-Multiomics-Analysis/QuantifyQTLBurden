# QTL Burden Workflows

The QTL burden workflows start from aFC weights and genotype dosages. They compute predicted genetically driven expression shifts for each individual and gene.

## Scientific Model

The primary burden statistic is:

```text
predicted_effect = sum(genotype_dosage_variant * log2_aFC_variant)
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

Inputs:

| Input | Type | Description |
| --- | --- | --- |
| `aFCWeights` | File | TSV or TSV.GZ containing fine-mapped QTL variants and allelic fold-change effects. |
| `VCF` | File | Genotype VCF, ideally subset to variants present in the aFC file. |
| `IndexVCF` | File | Index for `VCF`, typically `.tbi`. |
| `AlleleFrequencies` | File | Allele-frequency table for expected burden and variance. |
| `ExpressionZscores` | File | Matrix of observed expression z scores. |
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

The VCF should contain genotype fields for all samples to be scored. [`../scripts/QTLBurden.R`](../scripts/QTLBurden.R) parses GT values before the first `:` and supports:

```text
0/0, 0/1, 1/0, 1/1, 0|0, 0|1, 1|0, 1|1
```

Cleaning and annotation inputs should contain:

| File | Required columns or structure |
| --- | --- |
| `AlleleFrequencies` | `ID` plus one column per ancestry group. `ID` should match `sid`. |
| `ExpressionZscores` | `sample_id` plus one column per `pid`. |
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
| `n_contributing_variants` | Number of variants with nonzero contribution. |
| `dominant_variant_fraction_abs` | Fraction of absolute burden explained by the largest contribution. |
| `N_eff_abs` | Effective number of contributing variants based on absolute effects. |
| `has_missing_genotype` | Whether any genotypes were missing for that individual/gene. |

When `AnnotateBurden = true`, the workflow can also emit:

| Output | Description |
| --- | --- |
| `QTLBurdenSummary.cleaned.tsv.gz` | Annotated burden table. |
| `QTLGeneBurdenCounts.tsv.gz` | Gene-level summary of large predicted burden events. |
| `QTLGeneBurdenQC.tsv.gz` | Gene-level QC table. |
| `QTLBurdenMedianGenesPerBin.tsv` | Median genes per individual in each percent-change bin. |
| `QTLBurdenOutlierEnrichment.tsv` | Observed expression outlier enrichment by burden bin. |

## Scripts

[`../scripts/QTLBurden.R`](../scripts/QTLBurden.R) computes raw QTL burden from aFC weights and genotype dosages.

[`../scripts/CleanQTLBurden.R`](../scripts/CleanQTLBurden.R) joins annotations, ancestry assignments, expression z scores, allele frequencies, and coding-variant summaries.
