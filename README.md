# QuantifyQTLBurden

QuantifyQTLBurden is a WDL/R pipeline for estimating individual-level QTL burden from fine-mapped QTL effect sizes and genotypes. The repository also includes independent WDL workflows for computing allelic fold change (aFC), because aFC estimation is the upstream step that produces the effect-size weights used by QTL burden quantification.

The workflows remain separate:

1. Run the aFC workflows when starting from QTL mapping results, genotypes, expression BED files, and covariates.
2. Run the QuantifyQTLBurden workflows when starting from aFC weights and genotypes.

The expected use case is a set of fine-mapped QTL variants, usually filtered to high-confidence causal variants such as `PIP > 0.9`, where each variant has a `log2_aFC` effect size. For each gene, the pipeline computes the aggregate predicted expression effect carried by each individual. Optional annotation steps compare those predicted effects to observed expression z scores, ancestry-specific allele frequencies, gene annotations, and coding-variant annotations.

## Pipeline Overview

This repository has two independent pipeline groups.

### aFC Pipeline

The upstream aFC workflows were merged from <https://github.com/AoU-Multiomics-Analysis/aFC> and are based on the aFC tool described at <https://github.com/secastel/aFC>.

The preprocessing workflow, [`workflows/preprocess_AFC_inputs.wdl`](workflows/preprocess_AFC_inputs.wdl), prepares genotype and expression inputs so they meet aFC requirements:

- `preprocess_expression_bed`: re-compresses the expression BED with `bgzip` and creates a tabix index.
- `annotate_vcf_ids`: annotates VCF variant IDs in `CHROM:POS_REF_ALT` format so they match the `sid` field used in QTL files.

The aFC workflow, [`workflows/aFC.wdl`](workflows/aFC.wdl), runs aFC by chromosome and then merges per-chromosome reports:

- `split_vcf_by_chr`: subsets and indexes the VCF for one chromosome.
- `run_afc`: runs `/opt/aFC/aFC.py` on the chromosome-specific VCF, expression BED, covariates, and QTL file.
- `merge_afc_reports`: concatenates chromosome-level aFC outputs into one gzipped report.

Recommended aFC execution order:

1. Prepare a QTL file with at least `pid`, `sid`, `sid_chr`, and `sid_pos`, where `sid` uses `CHROM:POS_REF_ALT`.
2. Run [`workflows/preprocess_AFC_inputs.wdl`](workflows/preprocess_AFC_inputs.wdl) unless the VCF and expression BED are already correctly formatted and indexed.
3. Run [`workflows/aFC.wdl`](workflows/aFC.wdl) to produce `<prefix>.aFC.txt.gz`.

### QTL Burden Pipeline

The main burden workflow is [`workflows/QuantifyQTLBurden.wdl`](workflows/QuantifyQTLBurden.wdl). It is designed to run on Terra/Cromwell and has four major stages:

1. **Shard aFC weights by gene**
   - Splits the aFC input table into smaller gene-based shards.
   - Keeps all variants for the same gene together.
   - Controlled by `GenesPerShard`, which defaults to `500`.

2. **Compute QTL burden per shard**
   - Runs [`scripts/QTLBurden.R`](scripts/QTLBurden.R) on each shard.
   - Extracts genotype dosages from the input VCF using genomic ranges from the aFC table.
   - Computes per-individual, per-gene predicted expression burden.

3. **Aggregate shard outputs**
   - Concatenates all shard-level burden summaries into one file:
     `QTLBurdenSummary.AllGenes.tsv.gz`.

4. **Optionally clean and annotate burden results**
   - Runs [`scripts/CleanQTLBurden.R`](scripts/CleanQTLBurden.R) when `AnnotateBurden = true`.
   - Adds gene annotations, ancestry assignments, expected population-level burden moments, observed expression z scores, burden bins, gene-level burden summaries, and outlier enrichment summaries.

There is also a standalone cleaning workflow, [`workflows/CleanQTLBurden.wdl`](workflows/CleanQTLBurden.wdl), for annotating a previously generated merged QTL burden file.

## Scientific Summary

For a gene with causal QTL variants, the primary predicted burden is:

```text
predicted_effect = sum(genotype_dosage_variant * log2_aFC_variant)
```

where genotype dosage is encoded as 0, 1, or 2 alternate alleles. Positive `log2_aFC` values increase predicted expression, and negative values decrease predicted expression.

The cleaning step also computes ancestry/population-centered values using allele frequencies:

```text
ExpectedShift_Population = sum(2 * AF * log2_aFC)
GeneVariance_Population = sum(2 * AF * (1 - AF) * log2_aFC^2)
CenteredEffectPopulation = predicted_effect - ExpectedShift_Population
CenteredEffectZPopulation = CenteredEffectPopulation / sqrt(GeneVariance_Population)
PercentChangeCenteredEffectPopulation = (2^CenteredEffectPopulation - 1) * 100
```

These centered values estimate whether an individual's predicted genetically driven expression shift is unusually high or low relative to their assigned ancestry group.

## Inputs

### aFC Workflow Inputs

[`workflows/preprocess_AFC_inputs.wdl`](workflows/preprocess_AFC_inputs.wdl):

| Input | Type | Description |
| --- | --- | --- |
| `vcf_file` | File | Input genotype VCF. |
| `expression_bed` | File | Expression BED file. |
| `prefix` | String | Output file prefix. |
| `memory` | Int | Memory in GB. Default: `16`. |
| `disk_space` | Int | Extra disk space in GB. Default: `50`. |
| `num_threads` | Int | Number of CPU threads. Default: `8`. |
| `num_preempt` | Int | Number of preemptible retries. Default: `0`. |

[`workflows/aFC.wdl`](workflows/aFC.wdl):

| Input | Type | Description |
| --- | --- | --- |
| `vcf_file` | File | Annotated, indexed genotype VCF. |
| `vcf_index` | File | Tabix index for `vcf_file`. |
| `expression_bed` | File | bgzip-compressed, tabix-indexed expression BED. |
| `expression_bed_index` | File | Tabix index for `expression_bed`. |
| `covariates_file` | File | Covariates file in the format expected by aFC. |
| `afc_qtl_file` | File | QTL file containing at least `pid`, `sid`, `sid_chr`, and `sid_pos`. |
| `prefix` | String | Output file prefix. |
| `chromosomes` | Array[String]? | Optional chromosome list. Defaults to `chr1`-`chr22`, `chrX`, and `chrY`. |
| `memory` | Int | Memory in GB. Default: `16`. |
| `disk_space` | Int | Extra disk space in GB. Default: `50`. |
| `num_threads` | Int | Number of CPU threads. Default: `8`. |
| `num_preempt` | Int | Number of preemptible retries. Default: `0`. |

### Main Workflow Inputs

| Input | Type | Description |
| --- | --- | --- |
| `aFCWeights` | File | TSV or TSV.GZ containing fine-mapped QTL variants and allelic fold-change effects. |
| `VCF` | File | Genotype VCF, ideally subset to variants present in the aFC file. Must be queryable by genomic region. |
| `IndexVCF` | File | Index for `VCF`, typically `.tbi`. This is localized with the VCF for tabix access. |
| `AlleleFrequencies` | File | Allele-frequency table used to compute ancestry/population expected burden and variance. |
| `ExpressionZscores` | File | Matrix of observed expression z scores with samples as rows and genes/transcripts as columns. |
| `AncestryAssignments` | File | Table mapping individuals to ancestry groups. |
| `GTF` | File | Gene annotation file used to attach `gene_type` and `gene_name`. |
| `eQTLSusie` | File | Fine-mapping/SuSiE annotation table used to identify genes with likely causal coding variants. |
| `GenesPerShard` | Int | Number of genes per aFC shard. Default: `500`. |
| `AnnotateBurden` | Boolean | Whether to run the cleaning/annotation step. Default: `true`. |

### Expected aFC Columns

The aFC input is used by both the sharding task and the R scripts. It should contain at least:

| Column | Description |
| --- | --- |
| `pid` | Gene or molecular trait ID. Used as the gene-level grouping column. |
| `sid` | Variant ID matching the VCF variant format used internally: `CHROM:POS_REF_ALT`. |
| `sid_chr` | Variant chromosome. |
| `sid_pos` | Variant position. |
| `log2_aFC` | Allelic fold-change effect size on the log2 scale. |

### Expected VCF Format

The VCF should contain genotype fields for all samples to be scored. `scripts/QTLBurden.R` parses GT values before the first `:` and currently supports:

```text
0/0, 0/1, 1/0, 1/1, 0|0, 0|1, 1|0, 1|1
```

These are converted to alternate allele dosages of 0, 1, or 2.

### Expected Cleaning/Annotation Columns

| File | Required columns or structure |
| --- | --- |
| `AlleleFrequencies` | `ID` plus one column per ancestry group. `ID` should match `sid` in the aFC file. |
| `ExpressionZscores` | `sample_id` plus one column per `pid`; values are observed expression z scores. |
| `AncestryAssignments` | `research_id`, `ancestry_pred_other`. |
| `GTF` | Must include gene-level metadata convertible to `gene_id`, `gene_type`, and `gene_name`. Version suffixes are removed from `gene_id`. |
| `eQTLSusie` | `molecular_trait_id`, `pip`, and either `consequence` or boolean `frameshift` / `stop_gained` columns. |

## Outputs

### aFC Outputs

[`workflows/preprocess_AFC_inputs.wdl`](workflows/preprocess_AFC_inputs.wdl) emits:

| Output | Description |
| --- | --- |
| `processed_expression_bed` | bgzip-compressed expression BED: `<prefix>.processed_expression.bed.gz`. |
| `processed_expression_bed_tbi` | Tabix index for the processed expression BED. |
| `annotated_vcf` | VCF with variant IDs annotated as `CHROM:POS_REF_ALT`: `<prefix>.annotated.vcf.gz`. |
| `annotated_vcf_tbi` | Tabix index for the annotated VCF. |

[`workflows/aFC.wdl`](workflows/aFC.wdl) emits:

| Output | Description |
| --- | --- |
| `per_chr_afc_reports` | Per-chromosome aFC reports: `<prefix>.<chr>.aFC.txt.gz`. |
| `final_afc_report` | Merged aFC report across chromosomes: `<prefix>.aFC.txt.gz`. |

### Main Burden Output

`QTLBurdenSummary.AllGenes.tsv.gz`

This is the aggregated per-individual, per-gene burden output from all shards. Important columns include:

| Column | Description |
| --- | --- |
| `pid` | Gene or molecular trait ID. |
| `sample` | Individual/sample ID from the VCF. |
| `predicted_effect` | Signed aggregate predicted expression effect: `sum(dosage * log2_aFC)`. |
| `up_load` | Sum of positive-effect burden. |
| `down_load` | Sum of absolute negative-effect burden. |
| `net_load` | `up_load - down_load`. |
| `total_load` | `up_load + down_load`. |
| `total_abs_effect` | Sum of absolute variant contributions. |
| `n_up_variants` | Number of carried variants with positive effects. |
| `n_down_variants` | Number of carried variants with negative effects. |
| `n_contributing_variants` | Number of variants with nonzero contribution. |
| `dominant_variant_fraction_abs` | Fraction of total absolute burden explained by the largest variant contribution. |
| `N_eff_abs` | Effective number of contributing variants based on absolute effects. |
| `n_missing_genotypes` | Number of missing genotypes for that individual/gene. |
| `has_missing_genotype` | Whether any genotypes were missing for that individual/gene. |

### Annotated Outputs

When `AnnotateBurden = true`, the workflow also emits:

| Output | Description |
| --- | --- |
| `QTLBurdenSummary.cleaned.tsv.gz` | Annotated burden table with ancestry, gene metadata, observed expression z scores, expected population shifts, centered z scores, percent-change estimates, and percent-change bins. |
| `QTLGeneBurdenCounts.tsv.gz` | Gene-level summary of large predicted burden events, including deletion-like and duplication-like burden counts and expression-supported events. |
| `QTLGeneBurdenQC.tsv.gz` | Gene-level QC table summarizing variant counts, effect-size distributions, missing genotype rates, nonzero burden frequency, burden magnitude, centered z-score magnitude, and observed expression outlier counts. |
| `QTLBurdenMedianGenesPerBin.tsv` | Median number of genes per individual in each percent-change bin, stratified by gene type and coding-variant filtering status. |
| `QTLBurdenOutlierEnrichment.tsv` | Enrichment of observed expression outliers across burden percent-change bins, currently focused on protein-coding genes. |

## Scripts

### `scripts/QTLBurden.R`

Computes raw QTL burden from aFC weights and genotype dosages.

Main steps:

- Reads a shard of the aFC table.
- Groups variants by `pid`.
- Uses `bedr::tabix()` to extract matching VCF records by genomic interval.
- Converts VCF genotypes to alternate allele dosages.
- Aligns genotype rows to the aFC variant order.
- Computes signed and directional load metrics per sample.
- Writes `<OutputPrefix>.QTLBurdenSummary.tsv.gz`.

Core function:

- `compute_load_metrics(beta, G)`: computes burden metrics from a vector of variant effects and a genotype dosage matrix.

### `scripts/CleanQTLBurden.R`

Annotates and summarizes the raw burden output.

Main steps:

- Loads GTF gene annotations and removes gene ID version suffixes.
- Identifies genes with likely causal coding variants from the SuSiE/eQTL annotation file.
- Joins ancestry assignments and observed expression z scores.
- Computes expected ancestry-specific burden mean and variance from allele frequencies.
- Computes centered burden effects, z scores, and percent-change bins.
- Writes a cleaned burden table.
- Creates gene-level burden count summaries.
- Summarizes burden bins per individual.
- Computes expression outlier enrichment by burden bin.

## Workflows

### `workflows/preprocess_AFC_inputs.wdl`

Preprocessing-only workflow for aFC input preparation.

Workflow:

- `preprocess_workflow_parallel`

Tasks:

- `preprocess_expression_bed`
- `annotate_vcf_ids`

### `workflows/aFC.wdl`

Runs aFC per chromosome and merges the results.

Workflow:

- `aFC_workflow_split_by_chr`

Tasks:

- `split_vcf_by_chr`
- `run_afc`
- `merge_afc_reports`

### `workflows/QuantifyQTLBurden.wdl`

End-to-end workflow. Use this when starting from aFC weights and a VCF.

Tasks:

- `shard_afc_by_gene`
- `QuantifyQTLBurden`
- `AggregateQTLBurden`
- `CleanBurdenData`

### `workflows/CleanQTLBurden.wdl`

Annotation-only workflow. Use this when `QTLBurdenSummary.AllGenes.tsv.gz` has already been generated and only the cleaning/annotation step needs to be rerun.

Task:

- `CleanBurdenData`

## Docker Images

Dockerfiles are provided under `envs/`:

| Dockerfile | Purpose |
| --- | --- |
| `envs/aFC/Dockerfile` | Runtime reference for aFC preprocessing and aFC calculation. |
| `envs/QuantifyQTLBurden/Dockerfile` | Runtime for raw burden calculation. |
| `envs/CleanQTLBurden/Dockerfile` | Runtime for annotation and summary generation. |

The WDLs currently reference GitHub Container Registry images:

```text
gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8
ghcr.io/aou-multiomics-analysis/quantifyqtlburden/quantifyqtlburden:main
ghcr.io/aou-multiomics-analysis/quantifyqtlburden/cleanqtlburden:main
```

## Terra Notes

- The main workflow is intended for Terra/Cromwell execution.
- `GenesPerShard` controls scatter size. Smaller values create more shards with lower per-task memory needs; larger values reduce task count but increase per-task work.
- The raw burden task requests 32 GB memory and a large local SSD.
- The cleaning task requests 256 GB memory because it joins and summarizes all genes, samples, ancestry assignments, allele frequencies, expression z scores, and annotations.
- For best performance, subset the VCF to variants present in the aFC input before running this workflow.

## Repository Layout

```text
.
├── workflows/
│   ├── preprocess_AFC_inputs.wdl
│   ├── aFC.wdl
│   ├── QuantifyQTLBurden.wdl
│   └── CleanQTLBurden.wdl
├── scripts/
│   ├── QTLBurden.R
│   └── CleanQTLBurden.R
├── envs/
│   ├── aFC/Dockerfile
│   ├── QuantifyQTLBurden/Dockerfile
│   └── CleanQTLBurden/Dockerfile
├── .dockstore.yml
└── README.md
```
