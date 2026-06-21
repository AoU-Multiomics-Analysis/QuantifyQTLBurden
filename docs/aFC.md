# aFC Workflows

The aFC workflows were merged from <https://github.com/AoU-Multiomics-Analysis/aFC> and are based on the aFC tool described at <https://github.com/secastel/aFC>.

These workflows are independent from the QTL burden workflows. Their main output, `<prefix>.aFC.txt.gz`, can be used as an upstream effect-size source for burden quantification.

## Execution Order

1. Prepare a QTL file with at least `pid`, `sid`, `sid_chr`, and `sid_pos`, where `sid` uses `CHROM:POS_REF_ALT`.
2. Run [`../workflows/preprocess_AFC_inputs.wdl`](../workflows/preprocess_AFC_inputs.wdl) unless the VCF and expression BED are already correctly formatted and indexed.
3. Run [`../workflows/aFC.wdl`](../workflows/aFC.wdl) to produce `<prefix>.aFC.txt.gz`.

## Preprocessing Workflow

[`../workflows/preprocess_AFC_inputs.wdl`](../workflows/preprocess_AFC_inputs.wdl) prepares genotype and expression inputs for aFC.

Workflow:

- `preprocess_workflow_parallel`

Tasks:

- `preprocess_expression_bed`: re-compresses the expression BED with `bgzip` and creates a tabix index.
- `annotate_vcf_ids`: annotates VCF variant IDs in `CHROM:POS_REF_ALT` format so they match the QTL file `sid` field.

Inputs:

| Input | Type | Description |
| --- | --- | --- |
| `vcf_file` | File | Input genotype VCF. |
| `expression_bed` | File | Expression BED file. |
| `prefix` | String | Output file prefix. |
| `memory` | Int | Memory in GB. Default: `16`. |
| `disk_space` | Int | Extra disk space in GB. Default: `50`. |
| `num_threads` | Int | Number of CPU threads. Default: `8`. |
| `num_preempt` | Int | Number of preemptible retries. Default: `0`. |

Outputs:

| Output | Description |
| --- | --- |
| `processed_expression_bed` | bgzip-compressed expression BED: `<prefix>.processed_expression.bed.gz`. |
| `processed_expression_bed_tbi` | Tabix index for the processed expression BED. |
| `annotated_vcf` | VCF with variant IDs annotated as `CHROM:POS_REF_ALT`: `<prefix>.annotated.vcf.gz`. |
| `annotated_vcf_tbi` | Tabix index for the annotated VCF. |

## aFC Workflow

[`../workflows/aFC.wdl`](../workflows/aFC.wdl) runs aFC by chromosome and merges per-chromosome reports.

Workflow:

- `aFC_workflow_split_by_chr`

Tasks:

- `split_vcf_by_chr`: subsets and indexes the VCF for one chromosome.
- `run_afc`: runs `/opt/aFC/aFC.py` on the chromosome-specific VCF, expression BED, covariates, and QTL file.
- `merge_afc_reports`: concatenates chromosome-level aFC outputs into one gzipped report.

Inputs:

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

Outputs:

| Output | Description |
| --- | --- |
| `per_chr_afc_reports` | Per-chromosome aFC reports: `<prefix>.<chr>.aFC.txt.gz`. |
| `final_afc_report` | Merged aFC report across chromosomes: `<prefix>.aFC.txt.gz`. |

## Data Requirements

VCF files should be compressed and indexed. Variant IDs must match the QTL file `sid` format, `CHROM:POS_REF_ALT`.

Expression BED files must be bgzip-compressed and tabix-indexed.

The QTL file must contain at least:

| Column | Description |
| --- | --- |
| `pid` | Phenotype or gene ID. |
| `sid` | Variant ID in `CHROM:POS_REF_ALT` format. |
| `sid_chr` | Variant chromosome. |
| `sid_pos` | Variant position. |
