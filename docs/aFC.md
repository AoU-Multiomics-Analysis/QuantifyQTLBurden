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

[`../workflows/aFC.wdl`](../workflows/aFC.wdl) splits the VCF by chromosome, prepares gene-level inputs within each chromosome, runs aFC by gene, and merges per-gene reports.

Workflow:

- `aFC_workflow_split_by_gene`

Tasks:

- `build_gene_afc_manifest`: joins the QTL file to a cis-window BED and creates a gene manifest plus the chromosome list to process.
- `split_vcf_by_chr`: subsets and indexes the input VCF once per chromosome.
- `prepare_chromosome_afc_inputs`: creates gene-specific QTL files and gene cis-window VCFs for all genes on one chromosome.
- `run_afc`: runs `/opt/aFC/aFC.py` on the gene-specific VCF, expression BED, covariates, and QTL file.
- `merge_afc_reports`: concatenates gene-level aFC outputs into one gzipped report.

Inputs:

| Input | Type | Description |
| --- | --- | --- |
| `vcf_file` | File | Annotated, indexed genotype VCF. |
| `vcf_index` | File | Tabix index for `vcf_file`. |
| `expression_bed` | File | bgzip-compressed, tabix-indexed expression BED. |
| `expression_bed_index` | File | Tabix index for `expression_bed`. |
| `covariates_file` | File | Covariates file in the format expected by aFC. |
| `afc_qtl_file` | File | QTL file containing at least `pid`, `sid`, `sid_chr`, and `sid_pos`. |
| `cis_window_bed` | File | BED-like gene cis-window table with columns `#chr`/`chr`, `start`, `end`, and `gene_id`. Only genes present in `afc_qtl_file` are scattered. |
| `prefix` | String | Output file prefix. |
| `memory` | Int | Memory in GB. Default: `16`. |
| `disk_space` | Int | Extra disk space in GB. Default: `50`. |
| `num_threads` | Int | Number of CPU threads. Default: `8`. |
| `num_preempt` | Int | Number of preemptible retries. Default: `0`. |

Outputs:

| Output | Description |
| --- | --- |
| `gene_manifest` | Gene scatter manifest with `chr`, `start`, `end`, `gene_id`, and safe file-name gene ID columns. |
| `per_chr_vcfs` | Per-chromosome VCFs generated before gene-level preprocessing. |
| `per_chr_vcf_indexes` | Tabix indexes for `per_chr_vcfs`. |
| `chromosome_gene_windows` | Per-chromosome gene-window manifests used to generate gene-specific VCFs and QTL files. |
| `per_gene_afc_reports` | Per-gene aFC reports: `<prefix>.<gene_id>.aFC.txt.gz` using a file-safe gene ID. |
| `final_afc_report` | Merged aFC report across genes: `<prefix>.aFC.txt.gz`. |

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

The cis-window BED must contain at least:

| Column | Description |
| --- | --- |
| `#chr` or `chr` | Chromosome for the gene cis window. |
| `start` | Start coordinate for the cis window. |
| `end` | End coordinate for the cis window. |
| `gene_id` | Gene ID matching `pid` in the QTL file. |
