# Proteomics validation VCF subsetting

This repository includes [`workflows/SubsetVCF.wdl`](../workflows/SubsetVCF.wdl), which is used to build a validation cohort for proteomics analyses.

Current use case:
- Create a VCF containing only the cohort used for validation (for example, participants with proteomics but without RNA-seq profiling), and
- Restrict variants to the same genetic signals being used for burden-based validation.

In this project, we used this to:
- Build a proteomics-only genotype subset.
- Focus on high-confidence eQTL signals (high `pip`) at genes that have both proteomic and transcriptomic measurements.
- Use that subset as an independent validation set for downstream burden interpretation.

## Inputs

[`workflows/SubsetVCF.wdl`](../workflows/SubsetVCF.wdl) expects:

- `vcf_file`: `bgzip`-compressed, tabix-indexed VCF.
- `vcf_index`: matching index file (`.vcf.gz.tbi` or `.vcf.gz.csi`).
- `variant_list`: bcftools `-R` region file (no header). This should be plain region coordinates.
- `sample_list`: plain sample ID list (no header), one sample per line.
- `OutputPrefix`: optional output prefix for the subset VCF.

## Notes on variant filtering

The `variant_list` input is consumed by `bcftools view -R`, so it expects regions (for example `chr:start-end`), not `sid` strings (`CHROM:POS_REF_ALT`).

- Example accepted formats:
  - `chr1:1000000-1000000`
  - `chr1    1000000    1000000`

This is useful when you already have a filtered list of loci to keep and want fast, deterministic subsetting for validation.

## Typical flow for the proteomics cohort

1. Identify high-pip eQTL variants at the target gene set (intersection of genes with proteomic and transcriptomic observations).
2. Convert that variant list to regions for bcftools `-R`.
3. Build a sample list of proteomics-only individuals (exclude RNA-seq profiled individuals if that is the validation design).
4. Run `SubsetVCF.wdl` to generate:
   - `subset_vcf` (`.vcf.gz`)
   - `subset_vcf_index` (`.vcf.gz.tbi`)
5. Pass the subset outputs into the downstream burden steps as needed.
