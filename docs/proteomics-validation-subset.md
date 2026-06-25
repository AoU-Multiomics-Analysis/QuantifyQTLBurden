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
- `vcf_index`: matching index file (`.vcf.gz.tbi`). This workflow no longer generates indexes.
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
   - `input_index` (pass-through): the input `.vcf.gz.tbi` supplied to `vcf_index`.
5. Pass the subset outputs into the downstream burden steps as needed.

## Orthogonal proteomics validation strategy

The downstream logic uses the existing burden outputs plus measured proteomics values to run an **orthogonal validation**:

- We score a subset of samples that did not have transcriptomic profiling (proteomics-only individuals).
- For those samples, we still have proteomic z scores at a set of genes/proteins.
- We test whether predicted high-burden bins are enriched for proteomic outliers compared with a near-null reference bin.

In practice, `CleanQTLBurden.R` marks proteomics outliers as:
- `UpProteomicsOutlier = ObservedProteomicsZ > 4`
- `DownProteomicsOutlier = ObservedProteomicsZ < -4`

`QTLBurdenOutlierEnrichment` then performs bin-wise enrichment tests for:
- burden `type = Up` (`>= 50%` predicted gain-like change), and
- burden `type = Down` (`<= -50%` predicted loss-like change),

against the reference bin `(-10,10)` in percent change.

For each burden model (`AllCalls`, `HighConfidence`, `HighKORemoved`), it reports:
- enrichment odds ratio/log-odds ratio,
- Fisher test p-values,
- and permutation-based empirical p-values where available.

Observed enrichment of proteomic outlier calls in the extreme burden bins supports that the genotype-derived burden signal is mirrored at the proteome level, which is useful as validation independent of transcriptomic measurements.
