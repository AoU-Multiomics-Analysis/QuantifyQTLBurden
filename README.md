# QuantifyQTLBurden

QuantifyQTLBurden estimates individual-level QTL burden from fine-mapped QTL effect sizes and genotypes. This repository also includes independent workflows for computing allelic fold change (aFC), because aFC estimation is the upstream step that produces the effect-size weights used by QTL burden quantification.

The workflows stay separate:

1. Run the aFC workflows when starting from QTL mapping results, genotypes, expression BED files, and covariates.
2. Run the QuantifyQTLBurden workflows when starting from aFC weights and genotypes.

## Workflows

| Stage | Workflow | Use when |
| --- | --- | --- |
| aFC preprocessing | [`workflows/preprocess_AFC_inputs.wdl`](workflows/preprocess_AFC_inputs.wdl) | VCF IDs or expression BED indexing need to be prepared for aFC. |
| aFC calculation | [`workflows/aFC.wdl`](workflows/aFC.wdl) | Computing aFC weights by chromosome and merging them. |
| QTL burden | [`workflows/QuantifyQTLBurden.wdl`](workflows/QuantifyQTLBurden.wdl) | Computing per-individual, per-gene burden from aFC weights and genotypes. |
| Burden cleaning | [`workflows/CleanQTLBurden.wdl`](workflows/CleanQTLBurden.wdl) | Annotating an existing merged burden file. |
| VCF subsetting | [`workflows/SubsetVCF.wdl`](workflows/SubsetVCF.wdl) | Subsetting cohorts or variants (commonly used for proteomics-only validation cohorts). |

## Quick Start

1. Prepare a QTL file with at least `pid`, `sid`, `sid_chr`, and `sid_pos`.
2. Run [`workflows/preprocess_AFC_inputs.wdl`](workflows/preprocess_AFC_inputs.wdl) unless your VCF and expression BED are already formatted and indexed for aFC.
3. Run [`workflows/aFC.wdl`](workflows/aFC.wdl) to produce `<prefix>.aFC.txt.gz`.
4. Filter or prepare the resulting aFC weights for burden analysis.
5. Run [`workflows/QuantifyQTLBurden.wdl`](workflows/QuantifyQTLBurden.wdl).
6. Re-run [`workflows/CleanQTLBurden.wdl`](workflows/CleanQTLBurden.wdl) only when cleaning or annotation needs to be repeated separately.
7. Optionally run [`workflows/SubsetVCF.wdl`](workflows/SubsetVCF.wdl) first when generating a proteomics-only validation cohort.

## Documentation

- [aFC workflows](docs/aFC.md)
- [QTL burden workflows](docs/qtl-burden.md)
- [Burden tail probabilities](docs/burden-probability.md)
- [Runtime, Terra, and repository layout](docs/runtime-and-layout.md)
- [VCF subset workflow (proteomics validation)](docs/proteomics-validation-subset.md)

## Scientific Summary

For a gene with causal QTL variants, the primary predicted burden is:

```text
predicted_effect = sum(genotype_dosage_variant * log2_aFC_variant)
```

Positive `log2_aFC` values increase predicted expression, and negative values decrease predicted expression. The cleaning workflow can also center predicted burden by ancestry-specific allele frequencies and compare predicted expression shifts to observed expression z scores.

## Burden probability and extreme-call modeling

The quantification model computes analytic, one-sided tail probabilities for each sample/gene from the per-individual total log2 burden:

- [`docs/burden-probability.md`](docs/burden-probability.md) describes the exact formulas (`burden_probability_loss50`, `burden_probability_gain50`) and interpretation.
- In practice:
  - deletion/loss direction => `P(L <= loss_threshold)`
  - duplication/gain direction => `P(L >= gain_threshold)`
  - for each sample, this is the chance the dosage-weighted perturbation pushes the gene at least to the configured threshold, given the effect-size uncertainty (`log2_aFC` means and SEs).
- Both directional per-sample burden probabilities are written as:
    - `burden_probability_loss50`
    - `burden_probability_gain50`.

If you also want the median-count strategy (all calls vs noisy-extreme filtered and weighted medians), that is documented under [QTL burden outputs](docs/qtl-burden.md).

## Repository Layout

```text
.
├── workflows/
├── scripts/
├── envs/
├── docs/
├── .dockstore.yml
└── README.md
```
