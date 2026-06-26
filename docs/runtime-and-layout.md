# Runtime And Layout

## Docker Images

Dockerfiles are provided under [`../envs`](../envs):

| Dockerfile | Purpose |
| --- | --- |
| [`../envs/aFC/Dockerfile`](../envs/aFC/Dockerfile) | Runtime reference for aFC preprocessing and aFC calculation. |
| [`../envs/QuantifyQTLBurden/Dockerfile`](../envs/QuantifyQTLBurden/Dockerfile) | Runtime for raw burden calculation. |
| [`../envs/CleanQTLBurden/Dockerfile`](../envs/CleanQTLBurden/Dockerfile) | Runtime for burden cleaning and annotation. |
| [`../envs/QTLBurdenQC/Dockerfile`](../envs/QTLBurdenQC/Dockerfile) | Runtime for QC metrics, outlier enrichment, median-bin summaries, and QC plots. |

The WDLs currently reference these images:

```text
gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8
ghcr.io/aou-multiomics-analysis/quantifyqtlburden/quantifyqtlburden:main
ghcr.io/aou-multiomics-analysis/quantifyqtlburden/cleanqtlburden:main
ghcr.io/aou-multiomics-analysis/quantifyqtlburden/qtlburdenqc:main
```

## Terra Notes

- The workflows are intended for Terra/Cromwell execution.
- `GenesPerShard` controls QTL burden scatter size.
- Smaller `GenesPerShard` values create more shards with lower per-task memory needs.
- Larger `GenesPerShard` values reduce task count but increase per-task work.
- The raw burden task requests 32 GB memory and a large local SSD.
- The cleaning task requests 256 GB memory because it joins all genes, samples, ancestry assignments, allele frequencies, expression z scores, and annotations.
- The QC tasks use a separate `qtlburdenqc` image so QC script changes can rebuild independently from cleaning/annotation changes.
- For best performance, subset the VCF to variants present in the aFC weights before running the QTL burden workflow.

## Dockstore

[`../.dockstore.yml`](../.dockstore.yml) registers four independent WDL entries:

- `PreprocessAFCInputs`
- `aFC`
- `QuantifyQTLBurden`
- `CleanQTLBurdenTask`

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
│   ├── CleanQTLBurden/Dockerfile
│   └── QTLBurdenQC/Dockerfile
├── docs/
│   ├── aFC.md
│   ├── qtl-burden.md
│   └── runtime-and-layout.md
├── .dockstore.yml
└── README.md
```
