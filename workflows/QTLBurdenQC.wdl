version 1.0

task ComputeQTLBurdenQC {
  input {
    File CleanedQTLBurden
    File aFC
    Float MissingnessFailThreshold = 0.10
    Float MissingnessWarnThreshold = 0.05
    Float ContextMissingnessFailThreshold = 0.05
    Float VarianceRatioLower = 0.20
    Float VarianceRatioUpper = 5
    Float VarianceRatioWarnLower = 0.35
    Float VarianceRatioWarnUpper = 3
    Float TailZFailThreshold = 25
    Float TailZWarnThreshold = 15
    Float TailZWarnSampleProportionThreshold = 0.05
    Float DominantVariantWarnThreshold = 0.98
  }

  command <<<
  Rscript /tmp/QTLBurdenQCMetrics.R \
    --CleanedQTLBurden ~{CleanedQTLBurden} \
    --aFC ~{aFC} \
    --MissingnessFailThreshold ~{MissingnessFailThreshold} \
    --MissingnessWarnThreshold ~{MissingnessWarnThreshold} \
    --ContextMissingnessFailThreshold ~{ContextMissingnessFailThreshold} \
    --VarianceRatioLower ~{VarianceRatioLower} \
    --VarianceRatioUpper ~{VarianceRatioUpper} \
    --VarianceRatioWarnLower ~{VarianceRatioWarnLower} \
    --VarianceRatioWarnUpper ~{VarianceRatioWarnUpper} \
    --TailZFailThreshold ~{TailZFailThreshold} \
    --TailZWarnThreshold ~{TailZWarnThreshold} \
    --TailZWarnSampleProportionThreshold ~{TailZWarnSampleProportionThreshold} \
    --DominantVariantWarnThreshold ~{DominantVariantWarnThreshold}
  >>>

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/cleanqtlburden:main"
    memory: "256G"
    cpu: 2
    disks: "local-disk 2500 SSD"
  }

  output {
    File qtl_gene_burden_qc = "QTLGeneBurdenQC.tsv.gz"
    File qtl_gene_burden_status_list = "QTLGeneBurdenStatusList.tsv.gz"
  }
}

task ComputeQTLBurdenOutlierEnrichment {
  input {
    File CleanedQTLBurden
    Int OutlierPermutationIterations = 200
  }

  command <<<
  Rscript /tmp/QTLBurdenOutlierEnrichment.R \
    --CleanedQTLBurden ~{CleanedQTLBurden} \
    --OutlierPermutationIterations ~{OutlierPermutationIterations}
  >>>

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/cleanqtlburden:main"
    memory: "256G"
    cpu: 2
    disks: "local-disk 2500 SSD"
  }

  output {
    File qtl_burden_outlier_enrichment = "QTLBurdenOutlierEnrichment.tsv"
    File qtl_burden_outlier_enrichment_permutation = "QTLBurdenOutlierEnrichmentPermutation.tsv"
    File qtl_burden_outlier_enrichment_noisy_filter_impact = "QTLBurdenOutlierEnrichmentNoisyFilterImpact.tsv"
  }
}

task ComputeQTLBurdenMedianGenesPerBin {
  input {
    File CleanedQTLBurden
    File aFC
    Float DominantVariantWarnThreshold = 0.98
  }

  command <<<
  Rscript /tmp/QTLBurdenMedianGeneSummaries.R \
    --CleanedQTLBurden ~{CleanedQTLBurden} \
    --aFC ~{aFC} \
    --DominantVariantWarnThreshold ~{DominantVariantWarnThreshold}
  >>>

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/cleanqtlburden:main"
    memory: "512G"
    cpu: 2
    disks: "local-disk 2500 SSD"
  }

  output {
    File qtl_burden_median_genes_per_bin_all_calls = "QTLBurdenMedianGenesPerBin_AllCalls.tsv"
    File qtl_burden_median_genes_per_bin_high_confidence = "QTLBurdenMedianGenesPerBin_HighConfidence.tsv"
    File qtl_burden_median_genes_per_bin_high_ko_removed = "QTLBurdenMedianGenesPerBin_HighKORemoved.tsv"
    File qtl_burden_median_genes_per_bin_by_gene_set_all_calls = "QTLBurdenMedianGenesPerBinByGeneSet_AllCalls.tsv"
    File qtl_burden_median_genes_per_bin_by_gene_set_high_confidence = "QTLBurdenMedianGenesPerBinByGeneSet_HighConfidence.tsv"
    File qtl_burden_median_genes_per_bin_by_gene_set_high_ko_removed = "QTLBurdenMedianGenesPerBinByGeneSet_HighKORemoved.tsv"
    File qtl_burden_median_genes_per_bin_dosage_all_calls = "QTLBurdenMedianGenesPerBin_AllCalls_DosageNoisyFiltered.tsv"
    File qtl_burden_median_genes_per_bin_dosage_high_confidence = "QTLBurdenMedianGenesPerBin_HighConfidence_DosageNoisyFiltered.tsv"
    File qtl_burden_median_genes_per_bin_dosage_high_ko_removed = "QTLBurdenMedianGenesPerBin_HighKORemoved_DosageNoisyFiltered.tsv"
    File qtl_burden_median_genes_per_bin_by_gene_set_dosage_all_calls = "QTLBurdenMedianGenesPerBinByGeneSet_AllCalls_DosageNoisyFiltered.tsv"
    File qtl_burden_median_genes_per_bin_by_gene_set_dosage_high_confidence = "QTLBurdenMedianGenesPerBinByGeneSet_HighConfidence_DosageNoisyFiltered.tsv"
    File qtl_burden_median_genes_per_bin_by_gene_set_dosage_high_ko_removed = "QTLBurdenMedianGenesPerBinByGeneSet_HighKORemoved_DosageNoisyFiltered.tsv"
  }
}

task PlotQTLBurdenQC {
  input {
    File QTLGeneBurdenQC
    File QTLBurdenMedianGenesPerBinByGeneSet
    String QCPlotPrefix = "QTLGeneBurdenQC"
  }

  command <<<
    if [ -f /tmp/PlotQTLBurdenQC.R ]; then
      Rscript /tmp/PlotQTLBurdenQC.R \
        --QCFile ~{QTLGeneBurdenQC} \
        --MedianGenesPerBinByGeneSet ~{QTLBurdenMedianGenesPerBinByGeneSet} \
        --Prefix ~{QCPlotPrefix}
    else
      echo "PlotQTLBurdenQC.R not found in /tmp; creating placeholder outputs."
      touch ~{QCPlotPrefix}_StatusByGeneType.pdf
      touch ~{QCPlotPrefix}_StatusOverall.pdf
      touch ~{QCPlotPrefix}_Missingness.pdf
      touch ~{QCPlotPrefix}_VarianceRatio.pdf
      touch ~{QCPlotPrefix}_TailZ.pdf
      touch ~{QCPlotPrefix}_DominantVariantFraction.pdf
      touch ~{QCPlotPrefix}_ReasonCounts.pdf
      touch ~{QCPlotPrefix}_GeneSetMedianImpact.pdf
    fi
  >>>

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/cleanqtlburden:main"
    memory: "256G"
    cpu: 1
    disks: "local-disk 2500 SSD"
  }

  output {
    File qc_plot_status_by_gene_type = "~{QCPlotPrefix}_StatusByGeneType.pdf"
    File qc_plot_status_overall = "~{QCPlotPrefix}_StatusOverall.pdf"
    File qc_plot_missingness = "~{QCPlotPrefix}_Missingness.pdf"
    File qc_plot_variance_ratio = "~{QCPlotPrefix}_VarianceRatio.pdf"
    File qc_plot_tail_z = "~{QCPlotPrefix}_TailZ.pdf"
    File qc_plot_dominant_variant_fraction = "~{QCPlotPrefix}_DominantVariantFraction.pdf"
    File qc_plot_reason_counts = "~{QCPlotPrefix}_ReasonCounts.pdf"
    File qc_plot_gene_set_median_impact = "~{QCPlotPrefix}_GeneSetMedianImpact.pdf"
  }
}

workflow QTLBurdenQC {
  input {
    File input_cleaned_qtl_burden
    File input_aFC_weights
    Int input_outlier_permutation_iterations = 200
    Float input_missingness_fail_threshold = 0.10
    Float input_missingness_warn_threshold = 0.05
    Float input_context_missingness_fail_threshold = 0.05
    Float input_variance_ratio_lower = 0.20
    Float input_variance_ratio_upper = 5
    Float input_variance_ratio_warn_lower = 0.35
    Float input_variance_ratio_warn_upper = 3
    Float input_tail_z_fail_threshold = 25
    Float input_tail_z_warn_threshold = 15
    Float input_tail_z_warn_sample_proportion_threshold = 0.05
    Float input_dominant_variant_warn_threshold = 0.98
    String input_qc_plot_prefix = "QTLGeneBurdenQC"
  }

    call ComputeQTLBurdenOutlierEnrichment {
      input:
        CleanedQTLBurden = input_cleaned_qtl_burden,
        OutlierPermutationIterations = input_outlier_permutation_iterations
    }

    call ComputeQTLBurdenMedianGenesPerBin {
      input:
        CleanedQTLBurden = input_cleaned_qtl_burden,
        aFC = input_aFC_weights,
        DominantVariantWarnThreshold = input_dominant_variant_warn_threshold
    }

    call ComputeQTLBurdenQC {
      input:
        CleanedQTLBurden = input_cleaned_qtl_burden,
        aFC = input_aFC_weights,
        MissingnessFailThreshold = input_missingness_fail_threshold,
        MissingnessWarnThreshold = input_missingness_warn_threshold,
        ContextMissingnessFailThreshold = input_context_missingness_fail_threshold,
        VarianceRatioLower = input_variance_ratio_lower,
        VarianceRatioUpper = input_variance_ratio_upper,
        VarianceRatioWarnLower = input_variance_ratio_warn_lower,
        VarianceRatioWarnUpper = input_variance_ratio_warn_upper,
        TailZFailThreshold = input_tail_z_fail_threshold,
        TailZWarnThreshold = input_tail_z_warn_threshold,
        TailZWarnSampleProportionThreshold = input_tail_z_warn_sample_proportion_threshold,
        DominantVariantWarnThreshold = input_dominant_variant_warn_threshold
    }

  call PlotQTLBurdenQC {
    input:
      QTLGeneBurdenQC = ComputeQTLBurdenQC.qtl_gene_burden_qc,
      QTLBurdenMedianGenesPerBinByGeneSet = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_by_gene_set_all_calls,
      QCPlotPrefix = input_qc_plot_prefix
  }

  output {
    File CleanedBurden = input_cleaned_qtl_burden
    File QTLGeneBurdenQC = ComputeQTLBurdenQC.qtl_gene_burden_qc
    File QTLGeneBurdenStatusList = ComputeQTLBurdenQC.qtl_gene_burden_status_list
    File QTLBurdenOutlierEnrichment = ComputeQTLBurdenOutlierEnrichment.qtl_burden_outlier_enrichment
    File QTLBurdenOutlierEnrichmentNoisyFilterImpact = ComputeQTLBurdenOutlierEnrichment.qtl_burden_outlier_enrichment_noisy_filter_impact
    File QTLBurdenMedianGenesPerBin_AllCalls = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_all_calls
    File QTLBurdenMedianGenesPerBin_HighConfidence = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_high_confidence
    File QTLBurdenMedianGenesPerBin_HighKORemoved = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_high_ko_removed
    File QTLBurdenMedianGenesPerBinByGeneSet_AllCalls = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_by_gene_set_all_calls
    File QTLBurdenMedianGenesPerBinByGeneSet_HighConfidence = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_by_gene_set_high_confidence
    File QTLBurdenMedianGenesPerBinByGeneSet_HighKORemoved = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_by_gene_set_high_ko_removed
    File QTLBurdenMedianGenesPerBin_AllCalls_DosageNoisyFiltered = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_dosage_all_calls
    File QTLBurdenMedianGenesPerBin_HighConfidence_DosageNoisyFiltered = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_dosage_high_confidence
    File QTLBurdenMedianGenesPerBin_HighKORemoved_DosageNoisyFiltered = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_dosage_high_ko_removed
    File QTLBurdenMedianGenesPerBinByGeneSet_AllCalls_DosageNoisyFiltered = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_by_gene_set_dosage_all_calls
    File QTLBurdenMedianGenesPerBinByGeneSet_HighConfidence_DosageNoisyFiltered = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_by_gene_set_dosage_high_confidence
    File QTLBurdenMedianGenesPerBinByGeneSet_HighKORemoved_DosageNoisyFiltered = ComputeQTLBurdenMedianGenesPerBin.qtl_burden_median_genes_per_bin_by_gene_set_dosage_high_ko_removed
    File QTLBurdenOutlierEnrichmentPermutation = ComputeQTLBurdenOutlierEnrichment.qtl_burden_outlier_enrichment_permutation
    File QTLGeneBurdenQC_StatusByGeneType = PlotQTLBurdenQC.qc_plot_status_by_gene_type
    File QTLGeneBurdenQC_StatusOverall = PlotQTLBurdenQC.qc_plot_status_overall
    File QTLGeneBurdenQC_Missingness = PlotQTLBurdenQC.qc_plot_missingness
    File QTLGeneBurdenQC_VarianceRatio = PlotQTLBurdenQC.qc_plot_variance_ratio
    File QTLGeneBurdenQC_TailZ = PlotQTLBurdenQC.qc_plot_tail_z
    File QTLGeneBurdenQC_DominantVariantFraction = PlotQTLBurdenQC.qc_plot_dominant_variant_fraction
    File QTLGeneBurdenQC_ReasonCounts = PlotQTLBurdenQC.qc_plot_reason_counts
    File QTLGeneBurdenQC_GeneSetMedianImpact = PlotQTLBurdenQC.qc_plot_gene_set_median_impact
  }
}
