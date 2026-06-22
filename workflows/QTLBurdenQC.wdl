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
    --DominantVariantWarnThreshold ~{DominantVariantWarnThreshold}
  >>>

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/cleanqtlburden:main"
    memory: "32G"
    cpu: 2
    disks: "local-disk 2500 SSD"
  }

  output {
    File QTLGeneBurdenQC = "QTLGeneBurdenQC.tsv.gz"
    File QTLGeneBurdenStatusList = "QTLGeneBurdenStatusList.tsv.gz"
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
    memory: "32G"
    cpu: 2
    disks: "local-disk 2500 SSD"
  }

  output {
    File QTLBurdenOutlierEnrichment = "QTLBurdenOutlierEnrichment.tsv"
    File QTLBurdenOutlierEnrichmentPermutation = "QTLBurdenOutlierEnrichmentPermutation.tsv"
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
    memory: "32G"
    cpu: 2
    disks: "local-disk 2500 SSD"
  }

  output {
    File QTLBurdenMedianGenesPerBin = "QTLBurdenMedianGenesPerBin.tsv"
    File QTLBurdenMedianGenesPerBinByGeneSet = "QTLBurdenMedianGenesPerBinByGeneSet.tsv"
    File QTLBurdenMedianGenesPerBinDosageNoisyFiltered = "QTLBurdenMedianGenesPerBin_DosageNoisyFiltered.tsv"
    File QTLBurdenMedianGenesPerBinByGeneSetDosageNoisyFiltered = "QTLBurdenMedianGenesPerBinByGeneSet_DosageNoisyFiltered.tsv"
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
    memory: "8G"
    cpu: 1
    disks: "local-disk 2500 SSD"
  }

  output {
    File QTLGeneBurdenQC_StatusByGeneType = "~{QCPlotPrefix}_StatusByGeneType.pdf"
    File QTLGeneBurdenQC_StatusOverall = "~{QCPlotPrefix}_StatusOverall.pdf"
    File QTLGeneBurdenQC_Missingness = "~{QCPlotPrefix}_Missingness.pdf"
    File QTLGeneBurdenQC_VarianceRatio = "~{QCPlotPrefix}_VarianceRatio.pdf"
    File QTLGeneBurdenQC_TailZ = "~{QCPlotPrefix}_TailZ.pdf"
    File QTLGeneBurdenQC_DominantVariantFraction = "~{QCPlotPrefix}_DominantVariantFraction.pdf"
    File QTLGeneBurdenQC_ReasonCounts = "~{QCPlotPrefix}_ReasonCounts.pdf"
    File QTLGeneBurdenQC_GeneSetMedianImpact = "~{QCPlotPrefix}_GeneSetMedianImpact.pdf"
  }
}

workflow QTLBurdenQC {
  input {
    File CleanedQTLBurden
    File aFCWeights
    String QCPlotPrefix = "QTLGeneBurdenQC"
    Int OutlierPermutationIterations = 200
    Float MissingnessFailThreshold = 0.10
    Float MissingnessWarnThreshold = 0.05
    Float ContextMissingnessFailThreshold = 0.05
    Float VarianceRatioLower = 0.20
    Float VarianceRatioUpper = 5
    Float VarianceRatioWarnLower = 0.35
    Float VarianceRatioWarnUpper = 3
    Float TailZFailThreshold = 25
    Float TailZWarnThreshold = 15
    Float DominantVariantWarnThreshold = 0.98
    File? QTLBurdenCounts
    File? QTLBurdenPerSampleGene
  }

    call ComputeQTLBurdenOutlierEnrichment {
      input:
        CleanedQTLBurden = CleanedQTLBurden,
        OutlierPermutationIterations = OutlierPermutationIterations
    }

    call ComputeQTLBurdenMedianGenesPerBin {
      input:
        CleanedQTLBurden = CleanedQTLBurden,
        aFC = aFCWeights,
        DominantVariantWarnThreshold = DominantVariantWarnThreshold
    }

    call ComputeQTLBurdenQC {
      input:
        CleanedQTLBurden = CleanedQTLBurden,
        aFC = aFCWeights,
        MissingnessFailThreshold = MissingnessFailThreshold,
        MissingnessWarnThreshold = MissingnessWarnThreshold,
        ContextMissingnessFailThreshold = ContextMissingnessFailThreshold,
        VarianceRatioLower = VarianceRatioLower,
        VarianceRatioUpper = VarianceRatioUpper,
        VarianceRatioWarnLower = VarianceRatioWarnLower,
        VarianceRatioWarnUpper = VarianceRatioWarnUpper,
        TailZFailThreshold = TailZFailThreshold,
        TailZWarnThreshold = TailZWarnThreshold,
        DominantVariantWarnThreshold = DominantVariantWarnThreshold
    }

  call PlotQTLBurdenQC {
    input:
      QTLGeneBurdenQC = ComputeQTLBurdenQC.QTLGeneBurdenQC,
      QTLBurdenMedianGenesPerBinByGeneSet = ComputeQTLBurdenMedianGenesPerBin.QTLBurdenMedianGenesPerBinByGeneSet,
      QCPlotPrefix = QCPlotPrefix
  }

  output {
    File CleanedBurden = CleanedQTLBurden
    File? QTLBurdenCounts = QTLBurdenCounts
    File QTLGeneBurdenQC = ComputeQTLBurdenQC.QTLGeneBurdenQC
    File QTLGeneBurdenStatusList = ComputeQTLBurdenQC.QTLGeneBurdenStatusList
    File QTLBurdenOutlierEnrichment = ComputeQTLBurdenOutlierEnrichment.QTLBurdenOutlierEnrichment
    File QTLBurdenMedianGenesPerBin = ComputeQTLBurdenMedianGenesPerBin.QTLBurdenMedianGenesPerBin 
    File QTLBurdenMedianGenesPerBinByGeneSet = ComputeQTLBurdenMedianGenesPerBin.QTLBurdenMedianGenesPerBinByGeneSet
    File QTLBurdenMedianGenesPerBinDosageNoisyFiltered = ComputeQTLBurdenMedianGenesPerBin.QTLBurdenMedianGenesPerBinDosageNoisyFiltered
    File QTLBurdenMedianGenesPerBinByGeneSetDosageNoisyFiltered = ComputeQTLBurdenMedianGenesPerBin.QTLBurdenMedianGenesPerBinByGeneSetDosageNoisyFiltered
    File QTLBurdenOutlierEnrichmentPermutation = ComputeQTLBurdenOutlierEnrichment.QTLBurdenOutlierEnrichmentPermutation
    File? QTLBurdenPerSampleGene = QTLBurdenPerSampleGene
    File QTLGeneBurdenQC_StatusByGeneType = PlotQTLBurdenQC.QTLGeneBurdenQC_StatusByGeneType
    File QTLGeneBurdenQC_StatusOverall = PlotQTLBurdenQC.QTLGeneBurdenQC_StatusOverall
    File QTLGeneBurdenQC_Missingness = PlotQTLBurdenQC.QTLGeneBurdenQC_Missingness
    File QTLGeneBurdenQC_VarianceRatio = PlotQTLBurdenQC.QTLGeneBurdenQC_VarianceRatio
    File QTLGeneBurdenQC_TailZ = PlotQTLBurdenQC.QTLGeneBurdenQC_TailZ
    File QTLGeneBurdenQC_DominantVariantFraction = PlotQTLBurdenQC.QTLGeneBurdenQC_DominantVariantFraction
    File QTLGeneBurdenQC_ReasonCounts = PlotQTLBurdenQC.QTLGeneBurdenQC_ReasonCounts
    File QTLGeneBurdenQC_GeneSetMedianImpact = PlotQTLBurdenQC.QTLGeneBurdenQC_GeneSetMedianImpact
  }
}
