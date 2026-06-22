version 1.0

import "CleanQTLBurden.wdl" as clean

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
    File MergedQTLBurden
    File AlleleFrequencies
    File ExpressionZscores
    File aFCWeights
    File AncestryAssignments
    File eQTLSusie
    File GTF
    String QCPlotPrefix = "QTLGeneBurdenQC"
    Float BurdenTailProbability = 0.9
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
  }

    call clean.CleanBurdenData as CleanBurdenData {
      input:
        MergedQTLBurden = MergedQTLBurden,
      AlleleFrequencies = AlleleFrequencies,
      ExpressionZscores = ExpressionZscores,
      aFC = aFCWeights,
      AncestryAssignments = AncestryAssignments,
        eQTLSusie = eQTLSusie,
        GTF = GTF,
        BurdenTailProbability = BurdenTailProbability
    }

    call ComputeQTLBurdenOutlierEnrichment {
      input:
        CleanedQTLBurden = CleanBurdenData.QTLBurdenSummaryCleaned,
        OutlierPermutationIterations = OutlierPermutationIterations
    }

    call ComputeQTLBurdenQC {
      input:
        CleanedQTLBurden = CleanBurdenData.QTLBurdenSummaryCleaned,
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
      QTLBurdenMedianGenesPerBinByGeneSet = CleanBurdenData.QTLBurdenMedianGenesPerBinByGeneSet,
      QCPlotPrefix = QCPlotPrefix
  }

  output {
    File CleanedBurden = CleanBurdenData.QTLBurdenSummaryCleaned
    File QTLBurdenCounts = CleanBurdenData.QTLBurdenCounts
    File QTLGeneBurdenQC = ComputeQTLBurdenQC.QTLGeneBurdenQC
    File QTLGeneBurdenStatusList = ComputeQTLBurdenQC.QTLGeneBurdenStatusList
    File QTLBurdenOutlierEnrichment = ComputeQTLBurdenOutlierEnrichment.QTLBurdenOutlierEnrichment
    File QTLBurdenMedianGenesPerBin = CleanBurdenData.QTLBurdenMedianGenesPerBin 
    File QTLBurdenMedianGenesPerBinByGeneSet = CleanBurdenData.QTLBurdenMedianGenesPerBinByGeneSet
    File QTLBurdenMedianGenesPerBinDosageNoisyFiltered = CleanBurdenData.QTLBurdenMedianGenesPerBinDosageNoisyFiltered
    File QTLBurdenMedianGenesPerBinByGeneSetDosageNoisyFiltered = CleanBurdenData.QTLBurdenMedianGenesPerBinByGeneSetDosageNoisyFiltered
    File QTLBurdenOutlierEnrichmentPermutation = ComputeQTLBurdenOutlierEnrichment.QTLBurdenOutlierEnrichmentPermutation
    File QTLBurdenPerSampleGene = CleanBurdenData.QTLBurdenPerSampleGene
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
