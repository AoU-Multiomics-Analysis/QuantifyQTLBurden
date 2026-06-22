version 1.0 

import "QuantifyQTLBurdenCore.wdl" as core
import "QTLBurdenQC.wdl" as qctools
import "CleanQTLBurden.wdl" as clean

  workflow qtl_burden_workflow {

  input {
    File aFCWeights
    File VCF
    File IndexVCF
    File AlleleFrequencies
    File ExpressionZscores
    File AncestryAssignments
    Int GenesPerShard = 500
    Boolean AnnotateBurden = true
    File GTF 
    File eQTLSusie
    Float LossThreshold = -0.5849625007
    Float GainThreshold = 0.5849625007
    String BurdenDirection = "both"
    String VariantEffectSEColumn = "auto"
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

  call core.QuantifyQTLBurdenCore as QuantifyCoreBurden {
    input:
      aFCWeights = aFCWeights,
      VCF = VCF,
      IndexVCF = IndexVCF,
      GenesPerShard = GenesPerShard,
      LossThreshold = LossThreshold,
      GainThreshold = GainThreshold,
      BurdenDirection = BurdenDirection,
      VariantEffectSEColumn = VariantEffectSEColumn
  }

  if (AnnotateBurden) {
    call clean.CleanBurdenData as CleanBurdenData {
      input:
        MergedQTLBurden = QuantifyCoreBurden.AggregatedQTLBurden,
        AlleleFrequencies = AlleleFrequencies,
        ExpressionZscores = ExpressionZscores,
        aFC = aFCWeights,
        AncestryAssignments = AncestryAssignments,
        eQTLSusie = eQTLSusie,
        GTF = GTF,
        BurdenTailProbability = BurdenTailProbability
    }

    call qctools.QTLBurdenQC as QTLBurdenQC {
      input:
        CleanedQTLBurden = CleanBurdenData.QTLBurdenSummaryCleaned,
        aFCWeights = aFCWeights,
        QTLBurdenCounts = CleanBurdenData.QTLBurdenCounts,
        QTLBurdenPerSampleGene = CleanBurdenData.QTLBurdenPerSampleGene,
        OutlierPermutationIterations = OutlierPermutationIterations,
        MissingnessFailThreshold = MissingnessFailThreshold,
        MissingnessWarnThreshold = MissingnessWarnThreshold,
        ContextMissingnessFailThreshold = ContextMissingnessFailThreshold,
        VarianceRatioLower = VarianceRatioLower,
        VarianceRatioUpper = VarianceRatioUpper,
        VarianceRatioWarnLower = VarianceRatioWarnLower,
        VarianceRatioWarnUpper = VarianceRatioWarnUpper,
        TailZFailThreshold = TailZFailThreshold,
        TailZWarnThreshold = TailZWarnThreshold,
        DominantVariantWarnThreshold = DominantVariantWarnThreshold,
        QCPlotPrefix = "QTLGeneBurdenQC"
    }
  }

  output {
    File AggregatedBurden = QuantifyCoreBurden.AggregatedQTLBurden
    File FinalBurden = select_first([
      QTLBurdenQC.CleanedBurden,
      QuantifyCoreBurden.AggregatedQTLBurden
    ])
    File? CleanedBurden = QTLBurdenQC.CleanedBurden
    File? QTLBurdenCounts = QTLBurdenQC.QTLBurdenCounts
    File? QTLGeneBurdenQC = QTLBurdenQC.QTLGeneBurdenQC
    File? QTLGeneBurdenStatusList = QTLBurdenQC.QTLGeneBurdenStatusList
    File? QTLGeneBurdenQC_StatusByGeneType = QTLBurdenQC.QTLGeneBurdenQC_StatusByGeneType
    File? QTLGeneBurdenQC_StatusOverall = QTLBurdenQC.QTLGeneBurdenQC_StatusOverall
    File? QTLGeneBurdenQC_Missingness = QTLBurdenQC.QTLGeneBurdenQC_Missingness
    File? QTLGeneBurdenQC_VarianceRatio = QTLBurdenQC.QTLGeneBurdenQC_VarianceRatio
    File? QTLGeneBurdenQC_TailZ = QTLBurdenQC.QTLGeneBurdenQC_TailZ
    File? QTLGeneBurdenQC_DominantVariantFraction = QTLBurdenQC.QTLGeneBurdenQC_DominantVariantFraction
    File? QTLGeneBurdenQC_ReasonCounts = QTLBurdenQC.QTLGeneBurdenQC_ReasonCounts
    File? QTLGeneBurdenQC_GeneSetMedianImpact = QTLBurdenQC.QTLGeneBurdenQC_GeneSetMedianImpact
    File? QTLBurdenOutlierEnrichment = QTLBurdenQC.QTLBurdenOutlierEnrichment
    File? QTLBurdenMedianGenesPerBin = QTLBurdenQC.QTLBurdenMedianGenesPerBin 
    File? QTLBurdenMedianGenesPerBinByGeneSet = QTLBurdenQC.QTLBurdenMedianGenesPerBinByGeneSet
    File? QTLBurdenMedianGenesPerBinDosageNoisyFiltered = QTLBurdenQC.QTLBurdenMedianGenesPerBinDosageNoisyFiltered
    File? QTLBurdenMedianGenesPerBinByGeneSetDosageNoisyFiltered = QTLBurdenQC.QTLBurdenMedianGenesPerBinByGeneSetDosageNoisyFiltered
    File? QTLBurdenOutlierEnrichmentPermutation = QTLBurdenQC.QTLBurdenOutlierEnrichmentPermutation
    File? QTLBurdenPerSampleGene = QTLBurdenQC.QTLBurdenPerSampleGene

  }
}
