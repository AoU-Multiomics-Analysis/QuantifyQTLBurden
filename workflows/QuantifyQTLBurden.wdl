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
    File? ExpressionZscores
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
        input_cleaned_qtl_burden = CleanBurdenData.QTLBurdenSummaryCleaned,
        input_aFC_weights = aFCWeights,
        input_outlier_permutation_iterations = OutlierPermutationIterations,
        input_missingness_fail_threshold = MissingnessFailThreshold,
        input_missingness_warn_threshold = MissingnessWarnThreshold,
        input_context_missingness_fail_threshold = ContextMissingnessFailThreshold,
        input_variance_ratio_lower = VarianceRatioLower,
        input_variance_ratio_upper = VarianceRatioUpper,
        input_variance_ratio_warn_lower = VarianceRatioWarnLower,
        input_variance_ratio_warn_upper = VarianceRatioWarnUpper,
        input_tail_z_fail_threshold = TailZFailThreshold,
        input_tail_z_warn_threshold = TailZWarnThreshold,
        input_dominant_variant_warn_threshold = DominantVariantWarnThreshold,
        input_qc_plot_prefix = "QTLGeneBurdenQC"
    }
  }

  output {
    File AggregatedBurden = QuantifyCoreBurden.AggregatedQTLBurden
    File FinalBurden = select_first([
      QTLBurdenQC.CleanedBurden,
      QuantifyCoreBurden.AggregatedQTLBurden
    ])
    File? CleanedBurden = QTLBurdenQC.CleanedBurden
    File? QTLBurdenCounts = CleanBurdenData.QTLBurdenCounts
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
    File? QTLBurdenOutlierEnrichmentNoisyFilterImpact = QTLBurdenQC.QTLBurdenOutlierEnrichmentNoisyFilterImpact
    File? QTLBurdenMedianGenesPerBin_AllCalls = QTLBurdenQC.QTLBurdenMedianGenesPerBin_AllCalls
    File? QTLBurdenMedianGenesPerBin_HighConfidence = QTLBurdenQC.QTLBurdenMedianGenesPerBin_HighConfidence
    File? QTLBurdenMedianGenesPerBin_HighKORemoved = QTLBurdenQC.QTLBurdenMedianGenesPerBin_HighKORemoved
    File? QTLBurdenMedianGenesPerBinByGeneSet_AllCalls = QTLBurdenQC.QTLBurdenMedianGenesPerBinByGeneSet_AllCalls
    File? QTLBurdenMedianGenesPerBinByGeneSet_HighConfidence = QTLBurdenQC.QTLBurdenMedianGenesPerBinByGeneSet_HighConfidence
    File? QTLBurdenMedianGenesPerBinByGeneSet_HighKORemoved = QTLBurdenQC.QTLBurdenMedianGenesPerBinByGeneSet_HighKORemoved
    File? QTLBurdenMedianGenesPerBin_AllCalls_DosageNoisyFiltered = QTLBurdenQC.QTLBurdenMedianGenesPerBin_AllCalls_DosageNoisyFiltered
    File? QTLBurdenMedianGenesPerBin_HighConfidence_DosageNoisyFiltered = QTLBurdenQC.QTLBurdenMedianGenesPerBin_HighConfidence_DosageNoisyFiltered
    File? QTLBurdenMedianGenesPerBin_HighKORemoved_DosageNoisyFiltered = QTLBurdenQC.QTLBurdenMedianGenesPerBin_HighKORemoved_DosageNoisyFiltered
    File? QTLBurdenMedianGenesPerBinByGeneSet_AllCalls_DosageNoisyFiltered = QTLBurdenQC.QTLBurdenMedianGenesPerBinByGeneSet_AllCalls_DosageNoisyFiltered
    File? QTLBurdenMedianGenesPerBinByGeneSet_HighConfidence_DosageNoisyFiltered = QTLBurdenQC.QTLBurdenMedianGenesPerBinByGeneSet_HighConfidence_DosageNoisyFiltered
    File? QTLBurdenMedianGenesPerBinByGeneSet_HighKORemoved_DosageNoisyFiltered = QTLBurdenQC.QTLBurdenMedianGenesPerBinByGeneSet_HighKORemoved_DosageNoisyFiltered
    File? QTLBurdenOutlierEnrichmentPermutation = QTLBurdenQC.QTLBurdenOutlierEnrichmentPermutation
    File? QTLBurdenPerSampleGene = CleanBurdenData.QTLBurdenPerSampleGene

  }
}
