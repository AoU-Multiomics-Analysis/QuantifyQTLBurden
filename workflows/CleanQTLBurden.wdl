version 1.0

task CleanBurdenData {
    input {
        File MergedQTLBurden
        File AlleleFrequencies
        File ExpressionZscores 
        File aFC
        File AncestryAssignments
        File eQTLSusie
        File GTF
        Float BurdenTailProbability = 0.9
        Int OutlierPermutationIterations = 200
    }
    
    command <<<
    Rscript /tmp/CleanQTLBurden.R \
        --QTLBurden ~{MergedQTLBurden} \
        --AlleleFrequencies ~{AlleleFrequencies} \
        --ExpressionZscores ~{ExpressionZscores} \
        --eQTLSusie ~{eQTLSusie} \
        --aFC ~{aFC} \
        --GTF ~{GTF} \
        --AncestryAssignments ~{AncestryAssignments} \
        --BurdenTailProbability ~{BurdenTailProbability} \
        --OutlierPermutationIterations ~{OutlierPermutationIterations}

    if [ -f /tmp/PlotQTLBurdenQC.R ]; then
      Rscript /tmp/PlotQTLBurdenQC.R \
        --QCFile QTLGeneBurdenQC.tsv.gz \
        --MedianGenesPerBinByGeneSet QTLBurdenMedianGenesPerBinByGeneSet.tsv \
        --Prefix QTLGeneBurdenQC
    else
      echo "PlotQTLBurdenQC.R not found in /tmp; creating placeholder outputs."
      touch QTLGeneBurdenQC_StatusByGeneType.pdf
      touch QTLGeneBurdenQC_StatusOverall.pdf
      touch QTLGeneBurdenQC_Missingness.pdf
      touch QTLGeneBurdenQC_VarianceRatio.pdf
      touch QTLGeneBurdenQC_TailZ.pdf
      touch QTLGeneBurdenQC_DominantVariantFraction.pdf
      touch QTLGeneBurdenQC_ReasonCounts.pdf
      touch QTLGeneBurdenQC_GeneSetMedianImpact.pdf
    fi

    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/cleanqtlburden:main"
        memory: "512G"
        cpu: 2
        disks: "local-disk 2500 SSD"
    }
 
    output {
        File QTLBurdenSummaryCleaned = "QTLBurdenSummary.cleaned.tsv.gz"
        File QTLBurdenCounts = "QTLGeneBurdenCounts.tsv.gz"
        File QTLGeneBurdenQC = "QTLGeneBurdenQC.tsv.gz"
        File QTLGeneBurdenStatusList = "QTLGeneBurdenStatusList.tsv.gz"
        File QTLGeneBurdenQC_StatusByGeneType = "QTLGeneBurdenQC_StatusByGeneType.pdf"
        File QTLGeneBurdenQC_StatusOverall = "QTLGeneBurdenQC_StatusOverall.pdf"
        File QTLGeneBurdenQC_Missingness = "QTLGeneBurdenQC_Missingness.pdf"
        File QTLGeneBurdenQC_VarianceRatio = "QTLGeneBurdenQC_VarianceRatio.pdf"
        File QTLGeneBurdenQC_TailZ = "QTLGeneBurdenQC_TailZ.pdf"
        File QTLGeneBurdenQC_DominantVariantFraction = "QTLGeneBurdenQC_DominantVariantFraction.pdf"
        File QTLGeneBurdenQC_ReasonCounts = "QTLGeneBurdenQC_ReasonCounts.pdf"
        File QTLGeneBurdenQC_GeneSetMedianImpact = "QTLGeneBurdenQC_GeneSetMedianImpact.pdf"
        File QTLBurdenOutlierEnrichment = "QTLBurdenOutlierEnrichment.tsv"
        File QTLBurdenMedianGenesPerBin = "QTLBurdenMedianGenesPerBin.tsv"
        File QTLBurdenMedianGenesPerBinByGeneSet = "QTLBurdenMedianGenesPerBinByGeneSet.tsv"
        File QTLBurdenMedianGenesPerBinDosageNoisyFiltered = "QTLBurdenMedianGenesPerBin_DosageNoisyFiltered.tsv"
        File QTLBurdenMedianGenesPerBinByGeneSetDosageNoisyFiltered = "QTLBurdenMedianGenesPerBinByGeneSet_DosageNoisyFiltered.tsv"
        File QTLBurdenOutlierEnrichmentPermutation = "QTLBurdenOutlierEnrichmentPermutation.tsv"
        File QTLBurdenPerSampleGene = "QTLBurdenPerSampleGene.parquet"
    }
}

workflow CleanQTLBurden {
    input {
        File AlleleFrequencies
        File ExpressionZscores 
        File aFCWeights 
        File AncestryAssignments
        File MergedQTLBurden
        File eQTLSusie
        File GTF
        Float BurdenTailProbability = 0.9
        Int OutlierPermutationIterations = 200
    }
    call CleanBurdenData {
      input:
        MergedQTLBurden = MergedQTLBurden,
        AlleleFrequencies = AlleleFrequencies,
        ExpressionZscores = ExpressionZscores,
        aFC = aFCWeights,
        AncestryAssignments = AncestryAssignments,
        GTF = GTF,
        eQTLSusie = eQTLSusie,
        BurdenTailProbability = BurdenTailProbability,
        OutlierPermutationIterations = OutlierPermutationIterations
    }
    output {
        File CleanedBurden = CleanBurdenData.QTLBurdenSummaryCleaned
        File QTLBurdenCounts = CleanBurdenData.QTLBurdenCounts
        File QTLGeneBurdenQC = CleanBurdenData.QTLGeneBurdenQC
        File QTLGeneBurdenStatusList = CleanBurdenData.QTLGeneBurdenStatusList
        File QTLGeneBurdenQC_StatusByGeneType = CleanBurdenData.QTLGeneBurdenQC_StatusByGeneType
        File QTLGeneBurdenQC_StatusOverall = CleanBurdenData.QTLGeneBurdenQC_StatusOverall
        File QTLGeneBurdenQC_Missingness = CleanBurdenData.QTLGeneBurdenQC_Missingness
        File QTLGeneBurdenQC_VarianceRatio = CleanBurdenData.QTLGeneBurdenQC_VarianceRatio
        File QTLGeneBurdenQC_TailZ = CleanBurdenData.QTLGeneBurdenQC_TailZ
        File QTLGeneBurdenQC_DominantVariantFraction = CleanBurdenData.QTLGeneBurdenQC_DominantVariantFraction
        File QTLGeneBurdenQC_ReasonCounts = CleanBurdenData.QTLGeneBurdenQC_ReasonCounts
        File QTLGeneBurdenQC_GeneSetMedianImpact = CleanBurdenData.QTLGeneBurdenQC_GeneSetMedianImpact
        File QTLBurdenOutlierEnrichment = CleanBurdenData.QTLBurdenOutlierEnrichment
        File QTLBurdenMedianGenesPerBin = CleanBurdenData.QTLBurdenMedianGenesPerBin 
        File QTLBurdenMedianGenesPerBinByGeneSet = CleanBurdenData.QTLBurdenMedianGenesPerBinByGeneSet
        File QTLBurdenMedianGenesPerBinDosageNoisyFiltered = CleanBurdenData.QTLBurdenMedianGenesPerBinDosageNoisyFiltered
        File QTLBurdenMedianGenesPerBinByGeneSetDosageNoisyFiltered = CleanBurdenData.QTLBurdenMedianGenesPerBinByGeneSetDosageNoisyFiltered
        File QTLBurdenOutlierEnrichmentPermutation = CleanBurdenData.QTLBurdenOutlierEnrichmentPermutation
        File QTLBurdenPerSampleGene = CleanBurdenData.QTLBurdenPerSampleGene
    }


}
