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
        File QTLBurdenOutlierEnrichment = "QTLBurdenOutlierEnrichment.tsv"
        File QTLBurdenMedianGenesPerBin = "QTLBurdenMedianGenesPerBin.tsv"
        File QTLBurdenMedianGenesPerBinByGeneSet = "QTLBurdenMedianGenesPerBinByGeneSet.tsv"
        File QTLBurdenMedianGenesPerBinDosageNoisyFiltered = "QTLBurdenMedianGenesPerBin_DosageNoisyFiltered.tsv"
        File QTLBurdenMedianGenesPerBinByGeneSetDosageNoisyFiltered = "QTLBurdenMedianGenesPerBinByGeneSet_DosageNoisyFiltered.tsv"
        File QTLBurdenOutlierEnrichmentPermutation = "QTLBurdenOutlierEnrichmentPermutation.tsv"
        File QTLBurdenPerSampleGene = "QTLBurdenPerSampleGene.parquet"
    }
}
