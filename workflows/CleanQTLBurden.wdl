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
    }
    
    command <<<
    Rscript /tmp/CleanQTLBurden.R \
        --QTLBurden ~{MergedQTLBurden} \
        --AlleleFrequencies ~{AlleleFrequencies} \
        --ExpressionZscores ~{ExpressionZscores} \
        --eQTLSusie ~{eQTLSusie} \
        --aFC ~{aFC} \
        --GTF ~{GTF} \
        --AncestryAssignments ~{AncestryAssignments}

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
        File QTLBurdenOutlierEnrichment = "QTLBurdenOutlierEnrichment.tsv"
        File QTLBurdenMedianGenesPerBin = "QTLBurdenMedianGenesPerBin.tsv"
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
    }
    call CleanBurdenData {
      input:
        MergedQTLBurden = MergedQTLBurden,
        AlleleFrequencies = AlleleFrequencies,
        ExpressionZscores = ExpressionZscores,
        aFC = aFCWeights,
        AncestryAssignments = AncestryAssignments,
        GTF = GTF,
        eQTLSusie = eQTLSusie
    }
    output {
        File CleanedBurden = CleanBurdenData.QTLBurdenSummaryCleaned
        File QTLBurdenCounts = CleanBurdenData.QTLBurdenCounts
        File QTLGeneBurdenQC = CleanBurdenData.QTLGeneBurdenQC
        File QTLBurdenOutlierEnrichment = CleanBurdenData.QTLBurdenOutlierEnrichment
        File QTLBurdenMedianGenesPerBin = CleanBurdenData.QTLBurdenMedianGenesPerBin 
    }


}
