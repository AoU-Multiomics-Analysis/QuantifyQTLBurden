version 1.0

task CleanBurdenData {
    input {
        File MergedQTLBurden
        File AlleleFrequencies
        File ExpressionZscores 
        File aFC 
        File AncestryAssignments
    }
    
    command <<<
    Rscript /tmp/CleanQTLBurden.R \
        --QTLBurden ~{MergedQTLBurden} \
        --AlleleFrequencies ~{AlleleFrequencies} \
        --ExpressionZscores ~{ExpressionZscores} \
        --aFC ~{aFC} \
        --AncestryAssignments ~{AncestryAssignments}

    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/cleanqtlburden:main"
        memory: "256G"
        cpu: 2
        disks: "local-disk 2500 SSD"
    }
 
    output {
        File QTLBurdenSummaryCleaned = "QTLBurdenSummary.cleaned.tsv.gz"
        File QTLBurdenCounts = "QTLGeneBurdenCounts.tsv.gz"
    }
}

workflow CleanQTLBurden {
    input {
        File MergedQTLBurden
        File AlleleFrequencies
        File ExpressionZscores 
        File aFCWeights 
        File AncestryAssignments
        File AggregatedQTLBurden
    }
    call CleanBurdenData {
      input:
        MergedQTLBurden = AggregatedQTLBurden,
        AlleleFrequencies = AlleleFrequencies,
        ExpressionZscores = ExpressionZscores,
        aFC = aFCWeights,
        AncestryAssignments = AncestryAssignments
    }
    output {
        File CleanedBurden = CleanBurdenData.QTLBurdenSummaryCleaned
        File QTLBurdenCounts = CleanBurdenData.QTLBurdenCounts
    }


}
