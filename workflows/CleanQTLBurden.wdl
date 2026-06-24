version 1.0

task CleanBurdenData {
    input {
        File MergedQTLBurden
        File AlleleFrequencies
        File? ExpressionZscores 
        File aFC
        File AncestryAssignments
        File eQTLSusie
        File GTF
        Float BurdenTailProbability = 0.9
    }
    
    command <<<
    Rscript /tmp/CleanQTLBurden.R \
        --QTLBurden ~{MergedQTLBurden} \
        --AlleleFrequencies ~{AlleleFrequencies} \
        ~{if defined(ExpressionZscores) then "--ExpressionZscores " + ExpressionZscores else ""} \
        --eQTLSusie ~{eQTLSusie} \
        --aFC ~{aFC} \
        --GTF ~{GTF} \
        --AncestryAssignments ~{AncestryAssignments} \
        --BurdenTailProbability ~{BurdenTailProbability}
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
        File QTLBurdenPerSampleGene = "QTLBurdenPerSampleGene.parquet"
    }
}

workflow CleanQTLBurdenTask {
  input {
    File MergedQTLBurden
    File AlleleFrequencies
    File? ExpressionZscores
    File aFC
    File AncestryAssignments
    File eQTLSusie
    File GTF
    Float BurdenTailProbability = 0.9
  }

  call CleanBurdenData {
    input:
      MergedQTLBurden = MergedQTLBurden,
      AlleleFrequencies = AlleleFrequencies,
      ExpressionZscores = ExpressionZscores,
      aFC = aFC,
      AncestryAssignments = AncestryAssignments,
      eQTLSusie = eQTLSusie,
      GTF = GTF,
      BurdenTailProbability = BurdenTailProbability
  }

  output {
    File QTLBurdenSummaryCleaned = CleanBurdenData.QTLBurdenSummaryCleaned
    File QTLBurdenCounts = CleanBurdenData.QTLBurdenCounts
    File QTLBurdenPerSampleGene = CleanBurdenData.QTLBurdenPerSampleGene
  }
}
