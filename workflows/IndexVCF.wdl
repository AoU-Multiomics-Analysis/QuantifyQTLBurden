version 1.0

task IndexVCFTask {
    input {
        File vcf_file
    }

    command <<<
    set -euo pipefail

    bcftools index -t ~{vcf_file}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/subsetvcf:main"
        memory: "128GB"
        cpu: "1"
        disks: "local-disk 200 HDD"
        preemptible: "0"
    }

    output {
        File vcf_tbi = "~{vcf_file}.tbi"
    }
}

workflow IndexVCF {
    input {
        File vcf_file
    }

    call IndexVCFTask {
        input:
            vcf_file = vcf_file
    }

    output {
        File tbi = IndexVCFTask.vcf_tbi
    }
}
