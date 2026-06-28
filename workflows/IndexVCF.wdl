version 1.0

task IndexVCFTask {
    input {
        File vcf_file
        Int memory
    }

    command <<<
    set -euo pipefail

    bcftools index --tbi --force \
        --output "~{basename(vcf_file)}.tbi" \
        "~{vcf_file}"
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/subsetvcf:main"
        memory: "${memory}GB"
        cpu: "1"
        disks: "local-disk 1000 HDD"
        preemptible: "0"
    }

    output {
        File vcf_tbi = "~{basename(vcf_file)}.tbi"
    }
}

workflow IndexVCF {
    input {
        File vcf_file
        Int memory = 128
    }

    call IndexVCFTask {
        input:
            vcf_file = vcf_file,
            memory = memory
    }

    output {
        File tbi = IndexVCFTask.vcf_tbi
    }
}
