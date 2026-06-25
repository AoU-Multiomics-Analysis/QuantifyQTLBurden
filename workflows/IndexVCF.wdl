version 1.0

task IndexVCFTask {
    input {
        File vcf_file
    }

    command <<<
    set -euo pipefail

    INPUT_VCF="~{vcf_file}"
    WORK_VCF="${INPUT_VCF}"

    # bcftools index expects common VCF bgzip suffixes; handle .vcf.bgz inputs explicitly.
    if [[ "${INPUT_VCF}" == *.vcf.bgz ]]; then
        WORK_VCF="${INPUT_VCF}.tmp.vcf.gz"
        ln -s "${INPUT_VCF}" "${WORK_VCF}"
    fi

    bcftools index -t "${WORK_VCF}"

    if [[ ! -f "${WORK_VCF}.tbi" ]]; then
        echo "failed to create index for ${INPUT_VCF}" >&2
        exit 1
    fi

    mv "${WORK_VCF}.tbi" index.tbi
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/subsetvcf:main"
        memory: "128GB"
        cpu: "1"
        disks: "local-disk 200 HDD"
        preemptible: "0"
    }

    output {
        File vcf_tbi = "index.tbi"
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
