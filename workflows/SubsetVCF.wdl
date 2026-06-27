version 1.0

task SubsetVCFTask {
    input {
        File vcf_file
        File vcf_index
        File variant_list
        File? sample_list
        String OutputPrefix
    }

    command <<<
    set -euo pipefail

    INPUT_VCF="~{vcf_file}"
    SUBSET_VCF=~{OutputPrefix}.vcf.gz
    SAMPLE_LIST="~{sep=" " select_all([sample_list])}"
    SAMPLE_LIST_ARGS=()
    if [[ -n "${SAMPLE_LIST}" ]]; then
        SAMPLE_LIST_ARGS=(-S "${SAMPLE_LIST}")
    fi

    bcftools view \
        -R ~{variant_list} \
        "${SAMPLE_LIST_ARGS[@]}" \
        -Oz \
        -o "${SUBSET_VCF}" \
        "${INPUT_VCF}"

    bcftools index -t "${SUBSET_VCF}"
    >>>
    
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/subsetvcf:main"
        memory: "256GB"
        disks: "local-disk 500 HDD"
        cpu: "1" 
        preemptible: "0" 
    }
    
    output {
        File subset_vcf = "~{OutputPrefix}.vcf.gz"
        File subset_vcf_index = "~{OutputPrefix}.vcf.gz.tbi"
    }
}


workflow SubsetVCF {
    input {
        File vcf_file
        File vcf_index
        File variant_list
        File? sample_list
        String OutputPrefix = "subset_vcf"
    }
    
    call SubsetVCFTask {
        input:
            vcf_file = vcf_file,
            vcf_index = vcf_index,
            variant_list = variant_list,
            sample_list = sample_list,
            OutputPrefix = OutputPrefix

    }
    output {
        File vcf = SubsetVCFTask.subset_vcf
        File index = SubsetVCFTask.subset_vcf_index
    }
}
