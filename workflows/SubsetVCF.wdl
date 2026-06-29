version 1.0

task SubsetVCFTask {
    input {
        File vcf_file
        File vcf_index
        File variant_list
        File? sample_list
        String OutputPrefix
        Int num_threads
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
        -Ou \
        "${INPUT_VCF}" | \
    bcftools annotate --threads ~{num_threads} \
        --set-id '%CHROM\:%POS\_%REF\_%FIRST_ALT' \
        -Oz \
        -o "${SUBSET_VCF}" \
        -

    bcftools index -t "${SUBSET_VCF}"
    >>>
    
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/subsetvcf:main"
        memory: "256GB"
        disks: "local-disk 500 HDD"
        cpu: num_threads
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
        Int num_threads = 1
    }
    
    call SubsetVCFTask {
        input:
            vcf_file = vcf_file,
            vcf_index = vcf_index,
            variant_list = variant_list,
            sample_list = sample_list,
            OutputPrefix = OutputPrefix,
            num_threads = num_threads

    }
    output {
        File vcf = SubsetVCFTask.subset_vcf
        File index = SubsetVCFTask.subset_vcf_index
    }
}
