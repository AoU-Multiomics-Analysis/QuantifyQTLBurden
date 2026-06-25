version 1.0

task SubsetVCFTask {
    input {
        File vcf_file
        File? vcf_index
        File variant_list
        File sample_list
        String OutputPrefix
    }

    command <<<
    set -euo pipefail

    INPUT_VCF=~{vcf_file}
    INPUT_VCF_INDEX=${INPUT_VCF}.tbi
    TASK_INPUT_INDEX=~{OutputPrefix}.input.vcf.gz.tbi
    SUBSET_VCF=~{OutputPrefix}.vcf.gz

    if [[ "~{vcf_index}" == "null" ]] || [[ "~{vcf_index}" == "" ]] || [[ "~{vcf_index}" == "\${}" ]] || [[ "~{vcf_index}" == *.csi ]]; then
      # Index the original input VCF path.
      bcftools index -t "${INPUT_VCF}"
    else
      # Align provided index to the basename expected by bcftools.
      cp ~{vcf_index} "${INPUT_VCF_INDEX}"
    fi

    # Keep a reusable tabix path for downstream tasks.
    cp "${INPUT_VCF_INDEX}" "${TASK_INPUT_INDEX}"

    bcftools view \
        -R ~{variant_list} \
        -S  ~{sample_list} \
        -Oz \
        -o "${SUBSET_VCF}" \
        "${INPUT_VCF}"

    # Always emit a tabix index for downstream compatibility.
    bcftools index -t ~{OutputPrefix}.vcf.gz
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
        File? input_vcf_index = "~{OutputPrefix}.input.vcf.gz.tbi"
    }
}


workflow SubsetVCF {
    input {
        File vcf_file
        File? vcf_index
        File variant_list
        File sample_list
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
        File? input_index = SubsetVCFTask.input_vcf_index
    }
}
