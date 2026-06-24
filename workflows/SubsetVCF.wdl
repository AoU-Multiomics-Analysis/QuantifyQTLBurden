version 1.0

task SubsetVCFTask {
    input {
        File vcf_file
        File vcf_index
        File variant_list
        File sample_list
        String OutputPrefix
    }

    command <<<
    bcftools view \
        -R ~{variant_list} \
        -S  ~{sample_list} \
        -Oz \
        -o ~{OutputPrefix}.vcf.gz \
        ~{vcf_file} 

    bcftools index ~{OutputPrefix}.vcf.gz
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
        File sample_list
    }
    
    call SubsetVCFTask {
        input:
            vcf_file = vcf_file,
            vcf_index = vcf_index,
            variant_list = variant_list,
            sample_list = sample_list

    }
    output {
        File vcf = SubsetVCFTask.subset_vcf
        File index = SubsetVCFTask.subset_vcf_index
    }
}
