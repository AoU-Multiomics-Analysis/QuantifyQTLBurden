version 1.0

task preprocess_expression_bed {
  input {
    File expression_bed
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
  }

  Int actual_disk = ceil(size(expression_bed, "GB") * 2.5) + disk_space

  command <<<
    set -euo pipefail

    echo "Processing Expression BED..."
    gunzip -c ~{expression_bed} | bgzip -@ ~{num_threads} -c > ~{prefix}.processed_expression.bed.gz
    tabix -p bed ~{prefix}.processed_expression.bed.gz
  >>>

  runtime {
    docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
    memory: "${memory}GB"
    disks: "local-disk ${actual_disk} HDD"
    cpu: num_threads
    preemptible: num_preempt
  }

  output {
    File processed_expression_bed = "~{prefix}.processed_expression.bed.gz"
    File processed_expression_bed_tbi = "~{prefix}.processed_expression.bed.gz.tbi"
  }
}

task annotate_vcf_ids {
  input {
    File vcf_file
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
  }

  Int actual_disk = ceil(size(vcf_file, "GB") * 2.5) + disk_space

  command <<<
    set -euo pipefail

    echo "Annotating VCF IDs..."
    bcftools annotate --threads ~{num_threads} \
      --set-id '%CHROM\:%POS\_%REF\_%FIRST_ALT' \
      ~{vcf_file} -Oz -o ~{prefix}.annotated.vcf.gz
    tabix -p vcf ~{prefix}.annotated.vcf.gz
  >>>

  runtime {
    docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
    memory: "${memory}GB"
    disks: "local-disk ${actual_disk} HDD"
    cpu: num_threads
    preemptible: num_preempt
  }

  output {
    File annotated_vcf = "~{prefix}.annotated.vcf.gz"
    File annotated_vcf_tbi = "~{prefix}.annotated.vcf.gz.tbi"
  }
}

workflow preprocess_workflow_parallel {
  input {
    File vcf_file
    File expression_bed
    String prefix

    Int memory = 16
    Int disk_space = 50
    Int num_threads = 8
    Int num_preempt = 0
  }

  call preprocess_expression_bed {
    input:
      expression_bed = expression_bed,
      prefix = prefix,
      memory = memory,
      disk_space = disk_space,
      num_threads = num_threads,
      num_preempt = num_preempt
  }

  call annotate_vcf_ids {
    input:
      vcf_file = vcf_file,
      prefix = prefix,
      memory = memory,
      disk_space = disk_space,
      num_threads = num_threads,
      num_preempt = num_preempt
  }

  output {
    File processed_expression_bed = preprocess_expression_bed.processed_expression_bed
    File processed_expression_bed_tbi = preprocess_expression_bed.processed_expression_bed_tbi

    File annotated_vcf = annotate_vcf_ids.annotated_vcf
    File annotated_vcf_tbi = annotate_vcf_ids.annotated_vcf_tbi
  }
}
