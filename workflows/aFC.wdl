version 1.0

task split_vcf_by_chr {
  input {
    File vcf_file
    File vcf_index
    String chr
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
  }

  Int actual_disk = ceil(size(vcf_file, "GB")) + disk_space

  command <<<
    set -euo pipefail

    echo "Subsetting VCF to ~{chr}..."
    out_vcf="~{prefix}.~{chr}.vcf.gz"

    bcftools view \
      --threads ~{num_threads} \
      -r "~{chr}" \
      -Oz \
      -o "${out_vcf}" \
      "~{vcf_file}"

    echo "Indexing per-chromosome VCF..."
    tabix -p vcf "${out_vcf}"
  >>>

  runtime {
    docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
    memory: "${memory}GB"
    disks: "local-disk ${actual_disk} HDD"
    cpu: num_threads
    preemptible: num_preempt
  }

  output {
    File chr_vcf = "${prefix}.${chr}.vcf.gz"
    File chr_vcf_index = "${prefix}.${chr}.vcf.gz.tbi"
  }
}

task run_afc {
  input {
    File vcf_file
    File vcf_index
    File expression_bed
    File expression_bed_index
    File covariates_file
    File afc_qtl_file
    String prefix
    String chr

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
  }

  Int actual_disk = ceil(size(vcf_file, "GB") + size(expression_bed, "GB")) + disk_space

  command <<<
    set -euo pipefail

    echo "Running aFC.py for ~{chr}..."
    python3 /opt/aFC/aFC.py \
      --vcf "~{vcf_file}" \
      --chr "~{chr}" \
      --pheno "~{expression_bed}" \
      --qtl "~{afc_qtl_file}" \
      --cov "~{covariates_file}" \
      --log_xform 1 \
      --log_base 2 \
      --o "~{prefix}.~{chr}.aFC.txt"

    gzip "~{prefix}.~{chr}.aFC.txt"
  >>>

  runtime {
    docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
    memory: "${memory}GB"
    disks: "local-disk ${actual_disk} HDD"
    cpu: num_threads
    preemptible: num_preempt
  }

  output {
    File afc_report = "${prefix}.${chr}.aFC.txt.gz"
  }
}

task merge_afc_reports {
  input {
    Array[File] afc_reports
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
  }

  command <<<
    set -euo pipefail

    out="~{prefix}.aFC.txt"
    rm -f "${out}"

    first=true
    for f in ~{sep=' ' afc_reports}; do
      if [ "${first}" = "true" ]; then
        gzip -cd "$f" >> "${out}"
        first=false
      else
        gzip -cd "$f" | tail -n +2 >> "${out}"
      fi
    done

    gzip "${out}"
  >>>

  runtime {
    docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
    memory: "${memory}GB"
    disks: "local-disk ${disk_space} HDD"
    cpu: num_threads
    preemptible: num_preempt
  }

  output {
    File merged_afc_report = "${prefix}.aFC.txt.gz"
  }
}

workflow aFC_workflow_split_by_chr {
  input {
    File vcf_file
    File vcf_index
    File expression_bed
    File expression_bed_index
    File covariates_file
    File afc_qtl_file
    String prefix

    # Optional override. If not provided, defaults to chr1..chr22, chrX, chrY.
    Array[String]? chromosomes

    Int memory = 16
    Int disk_space = 50
    Int num_threads = 8
    Int num_preempt = 0
  }

  # Hard-coded default chromosome list (chr1 convention), used unless overridden.
  Array[String] default_chromosomes = [
    "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
    "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
    "chr21","chr22","chrX","chrY"
  ]

  Array[String] chr_list = select_first([chromosomes, default_chromosomes])

  scatter (chr in chr_list) {
    call split_vcf_by_chr {
      input:
        vcf_file = vcf_file,
        vcf_index = vcf_index,
        chr = chr,
        prefix = prefix,
        memory = memory,
        disk_space = disk_space,
        num_threads = num_threads,
        num_preempt = num_preempt
    }

    call run_afc {
      input:
        vcf_file = split_vcf_by_chr.chr_vcf,
        vcf_index = split_vcf_by_chr.chr_vcf_index,
        expression_bed = expression_bed,
        expression_bed_index = expression_bed_index,
        covariates_file = covariates_file,
        afc_qtl_file = afc_qtl_file,
        prefix = prefix,
        chr = chr,
        memory = memory,
        disk_space = disk_space,
        num_threads = num_threads,
        num_preempt = num_preempt
    }
  }

  call merge_afc_reports {
    input:
      afc_reports = run_afc.afc_report,
      prefix = prefix,
      memory = memory,
      disk_space = disk_space,
      num_threads = num_threads,
      num_preempt = num_preempt
  }

  output {
    Array[File] per_chr_afc_reports = run_afc.afc_report
    File final_afc_report = merge_afc_reports.merged_afc_report
  }
}
