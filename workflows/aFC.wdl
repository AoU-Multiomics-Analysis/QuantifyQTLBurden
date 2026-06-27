version 1.0

task build_gene_afc_manifest {
  input {
    File cis_window_bed
    File afc_qtl_file

    Int memory
    Int disk_space
    Int num_preempt
  }

  Int actual_disk = ceil(size(cis_window_bed, "GB") + size(afc_qtl_file, "GB")) + disk_space

  command <<<
    set -euo pipefail

    python3 <<'PY'
    import csv
    import gzip
    import re

    cis_window_bed = "~{cis_window_bed}"
    afc_qtl_file = "~{afc_qtl_file}"

    def open_text(path):
        if path.endswith(".gz"):
            return gzip.open(path, "rt")
        return open(path, "r")

    genes_with_qtls = set()
    with open_text(afc_qtl_file) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("QTL file has no header")
        if "pid" not in reader.fieldnames:
            raise ValueError("QTL file must contain a 'pid' column")
        for row in reader:
            gene_id = row.get("pid", "")
            if not gene_id or gene_id.startswith("<"):
                continue
            genes_with_qtls.add(gene_id)

    if not genes_with_qtls:
        raise ValueError("No gene IDs were found in the QTL file")

    def safe_name(value):
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)

    windows_by_gene = {}
    chromosome_order = []
    chromosomes_seen = set()
    with open_text(cis_window_bed) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = None
        for raw in reader:
            if not raw or all(field == "" for field in raw):
                continue
            if raw[0].startswith("#") and raw[0] != "#chr":
                continue
            header = [field.lstrip("#") for field in raw]
            break

        if header is None:
            raise ValueError("Cis-window BED has no header")

        required = {"chr", "start", "end", "gene_id"}
        missing = required - set(header)
        if missing:
            raise ValueError("Cis-window BED is missing required columns: {}".format(sorted(missing)))

        column_index = {column: idx for idx, column in enumerate(header)}
        for raw in reader:
            if not raw or all(field == "" for field in raw):
                continue
            gene_id = raw[column_index["gene_id"]]
            if gene_id not in genes_with_qtls:
                continue
            try:
                start = int(raw[column_index["start"]])
                end = int(raw[column_index["end"]])
            except ValueError:
                continue
            if end < start:
                raise ValueError("Cis window end is before start for {}: {}-{}".format(gene_id, start, end))
            if gene_id in windows_by_gene:
                raise ValueError("Multiple cis-window rows found for {}".format(gene_id))
            chrom = raw[column_index["chr"]]
            if chrom not in chromosomes_seen:
                chromosome_order.append(chrom)
                chromosomes_seen.add(chrom)
            windows_by_gene[gene_id] = (
                chrom,
                start,
                end,
                safe_name(gene_id),
            )

    missing_windows = sorted(genes_with_qtls - set(windows_by_gene))
    if missing_windows:
        preview = ", ".join(missing_windows[:10])
        raise ValueError(
            "{} gene IDs in the QTL file are missing cis-window rows. "
            "First missing gene IDs: {}".format(len(missing_windows), preview)
        )

    with open("gene_manifest.tsv", "w") as out:
        for gene_id in sorted(windows_by_gene):
            chrom, start, end, safe_gene_id = windows_by_gene[gene_id]
            out.write("{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, gene_id, safe_gene_id))

    with open("chromosomes.txt", "w") as out:
        for chrom in chromosome_order:
            out.write("{}\n".format(chrom))
    PY
  >>>

  runtime {
    docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
    memory: "${memory}GB"
    disks: "local-disk ${actual_disk} HDD"
    cpu: 1
    preemptible: num_preempt
  }

  output {
    File gene_manifest = "gene_manifest.tsv"
    Array[String] gene_records = read_lines("gene_manifest.tsv")
    Array[String] chromosomes = read_lines("chromosomes.txt")
  }
}

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

task prepare_chromosome_afc_inputs {
  input {
    File chr_vcf
    File chr_vcf_index
    File afc_qtl_file
    File gene_manifest
    String chr

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
  }

  Int actual_disk = ceil(size(chr_vcf, "GB") + size(afc_qtl_file, "GB") + size(gene_manifest, "GB")) + disk_space

  command <<<
    set -euo pipefail

    mkdir -p gene_inputs

    python3 <<'PY'
    import csv
    import gzip
    import os

    chromosome = "~{chr}"
    afc_qtl_file = "~{afc_qtl_file}"
    gene_manifest = "~{gene_manifest}"

    def open_text(path):
        if path.endswith(".gz"):
            return gzip.open(path, "rt")
        return open(path, "r")

    genes = []
    with open(gene_manifest, "r") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise ValueError(
                    "Expected 5 tab-delimited gene manifest fields, got {}: {!r}".format(
                        len(fields),
                        line,
                    )
                )
            gene_chr, gene_start, gene_end, gene_id, safe_gene_id = fields
            if gene_chr == chromosome:
                genes.append(
                    {
                        "chr": gene_chr,
                        "start": gene_start,
                        "end": gene_end,
                        "gene_id": gene_id,
                        "safe_gene_id": safe_gene_id,
                    }
                )

    if not genes:
        raise ValueError("No genes found in manifest for chromosome {}".format(chromosome))

    genes.sort(key=lambda gene: gene["safe_gene_id"])
    gene_ids = set(gene["gene_id"] for gene in genes)
    qtl_rows_by_gene = dict((gene["gene_id"], []) for gene in genes)

    with open_text(afc_qtl_file) as source:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("QTL file has no header")
        required = {"pid", "sid", "sid_chr", "sid_pos"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError("QTL file is missing required columns: {}".format(sorted(missing)))

        for row in reader:
            gene_id = row.get("pid")
            if gene_id in gene_ids:
                qtl_rows_by_gene[gene_id].append(row)

    with open("gene_windows.tsv", "w") as windows, \
         open("gene_chrs.txt", "w") as gene_chrs, \
         open("gene_ids.txt", "w") as gene_ids_out, \
         open("safe_gene_ids.txt", "w") as safe_gene_ids:
        for gene in genes:
            rows = qtl_rows_by_gene[gene["gene_id"]]
            if not rows:
                raise ValueError("No QTL rows found for {}".format(gene["gene_id"]))

            qtl_path = os.path.join("gene_inputs", "{}.qtl.tsv".format(gene["safe_gene_id"]))
            with open(qtl_path, "w", newline="") as target:
                writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)

            windows.write(
                "{}\t{}\t{}\t{}\t{}\n".format(
                    gene["safe_gene_id"],
                    gene["gene_id"],
                    gene["chr"],
                    gene["start"],
                    gene["end"],
                )
            )
            gene_chrs.write("{}\n".format(gene["chr"]))
            gene_ids_out.write("{}\n".format(gene["gene_id"]))
            safe_gene_ids.write("{}\n".format(gene["safe_gene_id"]))
    PY

    while IFS=$'\t' read -r safe_gene_id gene_id gene_chr gene_start gene_end; do
      region="${gene_chr}:${gene_start}-${gene_end}"
      out_vcf="gene_inputs/${safe_gene_id}.vcf.gz"

      echo "Subsetting chromosome VCF to ${gene_id} cis window ${region}..."
      bcftools view \
        --threads ~{num_threads} \
        -r "${region}" \
        -Oz \
        -o "${out_vcf}" \
        "~{chr_vcf}"

      tabix -p vcf "${out_vcf}"
    done < gene_windows.tsv
  >>>

  runtime {
    docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
    memory: "${memory}GB"
    disks: "local-disk ${actual_disk} HDD"
    cpu: num_threads
    preemptible: num_preempt
  }

  output {
    File chromosome_gene_windows = "gene_windows.tsv"
    Array[File] gene_vcfs = glob("gene_inputs/*.vcf.gz")
    Array[File] gene_vcf_indexes = glob("gene_inputs/*.vcf.gz.tbi")
    Array[File] gene_qtl_files = glob("gene_inputs/*.qtl.tsv")
    Array[String] gene_chrs = read_lines("gene_chrs.txt")
    Array[String] gene_ids = read_lines("gene_ids.txt")
    Array[String] safe_gene_ids = read_lines("safe_gene_ids.txt")
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
    String gene_id
    String safe_gene_id

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
  }

  Int actual_disk = ceil(size(vcf_file, "GB") + size(expression_bed, "GB")) + disk_space

  command <<<
    set -euo pipefail

    echo "Running aFC.py for ~{gene_id} on ~{chr}..."
    python3 /opt/aFC/aFC.py \
      --vcf "~{vcf_file}" \
      --chr "~{chr}" \
      --pheno "~{expression_bed}" \
      --qtl "~{afc_qtl_file}" \
      --cov "~{covariates_file}" \
      --log_xform 1 \
      --log_base 2 \
      --o "~{prefix}.~{safe_gene_id}.aFC.txt"

    gzip "~{prefix}.~{safe_gene_id}.aFC.txt"
  >>>

  runtime {
    docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
    memory: "${memory}GB"
    disks: "local-disk ${actual_disk} HDD"
    cpu: num_threads
    preemptible: num_preempt
  }

  output {
    File afc_report = glob("*.aFC.txt.gz")[0]
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

workflow aFC_workflow_split_by_gene {
  input {
    File vcf_file
    File vcf_index
    File expression_bed
    File expression_bed_index
    File covariates_file
    File afc_qtl_file
    File cis_window_bed
    String prefix

    Int memory = 16
    Int disk_space = 50
    Int num_threads = 8
    Int num_preempt = 0
  }

  call build_gene_afc_manifest {
    input:
      cis_window_bed = cis_window_bed,
      afc_qtl_file = afc_qtl_file,
      memory = memory,
      disk_space = disk_space,
      num_preempt = num_preempt
  }

  scatter (chr in build_gene_afc_manifest.chromosomes) {
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

    call prepare_chromosome_afc_inputs {
      input:
        chr_vcf = split_vcf_by_chr.chr_vcf,
        chr_vcf_index = split_vcf_by_chr.chr_vcf_index,
        afc_qtl_file = afc_qtl_file,
        gene_manifest = build_gene_afc_manifest.gene_manifest,
        chr = chr,
        memory = memory,
        disk_space = disk_space,
        num_threads = num_threads,
        num_preempt = num_preempt
    }
  }

  Array[File] gene_vcfs = flatten(prepare_chromosome_afc_inputs.gene_vcfs)
  Array[File] gene_vcf_indexes = flatten(prepare_chromosome_afc_inputs.gene_vcf_indexes)
  Array[File] gene_qtl_files = flatten(prepare_chromosome_afc_inputs.gene_qtl_files)
  Array[String] gene_chrs = flatten(prepare_chromosome_afc_inputs.gene_chrs)
  Array[String] gene_ids = flatten(prepare_chromosome_afc_inputs.gene_ids)
  Array[String] safe_gene_ids = flatten(prepare_chromosome_afc_inputs.safe_gene_ids)

  scatter (gene_index in range(length(gene_vcfs))) {
    call run_afc {
      input:
        vcf_file = gene_vcfs[gene_index],
        vcf_index = gene_vcf_indexes[gene_index],
        expression_bed = expression_bed,
        expression_bed_index = expression_bed_index,
        covariates_file = covariates_file,
        afc_qtl_file = gene_qtl_files[gene_index],
        prefix = prefix,
        chr = gene_chrs[gene_index],
        gene_id = gene_ids[gene_index],
        safe_gene_id = safe_gene_ids[gene_index],
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
    File gene_manifest = build_gene_afc_manifest.gene_manifest
    Array[File] per_chr_vcfs = split_vcf_by_chr.chr_vcf
    Array[File] per_chr_vcf_indexes = split_vcf_by_chr.chr_vcf_index
    Array[File] chromosome_gene_windows = prepare_chromosome_afc_inputs.chromosome_gene_windows
    Array[File] per_gene_afc_reports = run_afc.afc_report
    File final_afc_report = merge_afc_reports.merged_afc_report
  }
}
