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
            raise ValueError(f"Cis-window BED is missing required columns: {sorted(missing)}")

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
                raise ValueError(f"Cis window end is before start for {gene_id}: {start}-{end}")
            if gene_id in windows_by_gene:
                raise ValueError(f"Multiple cis-window rows found for {gene_id}")
            windows_by_gene[gene_id] = (
                raw[column_index["chr"]],
                start,
                end,
                safe_name(gene_id),
            )

    missing_windows = sorted(genes_with_qtls - set(windows_by_gene))
    if missing_windows:
        preview = ", ".join(missing_windows[:10])
        raise ValueError(
            f"{len(missing_windows)} gene IDs in the QTL file are missing cis-window rows. "
            f"First missing gene IDs: {preview}"
        )

    with open("gene_manifest.tsv", "w") as out:
        for gene_id in sorted(windows_by_gene):
            chrom, start, end, safe_gene_id = windows_by_gene[gene_id]
            out.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t{safe_gene_id}\n")
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
  }
}

task prepare_gene_afc_inputs {
  input {
    File vcf_file
    File vcf_index
    File afc_qtl_file
    String gene_record

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt
  }

  Int actual_disk = ceil(size(vcf_file, "GB") + size(afc_qtl_file, "GB")) + disk_space

  command <<<
    set -euo pipefail

    python3 <<'PY'
    import csv
    import gzip

    gene_record = """~{gene_record}"""
    afc_qtl_file = "~{afc_qtl_file}"
    fields = gene_record.rstrip("\n").split("\t")
    if len(fields) != 5:
        raise ValueError(f"Expected 5 tab-delimited gene manifest fields, got {len(fields)}: {gene_record!r}")

    gene_chr, gene_start, gene_end, gene_id, safe_gene_id = fields

    def open_text(path):
        if path.endswith(".gz"):
            return gzip.open(path, "rt")
        return open(path, "r")

    with open_text(afc_qtl_file) as source, open("gene.qtl.tsv", "w", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("QTL file has no header")
        required = {"pid", "sid", "sid_chr", "sid_pos"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"QTL file is missing required columns: {sorted(missing)}")

        writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        rows_written = 0
        for row in reader:
            if row.get("pid") == gene_id:
                writer.writerow(row)
                rows_written += 1

    if rows_written == 0:
        raise ValueError(f"No QTL rows found for {gene_id}")

    with open("gene_chr.txt", "w") as out:
        out.write(gene_chr)
    with open("gene_start.txt", "w") as out:
        out.write(gene_start)
    with open("gene_end.txt", "w") as out:
        out.write(gene_end)
    with open("gene_id.txt", "w") as out:
        out.write(gene_id)
    with open("safe_gene_id.txt", "w") as out:
        out.write(safe_gene_id)
    PY

    gene_chr="$(cat gene_chr.txt)"
    gene_start="$(cat gene_start.txt)"
    gene_end="$(cat gene_end.txt)"
    gene_id="$(cat gene_id.txt)"
    region="${gene_chr}:${gene_start}-${gene_end}"

    echo "Subsetting VCF to ${gene_id} cis window ${region}..."
    bcftools view \
      --threads ~{num_threads} \
      -r "${region}" \
      -Oz \
      -o gene.vcf.gz \
      "~{vcf_file}"

    echo "Indexing per-gene VCF..."
    tabix -p vcf gene.vcf.gz
  >>>

  runtime {
    docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
    memory: "${memory}GB"
    disks: "local-disk ${actual_disk} HDD"
    cpu: num_threads
    preemptible: num_preempt
  }

  output {
    File gene_vcf = "gene.vcf.gz"
    File gene_vcf_index = "gene.vcf.gz.tbi"
    File gene_qtl_file = "gene.qtl.tsv"
    String gene_chr = read_string("gene_chr.txt")
    String gene_id = read_string("gene_id.txt")
    String safe_gene_id = read_string("safe_gene_id.txt")
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

  scatter (gene_record in build_gene_afc_manifest.gene_records) {
    call prepare_gene_afc_inputs {
      input:
        vcf_file = vcf_file,
        vcf_index = vcf_index,
        afc_qtl_file = afc_qtl_file,
        gene_record = gene_record,
        memory = memory,
        disk_space = disk_space,
        num_threads = num_threads,
        num_preempt = num_preempt
    }

    call run_afc {
      input:
        vcf_file = prepare_gene_afc_inputs.gene_vcf,
        vcf_index = prepare_gene_afc_inputs.gene_vcf_index,
        expression_bed = expression_bed,
        expression_bed_index = expression_bed_index,
        covariates_file = covariates_file,
        afc_qtl_file = prepare_gene_afc_inputs.gene_qtl_file,
        prefix = prefix,
        chr = prepare_gene_afc_inputs.gene_chr,
        gene_id = prepare_gene_afc_inputs.gene_id,
        safe_gene_id = prepare_gene_afc_inputs.safe_gene_id,
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
    Array[File] per_gene_afc_reports = run_afc.afc_report
    File final_afc_report = merge_afc_reports.merged_afc_report
  }
}
