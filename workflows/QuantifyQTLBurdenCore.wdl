version 1.0

task shard_afc_by_gene {
  input {
    File afc_tsv
    String gene_column
    Int genes_per_shard
    String out_prefix
  }

    command <<<
      python3 <<PY
    import csv
    import json
    import os
    import gzip

    infile = "~{afc_tsv}"
    gene_col = "~{gene_column}"
    genes_per_shard = ~{genes_per_shard}
    out_prefix = "~{out_prefix}"

    os.makedirs("shards", exist_ok=True)

    # Open plain TSV or TSV.GZ
    if infile.endswith(".gz"):
        fh = gzip.open(infile, "rt")
    else:
        fh = open(infile, "r")

    # Read all rows so we can sort by gene ID
    reader = csv.DictReader(fh, delimiter="\t")
    header = reader.fieldnames
    if header is None:
        raise ValueError("Input file has no header")

    if gene_col not in header:
        raise ValueError(f"Gene column '{gene_col}' not found in header: {header}")

    rows = list(reader)
    fh.close()

    # ABSOLUTELY CRITICAL: sort by gene ID before sharding
    rows.sort(key=lambda x: x[gene_col])

    shard_idx = 0
    genes_in_current_shard = 0
    current_gene = None
    out = None
    writer = None
    shard_paths = []

    def open_shard(idx):
        path = f"shards/{out_prefix}.shard_{idx:04d}.tsv"
        handle = open(path, "w", newline="")
        w = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        w.writeheader()
        return path, handle, w

    for row in rows:
        gene = row[gene_col]

        # First row overall
        if writer is None:
            shard_path, out, writer = open_shard(shard_idx)
            shard_paths.append(shard_path)
            current_gene = gene
            genes_in_current_shard = 1

        # New gene encountered
        elif gene != current_gene:
            current_gene = gene

            # If current shard already has enough genes, start a new one
            if genes_in_current_shard >= genes_per_shard:
                out.close()
                shard_idx += 1
                shard_path, out, writer = open_shard(shard_idx)
                shard_paths.append(shard_path)
                genes_in_current_shard = 1
            else:
                genes_in_current_shard += 1

        writer.writerow(row)

    if out is not None:
        out.close()

    with open("shard_manifest.json", "w") as f:
        json.dump(shard_paths, f)
    PY
    >>>  
    output {
      Array[File] shard_files = glob("shards/*.tsv")
      File manifest = "shard_manifest.json"
    }  
    
    runtime {
        docker: "python:3.10"
        memory: "64G"
        cpu: 2
        disks: "local-disk 2500 SSD"
  }
}


task QuantifyQTLBurden {
    input {
        File aFCWeights
        File VCF
        File IndexVCF
        Float LossThreshold = -0.5849625007
        Float GainThreshold = 0.5849625007
        String VariantEffectSEColumn = "auto"
    }

    String shard_base = basename(aFCWeights, ".tsv")
    command <<<
    Rscript /tmp/QTLBurden.R \
        --AllelicFoldChangeData ~{aFCWeights} \
        --VariantEffectSEColumn ~{VariantEffectSEColumn} \
        --LossThreshold ~{LossThreshold} \
        --GainThreshold ~{GainThreshold} \
        --GenotypeData ~{VCF} \
        --OutputPrefix ~{shard_base}
    >>>
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/quantifyqtlburden:main"
        memory: "32G"
        cpu: 2
        disks: "local-disk 2500 SSD"
    }
    
    output {
        File ShardBurden = shard_base + ".QTLBurdenSummary.tsv.gz"
    }
}

task AggregateQTLBurden {

    input {
        Array[File] shard_outputs
    }

    command <<<
    set -euo pipefail

    first=1

    for f in ~{sep=' ' shard_outputs}; do
        if [ $first -eq 1 ]; then
            zcat "$f"
            first=0
        else
            zcat "$f" | tail -n +2
        fi
    done | gzip > QTLBurdenSummary.AllGenes.tsv.gz
    >>>
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden/quantifyqtlburden:main"
        memory: "32G"
        cpu: 2
        disks: "local-disk 2500 SSD"
    }
    

    output {
        File QTLBurdenSummary = "QTLBurdenSummary.AllGenes.tsv.gz"
    }
}



workflow QuantifyQTLBurdenCore {

  input {
    File aFCWeights
    File VCF
    File IndexVCF
    Int GenesPerShard = 500
    Float LossThreshold = -0.5849625007
    Float GainThreshold = 0.5849625007
    String VariantEffectSEColumn = "auto"
  }

  call shard_afc_by_gene {
    input:
      afc_tsv = aFCWeights,
      gene_column = "pid",
      genes_per_shard = GenesPerShard,
      out_prefix = "afc"
  }

  scatter (shard in shard_afc_by_gene.shard_files) {
    call QuantifyQTLBurden {
      input:
        aFCWeights = shard,
        VCF = VCF,
        IndexVCF = IndexVCF,
        LossThreshold = LossThreshold,
        GainThreshold = GainThreshold,
        VariantEffectSEColumn = VariantEffectSEColumn
    }
  }

  call AggregateQTLBurden {
    input:
      shard_outputs = QuantifyQTLBurden.ShardBurden
  }

  output {
    File AggregatedQTLBurden = AggregateQTLBurden.QTLBurdenSummary
  }


}
