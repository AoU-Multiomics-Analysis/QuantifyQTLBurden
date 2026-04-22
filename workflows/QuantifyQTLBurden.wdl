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
    import csv, json, os

    infile = "~{afc_tsv}"
    gene_col = "~{gene_column}"
    genes_per_shard = ~{genes_per_shard}
    out_prefix = "~{out_prefix}"

    os.makedirs("shards", exist_ok=True)

    # read and sort rows by gene
    with open(infile, "r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames
        rows = list(reader)

    rows.sort(key=lambda x: x[gene_col])

    shard_idx = 0
    genes_seen = 0
    current_gene = None
    out = None
    writer = None
    shard_paths = []

    def open_shard(idx):
        path = f"shards/{out_prefix}.shard_{idx:04d}.tsv"
        fh = open(path, "w")
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        writer.writeheader()
        return path, fh, writer

    for row in rows:
        gene = row[gene_col]

        if current_gene != gene:
            current_gene = gene
            genes_seen += 1

            if genes_seen > genes_per_shard:
                out.close()
                shard_idx += 1
                genes_seen = 1

        if genes_seen == 1 and writer is None:
            shard_path, out, writer = open_shard(shard_idx)
            shard_paths.append(shard_path)

        writer.writerow(row)

    if out:
        out.close()

    with open("shard_manifest.json", "w") as f:
        json.dump(shard_paths, f)

    PY
    >>>
  output {
    Array[File] shard_files = read_json("shard_manifest.json")
  }

  runtime {
    docker: "python:3.10"
    memory: "96G"
    cpu: 2
    disks: "local-disk 2500 SSD"

  }
}


task QuantifyQTLBurden {
    input {
        File aFCWeights
        File VCF
        File IndexVCF
    }

    String shard_base = basename(aFCWeights, ".tsv")
    command <<<
    Rscript /tmp/QTLBurden.R \
        --AllelicFoldChangeData ~{aFCWeights} \
        --GenotypeData ~{VCF} \
        --OutputPrefix ~{shard_base}
    >>>
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden:main"
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
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden:main"
        memory: "32G"
        cpu: 2
        disks: "local-disk 2500 SSD"
    }
    

    output {
        File QTLBurdenSummary = "QTLBurdenSummary.AllGenes.tsv.gz"
    }
}



task CleanBurdenData {
    input {
        File MergedQTLBurden
        File AlleleFrequencies
        File ExpressionZscores 
        File aFC 
        File AncestryAssignments
    }
    
    command <<<
    Rscript /tmp/CleanBurdenData.R \
        --QTLBurden ~{MergedQTLBurden} \
        --AlleleFrequencies ~{AlleleFrequencies} \
        --ExpressionZscores ~{ExpressionZscores} \
        --aFC ~{aFC} \
        --AncestryAssignments ~{AncestryAssignments}

    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/quantifyqtlburden:main"
        memory: "96G"
        cpu: 2
        disks: "local-disk 2500 SSD"
    }
 
    output {
        File QTLBurdenSummaryCleaned = "QTLBurdenSummary.cleaned.tsv.gz"
    }
}

workflow qtl_burden_workflow {

  input {
    File aFCWeights
    File VCF 
    File IndexVCF
    File AlleleFrequencies
    File ExpressionZscores 
    File aFC 
    File AncestryAssignments
  }

  call shard_afc_by_gene {
    input:
      afc_tsv = aFCWeights,
      gene_column = "pid",
      genes_per_shard = 500,
      out_prefix = "afc"
  }

  scatter (shard in shard_afc_by_gene.shard_files) {
    call QuantifyQTLBurden {
      input:
        aFCWeights = shard,
        VCF = VCF,
        IndexVCF = IndexVCF
        }
    }
    call AggregateQTLBurden {
        input:
            shard_outputs = QuantifyQTLBurden.ShardBurden
    }
    
    call CleanBurdenData {
        input:
            MergedQTLBurden = AggregateQTLBurden.QTLBurdenSummary,
            AlleleFrequencies = AlleleFrequencies,
            ExpressionZscores = ExpressionZscores,
            aFC = aFC,
            AncestryAssignments = AncestryAssignments
    } 

    output {
        File MergedBurden = CleanBurdenData.QTLBurdenSummaryCleaned 
    }
}


