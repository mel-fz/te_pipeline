version 1.0

# PrepareSample (work in progress): BAM -> barcodes -> cell-type labels.
# Notebook equivalents: ExtractBarcodes = Cell 4; MakeLabels = Cells 1 and 5.

task ExtractBarcodes {
  input {
    File bam
    String sample_id
  }

  command <<<
    set -euo pipefail
    samtools view ~{bam} \
      | grep -o 'CB:Z:[^[:space:]]*' \
      | sed 's/^CB:Z://' \
      | sort -u > ~{sample_id}_bam_CBs.txt
  >>>

  output {
    File barcodes = "~{sample_id}_bam_CBs.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
    cpu: 2
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
}

task MakeLabels {
  input {
    File script            # scripts/make_labels.py
    File bam_cbs
    File annotation        # cell_annotations.tsv
    File valid_files       # valid_files_filtered.txt
    String sample_id
  }

  command <<<
    set -euo pipefail
    pip install --quiet pandas==2.1.4
    python3 ~{script} \
      --bam_cbs ~{bam_cbs} \
      --annotation ~{annotation} \
      --sample_id ~{sample_id} \
      --valid_files ~{valid_files} \
      --out_labels ~{sample_id}_labels.tsv \
      --out_barcodes ~{sample_id}_labeled_barcodes.txt \
      --out_batch ~{sample_id}_batch.txt
  >>>

  output {
    File labels = "~{sample_id}_labels.tsv"
    File labeled_barcodes = "~{sample_id}_labeled_barcodes.txt"
    File batch = "~{sample_id}_batch.txt"
  }

  runtime {
    docker: "python:3.11-slim"
    cpu: 1
    memory: "8 GB"
    disks: "local-disk 20 HDD"
  }
}

workflow PrepareSample {
  input {
    File bam
    String sample_id
    File annotation
    File valid_files
    File make_labels_script
  }

  call ExtractBarcodes { input: bam = bam, sample_id = sample_id }

  call MakeLabels {
    input:
      script = make_labels_script,
      bam_cbs = ExtractBarcodes.barcodes,
      annotation = annotation,
      valid_files = valid_files,
      sample_id = sample_id
  }

  output {
    File bam_barcodes = ExtractBarcodes.barcodes
    File labels = MakeLabels.labels
    File labeled_barcodes = MakeLabels.labeled_barcodes
    File batch = MakeLabels.batch
  }
}
