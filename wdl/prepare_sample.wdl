version 1.0

# PrepareSample (work in progress): BAM -> barcodes -> cell-type labels.
# Notebook equivalents: ExtractBarcodes = Cell 4; MakeLabels = Cells 1 and 5;
# SubsetAndUnalign = Cells 6-7.

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

task SubsetAndUnalign {
  input {
    File bam
    File labeled_barcodes
    String sample_id
  }

  # Rough disk estimate from the notebook run: 2.4 GB BAM -> 17 GB FASTQ -> intermediate BAMs.
  # 12x input + 20 GB is a guess for larger BAMs; adjust after measuring.
  Int disk_gb = ceil(size(bam, "GB") * 12) + 20

  command <<<
    set -euo pipefail

    # keep only reads from labeled barcodes
    samtools view -@ 4 -b -D CB:~{labeled_barcodes} ~{bam} -o labeled.bam

    # BAM -> FASTQ -> unaligned BAM, keeping CB and UM tags (notebook Cell 7)
    samtools fastq -@ 4 -T CB,UM -0 labeled.fastq labeled.bam
    samtools import -0 labeled.fastq -T 'CB,UM' -o ~{sample_id}_labeled_unaligned.bam
    rm labeled.fastq

    # sanity check: unaligned reads must equal primary reads in the subset BAM
    n_unal=$(samtools view -@ 4 -c ~{sample_id}_labeled_unaligned.bam)
    n_prim=$(samtools view -@ 4 -c -F 0x900 labeled.bam)
    echo -e "unaligned_reads\t${n_unal}\nprimary_reads_in_subset\t${n_prim}" > ~{sample_id}_subset_counts.txt
    cat ~{sample_id}_subset_counts.txt
    if [ "$n_unal" != "$n_prim" ]; then
      echo "ERROR: unaligned read count ($n_unal) != primary read count ($n_prim)" >&2
      exit 1
    fi
  >>>

  output {
    File unaligned_bam = "~{sample_id}_labeled_unaligned.bam"
    File counts = "~{sample_id}_subset_counts.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
    cpu: 4
    memory: "8 GB"
    disks: "local-disk ~{disk_gb} HDD"
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

  call SubsetAndUnalign {
    input:
      bam = bam,
      labeled_barcodes = MakeLabels.labeled_barcodes,
      sample_id = sample_id
  }

  output {
    File bam_barcodes = ExtractBarcodes.barcodes
    File labels = MakeLabels.labels
    File labeled_barcodes = MakeLabels.labeled_barcodes
    File batch = MakeLabels.batch
    File unaligned_bam = SubsetAndUnalign.unaligned_bam
    File subset_counts = SubsetAndUnalign.counts
  }
}

