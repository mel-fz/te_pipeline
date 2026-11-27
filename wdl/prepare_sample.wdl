version 1.0

# PrepareSample (work in progress): BAM -> barcodes -> cell-type labels.
# Sample validity (valid-files list) is checked BEFORE the run, when building the sample table.
# Notebook equivalents: ExtractBarcodes = Cell 4; MakeLabels = Cell 5;
# SubsetAndUnalign = Cells 6-7; Align = Cell 9 (STAR) and Cell 11 (checks).

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
    maxRetries: 1
  }
}

task MakeLabels {
  input {
    File script            # scripts/make_labels.py
    File bam_cbs
    File annotation        # cell_annotations.tsv
    String sample_id
    String tools_image
  }

  command <<<
    set -euo pipefail
    python3 ~{script} \
      --bam_cbs ~{bam_cbs} \
      --annotation ~{annotation} \
      --sample_id ~{sample_id} \
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
    docker: tools_image
    cpu: 1
    memory: "8 GB"
    disks: "local-disk 20 HDD"
    maxRetries: 1
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

    # keep only reads from labeled barcodes (notebook Cell 6)
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
    maxRetries: 1
  }
}

task Align {
  input {
    File unaligned_bam
    File star_index_tar        # tar of the primary-contig STAR index (folder star_index_primary)
    String sample_id
    String tools_image         # e.g. <dockerhub-user>/te-pipeline:0.1.0
    Int threads = 12
    Int expected_contigs = 25  # primary-only index; stops the run if the wrong index is supplied
  }

  # Peak = tar + unpacked index (the tar is deleted after unpacking) + unaligned BAM + STAR sort scratch (>=24 GB).
  # Rough estimate from the notebook run, not measured on larger BAMs.
  Int disk_gb = ceil(size(star_index_tar, "GB") * 2 + size(unaligned_bam, "GB") * 6) + 40

  command <<<
    set -euo pipefail

    mkdir -p index
    tar -xf ~{star_index_tar} -C index
    rm -f ~{star_index_tar} || echo "could not delete the index tar; continuing (needs more disk)"
    IDX=$(find index -name genomeParameters.txt -printf '%h\n' | head -1)
    echo "STAR index: $IDX"

    # same settings as the notebook (Cell 9)
    STAR-avx2 --runThreadN ~{threads} --genomeDir "$IDX" \
      --readFilesType SAM SE --readFilesIn ~{unaligned_bam} \
      --readFilesCommand samtools view -h --readFilesSAMattrKeep CB UM \
      --outFilterMultimapNmax 500 --outFilterMultimapScoreRange 5 \
      --limitOutSJcollapsed 5000000 \
      --outSAMattributes NH HI AS nM NM MD --outSAMunmapped Within \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix ~{sample_id}_labeled_

    # checks (notebook Cell 11)
    test -s ~{sample_id}_labeled_Log.final.out
    n=$(samtools view -H ~{sample_id}_labeled_Aligned.sortedByCoord.out.bam | grep -c '^@SQ')
    echo "contigs in header: $n (expected ~{expected_contigs})"
    if [ "$n" != "~{expected_contigs}" ]; then
      echo "ERROR: expected ~{expected_contigs} contigs, found $n" >&2
      exit 1
    fi
    grep -E "Number of input reads|Uniquely mapped reads %|% of reads mapped to multiple loci|% of reads mapped to too many loci" ~{sample_id}_labeled_Log.final.out
  >>>

  output {
    File aligned_bam = "~{sample_id}_labeled_Aligned.sortedByCoord.out.bam"
    File star_log = "~{sample_id}_labeled_Log.final.out"
  }

  runtime {
    docker: tools_image
    cpu: threads
    memory: "40 GB"
    disks: "local-disk ~{disk_gb} HDD"
    maxRetries: 1
  }
}

workflow PrepareSample {
  input {
    File bam
    String sample_id
    File annotation
    File make_labels_script
    File star_index_tar
    String tools_image
  }

  call ExtractBarcodes { input: bam = bam, sample_id = sample_id }

  call MakeLabels {
    input:
      script = make_labels_script,
      bam_cbs = ExtractBarcodes.barcodes,
      annotation = annotation,
      sample_id = sample_id,
      tools_image = tools_image
  }

  call SubsetAndUnalign {
    input:
      bam = bam,
      labeled_barcodes = MakeLabels.labeled_barcodes,
      sample_id = sample_id
  }

  call Align {
    input:
      unaligned_bam = SubsetAndUnalign.unaligned_bam,
      star_index_tar = star_index_tar,
      sample_id = sample_id,
      tools_image = tools_image
  }

  output {
    File bam_barcodes = ExtractBarcodes.barcodes
    File labels = MakeLabels.labels
    File labeled_barcodes = MakeLabels.labeled_barcodes
    File batch = MakeLabels.batch
    File unaligned_bam = SubsetAndUnalign.unaligned_bam
    File subset_counts = SubsetAndUnalign.counts
    File aligned_bam = Align.aligned_bam
    File star_log = Align.star_log
  }
}
