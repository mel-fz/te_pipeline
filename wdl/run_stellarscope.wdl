version 1.0

# RunStellarscope: aligned BAM + cell-type labels -> TE counts per cell type with locations.
# Notebook equivalents: Cellsort = Cell 12; Assign = Cell 15; Resume = Cell 18; Summarize = Cell 21.
# Each task logs memory use every 30 s (mem_used_mb.txt) so peak memory can be measured
# and the memory settings below (guesses; the notebook machine had 50 GB) can be tuned.

task Cellsort {
  input {
    File aligned_bam
    File labeled_barcodes
    String sample_id
    String tools_image
    Int threads = 8
    Int memory_gb = 32
  }

  Int disk_gb = ceil(size(aligned_bam, "GB") * 4) + 20

  command <<<
    set -euo pipefail
    ( while true; do
        echo "$(date +%T) $(awk '/MemTotal/{t=$2}/MemAvailable/{a=$2}END{print int((t-a)/1024)}' /proc/meminfo)" >> mem_used_mb.txt
        sleep 30
      done ) &
    MON=$!
    trap 'kill $MON 2>/dev/null || true' EXIT

    mkdir -p tmp
    stellarscope cellsort --nproc ~{threads} --tempdir tmp \
      --outfile ~{sample_id}_labeled_sortedByCB.bam \
      ~{aligned_bam} ~{labeled_barcodes}

    echo "peak memory used (MB): $(sort -k2 -n mem_used_mb.txt | tail -1 | cut -d' ' -f2)"
  >>>

  output {
    File sorted_bam = "~{sample_id}_labeled_sortedByCB.bam"
    File mem_log = "mem_used_mb.txt"
  }

  runtime {
    docker: tools_image
    cpu: threads
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    maxRetries: 1
  }
}

task Assign {
  input {
    File sorted_bam
    File te_gtf
    File labeled_barcodes
    String sample_id
    String tools_image
    String umi_tag = "UM"
    String stranded_mode = "F"
    Int memory_gb = 48
  }

  Int disk_gb = ceil(size(sorted_bam, "GB") * 2 + size(te_gtf, "GB")) + 30

  command <<<
    set -euo pipefail
    ( while true; do
        echo "$(date +%T) $(awk '/MemTotal/{t=$2}/MemAvailable/{a=$2}END{print int((t-a)/1024)}' /proc/meminfo)" >> mem_used_mb.txt
        sleep 30
      done ) &
    MON=$!
    trap 'kill $MON 2>/dev/null || true' EXIT

    mkdir -p out
    stellarscope assign --exp_tag ~{sample_id}_ct --outdir out --skip_em \
      --umi_tag ~{umi_tag} --stranded_mode ~{stranded_mode} \
      --whitelist ~{labeled_barcodes} --logfile assign.log \
      ~{sorted_bam} ~{te_gtf}

    test -s out/~{sample_id}_ct-checkpoint.dedup_umi.pickle
    tail -3 assign.log
    echo "peak memory used (MB): $(sort -k2 -n mem_used_mb.txt | tail -1 | cut -d' ' -f2)"
  >>>

  output {
    File checkpoint = "out/~{sample_id}_ct-checkpoint.dedup_umi.pickle"
    File stats = "out/~{sample_id}_ct-stats.final.tsv"
    File log = "assign.log"
    File mem_log = "mem_used_mb.txt"
  }

  runtime {
    docker: tools_image
    cpu: 2
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    maxRetries: 1
  }
}

task Resume {
  input {
    File checkpoint
    File labels                # barcode<TAB>celltype, from PrepareSample
    String sample_id
    String tools_image
    Int nproc = 4
    Int max_iter = 500
    Int memory_gb = 48
  }

  Int disk_gb = ceil(size(checkpoint, "GB") * 6) + 30
  String prefix = "~{sample_id}_ct_celltype"

  command <<<
    set -euo pipefail
    ( while true; do
        echo "$(date +%T) $(awk '/MemTotal/{t=$2}/MemAvailable/{a=$2}END{print int((t-a)/1024)}' /proc/meminfo)" >> mem_used_mb.txt
        sleep 30
      done ) &
    MON=$!
    trap 'kill $MON 2>/dev/null || true' EXIT

    mkdir -p out
    stellarscope resume --exp_tag ~{prefix} --outdir out \
      --pooling_mode celltype --celltype_tsv ~{labels} --nproc ~{nproc} \
      --use_every_reassign_mode --max_iter ~{max_iter} \
      --logfile resume.log ~{checkpoint}

    test -s out/~{prefix}-TE_counts.mtx
    ls -lh out
    tail -3 resume.log
    echo "peak memory used (MB): $(sort -k2 -n mem_used_mb.txt | tail -1 | cut -d' ' -f2)"
  >>>

  output {
    File counts_mtx = "out/~{prefix}-TE_counts.mtx"          # default matrix (best_exclude)
    Array[File] all_mtx = glob("out/~{prefix}-TE_counts*.mtx") # every reassign mode
    File features = "out/~{prefix}-features.tsv"
    File barcodes = "out/~{prefix}-barcodes.tsv"
    File stats = "out/~{prefix}-stats.final.tsv"
    File log = "resume.log"
    File mem_log = "mem_used_mb.txt"
  }

  runtime {
    docker: tools_image
    cpu: nproc
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    maxRetries: 1
  }
}

task Summarize {
  input {
    File script                # scripts/summarize_by_celltype.py
    File counts_mtx
    File features
    File barcodes
    File labels
    File locus_table
    String sample_id
    String tools_image
  }

  command <<<
    set -euo pipefail
    python3 ~{script} \
      --mtx ~{counts_mtx} --features ~{features} --barcodes ~{barcodes} \
      --labels ~{labels} --locus_table ~{locus_table} \
      --out ~{sample_id}_counts_by_celltype_with_location.tsv \
      --out_summary ~{sample_id}_celltype_summary.tsv
  >>>

  output {
    File counts_by_celltype = "~{sample_id}_counts_by_celltype_with_location.tsv"
    File celltype_summary = "~{sample_id}_celltype_summary.tsv"
  }

  runtime {
    docker: tools_image
    cpu: 2
    memory: "16 GB"
    disks: "local-disk 30 HDD"
    maxRetries: 1
  }
}

workflow RunStellarscope {
  input {
    # per sample (outputs of PrepareSample)
    File aligned_bam
    File labels
    File labeled_barcodes
    String sample_id
    # shared
    File te_gtf                # retro_extended.hg38.gtf
    File locus_table
    File summarize_script
    String tools_image         # e.g. melfz/te-pipeline:0.1.0
  }

  call Cellsort {
    input: aligned_bam = aligned_bam, labeled_barcodes = labeled_barcodes,
           sample_id = sample_id, tools_image = tools_image
  }

  call Assign {
    input: sorted_bam = Cellsort.sorted_bam, te_gtf = te_gtf, labeled_barcodes = labeled_barcodes,
           sample_id = sample_id, tools_image = tools_image
  }

  call Resume {
    input: checkpoint = Assign.checkpoint, labels = labels,
           sample_id = sample_id, tools_image = tools_image
  }

  call Summarize {
    input: script = summarize_script, counts_mtx = Resume.counts_mtx,
           features = Resume.features, barcodes = Resume.barcodes,
           labels = labels, locus_table = locus_table,
           sample_id = sample_id, tools_image = tools_image
  }

  output {
    File counts_by_celltype = Summarize.counts_by_celltype
    File celltype_summary = Summarize.celltype_summary
    File counts_mtx = Resume.counts_mtx
    Array[File] all_mtx = Resume.all_mtx
    File features = Resume.features
    File barcodes = Resume.barcodes
    File assign_stats = Assign.stats
    File resume_stats = Resume.stats
    File assign_log = Assign.log
    File resume_log = Resume.log
    Array[File] mem_logs = [Cellsort.mem_log, Assign.mem_log, Resume.mem_log]
  }
}
