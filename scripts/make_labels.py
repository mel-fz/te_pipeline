#!/usr/bin/env python3
"""Build per-sample cell-type labels for Stellarscope --celltype_tsv.

Given the distinct cell barcodes (CB) found in one sample's BAM and the shared
cell-annotation table, this script
  1. (optional) checks the sample is listed in the valid-files list,
  2. works out which annotation batch the BAM belongs to (annotation barcodes
     look like "batchN_celltypist_<barcode>"),
  3. keeps only barcodes that are in that batch AND in the BAM,
  4. writes the labels file (barcode<TAB>celltype, no header) and the list of
     labeled barcodes.

Any failed check exits with a non-zero status and a message, so a workflow
stops instead of silently producing bad labels.
"""
import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def fail(msg):
    sys.exit(f"ERROR: {msg}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam_cbs", required=True, help="text file, one distinct BAM barcode per line")
    ap.add_argument("--annotation", required=True, help="cell_annotations.tsv (tab-delimited)")
    ap.add_argument("--sample_id", required=True, help="e.g. TOS019341")
    ap.add_argument("--valid_files", default=None,
                    help="optional valid-files list; sample_id must appear in at least one line")
    ap.add_argument("--celltype_col", default="l2_celltype")
    ap.add_argument("--min_frac", type=float, default=0.5,
                    help="minimum fraction of BAM barcodes that must receive a label")
    ap.add_argument("--min_ratio", type=float, default=2.0,
                    help="best batch overlap must be at least this many times the runner-up")
    ap.add_argument("--out_labels", required=True, help="output: barcode<TAB>celltype, no header")
    ap.add_argument("--out_barcodes", required=True, help="output: labeled barcodes, one per line")
    ap.add_argument("--out_batch", default=None, help="optional output: name of the matched batch")
    args = ap.parse_args()

    # 1) valid-files check
    if args.valid_files:
        with open(args.valid_files) as fh:
            hits = [l for l in fh if args.sample_id in l]
        if not hits:
            fail(f"{args.sample_id} is not in {args.valid_files}")
        print(f"valid-files: {len(hits)} line(s) mention {args.sample_id}")

    # 2) read inputs
    bam_cb = set(Path(args.bam_cbs).read_text().split())
    if not bam_cb:
        fail(f"no barcodes found in {args.bam_cbs}")
    ann = pd.read_csv(args.annotation, sep="\t", usecols=["barcode", args.celltype_col])
    parts = ann["barcode"].str.extract(r"^(batch\d+)_celltypist_(.+)$")
    ann["batch"], ann["cb"] = parts[0], parts[1]
    n_bad = int(ann["batch"].isna().sum())
    print(f"annotation rows not matching the batchN_celltypist_ pattern: {n_bad}")

    # 3) which batch does this BAM belong to?
    overlap = ann[ann["cb"].isin(bam_cb)].groupby("batch").size().sort_values(ascending=False)
    if overlap.empty:
        fail("no BAM barcode matches any annotation barcode")
    print(overlap.to_string())
    batch = overlap.index[0]
    runner_up = overlap.iloc[1] if len(overlap) > 1 else 0
    if overlap.iloc[0] < args.min_ratio * runner_up:
        fail(f"no clear batch match (best {overlap.iloc[0]} vs runner-up {runner_up})")

    in_batch = ann[ann["batch"] == batch]
    dups = int(in_batch["cb"].duplicated().sum())
    if dups:
        fail(f"{dups} duplicate barcodes within {batch}; barcode-only labeling is unsafe")

    b = in_batch[in_batch["cb"].isin(bam_cb)]
    frac = len(b) / len(bam_cb)
    print(f"batch {batch}: {len(b)} of {len(bam_cb)} BAM barcodes labeled ({frac:.1%})")
    if frac < args.min_frac:
        fail(f"only {frac:.1%} of BAM barcodes labeled (minimum {args.min_frac:.0%})")

    # 4) write outputs
    b[["cb", args.celltype_col]].to_csv(args.out_labels, sep="\t", header=False, index=False)
    b["cb"].to_csv(args.out_barcodes, header=False, index=False)
    if args.out_batch:
        Path(args.out_batch).write_text(batch + "\n")
    print(b[args.celltype_col].value_counts().to_string())


if __name__ == "__main__":
    main()
