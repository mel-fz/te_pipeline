#!/usr/bin/env python3
"""Sum Stellarscope TE counts by cell type and attach genomic locations.

Inputs are the outputs of `stellarscope resume --pooling_mode celltype`:
one count matrix (.mtx), its features.tsv and barcodes.tsv, plus the
labels file given to --celltype_tsv and the locus table (column "locus").

Output (one row per locus x cell type with a nonzero count):
  locus, celltype, total_count, n_cells, cells_in_type, <locus table columns>
The "__no_feature" row is reported on screen and left out of the table.
Notebook equivalent: Stage 8 (Cell 21).
"""
import argparse
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread


def fail(msg):
    sys.exit(f"ERROR: {msg}")


def first_col(path):
    with open(path) as fh:
        return [l.rstrip("\n").split("\t")[0] for l in fh]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mtx", required=True, help="TE count matrix (e.g. <prefix>-TE_counts.mtx)")
    ap.add_argument("--features", required=True, help="<prefix>-features.tsv")
    ap.add_argument("--barcodes", required=True, help="<prefix>-barcodes.tsv")
    ap.add_argument("--labels", required=True, help="barcode<TAB>celltype, no header")
    ap.add_argument("--locus_table", required=True, help="tab-delimited, must have a 'locus' column")
    ap.add_argument("--out", required=True, help="output table")
    ap.add_argument("--out_summary", default=None, help="optional: loci and total counts per cell type")
    args = ap.parse_args()

    m = mmread(args.mtx).tocsr()
    feat = first_col(args.features)
    bc = first_col(args.barcodes)
    print(f"matrix {m.shape}, features {len(feat)}, barcodes {len(bc)}")
    if m.shape == (len(bc), len(feat)) and len(bc) != len(feat):
        m = m.T.tocsr()                                  # make it features x barcodes
    if m.shape != (len(feat), len(bc)):
        fail("matrix shape does not match features/barcodes")

    lab = pd.read_csv(args.labels, sep="\t", header=None, names=["cb", "celltype"]).set_index("cb")["celltype"]
    ct = pd.Series(bc).map(lab)
    print("barcodes without a label:", int(ct.isna().sum()))

    feat_arr = np.asarray(feat)
    rows = []
    for c in sorted(ct.dropna().unique()):
        idx = np.where(ct.values == c)[0]
        sub = m[:, idx]
        tot = np.asarray(sub.sum(axis=1)).ravel()
        nc = np.asarray((sub > 0).sum(axis=1)).ravel()
        keep = tot > 0
        rows.append(pd.DataFrame({"locus": feat_arr[keep], "celltype": c,
                                  "total_count": tot[keep], "n_cells": nc[keep],
                                  "cells_in_type": len(idx)}))
    if not rows:
        fail("no labeled barcodes in the matrix")
    out = pd.concat(rows, ignore_index=True)

    nf = out[out.locus == "__no_feature"][["celltype", "total_count"]]
    print("__no_feature counts by cell type:")
    print(nf.to_string(index=False) if len(nf) else "  (none)")
    out = out[out.locus != "__no_feature"]

    loci = pd.read_csv(args.locus_table, sep="\t")
    if "locus" not in loci.columns:
        fail("locus table has no 'locus' column")
    out = out.merge(loci, on="locus", how="left")
    loc_cols = [c for c in loci.columns if c != "locus"]
    if loc_cols:
        print("loci without coordinates:", int(out[loc_cols[0]].isna().sum()))

    out = out.sort_values(["celltype", "total_count"], ascending=[True, False])
    out.to_csv(args.out, sep="\t", index=False)

    summ = out.groupby("celltype").agg(loci=("locus", "size"), total=("total_count", "sum"))
    print(summ.to_string())
    if args.out_summary:
        summ.to_csv(args.out_summary, sep="\t")


if __name__ == "__main__":
    main()
