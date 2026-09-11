#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cooler
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Overlay submatrix regions on Hi-C heatmap"
    )
    parser.add_argument("--mcool", required=True, help="Input .mcool file")
    parser.add_argument("--regions", required=True, help="Submatrix regions file (TSV)")
    parser.add_argument("--resolution", type=int, default=100000, help="Hi-C resolution")
    parser.add_argument("--outdir", default="hic_with_submatrices", help="Output directory")
    return parser.parse_args()


def plot_one_chrom(chrom, clr, regions_df, outdir):
    mat = clr.matrix(balance=False).fetch(chrom)
    bins = clr.bins().fetch(chrom)
    bin_starts = bins["start"].values

    def bp_to_bin(pos):
        return np.searchsorted(bin_starts, pos, side="right") - 1

    plt.figure(figsize=(12, 10))
    plt.imshow(np.log1p(mat), cmap="Reds", origin="upper", aspect="auto")
    plt.colorbar(label="log1p(contact)")
    plt.title(f"{chrom} Hi-C with submatrices")

    ax = plt.gca()

    colors = ["cyan", "lime", "magenta", "yellow", "orange"]

    for i, r in regions_df.iterrows():
        r0 = bp_to_bin(r.row_start)
        r1 = bp_to_bin(r.row_end)
        c0 = bp_to_bin(r.col_start)
        c1 = bp_to_bin(r.col_end)

        if r1 <= r0 or c1 <= c0:
            continue

        rect = Rectangle(
            (c0, r0),
            c1 - c0,
            r1 - r0,
            linewidth=1.5,
            edgecolor=colors[i % len(colors)],
            facecolor=colors[i % len(colors)],
            alpha=0.25
        )
        ax.add_patch(rect)

    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"{chrom}_hic_submatrices.png")
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()

    print(f"[INFO] Saved {out}")


def main():
    args = parse_args()
    clr = cooler.Cooler(f"{args.mcool}::/resolutions/{args.resolution}")

    regions = pd.read_csv(args.regions, sep="\t")

    for chrom in regions["chromosome"].unique():
        plot_one_chrom(
            chrom,
            clr,
            regions[regions.chromosome == chrom],
            args.outdir
        )


if __name__ == "__main__":
    main()
