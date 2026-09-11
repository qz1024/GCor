import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from multiprocessing import Pool

def parse_args():
    parser = argparse.ArgumentParser(
        description="Count Hi-C reads crossing structural variant breakpoints and output TSV + line plots."
    )

    parser.add_argument("-i", "--input", required=True, help="Input Hi-C pairs file (.pairs or .pairs.gz)")
    parser.add_argument("-r", "--regions", required=True, help="TSV file with columns: chromosome, cluster_id, start, end")
    parser.add_argument("-b", "--bin", type=int, default=100, help="Bin size for counting reads")
    parser.add_argument("-d", "--distance", type=int, default=100000, help="Maximum distance threshold for reads")
    parser.add_argument("-e", "--extend_len", type=int, default=250000, help="Extension length around each breakpoint")
    parser.add_argument("-o", "--output-prefix", default="hic_crossing_reads", help="Output file prefix")
    parser.add_argument("-n", "--nprocess", type=int, default=4, help="Number of parallel processes")
    return parser.parse_args()


def process_single_breakpoint(args):
    """
    Process a single breakpoint, return curve and minimum point info
    args: (df_pairs, chrom, breakPoint, extend_len, bin_size, max_distance)
    """
    df_pairs, chrom, breakPoint, extend_len, bin_size, max_distance = args

    start = max(0, breakPoint - extend_len)
    end = breakPoint + extend_len

    # Filter reads on the same chromosome
    df_chr = df_pairs[(df_pairs["chrom1"] == chrom) & (df_pairs["chrom2"] == chrom)]
    if df_chr.empty:
        return None

    # Filter reads in the extended region and distance threshold
    df_chr = df_chr[
        df_chr["start1"].between(start, end) &
        df_chr["end2"].between(start, end) &
        (df_chr["start1"] - df_chr["end2"]).abs() < max_distance
    ]
    if df_chr.empty:
        return None

    # Set bins
    bins = np.arange(start, end + bin_size, bin_size)
    bin_counts = np.zeros(len(bins)-1, dtype=int)

    # Count crossing reads
    for _, row in df_chr.iterrows():
        s = max(start, min(row["start1"], row["end2"]))
        e = min(end, max(row["start1"], row["end2"]))
        if s >= e:
            continue
        start_bin = np.searchsorted(bins, s, side='right') - 1
        end_bin = np.searchsorted(bins, e, side='right') - 1
        start_bin = max(0, start_bin)
        end_bin = min(len(bin_counts)-1, end_bin)
        bin_counts[start_bin:end_bin+1] += 1

    # Find minimum point
    n_bins = len(bin_counts)
    start_idx = n_bins // 5        # 前 1/5 开始
    end_idx = n_bins * 4 // 5      # 到后 4/5
    middle_bin_counts = bin_counts[start_idx:end_idx+1]
    min_idx_local = np.argmin(middle_bin_counts)
    min_idx = start_idx + min_idx_local
    min_coord = int((bins[min_idx] + bins[min_idx+1]) // 2)
    min_value = int(bin_counts[min_idx])


    return {"breakPoint": breakPoint, "bins": bins, "counts": bin_counts, "min_coord": min_coord, "min_value": min_value}


def main():
    args = parse_args()

    # Read pairs file
    compression = "gzip" if args.input.endswith(".gz") else None
    cols = ['readID', 'chrom1', 'start1', 'chrom2', 'end2', 'strand1', 'strand2', 'type']
    df_pairs = pd.read_csv(args.input, sep='\s+', names=cols, header=None, compression=compression)

    # Read regions TSV
    df_regions = pd.read_csv(args.regions, sep="\t")

    results_summary = []

    # Process each breakpoint separately
    for row in df_regions.itertuples(index=False):
        chrom = row.chromosome
        cluster_id = row.cluster_id
        breakPoint_list = [int(row.start), int(row.end)]

        # Prepare tasks for multiprocessing
        tasks = [(df_pairs, chrom, breakPoint, args.extend_len, args.bin, args.distance) for breakPoint in breakPoint_list]
        with Pool(args.nprocess) as pool:
            breakPoint_results = pool.map(process_single_breakpoint, tasks)

        # Plot each breakpoint separately
        for r in breakPoint_results:
            if r is None:
                continue

            plt.figure(figsize=(18,5))
            x = (r["bins"][:-1] + r["bins"][1:]) / 2
            y = r["counts"]
            plt.plot(x, y, linewidth=1, label="Crossing reads")

            # Draw vertical dashed line and mark the minimum point
            plt.axvline(r["min_coord"], color='red', linestyle='--', label=f"Min at {r['min_coord']} bp")
            plt.scatter(r["min_coord"], r["min_value"], color='red')

            # Annotate the minimum bp above the point
            plt.text(r["min_coord"], r["min_value"] + max(y)*0.02, f"{r['min_coord']}", color='red', ha='center', va='bottom', fontsize=10)

            plt.xlabel("Genomic Position (bp)")
            plt.ylabel("Crossing Reads Count")
            plt.title(f"{chrom} cluster {cluster_id} breakPoint:{r['breakPoint']}")
            plt.legend()
            plt.tight_layout()
            fig_out = f"{args.output_prefix}_{chrom}_{cluster_id}_breakPoint{r['breakPoint']}.png"
            plt.savefig(fig_out, dpi=300)
            plt.close()
            print(f"Saved figure: {fig_out}")


            # Save minimum point
            results_summary.append({
                "chrom": chrom,
                "cluster_id": cluster_id,
                "breakPoint": r["breakPoint"],
                "min_coord": r["min_coord"],
                "min_value": r["min_value"]
            })

    # Save all minimum points
    df_out = pd.DataFrame(results_summary)
    df_out.to_csv(f"{args.output_prefix}_min_points.tsv", sep="\t", index=False)
    print(f"Saved minimum points table: {args.output_prefix}_min_points.tsv")


if __name__ == "__main__":
    main()
