import argparse
import pandas as pd
import numpy as np
import networkx as nx
from collections import Counter

# -------------------- 0. Argument Parser --------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter Hi-C submatrices, cluster overlapping rectangles, and extract merged coordinates."
    )

    parser.add_argument("--input", type=str, required=True,
                        help="Path to all_submatrices_info.tsv input file")
    parser.add_argument("--out_prefix", type=str, default="result",
                        help="Prefix for output files (default: abnoral_trend_analyse_result)")
    parser.add_argument("--trend_p", type=float, default=0.05,
                        help="Maximum permutation p-value for row/column trends (default: 0.05)")
    parser.add_argument("--structure_p", type=float, default=0.05,
                        help="Maximum empirical p-value for diagonal/quadrant support (default: 0.05)")
    parser.add_argument("--structure_z", type=float, default=1.96,
                        help="Minimum absolute Z-score for diagonal/quadrant support (default: 1.96)")
    return parser.parse_args()


# -------------------- 1. 读取文件 --------------------
def main():
    args = parse_args()

    file_in = args.input
    prefix = args.out_prefix

    df = pd.read_csv(file_in, sep="\t")

    # -------------------- 2. 过滤数据 --------------------
    df['len_row'] = df['row_end'] - df['row_start']
    df['len_col'] = df['col_end'] - df['col_start']

    max_len = np.maximum(df['len_row'], df['len_col'])
    min_len = np.minimum(df['len_row'], df['len_col'])
    cond_length_ratio = max_len / min_len <= 5

    required_structure_columns = {
        'diag_p_value', 'diag_z_score', 'quad_p_value', 'quad_z_score'
    }
    missing = required_structure_columns.difference(df.columns)
    if missing:
        raise ValueError(
            "Input is missing structural feature columns: " + ", ".join(sorted(missing)) +
            ". Re-run find_abnormal_trend_from_matrix.py with the updated code."
        )

    # 二维结构作为原有趋势条件的支持证据。对角线或四象限任一显著即可，
    # 避免要求两个相关特征同时显著而过度损失召回率。
    diag_support = (df['diag_p_value'] <= args.structure_p) & \
                   (df['diag_z_score'].abs() >= args.structure_z)
    quad_support = (df['quad_p_value'] <= args.structure_p) & \
                   (df['quad_z_score'].abs() >= args.structure_z)
    structure_support = diag_support | quad_support

    cond1 = (df['linear_slope_row'] * df['linear_slope_col'] > 0) & \
            (df['spearman_rho_row'] * df['linear_slope_row'] > 0) & \
            (df['spearman_rho_col'] * df['linear_slope_col'] > 0) & \
            (df['linear_slope_row'] < 0) & (df['linear_p_row'] < args.trend_p) & \
            (df['perm_p_row'] <= args.trend_p) & (~df['unvaild_flag_row']) & \
            (~df['unvaild_flag_col']) & cond_length_ratio & structure_support


    cond2 = (df['linear_slope_row'] * df['linear_slope_col'] > 0) & \
            (df['spearman_rho_row'] * df['linear_slope_row'] > 0) & \
            (df['spearman_rho_col'] * df['linear_slope_col'] > 0) & \
            (df['linear_slope_col'] > 0) & (df['linear_p_col'] < args.trend_p) & \
            (df['perm_p_col'] <= args.trend_p) & (~df['unvaild_flag_row']) & \
            (~df['unvaild_flag_col']) & cond_length_ratio & structure_support

    df_filtered = df[cond1 | cond2].copy()
    df_filtered.reset_index(drop=True, inplace=True)

    

    out_filtered = f"{prefix}_submatrices_filtered.tsv"
    df_filtered.to_csv(out_filtered, sep="\t", index=False)
    print(f"Filtered data saved to {out_filtered}, total {len(df_filtered)} rows")

    output_columns = ['chromosome', 'cluster_id', 'start', 'end']
    if df_filtered.empty:
        pd.DataFrame(columns=df_filtered.columns.tolist() + ['cluster_id']).to_csv(
            f"{prefix}_submatrices_clustered.tsv", sep="\t", index=False
        )
        pd.DataFrame(columns=df_filtered.columns.tolist() + ['cluster_id']).to_csv(
            f"{prefix}_submatrices_clustered_filtered.tsv", sep="\t", index=False
        )
        pd.DataFrame(columns=output_columns).to_csv(
            f"{prefix}_cluster_merged_coords.tsv", sep="\t", index=False
        )
        print("No submatrices passed the combined trend and structural filters")
        return

    # -------------------- 3. 构建重叠矩形聚类 --------------------
    def rectangles_connect(r1, r2):
        # r1, r2: [row_start, row_end, col_start, col_end]
        r1_rs, r1_re, r1_cs, r1_ce = r1
        r2_rs, r2_re, r2_cs, r2_ce = r2
        
        # include boundaries
        row_overlap = r1_rs <= r2_re and r2_rs <= r1_re
        col_overlap = r1_cs <= r2_ce and r2_cs <= r1_ce
        
        return row_overlap or col_overlap

    G = nx.Graph()
    for i, row in df_filtered.iterrows():
        G.add_node(i)
    for i in range(len(df_filtered)):
        r1 = df_filtered.loc[i, ['row_start','row_end','col_start','col_end']].values
        for j in range(i+1, len(df_filtered)):
            r2 = df_filtered.loc[j, ['row_start','row_end','col_start','col_end']].values
            if rectangles_connect(r1, r2):
                G.add_edge(i, j)

    clusters = list(nx.connected_components(G))
    print(f"Detected {len(clusters)} rectangle clusters")

    # -------------------- 4. 保存聚类结果 --------------------
    cluster_result = []
    for cid, nodes in enumerate(clusters):
        for node in nodes:
            row_data = {'cluster_id': cid, **df_filtered.loc[node].to_dict()}
            for col in ['idx', 'row_start', 'row_end', 'col_start', 'col_end']:
                row_data[col] = int(row_data[col])
            cluster_result.append(row_data)

    df_clustered = pd.DataFrame(cluster_result)
    out_clustered = f"{prefix}_submatrices_clustered.tsv"
    df_clustered.to_csv(out_clustered, sep="\t", index=False)
    print(f"Clustered results saved to {out_clustered}")

    # -------------------- 5. 聚类过滤逻辑 --------------------
    def process_clusters(df):
        final_rows = []
        for cid, group in df.groupby('cluster_id'):
            if (group['linear_slope_row'].gt(0).all() or group['linear_slope_row'].lt(0).all() or
                group['linear_slope_col'].gt(0).all() or group['linear_slope_col'].lt(0).all()):
                continue
                
            # group = group.sort_values('row_start').copy()
            # if group.iloc[0]['linear_slope_row'] < 0:
            #     continue

            neg_rows = group[group['linear_slope_row'] < 0][['row_start', 'row_end']]
            pos_cols = group[group['linear_slope_col'] > 0][['col_start', 'col_end']]
            if not neg_rows.empty and not pos_cols.empty:
                neg_min, neg_max = neg_rows['row_start'].min(), neg_rows['row_end'].max()
                pos_min, pos_max = pos_cols['col_start'].min(), pos_cols['col_end'].max()
                if neg_max <= pos_min:
                    continue

            final_rows.append(group)
        if not final_rows:
            return df.iloc[0:0].copy()
        return pd.concat(final_rows, ignore_index=True)

    df_processed = process_clusters(df_clustered)
    out_processed = f"{prefix}_submatrices_clustered_filtered.tsv"
    df_processed.to_csv(out_processed, sep="\t", index=False)
    print(f"Processed clusters saved to {out_processed}")

    # -------------------- 6. 统计每簇最常出现坐标 --------------------
    def get_most_common_coords(df_clustered):
        result = []
        for cid, group in df_clustered.groupby('cluster_id'):
            chrom = group['chromosome'].iloc[0]
            pos_col_starts = group[group['linear_slope_col'] > 0]['col_start'].tolist()
            pos_col_ends   = group[group['linear_slope_col'] > 0]['col_end'].tolist()
            neg_row_starts = group[group['linear_slope_row'] < 0]['row_start'].tolist()
            neg_row_ends   = group[group['linear_slope_row'] < 0]['row_end'].tolist()
            all_starts = pos_col_starts + neg_row_starts
            all_ends   = pos_col_ends + neg_row_ends
            if not all_starts or not all_ends:
                continue
            result.append({
                'chromosome': chrom,
                'cluster_id': cid,
                'start': Counter(all_starts).most_common(1)[0][0],
                'end':   Counter(all_ends).most_common(1)[0][0]
            })
        return pd.DataFrame(result, columns=['chromosome', 'cluster_id', 'start', 'end'])

    df_coords = get_most_common_coords(df_processed)
    out_coords = f"{prefix}_cluster_merged_coords.tsv"
    df_coords.to_csv(out_coords, sep="\t", index=False)
    print(f"Merged coordinates saved to {out_coords}")


if __name__ == "__main__":
    main()
