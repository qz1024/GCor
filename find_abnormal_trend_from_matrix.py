#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import cooler
import math
import os
import argparse
import itertools
from multiprocessing import Pool
from scipy.stats import linregress, spearmanr
from scipy.ndimage import uniform_filter1d
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# ==================== 命令行参数 ====================
def parse_args():
    parser = argparse.ArgumentParser(
        description="Multi-chromosome Hi-C submatrix abnormal trend analysis with parallelization."
    )
    parser.add_argument("--mcool_file", type=str, required=True, help="Path to the input .mcool file")
    parser.add_argument("--boundary_file", type=str, required=True, help="Breakpoint GFF file (chrom in col1)")
    parser.add_argument("--gap_file", type=str, required=True, help="GAP file (chrom in col1)")
    parser.add_argument("--tr_file", type=str, required=False, help="TRF file (chrom in col1)")
    parser.add_argument("--global_res", type=int, default=100000, help="Global Hi-C resolution (default: 100kb)")
    parser.add_argument("--sub_res", type=int, default=50000, help="Submatrix resolution (default: 50kb)")
    parser.add_argument("--trend_smooth_window", type=int, default=5, help="Window size for smoothing trends")
    parser.add_argument("--diag_band_width", type=int, default=0,
                        help="Half-width in bins for diagonal feature bands (default: 0)")
    parser.add_argument("--perm_n", type=int, default=1000, help="Number of permutations for significance testing")
    parser.add_argument("--merge_dist", type=int, default=500000, help="Merge breakpoints closer than this distance (default: 500kb)")
    parser.add_argument("--nproces", type=int, default=12, help="Number of parallel processes (default: 12)")
    parser.add_argument("--boundary_perc", type=int, default=10, help="Percentage of top boundaries with smallest pvalues (0 for all, e.g., 10 for top 10%% )")
    return parser.parse_args()


# ==================== 工具函数 ====================
def load_breakpoints_for_chrom(chrom, args):
    """加载并合并指定染色体的断点/TR/GAP信息"""
    # --- boundary ---
    boundary_df = pd.read_csv(args.boundary_file, sep='\t', header=None)

    def extract_pvalue(attributes):
        for attr in attributes.split(';'):
            if attr.startswith('pvalue='):
                return float(attr.split('=')[1])
        return float('nan')  

    boundary_df['pvalue'] = boundary_df[8].apply(extract_pvalue)

    chrom_df = boundary_df[boundary_df[0] == chrom].copy()

    boundary_perc = max(0, min(args.boundary_perc, 100))

    if boundary_perc > 0 and not chrom_df.empty:
        chrom_df = chrom_df.sort_values('pvalue', ascending=True)
        num_to_keep = math.ceil(len(chrom_df) * (boundary_perc / 100.0))
        chrom_df = chrom_df.head(num_to_keep)

    boundary_bps = chrom_df[[3, 4]].astype(int).mean(axis=1).round().astype(int).tolist()

    # --- GAP ---
    gap_df = pd.read_csv(args.gap_file, sep='\t', header=None)
    gap_bps = gap_df[gap_df[0] == chrom][1].astype(int).tolist()
    
    # --- TRF (optional) ---
    tr_bps = []
    if args.tr_file:
        with open(args.tr_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split()
                if fields[0] == chrom:
                    tr_bps.append(int(fields[1]))

    # 合并去重并排序
    all_points = sorted(set(boundary_bps + gap_bps + tr_bps))

    # 合并距离小于 merge_dist 的点
    if len(all_points) == 0 or args.merge_dist <= 0:
        return all_points

    pts = np.array(all_points, dtype=int)
    diff = np.diff(pts)
    keep = np.concatenate(([True], diff >= args.merge_dist))
    bps_list = pts[keep].tolist()

    print(f"[INFO] BreakPoints counting: {len(bps_list)}", flush=True)

    # 保存为 txt 文件，每行一个断点
    output_file = "breakpoints.txt"
    with open(output_file, "w") as f:
        for bp in bps_list:
            f.write(f"{bp}\n")
    
    print(f"[INFO] BreakPoints saved to {output_file}", flush=True)
    return bps_list


def generate_two_intervals(positions, chrom, mcool_file,
                           min_len=2, max_len=3, max_gap=2, resolution=50000):
    clr = cooler.Cooler(f"{mcool_file}::/resolutions/{resolution}")
    chrom_length = int(clr.chromsizes.loc[chrom])
    positions = [0] + positions + [chrom_length]
    n = len(positions)
    intervals = []
    for start in range(n):
        for end in range(start + min_len - 1, min(start + max_len, n)):
            intervals.append((positions[start], positions[end]))
    valid_pairs = []
    for a, b in itertools.combinations(intervals, 2):
        first, second = sorted([a,b], key=lambda x: x[0])
        i1 = positions.index(first[1])
        i2 = positions.index(second[0])
        if i2 >= i1 and (i2-i1) <= max_gap:
            valid_pairs.append((*first,*second))
    return valid_pairs


def get_subM(row_start, row_end, col_start, col_end, cool_file, RES=50000, shrink=True):
    if shrink:
        row_start += int((row_end-row_start)/10)
        row_end   -= int((row_end-row_start)/10)
        col_start += int((col_end-col_start)/10)
        col_end   -= int((col_end-col_start)/10)
    origin_c = cooler.Cooler(f"{cool_file}::/resolutions/{RES}")
    bins = origin_c.bins()[:]
    row_idx = np.where((bins['start']<row_end)&(bins['end']>row_start))[0]
    col_idx = np.where((bins['start']<col_end)&(bins['end']>col_start))[0]
    mat = origin_c.matrix(balance=True)[:]
    sub_mat = mat[np.ix_(row_idx, col_idx)]
    return sub_mat, [row_start,row_end,col_start,col_end], RES


def _band_mean(matrix, anti=False, band_width=0):
    """Mean signal around a normalized diagonal of a rectangular matrix."""
    h, w = matrix.shape
    if h == 0 or w == 0:
        return np.nan
    values = []
    for i in range(h):
        center = 0.0 if h == 1 else i * (w - 1) / float(h - 1)
        if anti:
            center = (w - 1) - center
        left = max(0, int(math.floor(center - band_width)))
        right = min(w - 1, int(math.ceil(center + band_width)))
        values.extend(matrix[i, left:right + 1].tolist())
    return float(np.nanmean(values)) if values else np.nan


def calculate_structural_features(sub_mat, diag_band_width=0):
    """Calculate diagonal and quadrant features for one Hi-C submatrix."""
    matrix = np.asarray(sub_mat, dtype=float)
    empty = {
        'main_diag_signal': np.nan, 'anti_diag_signal': np.nan,
        'diag_contrast': np.nan, 'diag_ratio': np.nan,
        'q1_signal': np.nan, 'q2_signal': np.nan,
        'q3_signal': np.nan, 'q4_signal': np.nan,
        'quad_contrast': np.nan,
    }
    if matrix.ndim != 2 or matrix.size == 0:
        return empty

    main_signal = _band_mean(matrix, anti=False, band_width=diag_band_width)
    anti_signal = _band_mean(matrix, anti=True, band_width=diag_band_width)
    h, w = matrix.shape
    mid_h, mid_w = h // 2, w // 2
    if mid_h == 0 or mid_w == 0:
        q1 = q2 = q3 = q4 = np.nan
    else:
        q1 = float(np.nanmean(matrix[:mid_h, :mid_w]))
        q2 = float(np.nanmean(matrix[:mid_h, mid_w:]))
        q3 = float(np.nanmean(matrix[mid_h:, :mid_w]))
        q4 = float(np.nanmean(matrix[mid_h:, mid_w:]))

    return {
        'main_diag_signal': float(main_signal),
        'anti_diag_signal': float(anti_signal),
        'diag_contrast': float(anti_signal - main_signal),
        'diag_ratio': float(anti_signal / (main_signal + np.finfo(float).eps)),
        'q1_signal': q1, 'q2_signal': q2, 'q3_signal': q3, 'q4_signal': q4,
        'quad_contrast': float(((q2 + q3) - (q1 + q4)) / 2.0),
    }


def _empirical_two_sided_p(samples, observed):
    samples = np.asarray(samples, dtype=float)
    samples = samples[np.isfinite(samples)]
    if len(samples) == 0 or not np.isfinite(observed):
        return np.nan
    center = float(np.mean(samples))
    return float((np.sum(np.abs(samples - center) >= abs(observed - center)) + 1) / (len(samples) + 1))


def analyze_trend(idx, sub_mat, boundary, RES, outdir, mode='row', n_perm=1000, smooth_window=5, plot=True, remove_low_pct=0.05, remove_high_pct=0.05):

    sub_mat = np.nan_to_num(sub_mat, nan=0.0, posinf=0.0, neginf=0.0)

    # --- 去除每行或每列的最大值/最小值百分比 ---
    def trim_array(arr, low_pct, high_pct):
        arr_sorted = np.sort(arr)
        n = len(arr_sorted)
        low_k = int(n * low_pct)
        high_k = int(n * high_pct)
        if low_k + high_k >= n:
            return arr_sorted
        return arr_sorted[low_k : n - high_k]

    values = []
    values_median = []

    if mode=='row':
        # values = sub_mat.mean(axis=1)
        # values_median = np.median(sub_mat, axis=1)
        for row in sub_mat:
            trimmed = trim_array(row, remove_low_pct, remove_high_pct)
            values.append(np.mean(trimmed))
            values_median.append(np.median(trimmed))
        axis_name='row'
    else:
        # values = sub_mat.mean(axis=0)
        # values_median = np.median(sub_mat, axis=0)
        for col in sub_mat.T:
            trimmed = trim_array(col, remove_low_pct, remove_high_pct)
            values.append(np.mean(trimmed))
            values_median.append(np.median(trimmed))
        axis_name='column'

    values_valid = values
    values_median_valid = values_median

    if len(values_valid) > 1:
        values_smooth = uniform_filter1d(values_valid, size=smooth_window)
    else:
        values_smooth = values_valid

    x = np.arange(len(values_smooth))
    slope, intercept, r_value, linear_p, std_err = linregress(x, values_smooth)
    spearman_rho, spearman_p = spearmanr(x, values_smooth)

    # permutation test
    perm_rho = []
    for _ in range(n_perm):
        shuffled = np.random.permutation(values_smooth)
        rho,_ = spearmanr(x, shuffled)
        perm_rho.append(rho)
    perm_rho = np.array(perm_rho)
    # +1 修正避免置换 p 值为 0；后处理使用可配置的阈值判断
    perm_p = (np.sum(np.abs(perm_rho) >= np.abs(spearman_rho)) + 1) / (n_perm + 1)

    # filtrer unvaild submartix
    signal = np.array(values_smooth)
    threshold = np.mean(signal) / 2
    low_mask = signal < threshold
    padded = np.r_[False, low_mask, False]
    starts = np.where(~padded[:-1] & padded[1:])[0]
    ends = np.where(padded[:-1] & ~padded[1:])[0]
    lengths = ends - starts

    if len(lengths) > 0:
        max_idx = np.argmax(lengths)
        max_start, max_end = starts[max_idx], ends[max_idx]
    else:
        max_start, max_end = None, None

    max_len = lengths.max() if len(lengths) > 0 else 0
    if max_len >= len(signal) / 3:
        unvaild_flag = True  
    else:
        unvaild_flag = False


    if plot:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={'width_ratios': [1, 2]})
        im = ax1.imshow(sub_mat, cmap='Reds', origin='upper', aspect='auto')
        ax1.set_title(f'Submatrix Heatmap')
        fig.colorbar(im, ax=ax1)

        ax2.plot(x, values_valid, color='orange', label=f'{axis_name} mean')
        ax2.plot(x, values_smooth, color='blue', label=f'smoothed (w={smooth_window})')
        ax2.plot(x, values_median_valid, color='green', marker='x', linestyle='none', label=f'{axis_name} median')
        ax2.plot(x, intercept + slope * x, color='cyan', linestyle='--',
                 label=f'linear fit slope={slope:.5f} p={linear_p:.3e}')

        if max_start is not None and max_end is not None:
            ax2.axvspan(x[max_start], x[max_end-1], color='red', alpha=0.3, label='Longest low segment')

        ax2.set_xlabel(f'{axis_name} index (top -> bottom)' if mode == 'row' else f'{axis_name} index (left -> right)')
        ax2.set_ylabel('value')
        ax2.set_title(f'{axis_name.capitalize()} Trend\nSpearman={spearman_rho:.3f}, perm_p={perm_p:.3f}')
        ax2.legend()

        plt.tight_layout()
        plt.savefig(f"{outdir}/submatrix_trend_plots/submatrix_trend_{idx}_{mode}_{RES}_{boundary[0]}_{boundary[1]}_{boundary[2]}_{boundary[3]}.png",dpi=300)
        plt.close()

    return {'linear_slope':slope,'linear_p':linear_p,'spearman_rho':spearman_rho,'spearman_p':spearman_p,'perm_p':perm_p,'unvaild_flag':unvaild_flag}


def check_trnd(index, row_start, row_end, col_start, col_end, cool_file, outdir,
               RES=50000, smooth_window=5, perm_n=1000, diag_band_width=0):

    # ================== 子矩阵 ==================
    sub_mat, boundary, _ = get_subM(row_start, row_end, col_start, col_end, cool_file, RES)
    sub_mat_clean = np.nan_to_num(sub_mat, nan=0.0)

    observed_mean = np.mean(sub_mat_clean)
    h, w = sub_mat_clean.shape
    submatrix_shape = (h, w)

    genomic_dist_kb = abs(col_start - row_start) // 1000
    offset_bins = (col_start - row_start) // RES

    # 二维结构特征：和行/列趋势使用同一个子矩阵
    observed_structural = calculate_structural_features(sub_mat_clean, diag_band_width)
    background_diag_contrasts = []
    background_quad_contrasts = []

    # ================== 默认统计量（防炸） ==================
    background_means = []
    background_mean = np.nan
    background_std = np.nan
    perm_p_low = np.nan
    perm_p_high = np.nan
    p_value_signal = np.nan
    signal_effect_size = np.nan
    z_score = np.nan
    ci_low = np.nan
    ci_high = np.nan
    sig_depleted = False
    sig_enriched = False
    n_samples = 0

    # ================== 获取整条染色体矩阵 ==================
    clr = cooler.Cooler(f"{cool_file}::/resolutions/{RES}")

    try:
        bin_id = row_start // RES
        bins = clr.bins()[:]
        chrom = bins.iloc[bin_id]['chrom']

        full_mat = clr.matrix(balance=True).fetch(chrom)
        full_mat = np.nan_to_num(full_mat, nan=0.0)
        n_bins = full_mat.shape[0]

    except Exception as e:
        print(f"[WARN] Failed to fetch full matrix for bin {bin_id}: {e}")

    else:
        max_possible_start = n_bins - max(h, w) - abs(offset_bins)

        if max_possible_start > 0:
            attempts = 0
            max_attempts = perm_n * 100

            while len(background_means) < perm_n and attempts < max_attempts:
                attempts += 1

                row_bin_start = np.random.randint(0, max_possible_start + 1)
                col_bin_start = row_bin_start + offset_bins

                if (col_bin_start >= 0 and
                    row_bin_start + h <= n_bins and
                    col_bin_start + w <= n_bins):

                    rand_sub = full_mat[
                        row_bin_start:row_bin_start + h,
                        col_bin_start:col_bin_start + w
                    ]
                    background_means.append(np.mean(rand_sub))
                    bg_structural = calculate_structural_features(rand_sub, diag_band_width)
                    background_diag_contrasts.append(bg_structural['diag_contrast'])
                    background_quad_contrasts.append(bg_structural['quad_contrast'])

        # ================== 统计显著性 ==================
        if len(background_means) > 0:
            background_means = np.array(background_means)
            n_samples = len(background_means)

            background_mean = np.mean(background_means)
            background_std = np.std(background_means, ddof=1)

            # 经验置换 p 值（加 1 修正）
            perm_p_low = (np.sum(background_means <= observed_mean) + 1) / (n_samples + 1)
            perm_p_high = (np.sum(background_means >= observed_mean) + 1) / (n_samples + 1)

            # 双尾 p 值
            p_value_signal = min(1.0, 2 * min(perm_p_low, perm_p_high))

            # Z-score / 效应量
            if background_std > 0:
                z_score = (observed_mean - background_mean) / background_std
                signal_effect_size = z_score

            # Bootstrap CI
            ci_low, ci_high = np.percentile(background_means, [2.5, 97.5])
            sig_depleted = observed_mean < ci_low
            sig_enriched = observed_mean > ci_high

            # 对角线/四象限特征使用相同的随机子矩阵作为背景
            diag_bg = np.asarray(background_diag_contrasts, dtype=float)
            quad_bg = np.asarray(background_quad_contrasts, dtype=float)
            diag_bg = diag_bg[np.isfinite(diag_bg)]
            quad_bg = quad_bg[np.isfinite(quad_bg)]
            if len(diag_bg) > 0 and np.isfinite(observed_structural['diag_contrast']):
                diag_mean = float(np.mean(diag_bg))
                diag_std = float(np.std(diag_bg, ddof=1)) if len(diag_bg) > 1 else np.nan
                diag_z = ((observed_structural['diag_contrast'] - diag_mean) / diag_std
                          if np.isfinite(diag_std) and diag_std > 0 else np.nan)
                diag_p = _empirical_two_sided_p(diag_bg, observed_structural['diag_contrast'])
            else:
                diag_mean = diag_std = diag_z = diag_p = np.nan

            if len(quad_bg) > 0 and np.isfinite(observed_structural['quad_contrast']):
                quad_mean = float(np.mean(quad_bg))
                quad_std = float(np.std(quad_bg, ddof=1)) if len(quad_bg) > 1 else np.nan
                quad_z = ((observed_structural['quad_contrast'] - quad_mean) / quad_std
                          if np.isfinite(quad_std) and quad_std > 0 else np.nan)
                quad_p = _empirical_two_sided_p(quad_bg, observed_structural['quad_contrast'])
            else:
                quad_mean = quad_std = quad_z = quad_p = np.nan

        else:
            diag_mean = diag_std = diag_z = diag_p = np.nan
            quad_mean = quad_std = quad_z = quad_p = np.nan

    if len(background_means) == 0:
        diag_mean = diag_std = diag_z = diag_p = np.nan
        quad_mean = quad_std = quad_z = quad_p = np.nan

    # ================== 趋势分析（用 clean 数据） ==================
    res_row = analyze_trend(index, sub_mat_clean, boundary, RES, outdir, 'row',
                            n_perm=perm_n, smooth_window=smooth_window, plot=True)

    res_col = analyze_trend(index, sub_mat_clean, boundary, RES, outdir, 'col',
                            n_perm=perm_n, smooth_window=smooth_window, plot=True)

    # ================== 统一附加信息 ==================
    extra_info = {
        'genomic_dist_kb': float(genomic_dist_kb),
        'submatrix_shape': str(submatrix_shape),

        'background_mean_signal': float(background_mean),
        'background_std_signal': float(background_std),
        'n_background_samples': int(n_samples),

        'perm_p_low_signal': float(perm_p_low),
        # 'perm_p_high_signal': float(perm_p_high),
        'p_value_signal': float(p_value_signal),

        'z_score_signal': float(z_score),
        # 'signal_effect_size': float(signal_effect_size),

        'ci_low_95': float(ci_low),
        # 'ci_high_95': float(ci_high),
        'sig_depleted': bool(sig_depleted),
        # 'sig_enriched': bool(sig_enriched)

        # ===== 二维结构特征 =====
        'main_diag_signal': observed_structural['main_diag_signal'],
        'anti_diag_signal': observed_structural['anti_diag_signal'],
        'diag_contrast': observed_structural['diag_contrast'],
        'diag_ratio': observed_structural['diag_ratio'],
        'diag_background_mean': float(diag_mean),
        'diag_background_std': float(diag_std),
        'diag_z_score': float(diag_z),
        'diag_p_value': float(diag_p),
        'q1_signal': observed_structural['q1_signal'],
        'q2_signal': observed_structural['q2_signal'],
        'q3_signal': observed_structural['q3_signal'],
        'q4_signal': observed_structural['q4_signal'],
        'quad_contrast': observed_structural['quad_contrast'],
        'quad_background_mean': float(quad_mean),
        'quad_background_std': float(quad_std),
        'quad_z_score': float(quad_z),
        'quad_p_value': float(quad_p),
    }

    res_row.update(extra_info)
    res_col.update(extra_info)

    return (res_row, res_col), boundary


def save_all_submatrices_to_single_file(all_data, chrom, output_file=None):

    if output_file is None:
        output_file = f"{chrom}_submatrices_info.tsv"

    records = []

    for idx, d in enumerate(all_data):
        res_row = d['res_row']
        res_col = d['res_col']

        rec = {
            'chromosome': chrom,
            'idx': idx,
            'row_start': d['row_start'],
            'row_end': d['row_end'],
            'col_start': d['col_start'],
            'col_end': d['col_end'],

            # ===== 行 / 列趋势 =====
            'linear_slope_row': res_row.get('linear_slope', np.nan),
            'linear_p_row': res_row.get('linear_p', np.nan),
            'spearman_rho_row': res_row.get('spearman_rho', np.nan),
            'spearman_p_row': res_row.get('spearman_p', np.nan),
            'perm_p_row': res_row.get('perm_p', np.nan),
            'unvaild_flag_row': res_row.get('unvaild_flag', False),

            'linear_slope_col': res_col.get('linear_slope', np.nan),
            'linear_p_col': res_col.get('linear_p', np.nan),
            'spearman_rho_col': res_col.get('spearman_rho', np.nan),
            'spearman_p_col': res_col.get('spearman_p', np.nan),
            'perm_p_col': res_col.get('perm_p', np.nan),
            'unvaild_flag_col': res_col.get('unvaild_flag', False),

            # ===== 远距离信号显著性（核心）=====
            'genomic_dist_kb': res_col.get('genomic_dist_kb', np.nan),
            'submatrix_shape': res_col.get('submatrix_shape', None),

            'background_mean_signal': res_col.get('background_mean_signal', np.nan),
            'background_std_signal': res_col.get('background_std_signal', np.nan),
            'n_background_samples': res_col.get('n_background_samples', 0),

            'perm_p_low_signal': res_col.get('perm_p_low_signal', np.nan),
            # 'perm_p_high_signal': res_col.get('perm_p_high_signal', np.nan),
            'p_value_signal': res_col.get('p_value_signal', np.nan),

            'z_score_signal': res_col.get('z_score_signal', np.nan),
            # 'signal_effect_size': res_col.get('signal_effect_size', np.nan),

            'sig_depleted': res_col.get('sig_depleted', False),
            # 'sig_enriched': res_col.get('sig_enriched', False),

            # ===== 二维结构特征 =====
            'main_diag_signal': res_col.get('main_diag_signal', np.nan),
            'anti_diag_signal': res_col.get('anti_diag_signal', np.nan),
            'diag_contrast': res_col.get('diag_contrast', np.nan),
            'diag_ratio': res_col.get('diag_ratio', np.nan),
            'diag_background_mean': res_col.get('diag_background_mean', np.nan),
            'diag_background_std': res_col.get('diag_background_std', np.nan),
            'diag_z_score': res_col.get('diag_z_score', np.nan),
            'diag_p_value': res_col.get('diag_p_value', np.nan),
            'q1_signal': res_col.get('q1_signal', np.nan),
            'q2_signal': res_col.get('q2_signal', np.nan),
            'q3_signal': res_col.get('q3_signal', np.nan),
            'q4_signal': res_col.get('q4_signal', np.nan),
            'quad_contrast': res_col.get('quad_contrast', np.nan),
            'quad_background_mean': res_col.get('quad_background_mean', np.nan),
            'quad_background_std': res_col.get('quad_background_std', np.nan),
            'quad_z_score': res_col.get('quad_z_score', np.nan),
            'quad_p_value': res_col.get('quad_p_value', np.nan),
        }

        records.append(rec)

    df_all = pd.DataFrame(records)
    df_all.to_csv(output_file, sep='\t', index=False)
    print(f"[INFO] All submatrix info saved to {output_file}", flush=True)



def plot_chromosome_heatmap(chrom, mcool_file, analyzed_regions, all_bps, output_dir, global_res):

    c_global = cooler.Cooler(f"{mcool_file}::resolutions/{global_res}")
    full_mat = c_global.matrix(balance=False).fetch(chrom)
    bins_global = c_global.bins().fetch(chrom)
    bin_starts_global = bins_global['start'].values

    def pos_to_bin(pos):
        idx = np.searchsorted(bin_starts_global, pos, side='right') - 1
        return np.clip(idx, 0, len(bin_starts_global) - 1)

    plt.figure(figsize=(14, 12))
    im = plt.imshow(np.log1p(full_mat), cmap='Reds', origin='upper', aspect='auto')
    plt.colorbar(im, shrink=0.8, label='log1p(Raw Contact)')

    plt.title(f'Chromosome {chrom} Hi-C Matrix ({global_res//1000} kb)\n'
              f'{len(analyzed_regions)} regions | {len(all_bps)} breakpoints',
              fontsize=14, pad=25)
    plt.xlabel('Genomic Position (bin)')
    plt.ylabel('Genomic Position (bin)')
    ax = plt.gca()
    n_bins = full_mat.shape[0]

    colors = ['lime', 'cyan', 'yellow', 'magenta', 'orange',
              'deepskyblue', 'lightgreen', 'gold', 'pink', 'violet']
    for i, region in enumerate(analyzed_regions):
        r0, r1 = pos_to_bin(region['row_start']), pos_to_bin(region['row_end'])
        c0, c1 = pos_to_bin(region['col_start']), pos_to_bin(region['col_end'])
        if r1 <= r0 or c1 <= c0: 
            continue
        color = colors[i % len(colors)]
        rect = Rectangle((c0, r0), c1-c0, r1-r0, linewidth=1, edgecolor=color, facecolor='none', zorder=3)
        ax.add_patch(rect)

    marker_size = 20
    offset = 0.2
    for bp in all_bps:
        bin_idx = pos_to_bin(bp)
        if 0 <= bin_idx < n_bins:
            ax.scatter(bin_idx, bin_idx - offset, marker='^', s=marker_size,
                       color='red', edgecolor='darkred', linewidth=1, zorder=5)

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f'{chrom}_full_matrix_with_regions.png')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Chromosome {chrom} heatmap saved to {output_path}", flush=True)


# ==================== 单染色体处理 ====================
def process_one_chrom(chrom, args):
    try:
        print(f"[INFO] Processing chromosome {chrom}", flush=True)
        outdir=f"results/{chrom}"
        os.makedirs(f"{outdir}",exist_ok=True)
        os.makedirs(f"{outdir}/submatrix_trend_data",exist_ok=True)
        os.makedirs(f"{outdir}/submatrix_trend_plots",exist_ok=True)
        os.makedirs(f"{outdir}/chromsome_overview",exist_ok=True)

        all_bps = load_breakpoints_for_chrom(chrom,args)
        if len(all_bps)<2:
            print(f"[WARN] {chrom}: not enough breakpoints, skip", flush=True)
            return

        valid_pairs = generate_two_intervals(all_bps,chrom,args.mcool_file,
                                             min_len=2,max_len=3,max_gap=1,
                                             resolution=args.sub_res)
        all_submat_data=[]
        analyzed_regions=[]
        for idx,(rs,re,cs,ce) in enumerate(valid_pairs):
            try:
                (res_row,res_col),boundary=check_trnd(idx,
                                                     rs,re,cs,ce,
                                                     args.mcool_file,
                                                     RES=args.sub_res,
                                                     smooth_window=args.trend_smooth_window,
                                                     perm_n=args.perm_n,
                                                     diag_band_width=args.diag_band_width,
                                                     outdir=outdir)
                all_submat_data.append({'row_start':rs,'row_end':re,'col_start':cs,'col_end':ce,
                                        'res_row':res_row,'res_col':res_col})
                analyzed_regions.append({'row_start':rs,'row_end':re,'col_start':cs,'col_end':ce})
            except Exception as e:
                print(f"[ERROR] {chrom} submatrix failed: {e}", flush=True)

        save_all_submatrices_to_single_file(all_submat_data,chrom,
                                            output_file=f"{outdir}/submatrix_trend_data/{chrom}_submatrices_info.tsv")

        # 绘制全局热图
        plot_chromosome_heatmap(
            chrom=chrom,
            mcool_file=args.mcool_file,
            analyzed_regions=analyzed_regions,
            all_bps=all_bps,
            output_dir=f"{outdir}/chromsome_overview/",
            global_res=args.global_res
        )

        print(f"[INFO] Chromosome {chrom} finished", flush=True)
    except Exception as e:
        print(f"[FATAL] Chromosome {chrom} failed: {e}", flush=True)


# ==================== 多染色体并行 ====================
def run_parallel(args):
    clr = cooler.Cooler(f"{args.mcool_file}::/resolutions/{args.sub_res}")
    chroms = list(clr.chromnames)
    print(f"[INFO] Total chromosomes detected: {len(chroms)}")
    with Pool(processes=args.nproces) as pool:
        pool.starmap(process_one_chrom, [(chrom,args) for chrom in chroms])


# ==================== 主程序 ====================
if __name__=="__main__":
    args = parse_args()
    run_parallel(args)
