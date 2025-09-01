import os
import pandas as pd
import numpy as np
import sys
import yaml
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import poisson
from collections import defaultdict
from sklearn.metrics import roc_curve, auc

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'protinseq-main')))
import core.protinseq_tools as pst


# ─────────────────────────────
# 🔧 Config Loading
# ─────────────────────────────
def load_config(config_path="config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


# ─────────────────────────────
# 📥 Data Loaders
# ─────────────────────────────
def load_gold_from_excel(excel_path, sheet_name, model_name):
    df = pd.read_excel(excel_path, sheet_name=sheet_name)
    df = df[df['model'] == model_name].copy()
    # Sanity check for model match
    if df.empty:
        print(f"❌ No rows found for model: '{model_name}' in sheet '{sheet_name}' of {excel_path}")
        print("📋 Available model names in sheet:")
        print(df['model'].unique())
        raise ValueError(f"No data found for model '{model_name}' — check for typos.")

    df['identifier_clean'] = df['identifier'].astype(str).str.strip().str.split('.').str[0]
    return df
    return df


# ─────────────────────────────
# 📊 Summary Functions
# ─────────────────────────────
def summarize_insertions_by_class(all_reads, ncbi, gold_df, pst, out_path, selected_samples):
    i_metrics = []
    for sample_name, sample_data in all_reads.items():
        for gene_id, coords in ncbi.items():
            try:
                df = pst.get_f_metric(sample_data, coords)
                df_I = df[df['metric'] == 'I']
                for _, row in df_I.iterrows():
                    i_metrics.append({
                        'sample': sample_name,
                        'gene': gene_id,
                        'frame': row['frame'],
                        'I_value': row['value']
                    })
            except Exception:
                continue

    i_df = pd.DataFrame(i_metrics)
    i_df['gene_clean'] = i_df['gene'].astype(str).str.strip().str.split('.').str[0]
    i_df = i_df.merge(gold_df[['identifier_clean', 'class']],
                      left_on='gene_clean', right_on='identifier_clean', how='left')
    i_df.dropna(subset=['class'], inplace=True)
    i_df['sample_type'] = i_df['sample'].str.extract(r'(B_Cm\d+)')

    summary = i_df.groupby(['class', 'sample_type', 'frame'])['I_value'].sum().reset_index()
    total_per_sample = summary.groupby(['class', 'sample_type'])['I_value'].transform('sum')
    summary['percentage'] = 100 * summary['I_value'] / total_per_sample

    os.makedirs(out_path, exist_ok=True)
    out_file = os.path.join(out_path, "insertions_by_class_ncbi.csv")
    summary.to_csv(out_file, index=False)
    print("✅ Saved:", out_file)
    return summary


def summarize_insertions_per_frame(all_reads, gene_dict, pst, label, out_path, selected_samples):
    i_metrics = []
    for sample_name, sample_data in all_reads.items():
        for gene_name in gene_dict:
            try:
                df = pst.get_f_metric(sample_data, gene_dict[gene_name])
                df_I = df[df['metric'] == 'I']
                for _, row in df_I.iterrows():
                    i_metrics.append({
                        'sample': sample_name,
                        'gene': gene_name,
                        'frame': row['frame'],
                        'I_value': row['value']
                    })
            except Exception:
                continue

    df = pd.DataFrame(i_metrics)
    df['sample_type'] = df['sample'].str.extract(r'(B_Cm\d+)')
    frame_sums = df.groupby(['sample_type', 'frame'])['I_value'].sum().unstack(fill_value=0)
    frame_sums = frame_sums[[1, 2, 3]] if all(f in frame_sums.columns for f in [1, 2, 3]) else frame_sums
    frame_sums = frame_sums.reset_index().melt(id_vars='sample_type', var_name='frame', value_name='I_value')

    total_per_sample = frame_sums.groupby('sample_type')['I_value'].transform('sum')
    frame_sums['percentage'] = 100 * frame_sums['I_value'] / total_per_sample

    os.makedirs(out_path, exist_ok=True)
    out_file = os.path.join(out_path, f"{label}_insertions_per_frame.csv")
    frame_sums.to_csv(out_file, index=False)
    print("✅ Saved:", out_file)

    return frame_sums


# ─────────────────────────────
# 📈 Plotting
# ─────────────────────────────
frame_palette = {1: "#FFA500", 2: "#87CEFA", 3: "#9370DB"}


def plot_ncbi_by_class(csv_path, title, output_path):
    df = pd.read_csv(csv_path)
    df['frame'] = df['frame'].astype(int)

    g = sns.FacetGrid(df, col="class", col_order=["E", "F", "N"], sharey=False, height=5, aspect=1.2)
    g.map_dataframe(sns.barplot, x="sample_type", y="I_value", hue="frame", palette=frame_palette, dodge=True)

    for ax, class_label in zip(g.axes.flat, ["E", "F", "N"]):
        ax.grid(False)
        class_data = df[df["class"] == class_label]
        bars = ax.patches
        n_frames = class_data["frame"].nunique()
        bar_width = bars[0].get_width() if bars else 0.8 / n_frames
        xticks = ax.get_xticks()
        xticklabels = [t.get_text() for t in ax.get_xticklabels()]

        for bar in bars:
            center = bar.get_x() + bar.get_width() / 2
            height = bar.get_height()
            sample_idx = min(range(len(xticks)), key=lambda i: abs(center - xticks[i]))
            sample_type = xticklabels[sample_idx]
            offset = center - xticks[sample_idx]
            frame_order = sorted(class_data['frame'].unique())
            frame_index = round(offset / bar_width) + (n_frames // 2)
            if 0 <= frame_index < len(frame_order):
                frame = frame_order[frame_index]
            else:
                continue
            row = class_data[(class_data["sample_type"] == sample_type) & (class_data["frame"] == frame)]
            if not row.empty:
                I_value = int(row["I_value"].values[0])
                pct = row["percentage"].values[0]
                label = f"{I_value}\n({pct:.1f}%)"
                ax.text(center, height + 0.02 * height, label, ha="center", va="bottom", fontsize=9)

    g.set_axis_labels("Sample", "Insertion Count")
    g.add_legend(title="Frame")
    g.set_titles("Gene Class: {col_name}")
    g.fig.suptitle("Insertions per Frame by Gene Class (Annotated Genes)", fontsize=16)
    plt.subplots_adjust(top=0.85)
    os.makedirs(output_path, exist_ok=True)
    out_path = os.path.join(output_path, "ncbi_insertions_per_frame.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print("✅ Saved:", out_path)

def plot_insertions_per_frame(csv_path, title, save_path):
    df = pd.read_csv(csv_path)
    df['frame'] = df['frame'].astype(int)

    agg_df = df.groupby(['sample_type', 'frame'])['I_value'].sum().reset_index()
    sample_types = sorted(agg_df['sample_type'].unique())
    frames = sorted(agg_df['frame'].unique())

    bar_width = 0.25
    x = np.arange(len(sample_types))
    offsets = [-bar_width, 0, bar_width]

    fig, ax = plt.subplots(figsize=(10, 6))

    for i, frame in enumerate(frames):
        frame_data = agg_df[agg_df['frame'] == frame]
        means = frame_data.set_index('sample_type').reindex(sample_types)['I_value'].fillna(0).values
        bars = ax.bar(x + offsets[i], means, width=bar_width, label=f'Frame {frame}', color=frame_palette[frame])

        for j, bar in enumerate(bars):
            sample = sample_types[j]
            percent_row = df[(df['sample_type'] == sample) & (df['frame'] == frame)]
            if not percent_row.empty:
                percent_val = percent_row['percentage'].values[0]
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.02 * bar.get_height(),
                        f'{int(round(bar.get_height()))}\n({percent_val:.1f}%)',
                        ha='center', va='bottom', fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(sample_types)
    ax.set_xlabel('Sample')
    ax.set_ylabel('Insertion Count')
    ax.set_title(title)
    ax.legend(title='Frame')
    ax.grid(False)
    ax.set_ylim(0, agg_df['I_value'].max() * 1.15)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print("📊 Saved plot:", save_path)
    plt.close()



# -------------------------
# HELPER FUNCTIONS
# -------------------------
def get_inframe_insertions_from_df(df):
    df_i = df[(df['metric'] == 'I') & (df['frame'] == 1)]
    return df_i['value'].sum()

def get_inframe_reads_from_df(df):
    df_r = df[(df['metric'] == 'R') & (df['frame'] == 1)]
    if df_r.empty:
        return 0
    return df_r['value'].sum()

def get_inframe_coverage_from_df(df):
    df_cov = df[(df['metric'] == 'coverage') & (df['frame'] == 1)]
    if df_cov.empty:
        return 0.0
    return float(df_cov['value'].iloc[0])

def orf_poisson_pvalue(total_insertions, lambda_bg, bp_len, mean_intergenic_len):
    expected = lambda_bg * (bp_len / mean_intergenic_len)
    if expected <= 0:
        return 1.0 if total_insertions == 0 else 0.0
    return poisson.sf(int(total_insertions) - 1, expected)

def insertion_threshold(lambda_bg, bp_len, mean_intergenic_len, pval_thresh):
    expected = lambda_bg * (bp_len / mean_intergenic_len)
    if expected <= 0:
        return 0
    k = 0
    max_iter = 10000
    while poisson.sf(k - 1, expected) > pval_thresh:
        k += 1
        if k > max_iter:
            print(f"Warning: insertion_threshold loop reached max iterations (expected={expected})")
            break
    return k

def safe_mean_p(p_list):
    p_array = np.array(p_list, dtype=float)
    p_array[p_array <= 0] = 1e-300
    return np.mean(p_array)

# -------------------------
# BACKGROUND ESTIMATION
# -------------------------
def estimate_background_lambda(sample_reads, igan_dict, read_thr):
    counts = []
    lengths = []
    for ig_id, coords in igan_dict.items():
        df = pst.get_f_metric(sample_reads, coords, read_thr)
        cov = get_inframe_coverage_from_df(df)
        if np.log2(cov + 1e-9) < 2:
            ins = get_inframe_insertions_from_df(df)
            counts.append(ins)
            lengths.append(coords[1] - coords[0] + 1)
    lambda_bg = np.mean(counts) if counts else 0.0
    mean_len = np.mean(lengths) if lengths else 1.0
    return lambda_bg, mean_len

# -------------------------
# ORF PREDICTION
# -------------------------
def predict_orfs(all_reads=None, ncbi=None, igan=None, orfs=None, cfg=None, out_path=None):
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from collections import defaultdict
    from sklearn.metrics import roc_curve, auc

    if cfg is None:
        from core.core_p import load_config
        cfg = load_config()

    # Set output path
    if out_path is None:
        out_path = cfg['output']['protinseq']['predictions']
    os.makedirs(out_path, exist_ok=True)

    READ_THRESHOLD = cfg['settings']['read_threshold']
    PVALUE_THRESHOLD = cfg['settings']['pvalue_threshold']
    MIN_REPRO_SAMPLES = cfg['settings']['min_repro_samples']
    SAMPLES_TO_ANALYZE = cfg['settings']['predict_samples']

    import core.protinseq_tools as pst
    processed_dir = cfg['input']['processed']

    if all_reads is None:
        all_reads_path = cfg['output']['anubis']['all_reads']
        all_reads = pst.LoadPickle(all_reads_path)
    if ncbi is None or igan is None or orfs is None:
        annotations = pst.LoadPickle(os.path.join(processed_dir, 'annotations.dat'))
        orfs = orfs or annotations[0]
        ncbi = ncbi or annotations[1]
        igan = igan or pst.LoadPickle(os.path.join(processed_dir, 'intergenic.dat'))

    # -------------------------
    # CLIP REGIONS TO SAMPLE LENGTH (single sample)
    # -------------------------
    first_sample = next(iter(all_reads.values()))
    max_len = len(first_sample)
    print(f"Clipping all regions to sample length {max_len}")
    for region_dict_name, region_dict in zip(["ORF", "IGAN", "NCBI"], [orfs, igan, ncbi]):
        for rid, coords in region_dict.items():
            start, end, *rest = coords
            if end >= max_len:
                coords = (start, max_len - 1, *rest)
                region_dict[rid] = coords

    # -------------------------
    # BACKGROUND ESTIMATION
    # -------------------------
    sample_bg_cache = {s: estimate_background_lambda(all_reads[s], igan, READ_THRESHOLD)
                       for s in SAMPLES_TO_ANALYZE}

    # -------------------------
    # ORF P-VALUES AND OBSERVED COUNTS
    # -------------------------
    region_pvals = defaultdict(dict)
    region_obs_counts = defaultdict(dict)

    regions = {f"NCBI::{k}": v for k,v in ncbi.items()}
    regions.update({f"IGAN::{k}": v for k,v in igan.items()})
    regions.update({f"ORF::{k}": v for k,v in orfs.items()})

    for sample in SAMPLES_TO_ANALYZE:
        lambda_bg, mean_len = sample_bg_cache[sample]
        sample_reads = all_reads[sample]
        for rid, coords in regions.items():
            df = pst.get_f_metric(sample_reads, coords, READ_THRESHOLD)
            obs_ins = get_inframe_insertions_from_df(df)
            bp_len = coords[1] - coords[0] + 1
            p = orf_poisson_pvalue(obs_ins, lambda_bg, bp_len, mean_len)
            region_pvals[rid][sample] = p
            region_obs_counts[rid][sample] = obs_ins

    # -------------------------
    # ROC / Youden threshold
    # -------------------------
    ncbi_ids = [r for r in region_pvals.keys() if r.startswith("NCBI::")]
    igan_ids = [r for r in region_pvals.keys() if r.startswith("IGAN::")]

    pvals_for_roc = [-np.log10(safe_mean_p(list(region_pvals[r].values()))) for r in ncbi_ids + igan_ids]
    labels_for_roc = [1]*len(ncbi_ids) + [0]*len(igan_ids)

    fpr, tpr, thresholds = roc_curve(labels_for_roc, pvals_for_roc, pos_label=1)
    best_threshold_score = thresholds[np.argmax(tpr - fpr)]
    roc_auc_score = auc(fpr, tpr)

    avg_lambda_bg = np.mean([sample_bg_cache[s][0] for s in SAMPLES_TO_ANALYZE])
    avg_mean_len = np.mean([sample_bg_cache[s][1] for s in SAMPLES_TO_ANALYZE])

    # -------------------------
    # ORF-level predictions
    # -------------------------
    orf_results = []
    for rid in [r for r in region_pvals.keys() if r.startswith("ORF::")]:
        gene = rid.split("::",1)[1]
        p_list = list(region_pvals[rid].values())
        mean_p = safe_mean_p(p_list)
        avg_score = -np.log10(mean_p)
        sig_count = sum(1 for p in p_list if p < PVALUE_THRESHOLD)
        passes_repro = sig_count >= MIN_REPRO_SAMPLES
        passes_score = avg_score >= best_threshold_score
        total_insertions = sum(region_obs_counts[rid].values())
        bp_len = regions[rid][1] - regions[rid][0] + 1
        ins_threshold = insertion_threshold(avg_lambda_bg, bp_len, avg_mean_len, PVALUE_THRESHOLD)
        orf_results.append({
            "orf": gene,
            "mean_p": mean_p,
            "avg_score": avg_score,
            "sig_count": sig_count,
            "passes_repro": passes_repro,
            "passes_score": passes_score,
            "final_call": passes_repro and passes_score,
            "total_insertions": total_insertions,
            "insertion_threshold": ins_threshold
        })

    df_orfs = pd.DataFrame(orf_results).sort_values(['final_call','avg_score'], ascending=[False, False])
    df_orfs.to_csv(os.path.join(out_path, "all_orf_results.csv"), index=False)

    # Save only ORFs called as real
    df_orfs[df_orfs['final_call']].to_csv(os.path.join(out_path, "predicted_real_orfs.csv"), index=False)

    # -------------------------
    # Per-sample table
    # -------------------------
    per_sample_rows = []
    for rid, coords in regions.items():
        start, end, strand = coords
        ntlen = end - start + 1
        aalen = ntlen // 3
        ann_type = rid.split("::")[0]

        row_data = {
            "gene": rid.split("::")[1],
            "start": start,
            "end": end,
            "strand": strand,
            "ntlen": ntlen,
            "aalen": aalen,
            "ann_type": ann_type
        }

        sig_count = 0
        for sample in SAMPLES_TO_ANALYZE:
            sample_reads = all_reads[sample]
            lambda_bg, mean_len = sample_bg_cache[sample]

            df_metrics = pst.get_f_metric(sample_reads, coords, READ_THRESHOLD)
            I = get_inframe_insertions_from_df(df_metrics)
            R = get_inframe_reads_from_df(df_metrics)
            sfNC = orf_poisson_pvalue(I, lambda_bg, ntlen, mean_len)
            pred = 1 if sfNC <= PVALUE_THRESHOLD else 0
            if pred == 1:
                sig_count += 1

            row_data[f"{sample}_I"] = I
            row_data[f"{sample}_R"] = R
            row_data[f"{sample}_rNC"] = lambda_bg / mean_len if mean_len > 0 else 0
            row_data[f"{sample}_sfNC"] = sfNC
            row_data[f"{sample}_pred"] = pred

        row_data["samples_significant"] = sig_count
        per_sample_rows.append(row_data)

    df_per_sample = pd.DataFrame(per_sample_rows)
    df_per_sample.to_csv(os.path.join(out_path, "per_sample_region_stats.csv"), index=False)

    # -------------------------
    # Save ROC figure
    # -------------------------
    plt.figure(figsize=(6,6))
    plt.plot(fpr, tpr, label=f"AUC={roc_auc_score:.3f}")
    best_idx = np.argmax(tpr - fpr)
    plt.scatter(fpr[best_idx], tpr[best_idx], color='red', label=f"Youden (score={best_threshold_score:.3f})")
    plt.plot([0,1],[0,1],'k--')
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.title('ROC (NCBI positives vs IGAN negatives)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, "roc_ncb_vs_igan.png"), dpi=300)
    plt.close()

    # -------------------------
    # Log summary
    # -------------------------
    n_orfs_total = len(df_orfs)
    n_orfs_called = df_orfs['final_call'].sum()
    log_lines = [
        f"Total ORFs analyzed: {n_orfs_total}",
        f"ORFs called as real: {n_orfs_called}",
        "",
        f"ROC AUC: {roc_auc_score:.6f}",
        f"Best score threshold (Youden): {best_threshold_score:.6f}",
        f"TPR at best: {tpr[best_idx]:.6f}",
        f"FPR at best: {fpr[best_idx]:.6f}",
        "",
        "Top 10 ORFs by avg_score:",
        df_orfs.head(10).to_string(index=False)
    ]

    log_file = os.path.join(out_path, "log.txt")
    with open(log_file, "w") as f:
        f.write("\n".join(log_lines))

    print("\n".join(log_lines))
    print(f"\n✅ All outputs saved in: {out_path}")

    return df_orfs, df_per_sample, fpr, tpr, best_threshold_score

# -------------------------
# MAIN TABLE GENERATION
# -------------------------
def generate_all_tables(excel_path=None, sheet_name=None, model_name=None, out_path=None, config_path="config.yaml", predict=False):
    cfg = load_config(config_path)
    processed_dir = cfg['input']['processed']
    charts_dir = cfg['output']['protinseq']['charts']
    plots_dir = cfg['output']['protinseq']['plots']
    all_reads_path = cfg['output']['anubis']['all_reads']

    if not os.path.exists(all_reads_path):
        print(f"{all_reads_path} not found. Running ANUBIS pipeline...")
        from core.core_a import run_anubis_pipeline
        run_anubis_pipeline(cfg)
        if not os.path.exists(all_reads_path):
            raise FileNotFoundError(f"Failed to generate {all_reads_path}.")

    print("Loading required files...")
    all_reads = pst.LoadPickle(all_reads_path)
    orfs, ncbi, rna, myco2mg, mg2myco = pst.LoadPickle(os.path.join(processed_dir, 'annotations.dat'))
    igan = pst.LoadPickle(os.path.join(processed_dir, 'intergenic.dat'))
    aa_seqs, nt_seqs = pst.LoadPickle(os.path.join(processed_dir, 'sequences.dat'))
    myclone = pd.read_csv(os.path.join(processed_dir, 'myclone.tsv'), sep='\t', index_col=0)

    if predict:
        df_orfs, df_per_sample, fpr, tpr, best_threshold_score = predict_orfs(all_reads, ncbi, igan, orfs, cfg)
        os.makedirs(out_path, exist_ok=True)
        df_orfs.to_csv(os.path.join(out_path, "protinseq_orf_predictions.csv"), index=False)
        df_per_sample.to_csv(os.path.join(out_path, "per_sample_region_stats.csv"), index=False)

        # Save ROC figure
        plt.figure(figsize=(6,6))
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, label=f"AUC={roc_auc:.3f}")
        best_idx = np.argmax(tpr - fpr)
        plt.scatter(fpr[best_idx], tpr[best_idx], color='red', label=f"Youden (score={best_threshold_score:.3f})")
        plt.plot([0,1],[0,1],'k--')
        plt.xlabel('FPR')
        plt.ylabel('TPR')
        plt.title('ROC (NCBI positives vs IGAN negatives)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(out_path, "roc_ncb_vs_igan.png"), dpi=300)
        plt.close()
        print(f"✅ ORF predictions and per-sample stats saved to {out_path}")

    else:
        # Non-predict mode: generate the regular insertion tables and plots
        gold_df = load_gold_from_excel(excel_path, sheet_name, model_name)
        summarize_insertions_by_class(all_reads, ncbi, gold_df, pst, out_path=charts_dir, selected_samples=cfg['settings']['protinseq_samples'])
        summarize_insertions_per_frame(all_reads, igan, pst, label="igan", out_path=charts_dir, selected_samples=cfg['settings']['protinseq_samples'])
        summarize_insertions_per_frame(all_reads, orfs, pst, label="orfs", out_path=charts_dir, selected_samples=cfg['settings']['protinseq_samples'])

        plot_insertions_per_frame(os.path.join(charts_dir, "igan_insertions_per_frame.csv"),
                                  "Insertions per Frame (Intergenic Regions)",
                                  os.path.join(plots_dir, "igan_insertions_per_frame.png"))
        plot_insertions_per_frame(os.path.join(charts_dir, "orfs_insertions_per_frame.csv"),
                                  "Insertions per Frame (RanSEPs ORFs)",
                                  os.path.join(plots_dir, "orfs_insertions_per_frame.png"))
        plot_ncbi_by_class(os.path.join(charts_dir, "insertions_by_class_ncbi.csv"),
                           "Insertions by Class (NCBI Annotated Genes)",
                           plots_dir)