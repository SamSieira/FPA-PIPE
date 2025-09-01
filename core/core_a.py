import os
import logging
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import anubis
from tqdm import tqdm
from intervaltree import IntervalTree


# === Helper for saving plots ===
def saveplot(plot_dir, sample_name, model_name):
    os.makedirs(plot_dir, exist_ok=True)
    return os.path.join(plot_dir, f"{sample_name}_{model_name}.png")


def process_qins_file(sample_name, qins_path, genome_file, goldset_file,
                      filter_genes=False, filter_tails=False, filter_reads=False,
                      export_reads_dict=None, intermediate_dir=None, longform_csv_path=None):
    if not os.path.exists(qins_path):
        raise FileNotFoundError(f"File not found: {qins_path}")

    d = anubis.dataset(qins_path, genome=genome_file, goldset=goldset_file)
    d.annotation = {k: v for k, v in d.annotation.items() if len(k) == 10}

    if filter_tails:
        d = d.filter(filter_by='tails', lpercentile=15, rpercentile=90)
        logging.info(f"{sample_name}: Coverage after tails filtering: {np.sum(d.coverage)}")

    if filter_genes:
        d = d.filter(filter_by='genes', inplace=False)
        logging.info(f"{sample_name}: Coverage after gene filtering: {np.sum(d.coverage)}")

    if filter_reads:
        d = d.filter(filter_by='reads', min_reads=5, max_reads=1000)
        logging.info(f"{sample_name}: Coverage after reads filtering: {np.sum(d.coverage)}")

    d.cqn(condition='GC', bins=6, normalization='minmax', annotation='all')
    d.metric_by()

    if intermediate_dir:
        os.makedirs(intermediate_dir, exist_ok=True)
        metric_path = os.path.join(intermediate_dir, f"{sample_name}_filtered_metrics.csv")
        d.metric.to_csv(metric_path)
        logging.info(f"{sample_name}: Saved metrics to {metric_path}")

    if export_reads_dict is not None:
        export_reads_dict[sample_name] = d.zreads.astype(int)

    if longform_csv_path:
        df = pd.DataFrame({
            'position': range(1, d.genome_length + 1),
            'reads': d.zreads.astype(int),
            'sample': sample_name
        })
        df.to_csv(longform_csv_path, mode='a', index=False, header=not os.path.exists(longform_csv_path))
        logging.info(f"{sample_name}: Appended reads to longform csv")

    return d


def run_anubis_pipeline(config, filter_genes=False, filter_tails=False, filter_reads=False):
    qins_dir = config['input']['raw']
    output_dir = os.path.join(config['output']['anubis']['plots'], '..')
    os.makedirs(output_dir, exist_ok=True)

    output_excel = os.path.join(output_dir, "results.xlsx")
    intermediate_dir = os.path.join(output_dir, "metrics")
    plot_dir = config['output']['anubis']['plots']
    longform_csv_path = os.path.join(output_dir, "reads_longform.csv")
    os.makedirs(plot_dir, exist_ok=True)

    genome_files = [f for f in os.listdir(qins_dir) if f.endswith('.gb')]
    if not genome_files:
        raise FileNotFoundError("No .gb genome file found in raw input directory.")
    genome_file = os.path.join(qins_dir, genome_files[0])

    goldset_file = os.path.join(qins_dir, "goldset.csv")
    if not os.path.exists(goldset_file):
        raise FileNotFoundError("goldset.csv not found in input directory.")

    # Load sample list from config or auto-discover
    if config.get('settings', {}).get('use_sample_prefixes', False):
        base_files = config['settings']['sample_prefixes']
        logging.info(f"🧾 Using sample list from config.yaml: {base_files}")
    else:
        base_files = sorted([
            f.replace(".qins", "") for f in os.listdir(qins_dir)
            if f.endswith('.qins') and not f.endswith('_fw.qins') and not f.endswith('_rv.qins')
        ])
        logging.info(f"🔍 Auto-discovered samples: {base_files}")

    all_results = {}
    all_reads_dict = {}

    for sample_name in base_files:
        try:
            df = process_sample(
                sample_name=sample_name,
                qins_dir=qins_dir,
                genome_file=genome_file,
                goldset_file=goldset_file,
                plot_dir=plot_dir,
                filter_genes=filter_genes,
                filter_tails=filter_tails,
                filter_reads=filter_reads,
                export_reads_dict_dir=None,
                intermediate_dir=intermediate_dir,
                longform_csv_path=longform_csv_path,
                shared_reads_dict=all_reads_dict
            )
            all_results[sample_name] = df
        except Exception as e:
            logging.error(f"❌ Failed to process {sample_name}: {e}")

    if all_results:
        with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
            for sample_name, df in all_results.items():
                df.to_excel(writer, sheet_name=sample_name[:31], index=False)
        logging.info(f"✅ Results saved to {output_excel}")
    else:
        logging.warning("⚠️ No data processed; no Excel file created.")

    # Save all_reads_dict to all_reads.dat
    all_reads_dat_path = os.path.join(output_dir, "all_reads.dat")
    with open(all_reads_dat_path, "wb") as f:
        pickle.dump(all_reads_dict, f)
    logging.info(f"📦 All sample reads saved to {all_reads_dat_path}")


def process_sample(sample_name, qins_dir, genome_file, goldset_file, plot_dir,
                   filter_genes=False, filter_tails=False, filter_reads=False,
                   export_reads_dict_dir=None, intermediate_dir=None,
                   longform_csv_path=None, shared_reads_dict=None):
    base_path = os.path.join(qins_dir, f"{sample_name}.qins")
    fw_path = os.path.join(qins_dir, f"{sample_name}_fw.qins")
    rv_path = os.path.join(qins_dir, f"{sample_name}_rv.qins")

    logging.info(f"🔬 Processing sample: {sample_name}")

    d = anubis.dataset(base_path, genome=genome_file, goldset=goldset_file)
    d.annotation = {k: v for k, v in d.annotation.items() if len(k) == 10}

    total_reads_start = np.sum(d.zreads)
    total_insertions_start = np.count_nonzero(d.zreads)

    # Apply filtering and normalization to forward and reverse strands
    if shared_reads_dict is not None:
        # Process forward strand
        if os.path.exists(fw_path):
            fw_d = anubis.dataset(fw_path, genome=genome_file, goldset=goldset_file)
            if filter_tails:
                fw_d = fw_d.filter(filter_by='tails', lpercentile=15, rpercentile=90)
            if filter_reads:
                fw_d = fw_d.filter(filter_by='reads', min_reads=5, max_reads=1000)
            fw_d.cqn(condition='GC', bins=6, normalization='minmax', annotation='all', show_plot=False)
            shared_reads_dict[f"{sample_name}_pos"] = fw_d.zreads.astype(int)
            logging.info(f"{sample_name}_fw: Filtered and normalized reads stored as {sample_name}_pos")

        # Process reverse strand
        if os.path.exists(rv_path):
            rv_d = anubis.dataset(rv_path, genome=genome_file, goldset=goldset_file)
            if filter_tails:
                rv_d = rv_d.filter(filter_by='tails', lpercentile=15, rpercentile=90)
            if filter_reads:
                rv_d = rv_d.filter(filter_by='reads', min_reads=5, max_reads=1000)
            rv_d.cqn(condition='GC', bins=6, normalization='minmax', annotation='all', show_plot=False)
            shared_reads_dict[f"{sample_name}_neg"] = rv_d.zreads.astype(int)
            logging.info(f"{sample_name}_rv: Filtered and normalized reads stored as {sample_name}_neg")

    corr = d.correlation_by_distance(min_relative_distance=-20, max_relative_distance=20)
    positive_corr = {k: v for k, v in corr.items() if k > 0 and not np.isnan(v)}

    if positive_corr:
        best_tsd = max(positive_corr, key=positive_corr.get)
        max_corr_value = positive_corr[best_tsd]
        logging.info(f"{sample_name}: Detected TSD size: {best_tsd} with correlation {max_corr_value:.4f}")

        if max_corr_value >= 0.5:
            d = d.correct_TSD(fw=fw_path, rv=rv_path, tsd_size=best_tsd)
            logging.info(f"{sample_name}: Coverage after TSD correction: {np.sum(d.coverage)}")
        else:
            logging.info(f"{sample_name}: TSD correlation too low; skipping TSD correction")
    else:
        logging.info(f"{sample_name}: No positive correlations found; skipping TSD correction")

    d.set_goldset(goldset_file)
    d.annotation = {k: v for k, v in d.annotation.items() if len(k) == 10}

    def log_filter_stats(label):
        reads = np.sum(d.zreads)
        insertions = np.count_nonzero(d.zreads)
        logging.info(f"{sample_name}: {label} → Reads: {reads} | Insertions: {insertions}")
        return reads, insertions

    if filter_tails:
        d = d.filter(filter_by='tails', lpercentile=15, rpercentile=90)
        total_reads_start, total_insertions_start = log_filter_stats("After tails filtering")

    if filter_genes:
        d = d.filter(filter_by='genes', inplace=False)
        total_reads_start, total_insertions_start = log_filter_stats("After gene filtering")

    if filter_reads:
        d = d.filter(filter_by='reads', min_reads=5, max_reads=1000)
        total_reads_start, total_insertions_start = log_filter_stats("After reads filtering")

    d.cqn(condition='GC', bins=6, normalization='minmax', annotation='all', show_plot=False)
    logging.info(f"{sample_name}: Coverage after GC normalization: {np.sum(d.coverage)}")
    d.metric_by()

    if intermediate_dir:
        os.makedirs(intermediate_dir, exist_ok=True)
        metric_path = os.path.join(intermediate_dir, f"{sample_name}_filtered_metrics.csv")
        d.metric.to_csv(metric_path)
        logging.info(f"{sample_name}: Saved metrics to {metric_path}")

    if longform_csv_path:
        df = pd.DataFrame({
            'position': range(1, d.genome_length + 1),
            'reads': d.zreads.astype(int),
            'sample': sample_name
        })
        df.to_csv(longform_csv_path, mode='a', index=False, header=not os.path.exists(longform_csv_path))
        logging.info(f"{sample_name}: Appended reads to longform csv")

    return pd.concat([
        anubis.Poisson(d).predict(criteria='absolute', show_plot=False, save_plot=saveplot(plot_dir, sample_name, "Poisson_absolute")).assign(model='Poisson_absolute'),
        anubis.Gumbel(d).predict(criteria='threshold', thr=0.01, prior={'E': 'goldset', 'N': 'intergenic'}, show_plot=False, save_plot=saveplot(plot_dir, sample_name, "Gumbel_threshold")).assign(model='Gumbel_threshold'),
        anubis.Gumbel(d).predict(criteria='foldchange', foldchange=3, prior={'E': '0', 'N': 'intergenic'}, show_plot=False, save_plot=saveplot(plot_dir, sample_name, "Gumbel_foldchange")).assign(model='Gumbel_foldchange'),
        anubis.Gamma(d).predict(show_plot=False, save_plot=saveplot(plot_dir, sample_name, "Gamma_default")).assign(model='Gamma_default'),
        anubis.GaussianMixture(d).predict(show_plot=False, save_plot=saveplot(plot_dir, sample_name, "GMM_default")).assign(model='GMM_default'),
    ], ignore_index=True)


def annotate_genome_with_insertions(cfg):
    import yaml
    from intervaltree import IntervalTree

    if isinstance(cfg, str):
        with open(cfg, "r") as f:
            cfg = yaml.safe_load(f)

    input_file = cfg['input']['gene_regions']
    genome_size = cfg['settings'].get('genome_size', 580076)
    sample_prefixes = cfg['settings']['sample_prefixes']
    all_reads_path = cfg['output']['anubis']['all_reads']
    output_dir = cfg['output']['anubis']['charts']
    output_file = os.path.join(output_dir, "bp_annotations.csv")
    os.makedirs(output_dir, exist_ok=True)

    annotations = pd.read_csv(input_file, sep="\t", header=None, names=["gene", "start", "end", "strand"])
    annotations["start"] = pd.to_numeric(annotations["start"], errors="coerce")
    annotations["end"] = pd.to_numeric(annotations["end"], errors="coerce")
    annotations = annotations.dropna().copy()
    annotations["start"] = annotations["start"].astype(int)
    annotations["end"] = annotations["end"].astype(int)

    plus_tree = IntervalTree()
    minus_tree = IntervalTree()
    for _, row in annotations.iterrows():
        interval = (row["start"], row["end"] + 1)
        if row["strand"] == "+":
            plus_tree.addi(*interval, row["gene"])
        elif row["strand"] == "-":
            minus_tree.addi(*interval, row["gene"])

    gene_info = {row["gene"]: row for _, row in annotations.iterrows()}

    with open(all_reads_path, "rb") as f:
        all_reads = pickle.load(f)

    output = []
    for pos in tqdm(range(1, genome_size + 1), desc="Annotating genome"):
        row = {"pos": pos}

        # Positive strand
        hits_plus = plus_tree[pos]
        inframe_genes_pos = []
        labels_pos = []
        for hit in hits_plus:
            gene = hit.data
            gene_start = gene_info[gene]["start"]
            if (pos - gene_start) % 3 == 0:
                inframe_genes_pos.append(gene)
                label = "annotated" if gene.startswith("MG_RS") else "putative"
                labels_pos.append(label)
        row["in-frame positive"] = ";".join(inframe_genes_pos) if inframe_genes_pos else "NA"
        row["label_pos"] = ";".join(sorted(set(labels_pos))) if labels_pos else "non-coding"
        for sample in sample_prefixes:
            array_name = f"{sample}_pos"
            insert_array = all_reads.get(array_name, [])
            row[sample + "_pos"] = insert_array[pos - 1] if (pos - 1) < len(insert_array) else 0

        # Negative strand
        hits_minus = minus_tree[pos]
        inframe_genes_neg = []
        labels_neg = []
        for hit in hits_minus:
            gene = hit.data
            gene_end = gene_info[gene]["end"]
            if (gene_end - pos) % 3 == 0:
                inframe_genes_neg.append(gene)
                label = "annotated" if gene.startswith("MG_RS") else "putative"
                labels_neg.append(label)
        row["in-frame negative"] = ";".join(inframe_genes_neg) if inframe_genes_neg else "NA"
        row["label_neg"] = ";".join(sorted(set(labels_neg))) if labels_neg else "non-coding"
        for sample in sample_prefixes:
            array_name = f"{sample}_neg"
            insert_array = all_reads.get(array_name, [])
            row[sample + "_neg"] = insert_array[pos - 1] if (pos - 1) < len(insert_array) else 0

        output.append(row)

    pd.DataFrame(output).to_csv(output_file, index=False)
    print(f"✅ Annotated genome saved to: {output_file}")
