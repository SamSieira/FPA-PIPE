import os
import glob
import pickle
import logging
import pandas as pd
from Bio import SeqIO
import csv
import yaml

# ─────────────────────────────
# 📁 Config + Logging Setup
# ─────────────────────────────
def load_config(config_path="config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

config = load_config()
RAW_DIR = config["input"]["raw"]
PROCESSED_DIR = config["input"]["processed"]
os.makedirs(PROCESSED_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
)


# ─────────────────────────────
# 🔎 File Helper
# ─────────────────────────────
def find_file(ext, contains=None):
    files = glob.glob(os.path.join(RAW_DIR, f"*.{ext}"))
    if contains:
        files = [f for f in files if contains in os.path.basename(f)]
    if len(files) != 1:
        raise ValueError(f"Expected one .{ext} file containing '{contains}', found {len(files)}: {files}")
    return files[0]


# ─────────────────────────────
# 🧬 Sequence Data
# ─────────────────────────────
def generate_sequences():
    logging.info("📦 Generating sequences.dat...")
    csv_file = find_file("csv", contains="predicted")
    df = pd.read_csv(csv_file)

    aa_seqs = dict(zip(df["identifier"], df["aa_seq"]))
    nt_seqs = dict(zip(df["identifier"], df["nt_seq"]))

    out_file = os.path.join(PROCESSED_DIR, "sequences.dat")
    with open(out_file, "wb") as f:
        pickle.dump([aa_seqs, nt_seqs], f)
    logging.info("✅ Saved sequences.dat")


# ─────────────────────────────
# 🧬 Annotations
# ─────────────────────────────
def generate_annotations():
    logging.info("🧬 Generating annotations.dat...")
    csv_file = find_file("csv", contains="predicted")
    gb_file = find_file("gb")

    orfs, myco2mg, mg2myco = {}, {}, {}

    with open(csv_file, newline='') as f:
        for row in csv.DictReader(f):
            identifier = row["identifier"]
            start, stop = int(row["start"]), int(row["stop"])
            strand = row["strand"]
            annotation = row["annotation"]

            orfs[identifier] = [start, stop, strand]
            if annotation.startswith("MG"):
                myco2mg[identifier] = annotation
                mg2myco[annotation] = identifier

    ncbi, rna = {}, {}
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if not feature.qualifiers.get("locus_tag"):
                continue
            locus_tag = feature.qualifiers["locus_tag"][0]
            start = int(feature.location.start) + 1  # Convert to 1-based
            stop = int(feature.location.end)         # End is already exclusive, so no need to subtract
            strand = "+" if feature.location.strand == 1 else "-"

            if feature.type == "CDS":
                ncbi[locus_tag] = [start, stop, strand]
            elif feature.type in ["tRNA", "rRNA", "ncRNA", "RNA"]:
                rna[locus_tag] = [start, stop, strand]

    with open(os.path.join(PROCESSED_DIR, "annotations.dat"), "wb") as f:
        pickle.dump([orfs, ncbi, rna, myco2mg, mg2myco], f)

    logging.info("✅ Saved annotations.dat")


# ─────────────────────────────
# 🧬 Gene Regions
# ─────────────────────────────
def generate_gene_regions():
    logging.info("📊 Generating gene_regions.tsv...")
    with open(os.path.join(PROCESSED_DIR, "annotations.dat"), "rb") as f:
        data = pickle.load(f)

    rows = []
    for d in data[:2]:  # orfs and ncbi
        for gene, coords in sorted(d.items(), key=lambda x: x[1][0]):
            if isinstance(coords, list) and len(coords) == 3:
                start, end, strand = coords
                rows.append([gene, start, end, strand])

    out_path = os.path.join(PROCESSED_DIR, "gene_regions.tsv")
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gene", "start", "end", "strand"])
        writer.writerows(rows)
    logging.info("✅ Saved gene_regions.tsv")


# ─────────────────────────────
# 🧬 Intergenic Regions
# ─────────────────────────────
def generate_intergenic_regions():
    logging.info("🧬 Generating intergenic.dat...")
    regions = []
    with open(os.path.join(PROCESSED_DIR, "gene_regions.tsv"), "r") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            regions.append({
                "gene": row["gene"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": row["strand"]
            })

    positive = sorted([r for r in regions if r["strand"] == "+"], key=lambda x: x["start"])
    negative = sorted([r for r in regions if r["strand"] == "-"], key=lambda x: x["start"])

    def find_igs(regions, prefix):
        igs, prev_end, count = {}, 0, 1
        for region in regions:
            if region["start"] > prev_end + 1:
                ig_start = prev_end + 1
                ig_end = region["start"] - 1
                if (ig_end - ig_start) >= 5:
                    igs[f"{prefix}{count:04d}"] = [ig_start, ig_end, region["strand"]]
                    count += 1
            prev_end = max(prev_end, region["end"])
        return igs

    igs_all = {**find_igs(positive, "IGMP"), **find_igs(negative, "IGPM")}
    with open(os.path.join(PROCESSED_DIR, "intergenic.dat"), "wb") as f:
        pickle.dump(igs_all, f)
    logging.info("✅ Saved intergenic.dat")


# ─────────────────────────────
# 🧬 Myclone Table
# ─────────────────────────────
def generate_numbered_clone_table():
    logging.info("📋 Generating myclone.tsv...")
    rows = []

    with open(os.path.join(PROCESSED_DIR, "gene_regions.tsv"), newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append([row["gene"], row["start"], row["end"], row["strand"]])

    with open(os.path.join(PROCESSED_DIR, "intergenic.dat"), "rb") as f:
        igs = pickle.load(f)
        for gene, (start, end, strand) in igs.items():
            rows.append([gene, str(start), str(end), strand])

    rows = [[i+1] + row for i, row in enumerate(rows)]
    with open(os.path.join(PROCESSED_DIR, "myclone.tsv"), "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["row", "gene", "start", "end", "strand"])
        writer.writerows(rows)

    logging.info("✅ Saved myclone.tsv")


# ─────────────────────────────
# 🚀 Pipeline Orchestration
# ─────────────────────────────
def preprocess():
    logging.info("🚀 Starting data preprocessing pipeline...")
    generate_sequences()
    generate_annotations()
    generate_gene_regions()
    generate_intergenic_regions()
    generate_numbered_clone_table()
    logging.info("🎉 Preprocessing complete.")


# ─────────────────────────────
# 🖥️ CLI
# ─────────────────────────────
if __name__ == "__main__":
    preprocess()
