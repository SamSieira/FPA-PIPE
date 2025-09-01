import os
import argparse
import logging
import yaml

from core import data_processor as data_generation
from core.core_p import generate_all_tables, predict_orfs
from core.core_a import run_anubis_pipeline  
from core.core_a import annotate_genome_with_insertions

# ─────────────────────────────
# Config Loading
# ─────────────────────────────
def load_config(config_path="config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

# ─────────────────────────────
# Logging Setup
# ─────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    handlers=[
        logging.FileHandler("FPA_pipeline.log"),
        logging.StreamHandler()
    ]
)

# ─────────────────────────────
# Subcommand Handlers
# ─────────────────────────────
def run_anubis(args, config):
    run_anubis_pipeline(
        config,
        filter_genes=args.filter_genes,
        filter_tails=args.filter_tails,
        filter_reads=args.filter_reads
    )
    annotate_genome_with_insertions(config)


def run_protinseq(args, config):
    os.makedirs(args.output, exist_ok=True)
    if args.predict:
        logging.info("Running new Poisson ORF prediction (protinseq)...")
        # everything is now handled internally in predict_orfs
        predict_orfs(cfg=config, out_path=args.output)
    else:
        logging.info(f"Generating Protinseq tables and plots using sample {args.sample}, model {args.model}")
        generate_all_tables(
            excel_path=os.path.join("outputs", "anubis", "results.xlsx"),
            sheet_name=args.sample,
            model_name=args.model,
            out_path=args.output,
            config_path=None  # uses default config inside generate_all_tables
        )
        logging.info("Protinseq generation complete.")

def run_preprocessing():
    logging.info("Starting preprocessing...")
    data_generation.preprocess()
    logging.info("Preprocessing complete.")

# ─────────────────────────────
#  Main Entry Point
# ─────────────────────────────
def main():
    config = load_config()

    parser = argparse.ArgumentParser(description="FPA_pipeline CLI with subcommands")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ANUBIS
    anubis_parser = subparsers.add_parser("anubis", help="Run ANUBIS filtering and modeling pipeline")
    anubis_parser.add_argument("--filter_genes", action="store_true", help="Enable gene filtering")
    anubis_parser.add_argument("--filter_tails", action="store_true", help="Enable tails filtering")
    anubis_parser.add_argument("--filter_reads", action="store_true", help="Enable read count filtering")

    # Preprocessing
    subparsers.add_parser("preprocess", help="Run data generation preprocessing pipeline")

    # Protinseq
    protinseq_parser = subparsers.add_parser("protinseq", help="Run Protinseq table and plot generation")
    protinseq_parser.add_argument("--sample", required=False, help="Excel sheet name from results.xlsx of the sample from which you want to use essentiality predictions (e.g., A_Cm30)")
    protinseq_parser.add_argument("--model", required=False, help="Which model you want to use predictions from (e.g., Poisson_absolute)")
    protinseq_parser.add_argument("--output", default=os.path.join("outputs", "protinseq"), help="Output directory for Protinseq tables and plots")
    protinseq_parser.add_argument("--predict", action="store_true", help="Run new Poisson ORF prediction instead of old tables/plots")

    args = parser.parse_args()

    # Dispatch to correct handler
    if args.command == "anubis":
        run_anubis(args, config)
    elif args.command == "preprocess":
        run_preprocessing()
    elif args.command == "protinseq":
        run_protinseq(args, config)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
