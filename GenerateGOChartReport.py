#!/usr/bin/env python3
"""Generate DAVID GO reports for all CSV/XLSX gene-set files in a directory."""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from suds.client import Client

import ProcessGOresults as process

ROOT = Path(__file__).resolve().parent
MAPPINGS_PATH = ROOT / "Data" / "fb_synonym_fb_2024_06.tsv"
DAVID_WSDL = "https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl"
DAVID_ENDPOINT = (
    "https://davidbioinformatics.nih.gov/webservice/services/"
    "DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/"
)
DEFAULT_DAVID_EMAIL = "aadishms@umich.edu"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("indir", help="Directory with input gene-set CSV/XLSX files")
    parser.add_argument(
        "convert_to_name",
        type=int,
        choices=[0, 1],
        help="1: convert FBgn IDs to symbols in output, 0: keep FBgn IDs",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Optional output directory (default: <indir>/GO_Analysis)",
    )
    parser.add_argument(
        "--david-email",
        default=os.getenv("DAVID_EMAIL", DEFAULT_DAVID_EMAIL),
        help="Registered DAVID account email (default: DAVID_EMAIL env var or script default)",
    )
    return parser.parse_args()


def init_david_client(email):
    print(f"url={DAVID_WSDL}")
    client = Client(DAVID_WSDL)
    client.wsdl.services[0].setlocation(DAVID_ENDPOINT)
    client.service.authenticate(email)
    client.service.setCurrentSpecies("7227")  # Drosophila melanogaster
    return client


def load_input_table(path):
    if path.suffix.lower() == ".xlsx":
        df = pd.read_excel(path, dtype=str)
    else:
        df = pd.read_csv(path, dtype=str)

    if "flybase_gene_id" in df.columns:
        id_col = "flybase_gene_id"
    elif "gene_id" in df.columns:
        df = df.rename(columns={"gene_id": "flybase_gene_id"})
        id_col = "flybase_gene_id"
    else:
        return np.array([])

    values = (
        df[id_col]
        .dropna()
        .astype(str)
        .str.strip()
    )
    values = values[(values != "") & (values != "-")]
    return np.unique(values)


def fetch_go_report(client, flybase_ids, list_name):
    id_type = "FLYBASE_GENE_ID"
    list_type = 0
    client.service.addList(",".join(flybase_ids), id_type, list_name, list_type)
    client.service.setCategories("GOTERM_BP_DIRECT,GOTERM_CC_DIRECT,GOTERM_MF_DIRECT,KEGG_PATHWAY")
    chart_report = client.service.getChartReport(0.1, 1)

    go_results = []
    for record in chart_report:
        go_results.append(
            {
                "Category": record["categoryName"],
                "Term": record["termName"],
                "Count": int(record["listHits"]),
                "%": str(record["percent"]),
                "Pvalue": float(record["ease"]),
                "Genes": record["geneIds"],
                "List Total": int(record["listTotals"]),
                "Pop Hits": int(record["popHits"]),
                "Pop Total": int(record["popTotals"]),
                "Fold Enrichment": float(record["foldEnrichment"]),
                "Bonferroni": float(record["bonferroni"]),
                "Benjamini": float(record["benjamini"]),
                "FDR": float(record["afdr"]),
            }
        )

    go_df = pd.DataFrame(go_results)
    if go_df.empty:
        return go_df
    return go_df.sort_values("FDR").query("FDR <= 10")


def main():
    args = parse_args()
    input_dir = Path(args.indir).resolve()
    output_dir = args.outdir.resolve() if args.outdir else input_dir / "GO_Analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not MAPPINGS_PATH.exists():
        raise FileNotFoundError(f"Required mapping file not found: {MAPPINGS_PATH}")
    mappings_df = pd.read_csv(MAPPINGS_PATH, sep="\t", header=0, low_memory=False)

    data_files = sorted(input_dir.glob("*.csv")) + sorted(input_dir.glob("*.xlsx"))
    client = init_david_client(args.david_email)

    for data_file in data_files:
        if "FullDF" in data_file.name:
            continue

        print(data_file)
        flybase_ids = load_input_table(data_file)
        if len(flybase_ids) == 0:
            print(f"Skipping {data_file.name}: no flybase_gene_id/gene_id values found.")
            continue

        print(f"Before deduping: {len(flybase_ids)}")
        print(f"After deduping: {len(flybase_ids)}")
        go_df = fetch_go_report(client, flybase_ids, data_file.stem)
        if go_df.empty:
            print(f"No significant GO terms (FDR <= 10) for {data_file.name}. Skipping.")
            continue

        out_path = output_dir / f"{data_file.stem}_GO_Analysis.xlsx"
        go_df.to_excel(out_path, index=False)
        print(f"GO Analysis saved to: {out_path}")

    print("All files processed.")
    process.process_csv_files(str(output_dir), mappings_df, args.convert_to_name)


if __name__ == "__main__":
    main()
