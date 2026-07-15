from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
SOURCE_XLSX = REPO_ROOT / "GSE172268_merge_Annotation_Result.xlsx"
OUT_DIR = REPO_ROOT / "netact_input"
REF_DIR = OUT_DIR / "reference"

FPKM_COLS = ["WT-1_FPKM", "WT-2_FPKM", "Toluene-1_FPKM", "Toluene-2_FPKM"]
COUNT_COLS = [
    "WT-1_Read_Count",
    "WT-2_Read_Count",
    "Toluene-1_Read_Count",
    "Toluene-2_Read_Count",
]

TF_LIST_URL = "https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_tfs.txt"
GENE_PATTERN = re.compile(r"\(([A-Za-z0-9-]+)\)")


def extract_gene_symbol(nt_value: object) -> str | None:
    if not isinstance(nt_value, str) or nt_value == ".":
        return None
    match = GENE_PATTERN.search(nt_value)
    if not match:
        return None
    symbol = match.group(1).strip()
    if not symbol or symbol == ".":
        return None
    return symbol


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    REF_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_excel(SOURCE_XLSX, sheet_name="expression_profile")
    df["GeneSymbol"] = df["NT"].map(extract_gene_symbol)
    df = df[df["GeneSymbol"].notna()].copy()

    expr_fpkm = (
        df[["GeneSymbol", *FPKM_COLS]]
        .groupby("GeneSymbol", as_index=True)
        .mean()
        .sort_index()
    )
    expr_counts = (
        df[["GeneSymbol", *COUNT_COLS]]
        .groupby("GeneSymbol", as_index=True)
        .sum()
        .sort_index()
    )

    expr_fpkm.to_csv(OUT_DIR / "expression_fpkm.tsv", sep="\t", index=True, index_label="Gene")
    expr_counts.to_csv(
        OUT_DIR / "expression_counts.tsv", sep="\t", index=True, index_label="Gene"
    )

    sample_metadata = pd.DataFrame(
        {
            "sample": FPKM_COLS,
            "condition": ["WT", "WT", "Toluene", "Toluene"],
            "datatype": ["FPKM", "FPKM", "FPKM", "FPKM"],
        }
    )
    sample_metadata.to_csv(OUT_DIR / "sample_metadata.tsv", sep="\t", index=False)

    ranking = expr_fpkm.copy()
    ranking["WT_mean"] = ranking[["WT-1_FPKM", "WT-2_FPKM"]].mean(axis=1)
    ranking["Toluene_mean"] = ranking[["Toluene-1_FPKM", "Toluene-2_FPKM"]].mean(axis=1)
    ranking["log2FC_Toluene_vs_WT"] = np.log2(
        (ranking["Toluene_mean"] + 1e-6) / (ranking["WT_mean"] + 1e-6)
    )
    ranking = ranking[["WT_mean", "Toluene_mean", "log2FC_Toluene_vs_WT"]].sort_values(
        "log2FC_Toluene_vs_WT", ascending=False
    )
    ranking.to_csv(OUT_DIR / "gene_ranking_toluene_vs_wt.tsv", sep="\t", index=True, index_label="Gene")

    tf_reference = pd.read_csv(TF_LIST_URL, header=None, names=["TF"])
    tf_reference["TF"] = tf_reference["TF"].astype(str).str.strip()
    tf_reference = tf_reference[tf_reference["TF"] != ""].drop_duplicates().sort_values("TF")
    tf_reference.to_csv(REF_DIR / "hs_hgnc_tfs.txt", index=False, header=False)

    dataset_genes = set(expr_fpkm.index)
    tf_candidates = sorted(dataset_genes.intersection(set(tf_reference["TF"])))

    pd.DataFrame({"TF": tf_candidates}).to_csv(OUT_DIR / "tf_candidates_in_dataset.tsv", sep="\t", index=False)
    expr_fpkm.loc[expr_fpkm.index.intersection(tf_candidates)].sort_index().to_csv(
        OUT_DIR / "tf_expression_fpkm.tsv", sep="\t", index=True, index_label="TF"
    )

    summary_lines = [
        "Generated NetAct input files",
        "",
        "Files:",
        "- expression_fpkm.tsv: gene-by-sample FPKM matrix (required expression input).",
        "- expression_counts.tsv: gene-by-sample read count matrix.",
        "- sample_metadata.tsv: sample annotations (WT/Toluene groups).",
        "- gene_ranking_toluene_vs_wt.tsv: log2 fold-change ranking for enrichment workflows.",
        "- tf_candidates_in_dataset.tsv: TFs present in expression matrix, based on pySCENIC human TF list.",
        "- tf_expression_fpkm.tsv: expression matrix restricted to TF candidates.",
        "- reference/hs_hgnc_tfs.txt: downloaded human TF reference list.",
        "",
        "Gene symbols were extracted from the NT annotation column by parsing symbols inside parentheses, then aggregated by mean FPKM / summed read counts across contigs.",
    ]
    (OUT_DIR / "README.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    print(f"Done. Output directory: {OUT_DIR}")
    print(f"Genes in expression matrix: {expr_fpkm.shape[0]}")
    print(f"TF candidates in dataset: {len(tf_candidates)}")


if __name__ == "__main__":
    main()
