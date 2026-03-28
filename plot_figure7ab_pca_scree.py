#!/usr/bin/env python3
"""Generate manuscript Figure 7a/7b from center-level amino-acid vectors."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/CSAC_v2_primary_no_fallback_keep_nmr_firstmodel_2026-03-10_FE_DATABASE.csv"
OUTDIR = ROOT / "plots"
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


def parse_residues(raw: object) -> list[str]:
    if pd.isna(raw):
        return []
    toks = [x.strip().upper() for x in str(raw).split(",")]
    return [t for t in toks if t in AA_ORDER]


def build_composition_matrix(df: pd.DataFrame) -> np.ndarray:
    mat = np.zeros((len(df), len(AA_ORDER)), dtype=float)
    aa_idx = {aa: i for i, aa in enumerate(AA_ORDER)}
    for r, residues in enumerate(df["residues"]):
        if not residues:
            continue
        for aa in residues:
            mat[r, aa_idx[aa]] += 1.0
        mat[r, :] /= mat[r, :].sum()
    return mat


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    df = df[df["OxygenTOL"].isin(["aerobic", "anaerobic"])].copy()
    df["residues"] = df["Final Nearby Residues"].apply(parse_residues)
    df = df[df["residues"].map(len) > 0].copy().reset_index(drop=True)

    X = build_composition_matrix(df)

    pca2 = PCA(n_components=2)
    pcs = pca2.fit_transform(X)

    pca_all = PCA()
    pca_all.fit(X)
    evr = pca_all.explained_variance_ratio_

    # Save numeric outputs used in manuscript text/tables.
    scores = pd.DataFrame(
        {
            "PDB ID": df["PDB ID"].astype(str).str.upper().str.strip(),
            "OxygenTOL": df["OxygenTOL"].values,
            "PC1": pcs[:, 0],
            "PC2": pcs[:, 1],
        }
    )
    loadings = pd.DataFrame({"AA": AA_ORDER, "PC1_loading": pca2.components_[0], "PC2_loading": pca2.components_[1]})
    variance = pd.DataFrame(
        {
            "PC": np.arange(1, len(evr) + 1),
            "explained_variance_ratio": evr,
            "explained_variance_percent": evr * 100.0,
        }
    )

    scores.to_csv(OUTDIR / "figure7ab_pca_center_level_scores.csv", index=False)
    loadings.to_csv(OUTDIR / "figure7ab_pca_center_level_loadings.csv", index=False)
    variance.to_csv(OUTDIR / "figure7ab_pca_center_level_explained_variance.csv", index=False)

    sns.set_theme(style="whitegrid", context="talk")
    palette = {"anaerobic": "#1f77b4", "aerobic": "#d55e00"}

    # Figure 7a: PCA scatter
    fig, ax = plt.subplots(figsize=(8.6, 6.6))
    sns.scatterplot(
        data=scores,
        x="PC1",
        y="PC2",
        hue="OxygenTOL",
        hue_order=["anaerobic", "aerobic"],
        palette=palette,
        s=36,
        alpha=0.75,
        edgecolor="none",
        ax=ax,
    )
    pc1_pct = pca2.explained_variance_ratio_[0] * 100.0
    pc2_pct = pca2.explained_variance_ratio_[1] * 100.0
    ax.set_xlabel(f"PC1 ({pc1_pct:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pc2_pct:.1f}% variance)")
    ax.set_title("Figure 7a. PCA of Amino-Acid Composition Vectors")
    ax.legend(title="OxygenTOL", frameon=True)
    fig.tight_layout()
    fig.savefig(OUTDIR / "figure7a_pca_amino_acid_vectors_center_level.png", dpi=350, bbox_inches="tight")
    fig.savefig(OUTDIR / "figure7a_pca_amino_acid_vectors_center_level.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / "figure7a_pca_amino_acid_vectors_center_level.svg", bbox_inches="tight")
    plt.close(fig)

    # Figure 7b: Scree plot
    max_pc = min(10, len(evr))
    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    pcs_show = np.arange(1, max_pc + 1)
    evr_show = evr[:max_pc] * 100.0
    ax.plot(pcs_show, evr_show, marker="o", linewidth=1.8, markersize=5, color="#2c7fb8")
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Explained variance (%)")
    ax.set_title("Figure 7b. Scree Plot of Amino-Acid PCA")
    ax.set_xticks(pcs_show)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTDIR / "figure7b_scree_amino_acid_vectors_center_level.png", dpi=350, bbox_inches="tight")
    fig.savefig(OUTDIR / "figure7b_scree_amino_acid_vectors_center_level.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / "figure7b_scree_amino_acid_vectors_center_level.svg", bbox_inches="tight")
    plt.close(fig)

    print("Wrote Figure 7a/7b outputs to:", OUTDIR)


if __name__ == "__main__":
    main()
