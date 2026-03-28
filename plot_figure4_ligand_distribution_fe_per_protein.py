#!/usr/bin/env python3
"""Generate Figure 4 (per-protein): grouped Fe ligand prevalence."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/CSAC_v2_primary_no_fallback_keep_nmr_firstmodel_2026-03-10_FE_DATABASE.csv"
OUTDIR = ROOT / "plots"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    df = df[df["OxygenTOL"].isin(["aerobic", "anaerobic"])].copy()
    df["PDB ID"] = df["PDB ID"].astype(str).str.upper().str.strip()
    df["ligand type"] = df["ligand type"].astype(str).str.strip()

    # Grouping rule requested for Figure 4:
    # - keep all ligand types with count >= FE2 count
    # - always keep F3S, F4S, SRM, SF3
    raw_counts = df["ligand type"].value_counts()
    fe2_count = int(raw_counts.get("FE2", 0))
    keep_set = set(raw_counts[raw_counts >= fe2_count].index.tolist()) if fe2_count > 0 else set()
    keep_set.update({"F3S", "F4S", "SRM", "SF3"})

    df["ligand_group"] = df["ligand type"].where(df["ligand type"].isin(keep_set), "Other")

    # Per-protein prevalence: each protein counts once per ligand group.
    pp = df[["PDB ID", "OxygenTOL", "ligand_group"]].drop_duplicates()
    total_proteins = pp["PDB ID"].nunique()
    counts = pp["ligand_group"].value_counts().rename("n_proteins").reset_index()
    counts = counts.rename(columns={"index": "ligand_group", "ligand_group": "ligand_group"})
    counts["percent_proteins"] = counts["n_proteins"] / float(total_proteins) * 100.0

    ligand_order = counts["ligand_group"].tolist()
    palette = sns.color_palette("Spectral", n_colors=len(ligand_order))

    sns.set_theme(style="white", context="talk")
    fig, ax = plt.subplots(figsize=(14.0, 6.6))

    sns.barplot(
        data=counts,
        x="ligand_group",
        y="percent_proteins",
        order=ligand_order,
        palette=palette,
        edgecolor="black",
        linewidth=0.5,
        ax=ax,
    )

    ax.set_xlabel("Ligand group")
    ax.set_ylabel("Proteins containing ligand group (%)")
    ax.set_title("Figure 4. Fe Ligand Groups by Protein Prevalence")
    ax.tick_params(axis="x", rotation=30)
    ax.set_ylim(0, max(counts["percent_proteins"]) * 1.15)

    for i, row in counts.iterrows():
        ax.text(
            i,
            row["percent_proteins"] + 0.35,
            f"{row['percent_proteins']:.1f}%",
            ha="center",
            va="bottom",
            fontsize=9,
            rotation=90,
        )

    fig.tight_layout()

    out_png = OUTDIR / "figure4_ligand_distribution_fe_per_protein.png"
    out_pdf = OUTDIR / "figure4_ligand_distribution_fe_per_protein.pdf"
    out_svg = OUTDIR / "figure4_ligand_distribution_fe_per_protein.svg"
    out_csv = OUTDIR / "figure4_ligand_distribution_fe_per_protein_counts.csv"

    fig.savefig(out_png, dpi=350, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)

    by_oxy = (
        pp.groupby(["ligand_group", "OxygenTOL"], as_index=False)["PDB ID"]
        .nunique()
        .pivot(index="ligand_group", columns="OxygenTOL", values="PDB ID")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    out = counts.merge(by_oxy, on="ligand_group", how="left").fillna(0)
    for col in ["aerobic", "anaerobic"]:
        if col in out.columns:
            out[col] = out[col].astype(int)
    out["total_unique_proteins"] = int(total_proteins)
    out.to_csv(out_csv, index=False)

    print("Wrote:")
    print(out_png)
    print(out_pdf)
    print(out_svg)
    print(out_csv)


if __name__ == "__main__":
    main()
