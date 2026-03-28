#!/usr/bin/env python3
"""Generate manuscript Figure 4: grouped Fe ligand-type distribution bar chart."""

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
    df["ligand type"] = df["ligand type"].astype(str).str.strip()

    raw_counts = df["ligand type"].value_counts()

    # Keep ligand types at least as frequent as FE2, plus user-selected categories.
    fe2_count = int(raw_counts.get("FE2", 0))
    keep_set = set(raw_counts[raw_counts >= fe2_count].index.tolist()) if fe2_count > 0 else set()
    keep_set.update({"F3S", "F4S", "SRM", "SF3"})

    df["ligand_group"] = df["ligand type"].where(df["ligand type"].isin(keep_set), "Other")
    counts = df["ligand_group"].value_counts().rename("count").reset_index()
    counts = counts.rename(columns={"index": "ligand_group", "ligand_group": "ligand_group"})
    counts["percent"] = counts["count"] / counts["count"].sum() * 100.0

    # Keep bars ordered by abundance for readability.
    ligand_order = counts["ligand_group"].tolist()
    palette = sns.color_palette("Spectral", n_colors=len(ligand_order))

    sns.set_theme(style="white", context="talk")
    fig, ax = plt.subplots(figsize=(16.0, 6.8))

    sns.barplot(
        data=counts,
        x="ligand_group",
        y="percent",
        order=ligand_order,
        palette=palette,
        edgecolor="black",
        linewidth=0.5,
        ax=ax,
    )

    ax.set_xlabel("Ligand group")
    ax.set_ylabel("Percent of Fe centers (%)")
    ax.set_title("Figure 4. Distribution of Fe Ligand Groups (Aerobic + Anaerobic)")
    ax.tick_params(axis="x", rotation=35)
    ax.set_ylim(0, max(counts["percent"]) * 1.12)

    for i, row in counts.iterrows():
        ax.text(
            i,
            row["percent"] + 0.35,
            f"{row['percent']:.1f}%",
            ha="center",
            va="bottom",
            fontsize=9,
            rotation=90,
        )

    fig.tight_layout()

    out_png = OUTDIR / "figure4_ligand_distribution_fe.png"
    out_pdf = OUTDIR / "figure4_ligand_distribution_fe.pdf"
    out_svg = OUTDIR / "figure4_ligand_distribution_fe.svg"
    out_csv = OUTDIR / "figure4_ligand_distribution_fe_counts.csv"

    fig.savefig(out_png, dpi=350, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)

    # Add oxygen-specific counts for caption support.
    by_oxy = (
        df.groupby(["ligand_group", "OxygenTOL"], as_index=False)
        .size()
        .pivot(index="ligand_group", columns="OxygenTOL", values="size")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    out = counts.merge(by_oxy, on="ligand_group", how="left").fillna(0)
    for col in ["aerobic", "anaerobic"]:
        if col in out.columns:
            out[col] = out[col].astype(int)
    out.to_csv(out_csv, index=False)

    print("Wrote:")
    print(out_png)
    print(out_pdf)
    print(out_svg)
    print(out_csv)


if __name__ == "__main__":
    main()
