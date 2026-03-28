#!/usr/bin/env python3
"""Generate manuscript Figure 5: ligand-radius distribution by OxygenTOL."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/CSAC_v2_primary_no_fallback_keep_nmr_firstmodel_2026-03-10_FE_DATABASE.csv"
OUTDIR = ROOT / "plots"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    df = df[df["OxygenTOL"].isin(["aerobic", "anaerobic"])].copy()
    df["ligand radius"] = pd.to_numeric(df["ligand radius"], errors="coerce")
    df = df.dropna(subset=["ligand radius"]).copy()

    order = ["anaerobic", "aerobic"]
    palette = {"anaerobic": "#1f77b4", "aerobic": "#d55e00"}

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8.8, 6.4))

    sns.violinplot(
        data=df,
        x="OxygenTOL",
        y="ligand radius",
        order=order,
        palette=palette,
        inner="box",
        cut=0,
        linewidth=1.0,
        ax=ax,
    )

    # Light point overlay for transparency about sample density.
    sns.stripplot(
        data=df,
        x="OxygenTOL",
        y="ligand radius",
        order=order,
        color="black",
        alpha=0.16,
        size=2,
        jitter=0.22,
        ax=ax,
    )

    ax.set_xlabel("OxygenTOL")
    ax.set_ylabel("Ligand radius (Å)")
    ax.set_title("Figure 5. Distribution of Ligand Sizes by OxygenTOL", pad=16)

    aero = df.loc[df["OxygenTOL"] == "aerobic", "ligand radius"].to_numpy()
    ana = df.loc[df["OxygenTOL"] == "anaerobic", "ligand radius"].to_numpy()
    welch_p = float(stats.ttest_ind(aero, ana, equal_var=False, nan_policy="omit").pvalue)
    mwu_p = float(stats.mannwhitneyu(aero, ana, alternative="two-sided").pvalue)

    fig.tight_layout()

    out_png = OUTDIR / "figure5_ligand_size_violin.png"
    out_pdf = OUTDIR / "figure5_ligand_size_violin.pdf"
    out_svg = OUTDIR / "figure5_ligand_size_violin.svg"
    out_stats = OUTDIR / "figure5_ligand_size_stats.csv"

    fig.savefig(out_png, dpi=350, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)

    summary = (
        df.groupby("OxygenTOL", as_index=False)
        .agg(
            count=("ligand radius", "count"),
            mean=("ligand radius", "mean"),
            median=("ligand radius", "median"),
            std=("ligand radius", "std"),
            min=("ligand radius", "min"),
            max=("ligand radius", "max"),
        )
    )
    summary["welch_p"] = welch_p
    summary["mannwhitney_p"] = mwu_p
    summary.to_csv(out_stats, index=False)

    print("Wrote:")
    print(out_png)
    print(out_pdf)
    print(out_svg)
    print(out_stats)


if __name__ == "__main__":
    main()
