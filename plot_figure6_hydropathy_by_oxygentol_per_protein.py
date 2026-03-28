#!/usr/bin/env python3
"""Generate manuscript Figure 6b: per-protein Fe-site hydropathy by OxygenTOL."""

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
    df["PDB ID"] = df["PDB ID"].astype(str).str.upper().str.strip()
    df["Average Hydropathy"] = pd.to_numeric(df["Average Hydropathy"], errors="coerce")
    df = df.dropna(subset=["Average Hydropathy"]).copy()

    pp = (
        df.groupby(["PDB ID", "OxygenTOL"], as_index=False)["Average Hydropathy"]
        .mean()
        .rename(columns={"Average Hydropathy": "Protein Mean Hydropathy"})
    )

    order = ["anaerobic", "aerobic"]
    palette = {"anaerobic": "#1f77b4", "aerobic": "#d55e00"}

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8.8, 6.4))

    sns.violinplot(
        data=pp,
        x="OxygenTOL",
        y="Protein Mean Hydropathy",
        order=order,
        palette=palette,
        inner="box",
        cut=0,
        linewidth=1.0,
        ax=ax,
    )

    sns.stripplot(
        data=pp,
        x="OxygenTOL",
        y="Protein Mean Hydropathy",
        order=order,
        color="black",
        alpha=0.35,
        size=3,
        jitter=0.15,
        ax=ax,
    )

    ax.set_xlabel("OxygenTOL")
    ax.set_ylabel("Per-protein mean hydropathy")
    ax.set_title("Figure 6b. Per-Protein Fe-Site Hydropathy by OxygenTOL", pad=14)

    ax.annotate(
        "",
        xy=(1.02, 0.92),
        xytext=(1.02, 0.08),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="<->", lw=1.2, color="black"),
    )
    ax.text(1.045, 0.93, "More hydrophobic", transform=ax.transAxes, ha="left", va="bottom", fontsize=9)
    ax.text(1.045, 0.07, "More hydrophilic", transform=ax.transAxes, ha="left", va="top", fontsize=9)

    fig.tight_layout()

    out_png = OUTDIR / "figure6b_hydropathy_by_oxygentol_per_protein.png"
    out_pdf = OUTDIR / "figure6b_hydropathy_by_oxygentol_per_protein.pdf"
    out_svg = OUTDIR / "figure6b_hydropathy_by_oxygentol_per_protein.svg"
    out_stats = OUTDIR / "figure6b_hydropathy_by_oxygentol_per_protein_stats.csv"

    fig.savefig(out_png, dpi=350, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)

    aero = pp.loc[pp["OxygenTOL"] == "aerobic", "Protein Mean Hydropathy"].to_numpy()
    ana = pp.loc[pp["OxygenTOL"] == "anaerobic", "Protein Mean Hydropathy"].to_numpy()
    welch_p = float(stats.ttest_ind(aero, ana, equal_var=False, nan_policy="omit").pvalue)
    mwu_p = float(stats.mannwhitneyu(aero, ana, alternative="two-sided").pvalue)

    summary = (
        pp.groupby("OxygenTOL", as_index=False)
        .agg(
            n_proteins=("PDB ID", "nunique"),
            mean=("Protein Mean Hydropathy", "mean"),
            median=("Protein Mean Hydropathy", "median"),
            std=("Protein Mean Hydropathy", "std"),
            min=("Protein Mean Hydropathy", "min"),
            max=("Protein Mean Hydropathy", "max"),
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
