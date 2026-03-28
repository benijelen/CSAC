#!/usr/bin/env python3
"""Generate manuscript Figure 6a: Fe-center hydropathy by OxygenTOL."""

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
    df["Average Hydropathy"] = pd.to_numeric(df["Average Hydropathy"], errors="coerce")
    df = df.dropna(subset=["Average Hydropathy"]).copy()

    order = ["anaerobic", "aerobic"]
    palette = {"anaerobic": "#1f77b4", "aerobic": "#d55e00"}

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8.8, 6.4))

    sns.violinplot(
        data=df,
        x="OxygenTOL",
        y="Average Hydropathy",
        order=order,
        palette=palette,
        inner="box",
        cut=0,
        linewidth=1.0,
        ax=ax,
    )

    sns.stripplot(
        data=df,
        x="OxygenTOL",
        y="Average Hydropathy",
        order=order,
        color="black",
        alpha=0.16,
        size=2,
        jitter=0.22,
        ax=ax,
    )

    ax.set_xlabel("OxygenTOL")
    ax.set_ylabel("Average Hydropathy")
    ax.set_title("Figure 6a. Fe-Center Hydropathy by OxygenTOL", pad=14)

    # Reader orientation cue: higher hydropathy is more hydrophobic.
    ax.annotate(
        "",
        xy=(1.02, 0.92),
        xytext=(1.02, 0.08),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="<->", lw=1.2, color="black"),
    )
    ax.text(
        1.045,
        0.93,
        "More hydrophobic",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
    )
    ax.text(
        1.045,
        0.07,
        "More hydrophilic",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
    )

    fig.tight_layout()

    out_png = OUTDIR / "figure6a_hydropathy_by_oxygentol_center_level.png"
    out_pdf = OUTDIR / "figure6a_hydropathy_by_oxygentol_center_level.pdf"
    out_svg = OUTDIR / "figure6a_hydropathy_by_oxygentol_center_level.svg"
    out_stats = OUTDIR / "figure6a_hydropathy_by_oxygentol_center_level_stats.csv"

    fig.savefig(out_png, dpi=350, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)

    aero = df.loc[df["OxygenTOL"] == "aerobic", "Average Hydropathy"].to_numpy()
    ana = df.loc[df["OxygenTOL"] == "anaerobic", "Average Hydropathy"].to_numpy()
    welch_p = float(stats.ttest_ind(aero, ana, equal_var=False, nan_policy="omit").pvalue)
    mwu_p = float(stats.mannwhitneyu(aero, ana, alternative="two-sided").pvalue)

    summary = (
        df.groupby("OxygenTOL", as_index=False)
        .agg(
            count=("Average Hydropathy", "count"),
            mean=("Average Hydropathy", "mean"),
            median=("Average Hydropathy", "median"),
            std=("Average Hydropathy", "std"),
            min=("Average Hydropathy", "min"),
            max=("Average Hydropathy", "max"),
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
