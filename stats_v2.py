#!/usr/bin/env python3
"""v2 statistical analysis with pseudoreplication controls.

Runs two layers:
1) Per-center mixed-effects models (random effects for PDB + metabolism VC)
2) Per-protein aggregated robustness tests

Also supports NMR sensitivity by excluding NMR structures identified from
the combined metadata table.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf


def parse_residue_list(raw: object) -> List[str]:
    if pd.isna(raw):
        return []
    tokens = [x.strip().upper() for x in str(raw).split(",")]
    return [t for t in tokens if len(t) == 1 and t.isalpha()]


def compute_cys_fraction(residues: List[str]) -> float:
    if not residues:
        return np.nan
    return residues.count("C") / float(len(residues))


def cohen_d(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    pooled = ((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2)
    if pooled <= 0:
        return np.nan
    return (np.mean(x) - np.mean(y)) / np.sqrt(pooled)


def load_nmr_ids(metadata_csv: Path) -> Set[str]:
    nmr_ids: Set[str] = set()
    if not metadata_csv.exists():
        return nmr_ids
    with metadata_csv.open() as f:
        r = csv.DictReader(f)
        for row in r:
            pid = str(row.get("pdb_id", "")).strip().upper()
            method = str(row.get("method", "")).upper()
            if not pid or pid in {"NAN", "NA", "N/A", "UNKNOWN", "TBD_COMPUTED_MODEL"}:
                continue
            if "NMR" in method:
                nmr_ids.add(pid)
    return nmr_ids


def load_method_map(metadata_csv: Path) -> Dict[str, str]:
    method_map: Dict[str, str] = {}
    if not metadata_csv.exists():
        return method_map
    with metadata_csv.open() as f:
        r = csv.DictReader(f)
        for row in r:
            pid = str(row.get("pdb_id", "")).strip().upper()
            method = str(row.get("method", "")).strip()
            if not pid or pid in {"NAN", "NA", "N/A", "UNKNOWN", "TBD_COMPUTED_MODEL"}:
                continue
            if pid not in method_map and method:
                method_map[pid] = method
    return method_map


@dataclass
class AnalysisResult:
    dataset: str
    layer: str
    metric: str
    model: str
    n_rows: int
    n_proteins: int
    n_aerobic: int
    n_anaerobic: int
    estimate: float
    p_value: float
    effect_size: float
    notes: str


def run_mixed_model(
    df: pd.DataFrame,
    response: str,
) -> Dict[str, float]:
    """Fit mixed model response ~ oxy_bin with PDB random intercept and
    metabolism variance components. Falls back to PDB-only random intercept."""
    formula = f"{response} ~ oxy_bin"
    fit_note = "mixedlm_pdb_plus_metabolism_vc"
    try:
        model = smf.mixedlm(
            formula=formula,
            data=df,
            groups=df["PDB ID"],
            vc_formula={"Metabolism": "0 + C(Metabolism)"},
            re_formula="1",
        )
        fit = model.fit(reml=False, method="lbfgs", maxiter=500, disp=False)
    except Exception:
        fit_note = "mixedlm_pdb_only_random_intercept"
        model = smf.mixedlm(
            formula=formula,
            data=df,
            groups=df["PDB ID"],
            re_formula="1",
        )
        fit = model.fit(reml=False, method="lbfgs", maxiter=500, disp=False)

    estimate = float(fit.params.get("oxy_bin", np.nan))
    p_value = float(fit.pvalues.get("oxy_bin", np.nan))
    return {"estimate": estimate, "p_value": p_value, "model_note": fit_note}


def run_per_protein_tests(df: pd.DataFrame, metric: str) -> Dict[str, float]:
    agg = (
        df.groupby(["PDB ID", "OxygenTOL", "Metabolism"], as_index=False)[metric]
        .mean()
        .rename(columns={metric: f"{metric}_protein_mean"})
    )

    a = agg.loc[agg["OxygenTOL"] == "aerobic", f"{metric}_protein_mean"].dropna().values
    b = agg.loc[agg["OxygenTOL"] == "anaerobic", f"{metric}_protein_mean"].dropna().values
    if len(a) < 2 or len(b) < 2:
        return {"p_ttest": np.nan, "p_mwu": np.nan, "cohen_d": np.nan, "n_proteins": len(agg)}
    t = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
    mwu = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "p_ttest": float(t.pvalue),
        "p_mwu": float(mwu.pvalue),
        "cohen_d": float(cohen_d(a, b)),
        "n_proteins": int(len(agg)),
    }


def analyze_dataset(df: pd.DataFrame, name: str) -> List[AnalysisResult]:
    out: List[AnalysisResult] = []
    work = df.copy()
    work["oxy_bin"] = (work["OxygenTOL"] == "aerobic").astype(int)
    # Patsy formulas cannot handle spaces in column names.
    work["average_hydropathy"] = work["Average Hydropathy"]

    # per-center mixed effects for hydropathy
    metric_map = [("Average Hydropathy", "average_hydropathy"), ("cys_fraction", "cys_fraction")]
    for metric_label, metric_col in metric_map:
        mdf = work.dropna(subset=[metric_col]).copy()
        if len(mdf) < 10:
            continue
        mm = run_mixed_model(mdf, metric_col)
        eff = cohen_d(
            mdf.loc[mdf["OxygenTOL"] == "aerobic", metric_col].values,
            mdf.loc[mdf["OxygenTOL"] == "anaerobic", metric_col].values,
        )
        out.append(
            AnalysisResult(
                dataset=name,
                layer="per_center_mixed_effects",
                metric=metric_label,
                model=mm["model_note"],
                n_rows=int(len(mdf)),
                n_proteins=int(mdf["PDB ID"].nunique()),
                n_aerobic=int((mdf["OxygenTOL"] == "aerobic").sum()),
                n_anaerobic=int((mdf["OxygenTOL"] == "anaerobic").sum()),
                estimate=mm["estimate"],
                p_value=mm["p_value"],
                effect_size=float(eff),
                notes="fixed effect is aerobic-vs-anaerobic contrast (oxy_bin)",
            )
        )

        # per-protein robustness
        pt = run_per_protein_tests(mdf, metric_col)
        out.append(
            AnalysisResult(
                dataset=name,
                layer="per_protein_aggregated",
                metric=metric_label,
                model="welch_ttest_and_mannwhitney",
                n_rows=int(len(mdf)),
                n_proteins=int(pt["n_proteins"]),
                n_aerobic=int(mdf.loc[mdf["OxygenTOL"] == "aerobic", "PDB ID"].nunique()),
                n_anaerobic=int(mdf.loc[mdf["OxygenTOL"] == "anaerobic", "PDB ID"].nunique()),
                estimate=np.nan,
                p_value=float(pt["p_ttest"]),
                effect_size=float(pt["cohen_d"]),
                notes=f"mannwhitney_p={pt['p_mwu']}",
            )
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--database",
        default="results/CSAC_v2_primary_qc_2026-03-10_FE_DATABASE.csv",
        help="CSAC database output CSV",
    )
    ap.add_argument(
        "--metadata",
        default="combined_aerobic_anaerobic_master_review.csv",
        help="metadata CSV used to identify NMR PDB IDs",
    )
    ap.add_argument(
        "--out-prefix",
        default="results/stats_v2_primary_qc_2026-03-10",
        help="output prefix for generated files",
    )
    ap.add_argument(
        "--xray-only",
        action="store_true",
        help="keep only rows with X-RAY DIFFRACTION method",
    )
    ap.add_argument(
        "--exclude-pdb-list",
        default="",
        help="optional text file with one PDB ID per line to exclude before analysis",
    )
    args = ap.parse_args()

    db_path = Path(args.database)
    meta_path = Path(args.metadata)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(db_path)
    # keep only strict aerobic/anaerobic for hypothesis testing
    df = df[df["OxygenTOL"].isin(["aerobic", "anaerobic"])].copy()
    method_map = load_method_map(meta_path)
    df["method_meta"] = df["PDB ID"].astype(str).str.upper().map(method_map)
    if args.xray_only:
        df = df[df["method_meta"].astype(str).str.upper().str.contains("X-RAY DIFFRACTION", na=False)].copy()
    if args.exclude_pdb_list:
        ex_path = Path(args.exclude_pdb_list)
        if ex_path.exists():
            exclude_ids = {line.strip().upper() for line in ex_path.read_text().splitlines() if line.strip()}
            df = df[~df["PDB ID"].astype(str).str.upper().isin(exclude_ids)].copy()
    # drop rows without Fe neighborhood score
    df["residues_parsed"] = df["Final Nearby Residues"].apply(parse_residue_list)
    df["cys_fraction"] = df["residues_parsed"].apply(compute_cys_fraction)
    df["Average Hydropathy"] = pd.to_numeric(df["Average Hydropathy"], errors="coerce")

    nmr_ids = load_nmr_ids(meta_path)
    df_non_nmr = df[~df["PDB ID"].str.upper().isin(nmr_ids)].copy()

    all_results = []
    all_results.extend(analyze_dataset(df, "primary_all_structures"))
    all_results.extend(analyze_dataset(df_non_nmr, "sensitivity_exclude_nmr"))

    out_df = pd.DataFrame([r.__dict__ for r in all_results])
    out_csv = Path(f"{out_prefix}_summary.csv")
    out_df.to_csv(out_csv, index=False)

    # Also save per-protein aggregates for traceability
    protein_all = (
        df.groupby(["PDB ID", "OxygenTOL", "Metabolism"], as_index=False)
        .agg(
            hydropathy_mean=("Average Hydropathy", "mean"),
            cys_fraction_mean=("cys_fraction", "mean"),
            n_centers=("PDB ID", "size"),
        )
    )
    protein_nmr = (
        df_non_nmr.groupby(["PDB ID", "OxygenTOL", "Metabolism"], as_index=False)
        .agg(
            hydropathy_mean=("Average Hydropathy", "mean"),
            cys_fraction_mean=("cys_fraction", "mean"),
            n_centers=("PDB ID", "size"),
        )
    )
    protein_all.to_csv(f"{out_prefix}_per_protein_all.csv", index=False)
    protein_nmr.to_csv(f"{out_prefix}_per_protein_exclude_nmr.csv", index=False)

    with open(f"{out_prefix}_nmr_ids_used.txt", "w") as f:
        for pid in sorted(nmr_ids):
            f.write(pid + "\n")

    print(f"Rows (all): {len(df)}")
    print(f"Rows (exclude NMR): {len(df_non_nmr)}")
    print(f"NMR IDs: {sorted(nmr_ids)}")
    print(f"Wrote: {out_csv}")


if __name__ == "__main__":
    main()
