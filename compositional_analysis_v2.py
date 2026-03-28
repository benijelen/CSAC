#!/usr/bin/env python3
"""Compositional amino-acid analysis (CLR + model + permutation).

Primary target: reviewer-proof supplement for aerobic vs anaerobic differences.
Uses canonical strict FE database by default.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


def parse_residues(raw: object) -> List[str]:
    if pd.isna(raw):
        return []
    toks = [x.strip().upper() for x in str(raw).split(",")]
    return [t for t in toks if t in AA_ORDER]


def counts_from_residues(res: List[str]) -> np.ndarray:
    out = np.zeros(len(AA_ORDER), dtype=float)
    if not res:
        return out
    idx = {aa: i for i, aa in enumerate(AA_ORDER)}
    for r in res:
        out[idx[r]] += 1.0
    return out


def closure(x: np.ndarray) -> np.ndarray:
    s = x.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return x / s


def clr_transform(comp: np.ndarray, pseudocount: float = 1e-6) -> np.ndarray:
    x = comp + pseudocount
    x = closure(x)
    lx = np.log(x)
    gm = lx.mean(axis=1, keepdims=True)
    return lx - gm


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = min(prev, ranked[i] * n / rank)
        q[i] = val
        prev = val
    out = np.empty(n, dtype=float)
    out[order] = q
    return out


def perm_pvalue_delta(
    values: np.ndarray,
    groups: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
) -> float:
    # groups: 1=aerobic, 0=anaerobic
    obs = float(values[groups == 1].mean() - values[groups == 0].mean())
    hits = 0
    for _ in range(n_perm):
        perm = rng.permutation(groups)
        stat = float(values[perm == 1].mean() - values[perm == 0].mean())
        if abs(stat) >= abs(obs):
            hits += 1
    return (hits + 1) / (n_perm + 1)


def global_centroid_distance(clr: np.ndarray, groups: np.ndarray) -> float:
    a = clr[groups == 1]
    b = clr[groups == 0]
    return float(np.linalg.norm(a.mean(axis=0) - b.mean(axis=0)))


def global_perm_test(
    clr: np.ndarray,
    groups: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
) -> Tuple[float, float]:
    obs = global_centroid_distance(clr, groups)
    hits = 0
    for _ in range(n_perm):
        perm = rng.permutation(groups)
        stat = global_centroid_distance(clr, perm)
        if stat >= obs:
            hits += 1
    p = (hits + 1) / (n_perm + 1)
    return obs, p


def load_canonical_database(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[df["OxygenTOL"].isin(["aerobic", "anaerobic"])].copy()
    df["residues"] = df["Final Nearby Residues"].apply(parse_residues)
    df = df[df["residues"].map(len) > 0].copy()
    return df


def build_per_center(df: pd.DataFrame) -> pd.DataFrame:
    counts = np.vstack(df["residues"].apply(counts_from_residues).to_numpy())
    comp = closure(counts)
    clr = clr_transform(comp)
    out = df[["PDB ID", "OxygenTOL", "Metabolism"]].copy()
    out["oxy_bin"] = (out["OxygenTOL"] == "aerobic").astype(int)
    for i, aa in enumerate(AA_ORDER):
        out[f"prop_{aa}"] = comp[:, i]
        out[f"clr_{aa}"] = clr[:, i]
    return out


def build_per_protein(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pid, g in df.groupby("PDB ID"):
        counts = np.vstack(g["residues"].apply(counts_from_residues).to_numpy()).sum(axis=0)
        rows.append(
            {
                "PDB ID": pid,
                "OxygenTOL": g["OxygenTOL"].iloc[0],
                "Metabolism": g["Metabolism"].iloc[0] if "Metabolism" in g.columns else "Unknown",
                "counts": counts,
            }
        )
    out = pd.DataFrame(rows)
    comp = closure(np.vstack(out["counts"].to_numpy()))
    clr = clr_transform(comp)
    out["oxy_bin"] = (out["OxygenTOL"] == "aerobic").astype(int)
    for i, aa in enumerate(AA_ORDER):
        out[f"prop_{aa}"] = comp[:, i]
        out[f"clr_{aa}"] = clr[:, i]
    return out.drop(columns=["counts"])


def run_aa_tests(df: pd.DataFrame, layer: str, n_perm: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    groups = df["oxy_bin"].to_numpy()
    rows = []
    for aa in AA_ORDER:
        prop_col = f"prop_{aa}"
        clr_col = f"clr_{aa}"
        m_a = float(df.loc[df["oxy_bin"] == 1, prop_col].mean())
        m_b = float(df.loc[df["oxy_bin"] == 0, prop_col].mean())
        c_a = float(df.loc[df["oxy_bin"] == 1, clr_col].mean())
        c_b = float(df.loc[df["oxy_bin"] == 0, clr_col].mean())

        # model test (metabolism-adjusted OLS on CLR coordinate)
        mdf = df[[clr_col, "oxy_bin", "Metabolism"]].copy()
        fit = smf.ols(f"{clr_col} ~ oxy_bin + C(Metabolism)", data=mdf).fit()
        est = float(fit.params.get("oxy_bin", np.nan))
        p_model = float(fit.pvalues.get("oxy_bin", np.nan))

        p_perm = perm_pvalue_delta(df[clr_col].to_numpy(), groups, n_perm=n_perm, rng=rng)

        rows.append(
            {
                "layer": layer,
                "aa": aa,
                "n_rows": int(len(df)),
                "n_aerobic": int((groups == 1).sum()),
                "n_anaerobic": int((groups == 0).sum()),
                "mean_prop_aerobic": m_a,
                "mean_prop_anaerobic": m_b,
                "delta_prop_aero_minus_ana": m_a - m_b,
                "mean_clr_aerobic": c_a,
                "mean_clr_anaerobic": c_b,
                "delta_clr_aero_minus_ana": c_a - c_b,
                "model_estimate_oxy_bin": est,
                "p_model_metabolism_adjusted": p_model,
                "p_perm_delta_clr": float(p_perm),
            }
        )
    out = pd.DataFrame(rows)
    out["q_model_fdr"] = bh_fdr(out["p_model_metabolism_adjusted"].to_numpy())
    out["q_perm_fdr"] = bh_fdr(out["p_perm_delta_clr"].to_numpy())
    return out


def run_global_test(df: pd.DataFrame, layer: str, n_perm: int, seed: int) -> pd.DataFrame:
    clr_cols = [f"clr_{aa}" for aa in AA_ORDER]
    clr = df[clr_cols].to_numpy()
    groups = df["oxy_bin"].to_numpy()
    obs, p = global_perm_test(clr, groups, n_perm=n_perm, rng=np.random.default_rng(seed))
    return pd.DataFrame(
        [
            {
                "layer": layer,
                "n_rows": int(len(df)),
                "n_aerobic": int((groups == 1).sum()),
                "n_anaerobic": int((groups == 0).sum()),
                "global_stat_centroid_distance_clr": obs,
                "global_perm_p_value": p,
                "n_perm": n_perm,
            }
        ]
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--database",
        default="results/CSAC_v2_primary_no_fallback_keep_nmr_firstmodel_2026-03-10_FE_DATABASE.csv",
    )
    ap.add_argument("--out-prefix", default="results/compositional_v2_primary_final_2026-03-10")
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    db = load_canonical_database(Path(args.database))

    per_center = build_per_center(db)
    per_protein = build_per_protein(db)

    aa_center = run_aa_tests(per_center, "per_center", n_perm=args.n_perm, seed=args.seed)
    aa_protein = run_aa_tests(per_protein, "per_protein", n_perm=args.n_perm, seed=args.seed + 1)
    aa_all = pd.concat([aa_center, aa_protein], ignore_index=True)

    global_center = run_global_test(per_center, "per_center", n_perm=args.n_perm, seed=args.seed + 2)
    global_protein = run_global_test(per_protein, "per_protein", n_perm=args.n_perm, seed=args.seed + 3)
    global_all = pd.concat([global_center, global_protein], ignore_index=True)

    out_prefix = Path(args.out_prefix)
    aa_out = Path(f"{out_prefix}_aa_table.csv")
    global_out = Path(f"{out_prefix}_global_test.csv")
    supp_out = Path(f"{out_prefix}_supplement_table.csv")

    aa_all.to_csv(aa_out, index=False)
    global_all.to_csv(global_out, index=False)

    supp = aa_all[
        [
            "layer",
            "aa",
            "n_rows",
            "n_aerobic",
            "n_anaerobic",
            "mean_prop_aerobic",
            "mean_prop_anaerobic",
            "delta_prop_aero_minus_ana",
            "delta_clr_aero_minus_ana",
            "p_model_metabolism_adjusted",
            "q_model_fdr",
            "p_perm_delta_clr",
            "q_perm_fdr",
        ]
    ].copy()
    supp["significant_model_fdr_0.05"] = supp["q_model_fdr"] < 0.05
    supp["significant_perm_fdr_0.05"] = supp["q_perm_fdr"] < 0.05
    supp.to_csv(supp_out, index=False)

    print(f"Wrote: {aa_out}")
    print(f"Wrote: {global_out}")
    print(f"Wrote: {supp_out}")
    print("\\nTop per-protein hits by model FDR:")
    print(
        supp[(supp["layer"] == "per_protein")]
        .sort_values("q_model_fdr")
        .head(10)[["aa", "delta_prop_aero_minus_ana", "q_model_fdr", "q_perm_fdr"]]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()

