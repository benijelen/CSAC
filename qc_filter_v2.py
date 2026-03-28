#!/usr/bin/env python3
"""QC filtering pipeline for Fe-binding v2 master table.

This script applies metadata QC and deduplication to a combined aerobic/anaerobic CSV.
It writes pass/fail tables plus explicit notes for deduplicated rows.

Default input: combined_aerobic_anaerobic_master_review.csv
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


MISSING_TOKENS = {"", "na", "n/a", "nan", "none", "null", "unknown"}


@dataclass
class RowWrap:
    idx: int
    row: Dict[str, str]
    reasons: List[str]


def norm(s: Optional[str]) -> str:
    if s is None:
        return ""
    return str(s).strip()


def low(s: Optional[str]) -> str:
    return norm(s).lower()


def is_missing(s: Optional[str]) -> bool:
    return low(s) in MISSING_TOKENS


def parse_resolution(s: Optional[str]) -> Optional[float]:
    x = low(s)
    if is_missing(x):
        return None
    try:
        return float(x)
    except ValueError:
        return None


def method_bucket(method: str) -> str:
    m = low(method)
    if "x-ray" in m:
        return "xray"
    if "electron microscopy" in m:
        return "em"
    if "nmr" in m:
        return "nmr"
    if "neutron" in m:
        return "neutron"
    if "electron crystallography" in m:
        return "electron_crystallography"
    return "other"


def yes_no_unknown(s: Optional[str]) -> str:
    v = low(s)
    if v in {"yes", "y", "true", "1"}:
        return "yes"
    if v in {"no", "n", "false", "0"}:
        return "no"
    return "unknown"


def system_name(row: Dict[str, str]) -> str:
    for key in [
        "canonical_label",
        "canonical_gene_or_label",
        "enzyme_or_complex",
        "enzyme_class_or_family",
    ]:
        if key in row and not is_missing(row.get(key, "")):
            return norm(row[key])
    return "UNLABELED"


def normalize_key(s: str) -> str:
    return " ".join(low(s).replace("_", " ").split())


def method_rank(method: str) -> int:
    b = method_bucket(method)
    return {
        "xray": 0,
        "em": 1,
        "nmr": 2,
        "neutron": 3,
        "electron_crystallography": 4,
        "other": 5,
    }.get(b, 9)


def priority_rank(priority: str) -> int:
    p = low(priority)
    return {"high": 0, "medium": 1, "low": 2}.get(p, 9)


def find_structure_file(pdb_id: str, dirs: List[Path]) -> Optional[Path]:
    pid = norm(pdb_id).lower()
    if is_missing(pid) or pid == "tbd_computed_model":
        return None
    names = [f"{pid}.pdb", f"{pid}.cif", f"{pid}.mmcif"]
    for d in dirs:
        for n in names:
            p = d / n
            if p.exists():
                return p
    return None


def pdb_altloc_occupancy_qc(path: Path) -> Tuple[str, str]:
    """Return (status, note) for PDB text files only.

    status: pass | review | not_run_unsupported_format
    """
    suffix = path.suffix.lower()
    if suffix != ".pdb":
        return "not_run_unsupported_format", f"structure file found ({path.name}) but parser only checks .pdb altloc/occupancy"

    altloc_count = 0
    low_occ_count = 0
    total_atom = 0
    try:
        with path.open() as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                total_atom += 1
                if len(line) >= 17:
                    altloc = line[16].strip()
                    if altloc and altloc != "A":
                        altloc_count += 1
                if len(line) >= 60:
                    occ_raw = line[54:60].strip()
                    try:
                        occ = float(occ_raw)
                        if occ < 1.0:
                            low_occ_count += 1
                    except ValueError:
                        pass
    except OSError as exc:
        return "review", f"failed to parse {path.name}: {exc}"

    if total_atom == 0:
        return "review", f"no ATOM/HETATM records in {path.name}"

    if altloc_count == 0 and low_occ_count == 0:
        return "pass", f"no altloc/partial-occupancy flags detected in {path.name}"

    return "review", f"altloc_count={altloc_count}; occupancy_lt1_count={low_occ_count} in {path.name}"


def choose_best(rows: List[RowWrap]) -> RowWrap:
    def sort_key(w: RowWrap):
        r = w.row
        mut = yes_no_unknown(r.get("mutation_reported"))
        eng = yes_no_unknown(r.get("engineered_reported"))
        mut_rank = {"no": 0, "unknown": 1, "yes": 2}[mut]
        eng_rank = {"no": 0, "unknown": 1, "yes": 2}[eng]
        res = parse_resolution(r.get("resolution_a"))
        res_rank = res if res is not None else 999.0
        return (
            mut_rank,
            eng_rank,
            method_rank(r.get("method", "")),
            res_rank,
            priority_rank(r.get("priority", "")),
            w.idx,
        )

    return sorted(rows, key=sort_key)[0]


def main() -> None:
    ap = argparse.ArgumentParser(description="QC filter v2 for combined Fe master table")
    ap.add_argument("--input", default="combined_aerobic_anaerobic_master_review.csv")
    ap.add_argument("--out-pass", default="qc_pass.csv")
    ap.add_argument("--out-fail", default="qc_fail.csv")
    ap.add_argument("--out-summary", default="qc_summary_report.csv")
    ap.add_argument("--out-dedup-notes", default="qc_dedup_notes.csv")
    ap.add_argument("--xray-max", type=float, default=3.0)
    ap.add_argument("--em-max", type=float, default=3.5)
    ap.add_argument("--require-experimental", action="store_true", default=True)
    ap.add_argument("--no-require-experimental", dest="require_experimental", action="store_false")
    ap.add_argument("--exclude-mutant-engineered", action="store_true", default=True)
    ap.add_argument("--no-exclude-mutant-engineered", dest="exclude_mutant_engineered", action="store_false")
    ap.add_argument("--dedup-exact", action="store_true", default=True)
    ap.add_argument("--no-dedup-exact", dest="dedup_exact", action="store_false")
    ap.add_argument("--dedup-near", action="store_true", default=True)
    ap.add_argument("--no-dedup-near", dest="dedup_near", action="store_false")
    ap.add_argument(
        "--structure-dirs",
        default="baseline_v1_original,baseline_v1_frozen_2026-03-06",
        help="comma-separated dirs with local .pdb/.cif/.mmcif files for optional altloc/occupancy checks",
    )
    args = ap.parse_args()

    in_path = Path(args.input)
    with in_path.open() as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]

    output_fields = list(fieldnames)
    for k in [
        "system_name",
        "qc_drop_reasons",
        "dedup_stage",
        "dedup_key",
        "kept_row_idx",
        "structure_qc_status",
        "structure_qc_note",
    ]:
        if k not in output_fields:
            output_fields.append(k)

    struct_dirs = [Path(x.strip()) for x in args.structure_dirs.split(",") if x.strip()]

    # Stage 1: metadata filters
    candidates: List[RowWrap] = []
    fail_rows: List[RowWrap] = []

    for i, row in enumerate(rows, start=1):
        reasons: List[str] = []
        pdb_id = norm(row.get("pdb_id"))
        method = norm(row.get("method"))
        mut = yes_no_unknown(row.get("mutation_reported"))
        eng = yes_no_unknown(row.get("engineered_reported"))
        res = parse_resolution(row.get("resolution_a"))

        row["system_name"] = system_name(row)

        if is_missing(pdb_id) or low(pdb_id) == "tbd_computed_model":
            reasons.append("missing_or_computed_structure")

        if args.require_experimental:
            if is_missing(method) or "computed model" in low(method):
                reasons.append("non_experimental_method")

        if args.exclude_mutant_engineered:
            if mut == "yes":
                reasons.append("mutant_structure")
            if eng == "yes":
                reasons.append("engineered_structure")

        mb = method_bucket(method)
        if mb == "xray":
            if res is None:
                reasons.append("missing_resolution_xray")
            elif res > args.xray_max:
                reasons.append("resolution_above_xray_threshold")
        elif mb == "em":
            if res is None:
                reasons.append("missing_resolution_em")
            elif res > args.em_max:
                reasons.append("resolution_above_em_threshold")

        if reasons:
            fail_rows.append(RowWrap(i, row, reasons))
        else:
            candidates.append(RowWrap(i, row, []))

    dedup_notes: List[Dict[str, str]] = []

    # Stage 2: exact dedup (same oxygen + system + pdb)
    if args.dedup_exact:
        grouped: Dict[Tuple[str, str, str], List[RowWrap]] = defaultdict(list)
        for w in candidates:
            key = (
                normalize_key(w.row.get("combined_oxygen_group", "")),
                normalize_key(w.row.get("system_name", "")),
                normalize_key(w.row.get("pdb_id", "")),
            )
            grouped[key].append(w)

        next_candidates: List[RowWrap] = []
        for key, group in grouped.items():
            if len(group) == 1:
                group[0].row["dedup_stage"] = "exact_keep"
                group[0].row["dedup_key"] = "|".join(key)
                next_candidates.append(group[0])
                continue
            keep = choose_best(group)
            keep.row["dedup_stage"] = "exact_keep"
            keep.row["dedup_key"] = "|".join(key)
            next_candidates.append(keep)
            for w in group:
                if w is keep:
                    continue
                w.reasons.append("dedup_exact_system_pdb")
                w.row["dedup_stage"] = "exact_drop"
                w.row["dedup_key"] = "|".join(key)
                w.row["kept_row_idx"] = str(keep.idx)
                fail_rows.append(w)
                dedup_notes.append(
                    {
                        "dedup_type": "exact_system_pdb",
                        "dedup_key": "|".join(key),
                        "kept_row_idx": str(keep.idx),
                        "dropped_row_idx": str(w.idx),
                        "kept_pdb_id": keep.row.get("pdb_id", ""),
                        "dropped_pdb_id": w.row.get("pdb_id", ""),
                        "system_name": keep.row.get("system_name", ""),
                        "organism": keep.row.get("organism", ""),
                        "reason": "same oxygen_group+system_name+pdb_id",
                    }
                )
        candidates = next_candidates

    # Stage 3: near-identical dedup (same oxygen + system + organism)
    if args.dedup_near:
        grouped2: Dict[Tuple[str, str, str], List[RowWrap]] = defaultdict(list)
        for w in candidates:
            key = (
                normalize_key(w.row.get("combined_oxygen_group", "")),
                normalize_key(w.row.get("system_name", "")),
                normalize_key(w.row.get("organism", "")),
            )
            grouped2[key].append(w)

        next_candidates2: List[RowWrap] = []
        for key, group in grouped2.items():
            if len(group) == 1 or is_missing(key[2]):
                group[0].row["dedup_stage"] = group[0].row.get("dedup_stage", "near_keep")
                group[0].row["dedup_key"] = group[0].row.get("dedup_key", "|".join(key))
                next_candidates2.extend(group)
                continue

            keep = choose_best(group)
            keep.row["dedup_stage"] = "near_keep"
            keep.row["dedup_key"] = "|".join(key)
            next_candidates2.append(keep)
            for w in group:
                if w is keep:
                    continue
                w.reasons.append("dedup_near_identical_system_organism")
                w.row["dedup_stage"] = "near_drop"
                w.row["dedup_key"] = "|".join(key)
                w.row["kept_row_idx"] = str(keep.idx)
                fail_rows.append(w)
                dedup_notes.append(
                    {
                        "dedup_type": "near_identical_system_organism",
                        "dedup_key": "|".join(key),
                        "kept_row_idx": str(keep.idx),
                        "dropped_row_idx": str(w.idx),
                        "kept_pdb_id": keep.row.get("pdb_id", ""),
                        "dropped_pdb_id": w.row.get("pdb_id", ""),
                        "system_name": keep.row.get("system_name", ""),
                        "organism": keep.row.get("organism", ""),
                        "reason": "same oxygen_group+system_name+organism",
                    }
                )
        candidates = next_candidates2

    # Stage 4: optional structure altloc/occupancy screening by local files
    for w in candidates:
        pdb_id = w.row.get("pdb_id", "")
        fpath = find_structure_file(pdb_id, struct_dirs)
        if fpath is None:
            w.row["structure_qc_status"] = "pending_no_structure_file"
            w.row["structure_qc_note"] = "altloc/occupancy not checked (structure file not found locally)"
            continue
        status, note = pdb_altloc_occupancy_qc(fpath)
        w.row["structure_qc_status"] = status
        w.row["structure_qc_note"] = note

    # Write pass/fail
    pass_rows = [w.row for w in candidates]
    fail_out_rows = []
    for w in fail_rows:
        rr = dict(w.row)
        rr["qc_drop_reasons"] = ";".join(sorted(set(w.reasons)))
        if not rr.get("system_name"):
            rr["system_name"] = system_name(rr)
        if not rr.get("dedup_stage"):
            rr["dedup_stage"] = "base_filter_drop"
        fail_out_rows.append(rr)

    with Path(args.out_pass).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=output_fields)
        w.writeheader()
        for r in pass_rows:
            if not r.get("qc_drop_reasons"):
                r["qc_drop_reasons"] = ""
            w.writerow(r)

    with Path(args.out_fail).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=output_fields)
        w.writeheader()
        for r in fail_out_rows:
            w.writerow(r)

    # Write dedup notes
    dedup_fields = [
        "dedup_type",
        "dedup_key",
        "kept_row_idx",
        "dropped_row_idx",
        "kept_pdb_id",
        "dropped_pdb_id",
        "system_name",
        "organism",
        "reason",
    ]
    with Path(args.out_dedup_notes).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=dedup_fields)
        w.writeheader()
        for d in dedup_notes:
            w.writerow(d)

    # Write summary
    fail_reason_counter = Counter()
    for r in fail_out_rows:
        for reason in norm(r.get("qc_drop_reasons", "")).split(";"):
            if reason:
                fail_reason_counter[reason] += 1

    pass_by_group = Counter(low(r.get("combined_oxygen_group", "")) for r in pass_rows)
    fail_by_group = Counter(low(r.get("combined_oxygen_group", "")) for r in fail_out_rows)

    with Path(args.out_summary).open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["input_rows", len(rows)])
        w.writerow(["qc_pass_rows", len(pass_rows)])
        w.writerow(["qc_fail_rows", len(fail_out_rows)])
        w.writerow(["pass_aerobic", pass_by_group.get("aerobic", 0)])
        w.writerow(["pass_anaerobic", pass_by_group.get("anaerobic", 0)])
        w.writerow(["fail_aerobic", fail_by_group.get("aerobic", 0)])
        w.writerow(["fail_anaerobic", fail_by_group.get("anaerobic", 0)])
        w.writerow(["dedup_drops", sum(1 for x in fail_out_rows if "dedup_" in x.get("qc_drop_reasons", ""))])
        w.writerow([])
        w.writerow(["drop_reason", "count"])
        for reason, count in fail_reason_counter.most_common():
            w.writerow([reason, count])

    print(f"Input rows: {len(rows)}")
    print(f"Pass rows: {len(pass_rows)}")
    print(f"Fail rows: {len(fail_out_rows)}")
    print(f"Dedup note rows: {len(dedup_notes)}")
    print(f"Wrote: {args.out_pass}, {args.out_fail}, {args.out_summary}, {args.out_dedup_notes}")


if __name__ == "__main__":
    main()
