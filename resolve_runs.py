#!/usr/bin/env python3
"""
resolve_runs.py
===============

Turn a user‑supplied list of ENA accessions (study, sample or run)
into three files:

    resolved_accessions.tsv   sample_accession<TAB>run_accession
    available_runs.txt        every run with downloadable FASTQs
    missing_inputs.txt        original IDs that yielded no runs

Run once, then point the Snakefile to resolved_accessions.tsv.

Usage
-----
    python resolve_runs.py samples.txt
"""

import csv
import pathlib
import sys
from collections import defaultdict
from typing import List, Dict, Set

import requests

ENA_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={acc}&result=read_run"
    "&fields=run_accession,sample_accession,fastq_ftp"
    "&format=tsv"
)


def fetch_run_rows(acc: str) -> List[Dict[str, str]]:
    """Return a list of dicts (one per run) for any ENA accession."""
    try:
        txt = requests.get(ENA_URL.format(acc=acc), timeout=30).text
    except requests.exceptions.RequestException:
        return []
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    # some endpoints return header only when accession not found
    if rows and not rows[0].get("run_accession"):
        return []
    return rows


def main(samples_file: str):
    inputs = [l.strip() for l in open(samples_file) if l.strip()]

    mapping: Dict[str, Set[str]] = defaultdict(set)  # sample -> set(run)
    missing: List[str] = []

    for inp in inputs:
        rows = fetch_run_rows(inp)
        if not rows:
            missing.append(inp)
            continue

        for row in rows:
            run = row["run_accession"]
            sample = row["sample_accession"] or run   # fallback 1:1
            if row["fastq_ftp"]:
                mapping[sample].add(run)

    # ------------------------------------------------------------------
    # write outputs
    # ------------------------------------------------------------------
    pathlib.Path("resolved_accessions.tsv").write_text(
        "\n".join(f"{s}\t{r}" for s, runs in mapping.items() for r in sorted(runs))
        + "\n"
    )

    all_runs = sorted({r for runs in mapping.values() for r in runs})
    pathlib.Path("available_runs.txt").write_text("\n".join(all_runs) + "\n")
    pathlib.Path("missing_inputs.txt").write_text("\n".join(missing) + "\n")

    print(
        f"✓ {len(all_runs)} run(s) across {len(mapping)} sample(s) "
        f"→ resolved_accessions.tsv / available_runs.txt"
    )
    print(f"✗ {len(missing)} input IDs lacked downloadable reads "
          f"→ missing_inputs.txt")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: resolve_runs.py samples.txt")
    main(sys.argv[1])

