#!/usr/bin/env python
"""
Build the top-level summary tables from the per-species evaluation outputs.

Produces four tables in summary/:
  1. counts_summary.tsv   - species, gene_count, transcript_count
                            (from summary/counts/*_gc.txt and *_tc.txt)
  2. busco_summary.tsv     - species, lineage_used, lineage_completeness,
                            eukaryote_completeness  (BUSCO "Complete %" C, no % sign)
                            (from summary/busco_lineage/*_Lbusco.json and
                                  summary/busco_eukaryote/*_Ebusco.json)
  3. sum_stats.csv         - gffcompare metrics, one row per species
                            (from summary/gffcmp/*.stats via statsCompression.gffstats_2data)
  4. general_summary.tsv   - the three tables above merged on species

The shared key in every table is the bare species name (Genus_species[_extra]):
counts/busco stems drop the trailing _<taxonID>; gffcompare names drop the -Lycmp suffix.

Run from the repo root (needs pandas; run inside the buscomania conda env):
    python3 scripts/make_summary_tables.py
"""

import os
import re
import sys
import glob
import json
import pandas as pd

# reuse the existing gffcompare .stats parser
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from statsCompression import gffstats_2data

SUMMARY_DIR = "summary"
COUNTS_DIR = os.path.join(SUMMARY_DIR, "counts")
BUSCO_LINEAGE_DIR = os.path.join(SUMMARY_DIR, "busco_lineage")
BUSCO_EUK_DIR = os.path.join(SUMMARY_DIR, "busco_eukaryote")
GFFCMP_DIR = os.path.join(SUMMARY_DIR, "gffcmp")


def species_from_stem(stem):
    """<species_name>_<taxonID> -> <species_name> (taxon id is the trailing _<digits>)."""
    m = re.match(r"^(.+)_[0-9]+$", stem)
    return m.group(1) if m else stem


def read_value(path):
    try:
        with open(path) as fh:
            return fh.read().strip()
    except FileNotFoundError:
        return "NA"


def busco_completeness(path):
    """Return (Complete% without % sign, lineage dataset name) from a BUSCO short summary json."""
    with open(path) as fh:
        data = json.load(fh)
    results = data.get("results", {})
    completeness = results.get("Complete percentage")
    if completeness is None:  # fall back to parsing the one-line summary string
        m = re.search(r"C:\s*([0-9.]+)%", results.get("one_line_summary", ""))
        completeness = m.group(1) if m else "NA"
    lineage = data.get("lineage_dataset", {}).get("name", "NA")
    return completeness, lineage


# ---------------------------------------------------------------------------
# Table builders
# ---------------------------------------------------------------------------

def build_counts_table():
    rows = []
    for gc in sorted(glob.glob(os.path.join(COUNTS_DIR, "*_gc.txt"))):
        stem = os.path.basename(gc)[: -len("_gc.txt")]
        tc = os.path.join(COUNTS_DIR, f"{stem}_tc.txt")
        rows.append({
            "species": species_from_stem(stem),
            "gene_count": read_value(gc),
            "transcript_count": read_value(tc),
        })
    return pd.DataFrame(rows, columns=["species", "gene_count", "transcript_count"])


def build_busco_table():
    rows = {}
    for jp in sorted(glob.glob(os.path.join(BUSCO_LINEAGE_DIR, "*_Lbusco.json"))):
        stem = os.path.basename(jp)[: -len("_Lbusco.json")]
        sp = species_from_stem(stem)
        completeness, lineage = busco_completeness(jp)
        row = rows.setdefault(sp, {"species": sp})
        row["lineage_used"] = lineage
        row["lineage_completeness"] = completeness
    for jp in sorted(glob.glob(os.path.join(BUSCO_EUK_DIR, "*_Ebusco.json"))):
        stem = os.path.basename(jp)[: -len("_Ebusco.json")]
        sp = species_from_stem(stem)
        completeness, _ = busco_completeness(jp)
        row = rows.setdefault(sp, {"species": sp})
        row["eukaryote_completeness"] = completeness
    cols = ["species", "lineage_used", "lineage_completeness", "eukaryote_completeness"]
    return pd.DataFrame(list(rows.values()), columns=cols)


def build_gffcmp_table():
    stats_files = sorted(glob.glob(os.path.join(GFFCMP_DIR, "*.stats")))
    rows = [gffstats_2data(f) for f in stats_files]
    df = pd.DataFrame(rows).fillna("NA")
    if "Species" in df.columns:
        # gffcompare names the run <species_name>-Lycmp; normalise to the bare species key
        df.insert(0, "species", df["Species"].str.replace(r"-Lycmp.*$", "", regex=True))
        df = df.drop(columns=["Species"])
    elif "species" not in df.columns:
        df["species"] = pd.Series(dtype="object")
    return df


# ---------------------------------------------------------------------------

def main():
    os.makedirs(SUMMARY_DIR, exist_ok=True)

    counts_df = build_counts_table()
    busco_df = build_busco_table()
    gffcmp_df = build_gffcmp_table()

    counts_out = os.path.join(SUMMARY_DIR, "counts_summary.tsv")
    busco_out = os.path.join(SUMMARY_DIR, "busco_summary.tsv")
    gffcmp_out = os.path.join(SUMMARY_DIR, "sum_stats.csv")
    general_out = os.path.join(SUMMARY_DIR, "general_summary.tsv")

    counts_df.to_csv(counts_out, sep="\t", index=False)
    busco_df.to_csv(busco_out, sep="\t", index=False)
    gffcmp_df.to_csv(gffcmp_out, sep="\t", index=False)
    print(f"Wrote {counts_out} ({len(counts_df)} species)")
    print(f"Wrote {busco_out} ({len(busco_df)} species)")
    print(f"Wrote {gffcmp_out} ({len(gffcmp_df)} species)")

    # 4th table: merge the three on species (outer join keeps every species seen anywhere)
    general = counts_df.merge(busco_df, on="species", how="outer") \
                       .merge(gffcmp_df, on="species", how="outer")
    general = general.sort_values("species").fillna("NA")
    general.to_csv(general_out, sep="\t", index=False)
    print(f"Wrote {general_out} ({len(general)} species)")


if __name__ == "__main__":
    main()
