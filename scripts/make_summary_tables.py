#!/usr/bin/env python
"""
Build the top-level summary tables from the per-species evaluation outputs.

Modes (run after all species finish):
  regular (default): reads summary/{counts,busco_lineage,busco_eukaryote}
                     writes summary/{counts_summary,busco_summary,general_summary}.tsv
  --merged:          reads summary/merge/{counts,busco_lineage,busco_eukaryote}
                     writes summary/{merge_counts_summary,merge_busco_summary,merge_general_summary}.tsv
  --joint:           combines both general tables side by side (columns suffixed
                     _regular / _merged) into summary/joint_summary.tsv

Each mode produces three tables:
  1. counts_summary   - species + counts + derived metrics (gene/transcript
                        density, isoforms per gene, coding fraction)
                        (from <base>/counts/*_metrics.tsv, one row per species)
  2. busco_summary    - species, lineage_used, lineage_completeness,
                        eukaryote_completeness  (BUSCO "Complete %" C, no % sign)
                        (from <base>/busco_lineage/*_Lbusco.json and
                              <base>/busco_eukaryote/*_Ebusco.json)
  3. general_summary  - the two tables above merged on species

The shared key in every table is the bare species name (Genus_species[_extra]):
counts/busco stems drop the trailing _<taxonID>.

Run from the repo root (needs pandas; run inside the buscomania conda env):
    python3 scripts/make_summary_tables.py            # regular evaluation
    python3 scripts/make_summary_tables.py --merged   # merged (LyRic + reference) evaluation
    python3 scripts/make_summary_tables.py --joint    # regular + merged combined
"""

import os
import re
import glob
import json
import argparse
import pandas as pd

SUMMARY_DIR = "summary"

# counts + derived metrics, one row per species (written by evaluation.sh /
# merge_evaluation.sh). See isoquant_annotator/derived_metrics.md.
METRIC_COLUMNS = [
    "species",
    "gene_count",
    "transcript_count",
    "genome_size_bp",
    "coding_transcripts",
    "transcriptome_transcripts",
    "gene_density_per_mb",
    "transcript_density_per_mb",
    "isoforms_per_gene",
    "coding_fraction",
]


def species_from_stem(stem):
    """<species_name>_<taxonID> -> <species_name> (taxon id is the trailing _<digits>)."""
    m = re.match(r"^(.+)_[0-9]+$", stem)
    return m.group(1) if m else stem


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

def build_counts_table(counts_dir):
    """Concatenate every per-species _metrics.tsv; key on the bare species name."""
    frames = []
    for path in sorted(glob.glob(os.path.join(counts_dir, "*_metrics.tsv"))):
        stem = os.path.basename(path)[: -len("_metrics.tsv")]
        try:
            df = pd.read_csv(path, sep="\t", dtype=str)
        except (pd.errors.EmptyDataError, FileNotFoundError):
            continue
        if df.empty:
            continue
        df = df.reindex(columns=METRIC_COLUMNS)
        df["species"] = species_from_stem(stem)  # bare species key for the joins
        frames.append(df)
    if not frames:
        return pd.DataFrame(columns=METRIC_COLUMNS)
    return pd.concat(frames, ignore_index=True).reindex(columns=METRIC_COLUMNS)


def build_busco_table(busco_lineage_dir, busco_euk_dir):
    rows = {}
    for jp in sorted(glob.glob(os.path.join(busco_lineage_dir, "*_Lbusco.json"))):
        stem = os.path.basename(jp)[: -len("_Lbusco.json")]
        sp = species_from_stem(stem)
        completeness, lineage = busco_completeness(jp)
        row = rows.setdefault(sp, {"species": sp})
        row["lineage_used"] = lineage
        row["lineage_completeness"] = completeness
    for jp in sorted(glob.glob(os.path.join(busco_euk_dir, "*_Ebusco.json"))):
        stem = os.path.basename(jp)[: -len("_Ebusco.json")]
        sp = species_from_stem(stem)
        completeness, _ = busco_completeness(jp)
        row = rows.setdefault(sp, {"species": sp})
        row["eukaryote_completeness"] = completeness
    cols = ["species", "lineage_used", "lineage_completeness", "eukaryote_completeness"]
    return pd.DataFrame(list(rows.values()), columns=cols)


# ---------------------------------------------------------------------------

def build_general(base_dir):
    """Return the counts+busco tables (and their outer-join) for one evaluation base dir."""
    counts_df = build_counts_table(os.path.join(base_dir, "counts"))
    busco_df = build_busco_table(os.path.join(base_dir, "busco_lineage"),
                                 os.path.join(base_dir, "busco_eukaryote"))
    general = counts_df.merge(busco_df, on="species", how="outer")
    return counts_df, busco_df, general


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--merged", action="store_true",
                       help="Summarize the merged (LyRic + reference) evaluation instead of the regular one")
    group.add_argument("--joint", action="store_true",
                       help="Combine regular and merged evaluations into one side-by-side table (summary/joint_summary.tsv)")
    args = parser.parse_args()

    os.makedirs(SUMMARY_DIR, exist_ok=True)
    regular_base = SUMMARY_DIR
    merged_base = os.path.join(SUMMARY_DIR, "merge")

    if args.joint:
        #join the two general tables on species; suffix shared columns so both sides stay visible
        _, _, regular = build_general(regular_base)
        _, _, merged = build_general(merged_base)
        joint = regular.merge(merged, on="species", how="outer",
                              suffixes=("_regular", "_merged"))
        joint = joint.sort_values("species").fillna("NA")
        joint_out = os.path.join(SUMMARY_DIR, "joint_summary.tsv")
        joint.to_csv(joint_out, sep="\t", index=False)
        print(f"Wrote {joint_out} ({len(joint)} species)")
        return

    #regular reads summary/<sub>; merged reads summary/merge/<sub> and prefixes the outputs
    base_dir = merged_base if args.merged else regular_base
    out_prefix = "merge_" if args.merged else ""

    counts_df, busco_df, general = build_general(base_dir)

    counts_out = os.path.join(SUMMARY_DIR, f"{out_prefix}counts_summary.tsv")
    busco_out = os.path.join(SUMMARY_DIR, f"{out_prefix}busco_summary.tsv")
    general_out = os.path.join(SUMMARY_DIR, f"{out_prefix}general_summary.tsv")

    counts_df.to_csv(counts_out, sep="\t", index=False)
    busco_df.to_csv(busco_out, sep="\t", index=False)
    print(f"Wrote {counts_out} ({len(counts_df)} species)")
    print(f"Wrote {busco_out} ({len(busco_df)} species)")

    general = general.sort_values("species").fillna("NA")
    general.to_csv(general_out, sep="\t", index=False)
    print(f"Wrote {general_out} ({len(general)} species)")


if __name__ == "__main__":
    main()
