#!/usr/bin/env python3
"""
annotation_fallback.py

Finds the taxonomically closest species (by NCBI lineage overlap) that has a
GFF annotation file on disk, then copies it to the target species' data/input/
directory as Annotation.gff.
"""

import os
import sys
import csv
import gzip
import glob
import shutil
import argparse


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def folder_name(species_str: str) -> str:
    """'Plasmodium vivax' -> 'Plasmodium_vivax'  (first two words only)."""
    parts = species_str.strip().split()
    return "_".join(parts[:2])


def parse_lineage(lineage_str: str) -> list[str]:
    """'1;131567;2759;...' -> ['1','131567','2759',...]"""
    return lineage_str.strip().split(";")


def lineage_overlap(lin_a: list[str], lin_b: list[str]) -> int:
    """
    Count how many leading taxon IDs the two lineages share.
    A longer shared prefix means a closer taxonomic relationship.
    """
    score = 0
    for a, b in zip(lin_a, lin_b):
        if a == b:
            score += 1
        else:
            break
    return score


def find_gff(species_folder_name: str, species_root: str) -> str | None:
    """
    Searches for a GFF annotation file for the given species folder under species_root.
    Search order (stops at first hit):
      1. <species_root>/<species_folder>*/GCA*/*GCA*.gff
      2. <species_root>/<species_folder>*/GCA*/*GCA*.gff.gz
      3. <species_root>/<species_folder>*/GCF*/*GCF*.gff
      4. <species_root>/<species_folder>*/GCF*/*GCF*.gff.gz
    """
    folder_pattern = os.path.join(species_root, f"{species_folder_name}*")
    matched_folders = glob.glob(folder_pattern)

    if not matched_folders:
        return None

    assembly_prefixes = ["GCA", "GCF"]

    for species_dir in matched_folders:
        for prefix in assembly_prefixes:
            asm_dirs = glob.glob(os.path.join(species_dir, f"{prefix}*"))
            if not asm_dirs:
                continue

            for asm_dir in asm_dirs:
                for gff_glob in (f"*{prefix}*.gff", f"*{prefix}*.gff.gz"):
                    gff_matches = glob.glob(os.path.join(asm_dir, gff_glob))
                    if gff_matches:
                        return gff_matches[0]

    return None


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def load_candidates(tsv_path: str) -> dict[str, dict]:
    """
    Parse the TSV and return a dict keyed by normalised folder name.
    Keeps only the first occurrence of each species.
    """
    candidates: dict[str, dict] = {}
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 5:
                continue
            species_str = row[4].strip()
            lineage_str = row[3].strip()
            if not species_str or not lineage_str:
                continue
            key = folder_name(species_str)
            if key not in candidates:
                candidates[key] = {
                    "species": species_str,
                    "lineage": parse_lineage(lineage_str),
                }
    return candidates


def find_fallback_annotation(
    species_name: str,
    tsv_path: str,
    species_root: str,
    debug_top_n: int = 5,
) -> str | None:
    """
    Returns the GFF path of the closest related species that has one on disk,
    or None if nothing is found. Prints the top closest candidates.
    """
    if not os.path.isdir(species_root):
        return None

    candidates = load_candidates(tsv_path)
    target = candidates.get(species_name)
    target_lineage = target["lineage"] if target else []

    # Rank all candidates by lineage overlap
    ranked = sorted(
        [
            (key, meta)
            for key, meta in candidates.items()
            if key != species_name
        ],
        key=lambda kv: lineage_overlap(target_lineage, kv[1]["lineage"]),
        reverse=True,
    )

    # Print top N closest candidates
    print(f"Top {debug_top_n} closest candidates (by lineage overlap):")
    for i, (key, meta) in enumerate(ranked[:debug_top_n], 1):
        overlap = lineage_overlap(target_lineage, meta["lineage"])
        print(f"  #{i:02d}  overlap={overlap}  key='{key}'  species='{meta['species']}'")

    # Search for GFF
    for folder, meta in ranked:
        gff = find_gff(folder, species_root)
        if gff:
            return gff

    return None


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Copy the closest-relative annotation GFF when the target species lacks one."
    )
    parser.add_argument(
        "-s", "--species", required=True,
        help="Species folder name, e.g. Plasmodium_vivax"
    )
    parser.add_argument(
        "-d", "--database", required=True,
        help="Path to longread_protists.tsv"
    )
    parser.add_argument(
        "-r", "--species_root", default="../data/species",
        help="Root directory that contains per-species data folders (default: ../data/species)"
    )
    parser.add_argument(
        "-o", "--output", default=None,
        help="Destination path for the copied GFF (default: <species>/data/input/Annotation.gff)"
    )
    parser.add_argument(
        "-n", "--debug_n", type=int, default=5,
        help="Number of top candidates to show (default: 5)"
    )
    args = parser.parse_args()

    out_path = args.output or os.path.join(
        args.species, "data", "input", "Annotation.gff"
    )

    # Only run if the annotation is genuinely missing / empty
    if os.path.isfile(out_path) and os.path.getsize(out_path) > 0:
        sys.exit(0)

    gff = find_fallback_annotation(
        args.species, args.database, args.species_root, debug_top_n=args.debug_n
    )

    if gff is None:
        sys.exit(1)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if gff.endswith(".gz"):
        #decompress a gzipped source onto the plain Annotation.gff name
        tmp_out = f"{out_path}.alt"
        with gzip.open(gff, "rb") as f_in:
            with open(tmp_out, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.rename(tmp_out, out_path)
    else:
        shutil.copy(gff, out_path)
    
    # Print the copy operation details
    print(f"\nCopied annotation from: {gff}")
    print(f"Copied annotation to:   {out_path}")


if __name__ == "__main__":
    main()