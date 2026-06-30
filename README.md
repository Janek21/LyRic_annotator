# LyRic_annotator

LyRic wrapper for non-human genome annotation: sets up and runs [LyRic](https://github.com/guigolab/LyRic) per species, then evaluates the resulting annotation with gene/transcript counts and BUSCO. Built for an HPC/SLURM cluster.

Engine repository — long-read assembly + annotation, *unguided*. Used alongside [isoquant_annotator](https://github.com/Janek21/isoquant_annotator) (guided) and [geneid-training](https://github.com/Janek21/geneid-training) in a larger protist annotation pipeline.

## Overview

`lyric_execute.sh` chains the whole thing for one species, via SLURM job dependencies:

1. **`lyric_template.sh`** — clones the [LyRic_nonhuman](https://github.com/Janek21/LyRic_nonhuman) template, selects long-read SRA runs, writes the config, stages genome/annotation files, and submits the read-download array job.
2. **`scripts/runner.sh`** — runs LyRic itself (Snakemake), once downloads finish.
3. **`evaluation.sh`** — merges/cleans LyRic's output GFFs, infers CDS (TD2), counts gene/transcript models, and runs BUSCO (taxon-specific + eukaryote lineages).
4. **`scripts/merge_evaluation.sh`** — if a reference annotation exists for the species, merges it with the LyRic prediction (AGAT) and re-evaluates the merged result.

`turbo.sh` is a SLURM array job that runs `lyric_execute.sh` across a whole list of species.

## Repository structure

```
LyRic_annotator/
├── lyric_template.sh     # set up a species run + submit read downloads
├── lyric_execute.sh       # chain template → runner → evaluation → merge for one species
├── evaluation.sh          # merge/clean LyRic output, infer CDS, count models, run BUSCO
├── turbo.sh               # SLURM array: runs lyric_execute.sh over a species list
└── scripts/
    ├── LyRic_setup.py            # shortname / config / file_transfer / annotate_config
    ├── SRA_selector.py           # picks the best SRA runs for a species
    ├── annotation_fallback.py    # borrows a reference GFF from the closest related species
    ├── infer_cds.sh              # TD2 ORF prediction → longest-isoform proteome + CDS-augmented GFF
    ├── get_genetic_code.py       # resolves NCBI nuclear genetic code per taxon
    ├── get_busco_db.py           # resolves the BUSCO lineage per taxon
    ├── buscoPlot.py              # plots BUSCO completeness from result JSONs
    ├── make_summary_tables.py    # aggregates summary/ into final TSV reports
    ├── check_runs_status.sh      # reports each species' pipeline stage
    ├── runner.sh                 # SLURM job: runs the LyRic Snakemake pipeline
    ├── srr_dw.sh                 # SLURM array job: downloads SRA reads (ENA, falls back to SRA Toolkit)
    ├── merge_evaluation.sh       # merges LyRic + reference annotation, re-evaluates
    └── TD2/                      # vendored TransDecoder (TD2) fork, used by infer_cds.sh
```

## Requirements

- SLURM cluster
- conda env `buscomania`: AGAT, gffread, BUSCO, pigz
- a Snakemake environment (`~/bin/snakemake`) for running LyRic itself, plus `module load CMake`/`Python` as set up in `scripts/runner.sh`
- NCBI Entrez access (email / optional API key) for taxonomy lookups in `get_genetic_code.py` and `get_busco_db.py`
- `wget` / SRA Toolkit for read downloads
- Reference data one level up:
  `../data/species/<species_name>*/GC*/` (genome is mandatory, reference .gff is optional) and `../data/longread_protists.tsv`

## Usage

```bash
# run the full pipeline for one species (Genus_species[_extra] naming)
bash lyric_execute.sh <species_name> [longread_db] [busco_db]

# or run a batch of species as a SLURM array (one line per species in the list file)
sbatch turbo.sh <species_list.txt>

# check how far each species in a list got
bash scripts/check_runs_status.sh [species_list.txt] [base_dir]

# build the aggregate result tables once species have finished
python3 scripts/make_summary_tables.py            # regular evaluation
python3 scripts/make_summary_tables.py --merged   # LyRic + reference merge
python3 scripts/make_summary_tables.py --joint    # both, side by side
```

The individual stage scripts (`lyric_template.sh`, `evaluation.sh`, `scripts/merge_evaluation.sh`) can also be run/sbatch'd directly per species.

## Output

- `<species_name>/output/` — per-species working files (predicted proteome, CDS-augmented GFF, BUSCO run dirs)
- `summary/` — central, persistent results, independent of the per-species folders:
  - `counts/` — per-species gene/transcript counts + derived metrics (density, isoforms/gene, coding fraction)
  - `busco_lineage/`, `busco_eukaryote/` — BUSCO JSON results + plots
  - `pred/` — predicted (CDS-augmented) annotations
  - `merge/` — same structure, for species merged with a reference annotation
  - `counts_summary.tsv`, `busco_summary.tsv`, (+ `merge_*` / `joint_summary.tsv`) — particular summary tables from `make_summary_tables.py`
  - `general_summary.tsv` final aggregate table from `make_summary_tables.py`
