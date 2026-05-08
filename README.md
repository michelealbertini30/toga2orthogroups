# toga2orthogroups

Build gene orthogroups from [TOGA2](https://github.com/hillerlab/TOGA) pairwise orthology annotations.  

---

## Overview

toga2orthogroups constructs multi-species **orthogroups** (gene families) from TOGA2 pairwise orthology calls by modelling reference-to-query gene relationships as a graph and extracting connected components via a **Union-Find** algorithm. Each reference gene is registered as a node; whenever two reference genes share a query ortholog across any species, they are merged into the same component. Union-Find handles this incrementally through path compression and rank-balanced union operations, yielding an effectively linear time complexity per operation, making the approach scalable to genome-wide analyses across hundreds of species. Optionally, a gene family database such as [PANTHER](https://www.pantherdb.org/) can be provided to further consolidate orthogroups by adding family-level edges between reference genes prior to component extraction. The output is readily formatted for direct input into gene family evolution tools such as CAFE5.

---

## Requirements

- Python ≥ 3.9
- No external dependencies — standard library only

---

## Input files

| Flag | Description |
|------|-------------|
| `-t` | Directory containing one subdirectory per query species, each with TOGA2 output (`loss_summary.tsv`, `orthology_classification.tsv`) |
| `-s` | Plain-text file with one species name per line (names must match subdirectory names in `-t`) |
| `-b` | Reference transcript BED file (used to exclude sex-chromosome genes) |
| `-i` | TOGA2 isoforms TSV mapping gene IDs to transcript IDs |
| `-o` | Output directory (created automatically if absent) |

---

## Usage

```
usage: toga2orthogroups.py [-h] -t DIR -s FILE -b FILE -i FILE -o DIR [options]

Build orthogroups from TOGA2 pairwise orthology annotations.

required:
  -t DIR,  --toga-dir DIR          directory with per-species TOGA2 output subdirs
  -s FILE, --species-list FILE     newline separated list of species
  -b FILE, --transcripts-bed FILE  reference transcript BED file
  -i FILE, --isoforms FILE         TOGA2 isoforms file
  -o DIR,  --out-dir DIR           output directory

optional:
  -f,      --force                 overwrite output files if they already exist
  -v,      --verbose               print per-species processing stats
  -ul,     --include-ul            include UL (Uncertain Loss) transcripts
           --panther FILE          PANTHER database TSV (replaces TOGA-only run)
           --one-to-one            write one2one.lst of reference genes with exactly
                                   1 ortholog in every species; skips PANTHER and QC

QC:
  -z FLOAT,  --z-threshold FLOAT   z-score threshold for outlier detection  (default: 3.0)
             --no-qc               skip species QC diagnostics

example:
  toga2orthogroups.py \
    -t TOGA2 \
    -s species.lst \
    -b toga2.transcripts.bed \
    -i toga2.isoforms.tsv \
    -o orthogroups
```

---

## Output files

### TOGA mode (default)

| File | Description |
|------|-------------|
| `TOGA2.orthogroups.tsv` | Copy-number count table — one row per orthogroup, one column per species |
| `TOGA2.ortho_map.tsv` | Full membership map — family ID followed by all `species\|query_gene` entries |

### PANTHER mode (`--panther`)

| File | Description |
|------|-------------|
| `PANTHER.orthogroups.tsv` | Same format as above; family IDs are PANTHER IDs (e.g. `PTHR12371`) |
| `PANTHER.ortho_map.tsv` | Full membership map with PANTHER family IDs |

### One-to-one mode (`--one-to-one`)

| File | Description |
|------|-------------|
| `one2one.lst` | Reference gene names (one per line) for families where every species has exactly 1 ortholog |

This mode skips PANTHER merging, the count table, the membership map, and QC. It is useful for extracting a high-confidence set of strictly conserved single-copy orthologs.

---

## Species QC

After building orthogroups, a diagnostic report is printed to stderr. For each species, the mean number of reference genes per query gene (`Mean R/Q`) is computed from that species' orthology data alone, independently of any other species. A value of 1.0 means every query gene maps to exactly one reference gene; higher values indicate that query genes span multiple reference genes, causing those reference genes to be merged into the same family and inflating copy-number counts for all other species.

A species is flagged (`***`) if its `Mean R/Q` z-score exceeds the threshold (default 3.0, one-sided), which is indicative of assembly fragmentation or annotation artefacts.

```
=== Running toga2orthogroups ===
...
=== Species QC Diagnostics ===
  z-score threshold (one-sided): 3.0
  mean refs/query across species: 1.0031
  one2one orthologs: 14203
  total families: 17511

  Species                   Mean R/Q    Total       Z   Flag
  ------------------------- --------- ------- ------- ------
  HLmarVanc2                   1.0891    7200   4.231    ***

  WARNING: 1 species flagged as potential orthogroup inflators:
    - HLmarVanc2: mean 1.0891 ref genes per query gene (z=4.23). Consider re-running without this species.

Time elapsed: 23.4s
```

Use `--no-qc` to suppress this report entirely.

---

## Use with CAFE5

The `TOGA2.orthogroups.tsv` (or `PANTHER.orthogroups.tsv`) output is directly compatible with [CAFE5](https://github.com/hahnlab/CAFE5) as the gene family count table. Pass it with the `-i` flag:

```bash
cafe5 -i orthogroups/TOGA2.orthogroups.tsv -t species.tree -o cafe_out
```

It is recommended to run the species QC step first and exclude any flagged species from both the count table and the species tree before running CAFE5, as assembly artefacts can inflate apparent gene family expansions and produce spurious evolutionary rate estimates.

---

## Author

**Michele Albertini** — [@michelealbertini30](https://github.com/michelealbertini30)
