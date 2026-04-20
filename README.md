# toga2orthogroups

Build gene orthogroups from [TOGA2](https://github.com/hillerlab/TOGA) pairwise orthology annotations.  
Outputs a copy-number count table ready for downstream evolutionary analyses, including [CAFE5](https://github.com/hahnlab/CAFE5).

---

## Overview

TOGA2 produces pairwise orthology between a reference genome and each query species independently. This script integrates those pairwise mappings across all species into multi-species **orthogroups** (gene families), using a **Union-Find** graph algorithm.

Optionally, orthogroups can be extended with [PANTHER](https://www.pantherdb.org/) family assignments to capture paralogs that TOGA does not link directly.

---

## Algorithm

1. **Reference gene set** — Load reference transcript coordinates from a BED file and filter out sex-chromosome genes (chrX, chrY) to avoid copy-number biases.

2. **Per-species ortholog loading** — For each query species, read `loss_summary.tsv` to identify intact transcripts (status `I`, `FI`, `PI`; optionally `UL`) and `orthology_classification.tsv` to extract reference → query gene mappings.

3. **Union-Find grouping** — All reference genes start as singleton components. Within each species, an inverted index (query gene → reference genes) is built. Any two reference genes that share a query ortholog are merged into the same component. This is repeated across all species on the same Union-Find, so transitive merges propagate automatically across species.

4. **PANTHER merge (optional)** — If a PANTHER database is provided, additional union edges are added between reference genes in the same PANTHER protein family, capturing paralogs not linked by TOGA.

5. **Output** — Connected components are extracted as orthogroups. Each family is named by its lexicographically smallest member ID (ENSG in TOGA mode, PTHR in PANTHER mode). A copy-number count table and a full membership map are written to disk.

6. **Species QC** — Per-family z-scores of copy number are computed across species. Species with an anomalously high fraction of outlier families are flagged as potentially misassembled or mis-annotated.

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

QC:
  -z FLOAT,  --z-threshold FLOAT        z-score threshold for outlier detection  (default: 3.0)
  -of FLOAT, --outlier-fraction FLOAT   max outlier family fraction for flagging  (default: 0.05)
             --no-qc                    skip species QC diagnostics

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

---

## Loss status filtering

TOGA2 assigns each query transcript one of the following status codes:

| Status | Meaning | Included by default |
|--------|---------|-------------------|
| `I` | Intact | Yes |
| `FI` | Frameshift Intact | Yes |
| `PI` | Partial Intact | Yes |
| `UL` | Uncertain Loss | No (add `-ul` to include) |
| `L` | Lost | No |
| `M` | Missing | No |

Use `-ul` to be more permissive and recover genes that TOGA could not confidently call as lost.

---

## Species QC

After building orthogroups, a diagnostic report is printed to stderr. For each gene family, per-species copy numbers are z-scored relative to the cross-species distribution. A species is flagged (`***`) if it is an outlier in more than 5% of all families (adjustable with `-of`), which is indicative of assembly fragmentation or annotation artefacts.

```
=== Species QC Diagnostics ===
  z-score threshold: 3.0
  outlier family fraction threshold: 0.05
  one2one orthologs: 14203
  total families: 17511

  Species                   Outliers   Fraction   Flag
  ------------------------- -------- ---------- ------
  HLmarVanc2                     975     0.0557    ***

  WARNING: 1 species flagged with high orthogroup variance:
    - HLmarVanc2: 5.6% of families are outliers. Consider re-running without this species.
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
