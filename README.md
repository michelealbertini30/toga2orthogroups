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
| `-b` | Reference transcript BED file (chromosomes/scaffolds are filtered via `-bl`) |
| `-i` | TOGA2 isoforms TSV mapping gene IDs to transcript IDs |
| `-o` | Output directory (created automatically if absent) |

---

## Usage

```
usage: toga2orthogroups.py [-h] -t DIR -s FILE -b FILE -i FILE -o DIR [options]

Build orthogroups from TOGA2 pairwise orthology annotations.

required:
  -t DIR,  --toga-dir DIR          Directory with per-species TOGA2 output subdirs
  -s FILE, --species-list FILE     Newline separated list of species
  -b FILE, --transcripts-bed FILE  Reference transcript BED file
  -i FILE, --isoforms FILE         TOGA2 isoforms file
  -o DIR,  --out-dir DIR           Output directory

optional:
  -f,      --force                 Overwrite output files if they already exist
  -v,      --verbose               Print per-species processing stats
  -ul,     --include-ul            Include UL (Uncertain Loss) transcripts
  -sf INT, --size-filter INT       Remove gene families where any species has >= INT gene copies
  -bl LIST|FILE, --blacklist       Chromosomes/scaffolds to exclude: comma-separated list
                                   (e.g. chrX,chrY) or path to a file with one name per line
           --panther FILE          PANTHER database flat file; enables PANTHER-guided family merging
           --one-to-one            Per-reference-gene matrix with the query gene name if
                                   single-copy, dash otherwise

QC:
             --ortho-z FLOAT       OrthoZ threshold for orthogroup inflation detection  (default: 3.0)
             --fam-z FLOAT         FamZ threshold for family copy-number outlier detection  (default: 3.0)
             --no-qc               Skip species QC diagnostics

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
| `orthogroups_matrix.tsv` | Copy-number count table — `Desc`, family ID, then one copy-count column per species (CAFE5-compatible) |
| `orthogroups_map.tsv` | Full membership map — family ID followed by all `species\|query_gene` entries |
| `one2one_matrix.tsv` | Per-reference-gene matrix (`--one-to-one`) — query gene name when single-copy, `-` otherwise |

### PANTHER mode (`--panther`)

Same files as TOGA mode; family IDs are replaced by PANTHER IDs (e.g. `PTHR12371`).


---

## One-to-one mode (`--one-to-one`)

`--one-to-one` is a separate, simpler mode that bypasses Union-Find entirely and works **per reference gene**: for each reference gene it compiles how many query orthologs each species has.

**Output** — `one2one_matrix.tsv`:

```
Ref_ID           HLammLeuc1            HLcalLat1A            HLcynGun1             HLictTrid4A           …
ENSG00000256043  ENSG00000256043       ENSG00000256043       ENSG00000256043       ENSG00000256043       …
ENSG00000256053  ENSG00000256053       ENSG00000256053       ENSG00000256053       ENSG00000256053       …
ENSG00000256061  ENSG00000256061       ENSG00000256061       ENSG00000256061       ENSG00000256061       …
ENSG00000256087  -                     -                     -                     ENSG00000256087       …
ENSG00000256188  ENSG00000256188_1+    ENSG00000256188_1+    -                     ENSG00000256188_1+    …
```

- Each row is one reference gene.
- Each cell contains the query gene ID when that species has **exactly one** ortholog for that reference gene, or `-` when it has zero or more than one.
- Rows where fewer than two species have a single-copy ortholog are dropped.

**Use case**: derive a list of directly comparable one-to-one orthologous genes for comparative genomics analyses.

**Notes**:
- `--panther` and `--size-filter` are both ignored in this mode (a warning is printed for each).
- Species QC diagnostics are not run.
- The output is not intended as CAFE5 input.

---

## Species QC

After building orthogroups, a diagnostic report is printed to stderr. Two complementary signals are computed per species:

**OrthoZ — spanning-rate signal**  
For each qualifying family, the spanning rate is computed as the sum of query genes spanning multiple reference genes, divided by the total number of query genes. Rates are z-scored within each family; positive z-scores are accumulated across families and a final cross-species z-score gives OrthoZ. A high OrthoZ indicates orthogroup inflation.

**FamZ — copy-number signal**  
For each qualifying family, copy numbers are z-scored across species. FamZ is the cross-species z-score of the per-species count of families where that species was a copy-number outlier. A high FamZ indicates consistently abnormal gene copy numbers, possibly due to assembly issues (high proportion of inactivated and missing genes) or inflated one-to-many orthologs.

A species is flagged (`***`) if OrthoZ or FamZ exceeds the threshold.

```
=== Species QC Diagnostics ===
  ortho-z threshold: 3.0
  fam-z threshold:  3.0
  one-to-one orthologs: 11088
  total families: 17511
  scored families (OrthoZ): 288
  
  Species                        FlagFam     FamZ   OrthoZ   Flag
  ------------------------- ------------ -------- -------- ------
  HLammLeuc1                     49/5377   -0.821   -0.982
  HLcalLat1A                    199/5377   -0.190   -0.058
  HLcynGun1                     134/5377   -0.463   -0.237
  HLictTrid4A                   108/5377   -0.573   -0.679
  HLmarFlav2A                    66/5377   -0.749   -0.423
  HLmarHim1                     436/5377    0.808   -0.640
  HLmarMar1                      94/5377   -0.632    0.071
  HLmarMon3                      52/5377   -0.808   -0.480
  HLmarVanc2                    969/5377    3.051   -0.482    ***
  HLsciCar1                     179/5377   -0.274   -0.437
  HLsciNig1                     256/5377    0.050    0.627
  HLsciVul1                     167/5377   -0.324   -0.669
  HLspeDau1                     529/5377    1.199    0.012
  HLuroPar2A                    355/5377    0.467    0.128
  HLxerIna1                     115/5377   -0.543    1.034
  HLxerRut1                     197/5377   -0.198    3.216    ***

  WARNING: 1 species flagged as potential orthogroup inflators (ortho_z > 3.0):
    - HLxerRut1: ortho_z=3.22

  WARNING: 1 species flagged for abnormal copy-number estimates (fam_z > 3.0):
    - HLmarVanc2: fam_z=3.05 (969/5377 families)

  Consider re-running without these species.

Time elapsed: 2.9s
```

Use `--no-qc` to suppress this report entirely.

Use `--verbose` to print per-species orthogroups statistics

---

## Use with CAFE5

The `orthogroups_matrix.tsv` output is directly compatible with [CAFE5](https://github.com/hahnlab/CAFE5). Pass it with the `-i` flag:

```bash
cafe5 -i orthogroups/orthogroups_matrix.tsv -t species.tree -o cafe_out
```

>[!TIP]
It is recommended to run the species QC step first and exclude any flagged species before running CAFE5, as assembly artefacts can inflate apparent gene family expansions and produce spurious evolutionary rate estimates.

>[!TIP]
It is recommended to run with `--size-filter` as big gene families cause the variance of gene copy number to be too large and lead to noninformative parameter estimates.

---

## Author

**Michele Albertini** — [@michelealbertini30](https://github.com/michelealbertini30)
