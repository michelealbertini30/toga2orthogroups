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
  -t DIR,  --toga-dir DIR          Directory with per-species TOGA2 output subdirs
  -s FILE, --species-list FILE     Newline separated list of species
  -b FILE, --transcripts-bed FILE  Reference transcript BED file
  -i FILE, --isoforms FILE         TOGA2 isoforms file
  -o DIR,  --out-dir DIR           Output directory

optional:
  -f,      --force                 Overwrite output files if they already exist
  -v,      --verbose               Print per-species processing stats
  -ul,     --include-ul            Include UL (Uncertain Loss) transcripts
           --panther FILE          PANTHER database flat file; enables PANTHER-guided family merging
           --one-to-one            Only write list of one-to-one orthologs

QC:
             --span-z FLOAT        SpanZ threshold for spanning-rate outlier detection  (default: 3.0)
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
| `orthogroups_matrix.tsv` | Copy-number count table — one row per orthogroup, one column per species |
| `orthogroups_map.tsv` | Full membership map — family ID followed by all `species\|query_gene` entries |
| `one2one.lst` | Reference gene names (one per line) for families where every species has exactly 1 ortholog (optional)|

### PANTHER mode (`--panther`)

Same files as TOGA mode; family IDs are replaced by PANTHER IDs (e.g. `PTHR12371`).


---

## Species QC

After building orthogroups, a diagnostic report is printed to stderr. Two complementary signals are computed per species:

**SpanZ — spanning-rate signal**  
For each qualifying family, the spanning rate is computed as the sum of query genes spanning multiple reference genes, divided by the total number of query genes. Rates are z-scored within each family; positive z-scores are accumulated across families and a final cross-species z-score gives SpanZ. A high SpanZ indicates orthogroup inflation.

**FamZ — copy-number signal**  
For each qualifying family, copy numbers are z-scored across species. FamZ is the cross-species z-score of the per-species count of families where that species was a copy-number outlier. A high FamZ indicates consistently abnormal gene copy numbers, possibly due to assembly issues (high proportion of inactivated and missing genes) or inflated one-to-many orthologs.

A species is flagged (`***`) if SpanZ or FamZ exceeds the threshold.

```
=== Species QC Diagnostics ===
  span-z threshold: 3.0
  fam-z threshold:  3.0
  one-to-one orthologs: 11088
  total families: 17511
  scored families (spanning): 288
  
  Species                        FlagFam     FamZ    SpanZ   Flag
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

  WARNING: 1 species flagged as potential orthogroup inflators (span_z > 3.0):
    - HLxerRut1: span_z=3.22

  WARNING: 1 species flagged for abnormal copy-number estimates (fam_z > 3.0):
    - HLmarVanc2: fam_z=3.05 (969/5377 families)

  Consider re-running without these species.

Time elapsed: 2.9s
```

Use `--no-qc` to suppress this report entirely.

Use `--verbose` to print per-species orthogroups statistics

---

## Use with CAFE5

The `orthogroups_matrix.tsv` output is directly compatible with [CAFE5](https://github.com/hahnlab/CAFE5) as the gene family count table. Pass it with the `-i` flag:

```bash
cafe5 -i orthogroups/orthogroups_matrix.tsv -t species.tree -o cafe_out
```

>[!TIP]
It is recommended to run the species QC step first and exclude any flagged species before running CAFE5, as assembly artefacts can inflate apparent gene family expansions and produce spurious evolutionary rate estimates.

---

## Author

**Michele Albertini** — [@michelealbertini30](https://github.com/michelealbertini30)
