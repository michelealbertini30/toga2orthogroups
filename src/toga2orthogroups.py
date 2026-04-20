#!/usr/bin/env python3
"""
toga2orthogroups — Build orthogroups from TOGA2 pairwise orthology annotations.

Constructs orthogroups by modelling the (reference_gene - query_gene) orthology
relationships as a graph and extracting connected components with Union-Find.
Optionally merges with PANTHER family assignments.

Usable as:
    - CLI:  python toga2orthogroups.py -t DIR -s FILE -b FILE -i FILE -o DIR
    - Module:  from src.python.modules.toga2orthogroups import build_orthogroups, generate_count_table
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass, field
from math import inf, sqrt


__author__ = "Michele Albertini"
__email__ = "michelealbertini30@gmail.com"
__github__ = "https://github.com/michelealbertini30"


# Module-level logger
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Union-Find
# ---------------------------------------------------------------------------
# Union-Find tracks which elements belong to the same connected component.
#
# Core architecture:
#   - Every reference gene starts as its own singleton component
#   - For each query species, an inverted index is built: query gene -> ref genes
#   - If a query genes shares multiple ref genes, those get unioned in the same component
#   - This is repeated for all species, modifying the same components
#   - After all species are processed, connected components are extracted (orthogroups)

class UnionFind:
    """Weighted Union-Find with path compression"""

    __slots__ = ("_parent", "_rank")

    def __init__(self) -> None:
        # _parent[x] -> points towards the root of component x.
        self._parent: dict[str, str] = {}

        # _rank[x] -> height of the subtree rooted at x.
        self._rank: dict[str, int] = {}

    def find(self, x: str) -> str:
        """Find representative with path compression.

        If x has not been seen before, register it as a new singleton.
        Otherwise, walk up the parent chain to find the root, then
        compress the path so every node on the way points directly to the root.
        """
        # First encounter: initialise x as its own parent with rank 0.
        if x not in self._parent:
            self._parent[x] = x
            self._rank[x] = 0
            return x

        # Iterative root-finding.
        root = x
        while self._parent[root] != root:
            root = self._parent[root]

        # Path compression.
        while self._parent[x] != root:
            self._parent[x], x = root, self._parent[x]

        return root

    def union(self, a: str, b: str) -> None:
        """Merge the components containing a and b (union by rank).

        Otherwise, attach the root with the lower rank under the root with
        the higher rank. If ranks are equal, pick arbitrarily and increment
        the winner's rank.
        """
        ra, rb = self.find(a), self.find(b)

        # Already in the same component
        if ra == rb:
            return

        # Ensure ra is always the higher-rank root
        if self._rank[ra] < self._rank[rb]:
            ra, rb = rb, ra

        # Attach rb's tree under ra.
        self._parent[rb] = ra

        # Increment rank when the two trees have equal rank
        if self._rank[ra] == self._rank[rb]:
            self._rank[ra] += 1

    def components(self) -> dict[str, list[str]]:
        """Return {root: [members]} for every element currently tracked.

        Each key is the canonical root of one connected component.
        The associated list contains all members of that component,
        including the root itself.
        """
        groups: dict[str, list[str]] = defaultdict(list)
        for x in self._parent:
            # Add each member to its root's list
            # self.find(x) returns the root
            groups[self.find(x)].append(x)
        return dict(groups)


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class ReferenceGeneSet:
    """Autosomal reference genes and their transcripts.

    Built from the reference isoforms file after filtering out
    transcripts on sex chromosomes (chrX, chrY). Only autosomal genes are
    used downstream to avoid biases from copy-number differences on sex chroms.
    """
    genes: set[str]
    transcripts: set[str]
    gene_to_transcripts: dict[str, list[str]] = field(default_factory=dict)


@dataclass
class SpeciesOrthologs:
    """Per-species mapping of reference genes to query genes.

    Populated from two TOGA2 output files for each query species:
      - loss_summary.tsv -> which query transcripts have acceptable loss status.
      - orthology_classification.tsv → ref_gene + query_gene pairs
    """
    species: str
    # Set of query gene IDs with at least one intact transcript.
    ref_to_query: dict[str, set[str]]


@dataclass
class Orthogroups:
    """Connected-component orthogroups with per-species copy-number info.

    The families dict encodes the full membership needed to build count table.
    """
    # family_id -> set of "species|query_gene" strings, one per orthologous gene
    families: dict[str, set[str]]
    # Reference gene to family membership
    ref_gene_to_family: dict[str, str]
    # Full set of autosomal reference genes
    reference_genes: set[str]


# ---------------------------------------------------------------------------
# Main Functions
# ---------------------------------------------------------------------------

def load_reference_genes(
    isoforms_path: str | Path,
    transcripts_bed_path: str | Path,
) -> ReferenceGeneSet:
    """Load TOGA isoforms, filter to autosomal transcripts, return gene set.

    Parameters
    ----------
    isoforms_path : path to TOGA isoforms TSV (gene \\t transcript)
    transcripts_bed_path : path to transcript BED (chr in col1, name in col4)
    """
    # --- Collect autosomal transcript IDs from the BED file ---
    autosomal_transcripts: set[str] = set()
    sex_chroms = {"chrX", "chrY"}
    n_excluded = 0  # counter for sex-chrom transcripts dropped

    with open(transcripts_bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.split("\t", 5)
            chrom, name = parts[0], parts[3]
            if chrom in sex_chroms:
                n_excluded += 1
            else:
                autosomal_transcripts.add(name)

    log.info(
        "Transcripts BED: %d autosomal, %d chrX/Y excluded",
        len(autosomal_transcripts), n_excluded,
    )

    # --- Read isoforms TSV, keep only entries with at least one autosomal transcript---
    genes: set[str] = set()
    gene_to_tx: dict[str, list[str]] = defaultdict(list)
    n_iso_total = 0
    n_iso_kept = 0

    with open(isoforms_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            n_iso_total += 1
            gene, transcript = row[0], row[1]
            if transcript in autosomal_transcripts:
                n_iso_kept += 1
                genes.add(gene)
                gene_to_tx[gene].append(transcript)

    log.info(
        "TOGA isoforms: %d total, %d kept (autosomal), %d genes",
        n_iso_total, n_iso_kept, len(genes),
    )

    return ReferenceGeneSet(
        genes=genes,
        transcripts=autosomal_transcripts,
        gene_to_transcripts=dict(gene_to_tx),
    )


def load_species_orthologs(
    species: str,
    toga_dir: str | Path,
    reference_genes: set[str],
    include_ul: bool = False,
) -> SpeciesOrthologs:
    """Load orthology data filtered by loss status.

    Only keeps relationships where the query transcript loss status is
    FI/I/PI and UL (if -ul argument in provided).
    """
    toga_dir = Path(toga_dir)
    species_dir = toga_dir / species
    # The set of loss status codes we consider as gene presence (intact).
    intact_statuses = {"I", "FI", "PI"}
    if include_ul:
        intact_statuses.add("UL")

    # --- Parse loss_summary.tsv to find intact query transcripts ---
    intact_transcripts: set[str] = set()
    loss_path = species_dir / "loss_summary.tsv"

    with open(loss_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) >= 3 and row[2] in intact_statuses:
                intact_transcripts.add(row[1])

    log.debug(
        "  %s loss_summary: %d intact transcripts",
        species, len(intact_transcripts),
    )

    # --- Parse orthology_classification.tsv ---
    ref_to_query: dict[str, set[str]] = defaultdict(set)
    ortho_path = species_dir / "orthology_classification.tsv"
    n_rows = 0
    n_kept = 0   # not None, intact, autosomal

    with open(ortho_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 5:
                continue
            n_rows += 1
            t_gene, q_gene, q_transcript = row[0], row[2], row[3]

            if q_gene == "None":
                continue

            if q_transcript not in intact_transcripts:
                continue

            if t_gene not in reference_genes:
                continue

            n_kept += 1
            ref_to_query[t_gene].add(q_gene)

    log.debug(
        "  %s orthology: %d rows, %d kept, %d ref genes with orthologs",
        species, n_rows, n_kept, len(ref_to_query),
    )

    return SpeciesOrthologs(species=species, ref_to_query=dict(ref_to_query))


def build_orthogroups(
    ref_genes: ReferenceGeneSet,
    species_list: list[str],
    toga_dir: str | Path,
    panther_path: str | Path | None = None,
    isoforms_path: str | Path | None = None,
    include_ul: bool = False,
) -> Orthogroups:
    """Build gene families using Union-Find in two passes.

    Pass 1 — within-species orthogroups
    Pass 2 — across-species merging

    Optionally, PANTHER family assignments add additional ref->ref edges.

    Returns an Orthogroups object with family memberships and a mapping from
    reference genes to family IDs.
    """

    uf = UnionFind()

    # Pre-register reference gene so that genes with no orthologs in any
    # species still appear as singleton families in the final output.
    for g in ref_genes.genes:
        uf.find(g)

    # For each reference gene, the set of "species|query_gene"
    ref_gene_orthologs: dict[str, set[str]] = defaultdict(set)

    # --- Process each species ---
    for spe in species_list:
        log.debug("Processing species: %s", spe)

        # Load the filtered ref->query mappings for this species.
        sp_ortho = load_species_orthologs(spe, toga_dir, ref_genes.genes, include_ul=include_ul)

        # Build the inverted index: query_gene -> ref_genes
        query_to_refs: dict[str, list[str]] = defaultdict(list)

        for ref_gene, q_genes in sp_ortho.ref_to_query.items():
            for qg in q_genes:
                query_to_refs[qg].append(ref_gene)
                ref_gene_orthologs[ref_gene].add(f"{spe}|{qg}")

        # Union all reference genes that share any query gene.
        for qg, ref_list in query_to_refs.items():
            if len(ref_list) > 1:
                anchor = ref_list[0]
                for rg in ref_list[1:]:
                    uf.union(anchor, rg)

        n_merged = sum(1 for refs in query_to_refs.values() if len(refs) > 1)
        log.debug(
            "  %s: %d query genes shared by multiple ref genes",
            spe, n_merged,
        )

    # --- Optional PANTHER merge ---
    gene_to_panther: dict[str, str] = {}
    if panther_path is not None and isoforms_path is not None:
        log.info("Loading PANTHER families...")
        panther_families, gene_to_panther = _load_panther_families(panther_path, isoforms_path)
        for gene_list in panther_families.values():
            ref_in_family = [g for g in gene_list if g in ref_genes.genes]
            if len(ref_in_family) < 2:
                continue
            anchor = ref_in_family[0]
            for g in ref_in_family[1:]:
                uf.union(anchor, g)

    # --- Extract connected components and build family objects ---
    components = uf.components()

    families: dict[str, set[str]] = {}
    ref_gene_to_family: dict[str, str] = {}

    for _root, members in components.items():
        if not members:
            continue

        # Family ID: use PANTHER ID when available (PANTHER mode)
        if gene_to_panther:
            panther_ids = {gene_to_panther[m] for m in members if m in gene_to_panther}
            family_id = min(panther_ids) if panther_ids else min(members)
        else:
            family_id = min(members)

        # Aggregate all "species|query_gene" orthologs from every reference
        family_orthologs: set[str] = set()
        for rg in members:
            family_orthologs.update(ref_gene_orthologs.get(rg, set()))

        families[family_id] = family_orthologs

        # Record the family ID for every member reference gene.
        for rg in members:
            ref_gene_to_family[rg] = family_id

    log.info(
        "Orthogroups: %d families from %d reference genes",
        len(families), len(ref_gene_to_family),
    )

    return Orthogroups(
        families=families,
        ref_gene_to_family=ref_gene_to_family,
        reference_genes=ref_genes.genes,
    )


def _load_panther_families(
    panther_path: str | Path,
    isoforms_path: str | Path,
) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Load PANTHER DB and map gene symbols to reference gene IDs.

    Returns:
        family_to_genes : {panther_family_id: [ref_gene_ids]}
        gene_to_panther : {ref_gene_id: panther_family_id}

    The PANTHER flat file uses gene symbols (e.g. "BRCA2") rather than ENSG
    IDs. The isoforms TSV encodes the gene symbol in the transcript ID field
    after the last "#" character (e.g. "ENST00000380152.7#BRCA2").
    We exploit this to build a symbol -> ENSG mapping and then look up each
    reference gene in the PANTHER table.
    """
    # --- Build gene-symbol -> PANTHER family ID mapping ---
    symbol_to_family: dict[str, str] = {}

    with open(panther_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 5:
                continue
            symbol = row[2]
            # Trim to family-level ID (9 chars = "PTHRxxxxx")
            family_id = row[3][:9] if len(row[3]) >= 9 else row[3]
            symbol_to_family[symbol] = family_id

    # --- Build reference gene ID -> gene symbol mapping via isoforms ---
    gene_to_symbol: dict[str, str] = {}

    with open(isoforms_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 2:
                continue
            gene = row[0]
            symbol = row[1].rsplit("#", 1)[-1]
            if gene not in gene_to_symbol:
                gene_to_symbol[gene] = symbol

    # --- Group reference genes by their PANTHER family ---
    family_to_genes: dict[str, list[str]] = defaultdict(list)
    gene_to_panther: dict[str, str] = {}
    matched = 0  # ref genes matched to PANTHER

    for gene, symbol in gene_to_symbol.items():
        if symbol in symbol_to_family:
            pid = symbol_to_family[symbol]
            family_to_genes[pid].append(gene)
            gene_to_panther[gene] = pid
            matched += 1

    log.info(
        "PANTHER: %d symbols matched, %d families",
        matched, len(family_to_genes),
    )
    return dict(family_to_genes), gene_to_panther


def generate_count_table(
    orthogroups: Orthogroups,
    species_list: list[str],
) -> tuple[list[str], list[list]]:
    """Build the copy-number count table.

    CAFE5 expects a tab-separated file where:
      - Header: [Family ID, species1, species2, ...]
      - Rows: Family ID and copy number of that gene family in each species.
    """
    header = ["Family ID"] + species_list
    rows = []

    # Iterate over families in sorted order so the output is deterministic
    for family_id in sorted(orthogroups.families):
        members = orthogroups.families[family_id]

        # Count unique query genes per species.
        counts: dict[str, int] = {sp: 0 for sp in species_list}
        for member in members:
            if "|" in member:
                sp, _qg = member.split("|", 1)
                if sp in counts:
                    counts[sp] += 1

        rows.append([family_id] + [counts[sp] for sp in species_list])

    return header, rows


def write_count_table(
    header: list[str],
    rows: list[list],
    output_path: str | Path,
) -> None:
    """Write copy-number count table to disk."""

    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    log.info("Wrote count table: %s (%d families)", output_path, len(rows))


def write_orthogroup_membership(
    orthogroups: Orthogroups,
    output_path: str | Path,
) -> None:
    """Write orthogroup membership file.

    The first field is the family ID, remaining fields are the "species|query_gene"
    """
    with open(output_path, "w") as fh:
        for family_id in sorted(orthogroups.families):
            members = sorted(orthogroups.families[family_id])
            fh.write(family_id + "\t" + "\t".join(members) + "\n")
    log.info("Wrote membership: %s (%d families)", output_path, len(orthogroups.families))


# ---------------------------------------------------------------------------
# Species QC — inflation diagnostics
# ---------------------------------------------------------------------------

def species_qc_diagnostics(
    header: list[str],
    rows: list[list],
    species_list: list[str],
    z_threshold: float = 3.0,
    outlier_family_fraction: float = 0.05,
    output=None,
) -> dict[str, dict]:
    """Detect species with inflated orthogroup sizes.

    For each gene family, compute per-species z-scores of copy number
    relative to the cross-species distribution for that family.  A species
    is flagged if its proportion of z-score outlier families exceeds
    outlier_family_fraction.

    z_threshold : families where higher z are counted as outliers
    outlier_family_fraction : species are flagged if their outlier
        proportion exceeds this

    Returns
    -------
    dict mapping species -> {outlier_count, outlier_fraction, flagged}
    """
    # Default output target is stderr so QC text doesn't pollute stdout
    if output is None:
        output = sys.stderr

    # Map each species name to its column index in the count table.
    sp_indices = {sp: i for i, sp in enumerate(species_list)}
    n_sp = len(species_list)

    count_matrix: list[list[int]] = [
        [row[sp_indices[sp] + 1] for sp in species_list] for row in rows
    ]
    n_families = len(count_matrix)

    # Compute per-family z-scores across species.
    z_scores: list[list[float]] = []
    for fam_row in count_matrix:
        mean_val = sum(fam_row) / n_sp
        # Sample variance (divide by n-1) to be conservative.
        var_val = sum((x - mean_val) ** 2 for x in fam_row) / (n_sp - 1) if n_sp > 1 else 0.0
        std_val = sqrt(var_val) if var_val > 0 else inf
        z_scores.append([(x - mean_val) / std_val for x in fam_row])

    # For each species, count how many families it is an outlier in.
    results: dict[str, dict] = {}
    flagged_species: list[str] = []

    for sp in species_list:
        col = sp_indices[sp]
        outlier_count = sum(1 for zrow in z_scores if abs(zrow[col]) > z_threshold)
        outlier_frac = outlier_count / n_families if n_families > 0 else 0.0
        flagged = outlier_frac > outlier_family_fraction

        results[sp] = {
            "outlier_count": outlier_count,
            "outlier_fraction": round(outlier_frac, 4),
            "flagged": flagged,
        }
        if flagged:
            flagged_species.append(sp)

    # One-to-one orthologs: families where no species has more than 1 copy.
    n_one2one = sum(1 for fam_row in count_matrix if max(fam_row) <= 1 and sum(fam_row) > 0)

    # --- Print formatted report ---
    output.write("\n=== Species QC Diagnostics ===\n")
    output.write(f"  z-score threshold: {z_threshold}\n")
    output.write(f"  outlier family fraction threshold: {outlier_family_fraction}\n")
    output.write(f"  one2one orthologs: {n_one2one}\n")
    output.write(f"  total families: {n_families}\n\n")
    output.write(f"  {'Species':<25} {'Outliers':>8} {'Fraction':>10} {'Flag':>6}\n")
    output.write(f"  {'-'*25} {'-'*8} {'-'*10} {'-'*6}\n")

    for sp in species_list:
        r = results[sp]
        flag_str = " ***" if r["flagged"] else ""
        output.write(
            f"  {sp:<25} {r['outlier_count']:>8} {r['outlier_fraction']:>10.4f}"
            f" {flag_str:>6}\n"
        )

    if flagged_species:
        output.write(f"\n  WARNING: {len(flagged_species)} species flagged with high orthogroup variance:\n")
        for sp in flagged_species:
            output.write(
                f"    - {sp}: {results[sp]['outlier_fraction']:.1%} of families are outliers. "
                f"Consider re-running without this species.\n"
            )
    else:
        output.write("\n  No species flagged.\n")

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    app = argparse.ArgumentParser(
        description=(
            "Build orthogroups from TOGA2 pairwise orthology annotations.\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s [-h] -t DIR -s FILE -b FILE -i FILE -o DIR [options]",
        epilog=(
            "example:\n"
            "  toga2orthogroups.py \\\n"
            "    -t TOGA2 \\\n"
            "    -s species.lst \\\n"
            "    -b toga2.transcripts.bed \\\n"
            "    -i toga2.isoforms.tsv \\\n"
            "    -o orthogroups \\\n"
        ),
    )

    req = app.add_argument_group("required")
    req.add_argument(
        "-t", "--toga-dir",
        required=True,
        metavar="DIR",
        help="directory with per-species TOGA2 output subdirs",
    )
    req.add_argument(
        "-s", "--species-list",
        required=True,
        metavar="FILE",
        help="newline separated list of species",
    )
    req.add_argument(
        "-b", "--transcripts-bed",
        required=True,
        metavar="FILE",
        help="reference transcript BED file",
    )
    req.add_argument(
        "-i", "--isoforms",
        required=True,
        metavar="FILE",
        help="TOGA2 isoforms file",
    )
    req.add_argument(
        "-o", "--out-dir",
        required=True,
        metavar="DIR",
        help="output directory",
    )

    opt = app.add_argument_group("optional")
    opt.add_argument(
        "-f", "--force",
        action="store_true",
        help="overwrite output files if they already exist",
    )
    opt.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="print per-species processing stats",
    )
    opt.add_argument(
        "-ul", "--include-ul",
        action="store_true",
        help="include UL (Uncertain Loss) transcripts",
    )
    opt.add_argument(
        "--panther",
        metavar="FILE",
        default=None,
        help="panther database",
    )

    qc = app.add_argument_group("QC")
    qc.add_argument(
        "-z",
        "--z-threshold",
        type=float,
        metavar="FLOAT",
        default=3.0,
        help="z-score threshold for outlier detection  (default: 3.0)",
    )
    qc.add_argument(
        "-of",
        "--outlier-fraction",
        type=float,
        metavar="FLOAT",
        default=0.05,
        help="max outlier family fraction for flagging  (default: 0.05)",
    )
    qc.add_argument(
        "--no-qc",
        action="store_true",
        help="skip species QC diagnostics",
    )

    return app.parse_args(argv)


def run(
    toga_dir: str,
    species_list_path: str,
    transcripts_bed: str,
    isoforms: str,
    out_dir: str,
    force: bool = False,
    verbose: bool = False,
    include_ul: bool = False,
    panther: str | None = None,
    z_threshold: float = 3.0,
    outlier_fraction: float = 0.05,
    no_qc: bool = False,
) -> None:
    """Core runner — called by both the argparse main() and the Click subcommand.

    All filesystem I/O, orthogroup building, and QC happen here.
    Logging is configured by the caller before invoking this function.
    """
    # --- Resolve output directory and expected file paths ---
    _out_dir = Path(out_dir)
    prefix = "PANTHER" if panther else "TOGA2"
    out_tsv = _out_dir / f"{prefix}.orthogroups.tsv"
    out_map = _out_dir / f"{prefix}.ortho_map.tsv"

    # Create output directory if it doesn't exist.
    _out_dir.mkdir(parents=True, exist_ok=True)

    # Check for existing output files unless --force is set.
    if not force:
        existing = [f for f in (out_tsv, out_map) if f.exists()]
        if existing:
            log.error(
                "Output file(s) already exist: %s\n"
                "Use -f / --force to overwrite.",
                ", ".join(str(f) for f in existing),
            )
            sys.exit(1)

    # --- Load species list ---
    # One species name per line; tab-separated lines are also accepted
    species_list: list[str] = []
    with open(species_list_path) as fh:
        for line in fh:
            sp = line.strip().split("\t")[0]
            if sp:
                species_list.append(sp)
    log.info("Species: %d loaded", len(species_list))

    # --- Load reference gene set ---
    ref_genes = load_reference_genes(isoforms, transcripts_bed)

    # --- Build orthogroups ---
    if panther:
        log.info("Running in PANTHER mode...")
        orthogroups = build_orthogroups(
            ref_genes, species_list, toga_dir,
            panther_path=panther,
            isoforms_path=isoforms,
            include_ul=include_ul,
        )
    else:
        orthogroups = build_orthogroups(
            ref_genes, species_list, toga_dir,
            panther_path=None, isoforms_path=None,
            include_ul=include_ul,
        )

    header, rows = generate_count_table(orthogroups, species_list)
    write_count_table(header, rows, out_tsv)
    write_orthogroup_membership(orthogroups, out_map)

    if not no_qc:
        species_qc_diagnostics(
            header, rows, species_list,
            z_threshold=z_threshold,
            outlier_family_fraction=outlier_fraction,
        )


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)

    # Configure logging.
    # With --verbose: DEBUG level with per-species stats.
    logging.basicConfig(format="%(message)s", stream=sys.stderr)
    logging.getLogger().setLevel(logging.DEBUG if args.verbose else logging.INFO)

    run(
        toga_dir=args.toga_dir,
        species_list_path=args.species_list,
        transcripts_bed=args.transcripts_bed,
        isoforms=args.isoforms,
        out_dir=args.out_dir,
        force=args.force,
        verbose=args.verbose,
        include_ul=args.include_ul,
        panther=args.panther,
        z_threshold=args.z_threshold,
        outlier_fraction=args.outlier_fraction,
        no_qc=args.no_qc,
    )


if __name__ == "__main__":
    main()
