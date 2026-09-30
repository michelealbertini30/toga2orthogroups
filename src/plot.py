#!/usr/bin/env python3
"""
plot — draw pairwise orthology relationships for a single gene family.

Renders the (reference_gene - query_gene) relationships of one gene family in
up to three views: bipartite (two columns), projection (force-directed, query
loci drawn as diamond nodes) and collapsed (reference genes only). Pairwise
only: one reference annotation against one query species.

The orthology_classification.tsv file is taken as ground truth. Every kept row
becomes one edge; nothing is inferred, grouped, highlighted or colour-coded.

Usable as:
    - CLI:     run_toga2_orthogroups.py plot -t DIR -q SPECIES -g GENE -o OUTDIR
    - Module:  from src.plot import load_pairwise_edges, family_edges, render
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from pathlib import Path
from collections import Counter, defaultdict
from dataclasses import dataclass, field


log = logging.getLogger(__name__)

# Same guard as the main module: TOGA2 orthology files contain very long fields.
_csv_limit = sys.maxsize
while True:
    try:
        csv.field_size_limit(_csv_limit)
        break
    except OverflowError:
        _csv_limit //= 2


_DEPS_MISSING = """\
Plotting requires two extra packages:

    pip install networkx matplotlib

Everything else in toga2orthogroups runs on the standard library alone.
"""


def _load_deps():
    """Import networkx/matplotlib on demand so the core stays stdlib-only."""
    try:
        import networkx as nx
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patheffects as pe
    except ImportError:
        sys.stderr.write(_DEPS_MISSING)
        raise SystemExit(1)
    return nx, plt, pe


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
# Deliberately monochrome. The layout is the only thing that says anything
# about structure; the plotter asserts nothing on top of the data.

SURFACE     = "#fcfcfb"
REF_FILL    = "#ffffff"
REF_EDGE    = "#3d3c39"
QUERY_FILL  = "#6f6e69"
LINK        = "#c9c8c0"
INK         = "#0b0b0b"
INK_QUERY   = "#52514e"


# The three figures the module can draw, in the order they are written.
LAYOUTS = ("bipartite", "projection", "collapsed")


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class PairwiseOrthology:
    """Reference-to-query orthology edges for one query species."""

    species: str
    # (ref_gene, query_gene, orthology_class) — one per distinct relationship
    edges: list[tuple[str, str, str]] = field(default_factory=list)
    # ref_gene -> gene symbol, parsed from the reference transcript name
    symbols: dict[str, str] = field(default_factory=dict)
    # every reference gene in the file, including rows the filters dropped
    genes: set[str] = field(default_factory=set)
    # relationships rejected *only* because their projection is UL; empty when
    # include_ul is set. Used to explain what the default filter left out.
    ul_edges: list[tuple[str, str, str]] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_pairwise_edges(
    species: str,
    toga_dir: str | Path,
    reference_genes: set[str] | None = None,
    include_ul: bool = False,
) -> PairwiseOrthology:
    """Read one species' orthology relationships.

    Applies the same loss-status filter as the orthogroup builder: a projection
    counts only when its query transcript is FI/I/PI, plus UL if include_ul.

    Parameters
    ----------
    species : query species name (a subdirectory of toga_dir)
    toga_dir : directory holding one subdirectory per query species
    reference_genes : restrict to these reference genes (default: no restriction)
    include_ul : also accept UL (Uncertain Loss) projections
    """
    species_dir = Path(toga_dir) / species
    intact_statuses = {"I", "FI", "PI"}
    if include_ul:
        intact_statuses.add("UL")

    # --- Query transcripts passing the loss filter ---
    # `ul` collects what the default filter rejects for being UL and nothing
    # else, so the caller can say what a plain run leaves out. With include_ul
    # those transcripts land in `intact` and `ul` stays empty.
    intact: set[str] = set()
    ul: set[str] = set()
    with open(species_dir / "loss_summary.tsv") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 3:
                continue
            if row[2] in intact_statuses:
                intact.add(row[1])
            elif row[2] == "UL":
                ul.add(row[1])

    # --- Orthology relationships ---
    seen: set[tuple[str, str]] = set()
    edges: list[tuple[str, str, str]] = []
    symbols: dict[str, str] = {}
    genes: set[str] = set()
    ul_seen: set[tuple[str, str]] = set()
    ul_edges: list[tuple[str, str, str]] = []
    n_rows = 0

    with open(species_dir / "orthology_classification.tsv") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 5:
                continue
            n_rows += 1
            t_gene, t_transcript, q_gene, q_transcript, ortho_class = row[:5]

            # Every reference gene in the file is a legal seed, and symbols
            # come from the reference transcript name (ENST...#SYMBOL), which
            # is not always present. Record both even for rows the filters drop.
            genes.add(t_gene)
            if "#" in t_transcript:
                symbols.setdefault(t_gene, t_transcript.split("#")[1])

            if q_gene == "None":
                continue
            if reference_genes is not None and t_gene not in reference_genes:
                continue

            # Several transcripts of the same gene pair produce the same edge.
            if q_transcript in intact:
                if (t_gene, q_gene) in seen:
                    continue
                seen.add((t_gene, q_gene))
                edges.append((t_gene, q_gene, ortho_class))
            elif q_transcript in ul:
                if (t_gene, q_gene) in ul_seen:
                    continue
                ul_seen.add((t_gene, q_gene))
                ul_edges.append((t_gene, q_gene, ortho_class))

    # A pair reachable through an accepted projection is not "UL-only", even
    # if some other transcript of the same pair was UL.
    ul_edges = [e for e in ul_edges if (e[0], e[1]) not in seen]

    log.debug(
        "  %s: %d rows, %d distinct gene-level relationships (%d UL-only)",
        species, n_rows, len(edges), len(ul_edges),
    )
    return PairwiseOrthology(
        species=species, edges=edges, symbols=symbols, genes=genes,
        ul_edges=ul_edges,
    )


# ---------------------------------------------------------------------------
# Family selection
# ---------------------------------------------------------------------------

def resolve_seeds(
    seeds: list[str],
    symbols: dict[str, str],
    known_genes: set[str] | None = None,
) -> set[str]:
    """Map user-supplied gene IDs or gene symbols onto reference gene IDs.

    A reference gene ID (ENSG00000002586) matches directly, a gene symbol
    case-insensitively. IDs are checked against every gene in the annotation,
    not only those whose transcript name carried a symbol.
    """
    known = set(symbols) | (known_genes or set())
    by_symbol: dict[str, list[str]] = defaultdict(list)
    for gene, sym in symbols.items():
        by_symbol[sym.upper()].append(gene)

    resolved: set[str] = set()
    for seed in seeds:
        if seed in known:
            resolved.add(seed)
        elif seed.upper() in by_symbol:
            resolved.update(by_symbol[seed.upper()])
        else:
            log.warning("Seed not found in reference annotation: %s", seed)
    return resolved


def _name(gene: str, symbols: dict[str, str]) -> str:
    """Gene symbol when the annotation carries one, otherwise the ID."""
    sym = symbols.get(gene)
    return f"{sym} ({gene})" if sym else gene


def family_edges(
    ortho: PairwiseOrthology,
    seeds: set[str],
) -> list[tuple[str, str, str]]:
    """Return the edges of the connected component(s) containing the seeds.

    Reference and query genes form one bipartite graph; the family is whatever
    is reachable from a seed through shared query genes.
    """
    ref_to_q: dict[str, set[str]] = defaultdict(set)
    q_to_ref: dict[str, set[str]] = defaultdict(set)
    for ref, qry, _cls in ortho.edges:
        ref_to_q[ref].add(qry)
        q_to_ref[qry].add(ref)

    keep_ref: set[str] = set()
    stack = [s for s in seeds if s in ref_to_q]
    for s in seeds:
        if s not in ref_to_q:
            log.warning(
                "Seed has no orthologs in %s under the current filter: %s",
                ortho.species, s,
            )
    keep_ref.update(stack)
    keep_q: set[str] = set()

    while stack:
        ref = stack.pop()
        for qry in ref_to_q[ref]:
            if qry in keep_q:
                continue
            keep_q.add(qry)
            for other in q_to_ref[qry]:
                if other not in keep_ref:
                    keep_ref.add(other)
                    stack.append(other)

    return [e for e in ortho.edges if e[0] in keep_ref]


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def _collapsed_edges(
    edges: list[tuple[str, str, str]],
) -> list[tuple[str, str, set[str]]]:
    """Collapse the bipartite graph onto reference genes only.

    Two reference genes are joined when they share at least one query locus.
    One shared locus with n reference genes becomes an n-clique, so a single
    shared locus stops being distinguishable from several separate ones —
    that information survives only in the bipartite and projection views.
    """
    q_to_refs: dict[str, set[str]] = defaultdict(set)
    for ref, qry, _cls in edges:
        q_to_refs[qry].add(ref)

    pairs: dict[tuple[str, str], set[str]] = defaultdict(set)
    for qry, refs in q_to_refs.items():
        members = sorted(refs)
        for i, a in enumerate(members):
            for b in members[i + 1:]:
                pairs[(a, b)].add(qry)
    return [(a, b, q) for (a, b), q in pairs.items()]


def _bipartite_positions(graph, ref_nodes, query_nodes, sweeps: int = 12):
    """Two columns, reference left and query right, with crossing reduction.

    Orders each column by the barycentre of its neighbours in the other
    column, sweeping back and forth — the standard crossing-reduction
    heuristic for two-layer graphs.
    """
    left, right = list(ref_nodes), list(query_nodes)

    def sweep(movable, fixed):
        rank = {n: i for i, n in enumerate(fixed)}
        def bary(n):
            nb = [rank[m] for m in graph.neighbors(n) if m in rank]
            return sum(nb) / len(nb) if nb else 0.0
        movable.sort(key=bary)

    for _ in range(sweeps):
        sweep(right, left)
        sweep(left, right)

    pos = {}
    for column, x in ((left, -1.0), (right, 1.0)):
        span = max(len(column) - 1, 1)
        for i, node in enumerate(column):
            pos[node] = (x, 1.0 - 2.0 * i / span)
    return pos


def render(
    edges: list[tuple[str, str, str]],
    symbols: dict[str, str],
    output_path: str | Path,
    layout: str = "projection",
    labels: str = "all",
    hide_pendants: bool = False,
    layout_seed: int = 1,
    spread: float = 0.55,
    iterations: int = 800,
    figsize: tuple[float, float] = (12.0, 9.0),
    dpi: int = 150,
) -> None:
    """Draw one view of the orthology graph and write it to output_path.

    Reference genes are open circles, query genes small filled diamonds, and
    every kept relationship is one edge. A diamond carrying several edges is
    one query locus that several reference genes map onto. Query nodes are
    labelled with their TOGA ID verbatim.

    layout is one of LAYOUTS: "bipartite" puts the two species in their own
    column, "projection" lays the same graph out force-directed, "collapsed"
    drops the query nodes and joins reference genes sharing a locus. Output
    format follows the file suffix (.svg, .pdf, .png).
    """
    nx, plt, pe = _load_deps()

    graph = nx.Graph()
    if layout == "collapsed":
        # Reference genes only: query loci become edges rather than nodes.
        for ref, _qry, _cls in edges:
            graph.add_node(("R", ref))
        for a, b, _shared in _collapsed_edges(edges):
            graph.add_edge(("R", a), ("R", b))
    else:
        for ref, qry, _cls in edges:
            graph.add_edge(("R", ref), ("Q", qry))

    if hide_pendants:
        # A query gene joined to a single reference gene adds no connectivity.
        pendants = [n for n in graph if n[0] == "Q" and graph.degree(n) == 1]
        graph.remove_nodes_from(pendants)
        log.info("Hid %d pendant query genes", len(pendants))

    ref_nodes = [n for n in graph if n[0] == "R"]
    query_nodes = [n for n in graph if n[0] == "Q"]
    if not ref_nodes:
        log.error("Nothing to plot: no reference genes in the selected family.")
        sys.exit(1)

    # Shrink marks and text as the family grows so large graphs stay legible.
    n = len(ref_nodes)
    fit = 1.0 if n <= 20 else max(0.28, (20.0 / n) ** 0.5)
    ref_size = max(70.0, 330.0 * fit)
    query_size = max(16.0, 60.0 * fit)
    ref_font = max(4.5, 9.0 * max(0.55, fit))
    query_font = ref_font * 0.76

    if layout == "bipartite":
        pos = _bipartite_positions(graph, ref_nodes, query_nodes)
    else:
        pos = nx.spring_layout(
            graph, k=spread, iterations=iterations, seed=layout_seed,
        )

    fig, ax = plt.subplots(figsize=figsize, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    nx.draw_networkx_edges(graph, pos, edge_color=LINK, width=1.0, ax=ax)
    nx.draw_networkx_nodes(
        graph, pos, nodelist=query_nodes, node_color=QUERY_FILL,
        node_shape="d", node_size=query_size, linewidths=0, ax=ax,
    )
    nx.draw_networkx_nodes(
        graph, pos, nodelist=ref_nodes, node_color=REF_FILL, node_shape="o",
        node_size=ref_size, edgecolors=REF_EDGE, linewidths=1.3, ax=ax,
    )

    halo = lambda w: [pe.withStroke(linewidth=w, foreground=SURFACE)]
    offset = 0.055 * max(0.5, fit)
    # In the two-column layout each species owns a side, so labels sit
    # outward from their own column instead of above the node.
    side = layout == "bipartite"

    if labels in ("ref", "all"):
        for node in ref_nodes:
            x, y = pos[node]
            text = symbols.get(node[1], node[1])
            if side:
                ax.text(x - 0.06, y, text, fontsize=ref_font, color=INK,
                        ha="right", va="center", zorder=7,
                        path_effects=halo(2.6))
            else:
                ax.text(x, y + offset, text, fontsize=ref_font, color=INK,
                        ha="center", va="bottom", zorder=7,
                        path_effects=halo(2.6))
    if labels == "all":
        # Query nodes keep their TOGA ID verbatim.
        for node in query_nodes:
            x, y = pos[node]
            text = node[1]
            if side:
                ax.text(x + 0.06, y, text, fontsize=query_font,
                        color=INK_QUERY, ha="left", va="center", zorder=6,
                        style="italic", path_effects=halo(2.4))
            else:
                ax.text(x, y - offset * 0.8, text, fontsize=query_font,
                        color=INK_QUERY, ha="center", va="top", zorder=6,
                        style="italic", path_effects=halo(2.4))

    # The two-column layout puts long TOGA IDs outside the data range.
    ax.margins(x=0.55, y=0.05) if side else ax.margins(0.12)
    ax.axis("off")
    fig.tight_layout(pad=0.6)
    fig.savefig(output_path, dpi=dpi, facecolor=SURFACE)
    plt.close(fig)

    if layout == "collapsed":
        log.info(
            "Wrote %s (%d reference genes, %d pairs sharing a query locus)",
            output_path, len(ref_nodes), graph.number_of_edges(),
        )
    else:
        log.info(
            "Wrote %s (%d reference genes, %d query genes, %d relationships)",
            output_path, len(ref_nodes), len(query_nodes),
            graph.number_of_edges(),
        )


# ---------------------------------------------------------------------------
# Output naming
# ---------------------------------------------------------------------------

_UNSAFE = re.compile(r"[^A-Za-z0-9._-]+")


def parse_layouts(spec: str) -> list[str]:
    """Expand a --layout value into an ordered list of layout names."""
    if spec.strip().lower() == "all":
        return list(LAYOUTS)

    chosen: list[str] = []
    for item in spec.split(","):
        item = item.strip().lower()
        if not item:
            continue
        if item not in LAYOUTS:
            log.error("Unknown layout %r; choose from %s, or 'all'",
                      item, ", ".join(LAYOUTS))
            sys.exit(1)
        if item not in chosen:
            chosen.append(item)

    if not chosen:
        log.error("--layout selected nothing")
        sys.exit(1)
    # Keep the canonical order whatever order they were given in.
    return [lay for lay in LAYOUTS if lay in chosen]


def output_stem(gene: str) -> str:
    """Filename stem from the -g value, e.g. ENSG00000002586 or LILRB1.

    Several comma-separated seeds are joined with underscores; anything that
    does not belong in a filename is replaced.
    """
    parts = [_UNSAFE.sub("_", s.strip()) for s in gene.split(",") if s.strip()]
    return "_".join(parts) or "family"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    app = argparse.ArgumentParser(
        prog="toga2orthogroups.py plot",
        description="Plot pairwise orthology relationships for one gene family.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
example:
  run_toga2_orthogroups.py plot \\
    -t TOGA2 \\
    -q mm39 \\
    -g ENSG00000104974 \\
    -o figures/

  writes figures/ENSG00000104974_bipartite.png
         figures/ENSG00000104974_projection.png
         figures/ENSG00000104974_collapsed.png

notes:
  -g accepts a reference gene ID (ENSG00000104974), a gene symbol, or a family
  ID from orthogroups_map.tsv (family IDs are reference gene IDs). Several may
  be given comma-separated; the plot covers every family they fall into, and
  the file names join them with underscores.
""",
    )

    req = app.add_argument_group("required")
    req.add_argument(
        "-t", "--toga-dir", metavar="DIR", required=True,
        help="Directory with per-species TOGA2 output subdirs",
    )
    req.add_argument(
        "-q", "--query-species", metavar="NAME", required=True,
        help="Query species to plot against (a subdirectory of --toga-dir)",
    )
    req.add_argument(
        "-g", "--gene", metavar="LIST", required=True,
        help="Seed gene(s): symbol, reference gene ID, or family ID (comma-separated)",
    )
    req.add_argument(
        "-o", "--out-dir", metavar="DIR", required=True,
        help="Output directory; files are named <gene>_<layout>.<format>",
    )

    opt = app.add_argument_group("optional")
    opt.add_argument(
        "-i", "--isoforms", metavar="FILE",
        help="TOGA2 isoforms file; with -b, restricts to the same reference "
             "gene set the orthogroup builder uses",
    )
    opt.add_argument(
        "-b", "--transcripts-bed", metavar="FILE",
        help="Reference transcript BED file (needed for -bl)",
    )
    opt.add_argument(
        "-bl", "--blacklist", metavar="LIST|FILE",
        help="Chromosomes/scaffolds to exclude: comma-separated list or file",
    )
    opt.add_argument(
        "-ul", "--include-ul", action="store_true",
        help="Include UL (Uncertain Loss) transcripts",
    )
    opt.add_argument(
        "--layout", metavar="LIST", default="all",
        help="Comma-separated list of figures to plot: bipartite, projection, collapsed (default: all)",
    )
    opt.add_argument(
        "--format", choices=("png", "svg", "pdf"), default="png",
        help="Output file format (default: png)",
    )
    opt.add_argument(
        "--labels", choices=("all", "ref", "none"), default="all",
        help="Which nodes to label: all, ref, none (default: all)",
    )
    opt.add_argument(
        "--hide-pendants", action="store_true",
        help="Hide query genes joined to only one reference gene",
    )
    opt.add_argument(
        "--layout-seed", type=int, default=1, metavar="INT",
        help="Layout random seed; same seed gives the same figure (default: 1)",
    )
    opt.add_argument(
        "--spread", type=float, default=0.55, metavar="FLOAT",
        help="Force-directed repulsion constant (default: 0.55)",
    )
    opt.add_argument(
        "--figsize", default="12x9", metavar="WxH",
        help="Figure size in inches (default: 12x9)",
    )
    opt.add_argument(
        "--dpi", type=int, default=150, metavar="INT",
        help="Raster resolution (default: 150)",
    )
    opt.add_argument(
        "-v", "--verbose", action="store_true",
        help="Print per-species processing stats",
    )

    return app.parse_args(argv)


def run(
    toga_dir: str,
    query_species: str,
    gene: str,
    out_dir: str,
    isoforms: str | None = None,
    transcripts_bed: str | None = None,
    blacklist: str | None = None,
    include_ul: bool = False,
    layout: str = "all",
    fmt: str = "png",
    labels: str = "all",
    hide_pendants: bool = False,
    layout_seed: int = 1,
    spread: float = 0.55,
    figsize: str = "12x9",
    dpi: int = 150,
) -> None:
    """Core runner — load, select the family, render every requested view."""
    # Fail before doing any work: missing packages, bad flags, unwritable dir.
    _load_deps()
    layouts = parse_layouts(layout)
    try:
        w, h = (float(v) for v in figsize.lower().split("x"))
    except ValueError:
        log.error("--figsize must look like WxH, e.g. 12x9")
        sys.exit(1)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    log.info("\n=== Plotting orthology relationships ===")

    # --- Optional reference gene filtering, matching the orthogroup builder ---
    ref_genes: set[str] | None = None
    if isoforms and transcripts_bed:
        try:
            from src.toga2orthogroups import load_reference_genes, _parse_blacklist
        except ImportError:
            from toga2orthogroups import load_reference_genes, _parse_blacklist
        ref_genes = load_reference_genes(
            isoforms, transcripts_bed, blacklist=_parse_blacklist(blacklist),
        ).genes
    elif blacklist:
        log.warning("--blacklist needs both -i and -b; ignoring it")

    ortho = load_pairwise_edges(
        query_species, toga_dir, reference_genes=ref_genes, include_ul=include_ul,
    )
    log.info(
        "%s: %d relationships over %d reference genes",
        query_species, len(ortho.edges), len({e[0] for e in ortho.edges}),
    )

    seeds = resolve_seeds([s.strip() for s in gene.split(",") if s.strip()],
                          ortho.symbols, ortho.genes)
    if not seeds:
        log.error("None of the requested genes were found.")
        sys.exit(1)

    edges = family_edges(ortho, seeds)
    if not edges:
        # Distinguish "nothing in the file" from "everything filtered out".
        ul_only = sorted(s for s in seeds if any(e[0] == s for e in ortho.ul_edges))
        if ul_only:
            log.error(
                "No orthologs in %s under the default loss filter (FI/I/PI).",
                query_species,
            )
            log.error(
                "  Every projection of %s is UL (Uncertain Loss). Add -ul to keep them.",
                ", ".join(_name(g, ortho.symbols) for g in ul_only),
            )
        else:
            log.error("The requested gene(s) have no orthologs in %s.", query_species)
        sys.exit(1)

    members = sorted({e[0] for e in edges}, key=lambda g: ortho.symbols.get(g, g))
    classes = Counter(e[2] for e in edges)
    log.info(
        "Family: %d reference genes, %d query genes, %d relationships",
        len(members), len({e[1] for e in edges}), len(edges),
    )
    log.info("  %s", ", ".join(f"{k}={v}" for k, v in sorted(classes.items())))
    log.info("  %s", ", ".join(_name(g, ortho.symbols) for g in members))

    # The loss filter can hold genes out of a family that is otherwise intact,
    # which looks like a missing gene rather than a filtered one. Say so.
    if ortho.ul_edges:
        wider = family_edges(
            PairwiseOrthology(
                species=ortho.species,
                edges=ortho.edges + ortho.ul_edges,
                symbols=ortho.symbols,
                genes=ortho.genes,
            ),
            seeds,
        )
        extra = sorted({e[0] for e in wider} - set(members),
                       key=lambda g: _name(g, ortho.symbols))
        if extra:
            log.info(
                "  note: %d more reference gene(s) join this family with -ul: %s",
                len(extra), ", ".join(_name(g, ortho.symbols) for g in extra),
            )

    stem = output_stem(gene)
    for lay in layouts:
        render(
            edges, ortho.symbols, out_path / f"{stem}_{lay}.{fmt}",
            layout=lay,
            labels=labels, hide_pendants=hide_pendants, layout_seed=layout_seed,
            spread=spread, figsize=(w, h), dpi=dpi,
        )


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)

    logging.basicConfig(format="%(message)s", stream=sys.stderr)
    logging.getLogger().setLevel(logging.DEBUG if args.verbose else logging.INFO)

    run(
        toga_dir=args.toga_dir,
        query_species=args.query_species,
        gene=args.gene,
        out_dir=args.out_dir,
        isoforms=args.isoforms,
        transcripts_bed=args.transcripts_bed,
        blacklist=args.blacklist,
        include_ul=args.include_ul,
        layout=args.layout,
        fmt=args.format,
        labels=args.labels,
        hide_pendants=args.hide_pendants,
        layout_seed=args.layout_seed,
        spread=args.spread,
        figsize=args.figsize,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
