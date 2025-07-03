"""GO_assessment.py  ▸  Gene Ontology enrichment helper for BCTool
-------------------------------------------------------------------
Given a mapping of *algorithm → list[Bicluster]* produced by the GUI, run
BH‑corrected GO enrichment for every bicluster and return a tidy
pandas DataFrame that Streamlit can render or download.

External dependency:  ``pip install goatools``

A *Bicluster* must expose one of these attributes:
    • ``.genes``  – iterable of gene identifiers (preferred)
    • ``.rows``   – fallback name used by some adapters

The background (universe) set defaults to *all* uploaded gene IDs.
"""
from __future__ import annotations

import pandas as pd
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Mapping, Sequence, Any

from goatools.base import (
    download_go_basic_obo,
    download_ncbi_associations,
)
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

__all__ = ["go_assessment"]

# ────────────────────────────────────────────────────────────
# Internal: download + cache GOEA object
# ────────────────────────────────────────────────────────────

@lru_cache(maxsize=2)
def _init_goea(taxid: int, universe: tuple[int, ...]) -> GOEnrichmentStudyNS:  # type: ignore[name‑defined]
    """Return a cached GOEnrichmentStudyNS for the given species/universe."""
    obo_path: Path = Path(download_go_basic_obo())
    gene2go_path: Path = Path(download_ncbi_associations(taxid))

    # Build GO DAG
    obodag = GODag(obo_path)

    # Parse gene2go → {gene: {go1, go2, …}}
    gene2go: dict[int, set[str]] = {}
    with gene2go_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            _tax, gene, go_id, *_ = line.split("\t")
            gene = int(gene)
            gene2go.setdefault(gene, set()).add(go_id)

    return GOEnrichmentStudyNS(
        set(universe),
        gene2go,
        obodag,
        propagate_counts=False,
        methods=["fdr_bh"],  # Benjamini–Hochberg
    )

# ────────────────────────────────────────────────────────────
# Public API
# ────────────────────────────────────────────────────────────

def _extract_genes(bic: Any) -> Sequence[str]:
    """Return the list of gene IDs from a Bicluster object."""
    if hasattr(bic, "genes"):
        return bic.genes  # type: ignore[attr‑defined]
    if hasattr(bic, "rows"):
        return bic.rows   # type: ignore[attr‑defined]
    raise AttributeError("Bicluster lacks .genes / .rows attribute")

def go_assessment(
    results: Mapping[str, Sequence[Any]],
    gene_ids: Iterable[Any],
    *,
    taxid: int = 9606,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Compute GO enrichment for every bicluster.

    Parameters
    ----------
    results : dict
        Mapping *algorithm → list[Bicluster]*.
    gene_ids : iterable
        All gene identifiers present in the uploaded matrix (background).
    taxid : int, default 9606
        NCBI taxonomy ID (9606 = Homo sapiens).
    alpha : float, default 0.05
        FDR threshold for significance.

    Returns
    -------
    pd.DataFrame
        Columns: [algorithm, bicluster_id, GO_ID, name, p_FDR,
                  study_hits, pop_hits]
    """
    # Build universe of integer gene IDs (required by goatools)
    universe_int = [int(g) for g in gene_ids if str(g).isdigit()]
    goea = _init_goea(taxid, tuple(universe_int))

    rows: list[dict[str, Any]] = []
    for algo, bicls in results.items():
        for idx, bic in enumerate(bicls):
            genes = [int(g) for g in _extract_genes(bic) if str(g).isdigit()]
            if not genes:
                continue
            for r in goea.run_study(genes):
                if r.p_fdr_bh <= alpha:
                    rows.append({
                        "algorithm"   : algo,
                        "bicluster_id": idx,
                        "GO_ID"       : r.GO,
                        "name"        : r.name,
                        "p_FDR"       : r.p_fdr_bh,
                        "study_hits"  : f"{r.study_count}/{r.study_n}",
                        "pop_hits"    : f"{r.pop_count}/{r.pop_n}",
                    })

    if not rows:
        return pd.DataFrame(columns=[
            "algorithm","bicluster_id","GO_ID","name",
            "p_FDR","study_hits","pop_hits",
        ])

    return (pd.DataFrame(rows)
              .sort_values(["algorithm", "bicluster_id", "p_FDR"])
              .reset_index(drop=True))
