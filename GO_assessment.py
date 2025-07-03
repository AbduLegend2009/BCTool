from __future__ import annotations

"""Gene Ontology (GO) enrichment helper for BCTool.

This module exposes a single public function, ``go_assessment``.
Given a mapping of algorithm → list[ Bicluster ], it computes
Benjamini–Hochberg‑corrected GO‑term enrichment for every bicluster
and returns a tidy ``pandas.DataFrame`` that Streamlit can render
or the user can download as CSV.

Requirements
------------
$ pip install goatools pandas

Notes
-----
*   GO data (``go-basic.obo`` and NCBI ``gene2go``) are downloaded once
    and cached in the user’s home directory. The heavy parsing work is
    memoised via ``functools.lru_cache`` so repeated GUI calls are fast.
*   ``taxid`` defaults to **9606 (Homo sapiens)**. Override if you’re
    analysing a different organism.
"""

from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
from goatools.base import download_go_basic_obo, download_ncbi_associations
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.obo_parser import GODag

__all__ = ["go_assessment"]

# ──────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────

def _ensure_int_ids(ids: Iterable) -> set[int]:
    """Return a *set* of numeric gene IDs, discarding any junk strings."""
    return {int(g) for g in ids if str(g).isdigit()}


@lru_cache(maxsize=1)
def _build_goea(taxid: int, universe: Tuple[int, ...]) -> GOEnrichmentStudyNS:  # type: ignore[name-defined]
    """Download resources (if needed) and build a memoised GOEA object."""
    obo_path = Path(download_go_basic_obo())
    gene2go_path = Path(download_ncbi_associations(taxid))

    obodag = GODag(obo_path)

    gene2go: dict[int, set[str]] = {}
    with gene2go_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue  # skip comment/header
            _tax, gene, go_id, *_ = line.split("\t")
            gene2go.setdefault(int(gene), set()).add(go_id)

    return GOEnrichmentStudyNS(
        universe,
        gene2go,
        obodag,
        propagate_counts=False,
        methods=["fdr_bh"],  # Benjamini–Hochberg FDR
    )


# ──────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────

def go_assessment(
    results: Dict[str, Iterable],
    gene_ids: Iterable,
    *,
    taxid: int = 9606,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Compute GO enrichment for every bicluster in *results*.

    Parameters
    ----------
    results
        ``{"algo_name": list[bicluster, …], …}``.  Each bicluster must
        expose either a ``genes`` attribute **or** a ``rows`` attribute
        containing an iterable of gene IDs compatible with the chosen
        organism’s NCBI gene2go mapping.
    gene_ids
        Background universe — typically ``expr_df.index`` from the GUI.
    taxid
        NCBI Taxonomy ID.  Defaults to *9606* (human).
    alpha
        Significance threshold after BH‑FDR correction (default 0.05).

    Returns
    -------
    pandas.DataFrame
        Columns: *algorithm, bicluster_id, GO_ID, term, p_FDR,
        study_hits, pop_hits*.
    """
    universe_set = _ensure_int_ids(gene_ids)
    goea = _build_goea(taxid, tuple(universe_set))

    hits: List[dict] = []

    for algo, bicls in results.items():
        for idx, bc in enumerate(bicls):
            genes = _ensure_int_ids(
                getattr(bc, "genes", getattr(bc, "rows", []))
            )
            if not genes:
                continue  # skip empty or non‑numeric sets

            for res in goea.run_study(genes):  # type: ignore[attr-defined]
                if res.p_fdr_bh <= alpha:
                    hits.append(
                        {
                            "algorithm": algo,
                            "bicluster_id": idx,
                            "GO_ID": res.GO,
                            "term": res.name,
                            "p_FDR": res.p_fdr_bh,
                            "study_hits": f"{res.study_count}/{res.study_n}",
                            "pop_hits": f"{res.pop_count}/{res.pop_n}",
                        }
                    )

    if not hits:
        return pd.DataFrame(
            columns=[
                "algorithm",
                "bicluster_id",
                "GO_ID",
                "term",
                "p_FDR",
                "study_hits",
                "pop_hits",
            ]
        )

    return (
        pd.DataFrame(hits)
        .sort_values(["algorithm", "bicluster_id", "p_FDR"])
        .reset_index(drop=True)
    )
