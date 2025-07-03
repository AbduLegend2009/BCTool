"""GO_assessment.py – Gene Ontology enrichment for BCTool

This module exposes one function, ``go_assessment``.
Given a mapping ``{algo_name: list[Bicluster]}`` and the full set of
uploaded gene IDs, it returns a Benjamini–Hochberg‑corrected GO‑term
enrichment table (``pandas.DataFrame``).

Requirements
------------
$ pip install goatools pandas

A *Bicluster* must expose either ``.genes`` **or** ``.rows`` as an
iterable of gene identifiers (Entrez IDs preferred).
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Iterable, Dict, List

import pandas as pd
from goatools.base import download_go_basic_obo, download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

# ────────────────────────────────────────────────────────────
# Internal helper – build & cache GOEA object
# ────────────────────────────────────────────────────────────

@lru_cache(maxsize=4)
def _init_goea(taxid: int, gene_ids: tuple[str, ...]) -> GOEnrichmentStudyNS:  # type: ignore[name-defined]
    """Download GO files (once) and return a GOEnrichmentStudyNS."""
    obo_path = Path(download_go_basic_obo())
    gene2go_path = Path(download_ncbi_associations())

    obodag = GODag(obo_path)

    # Build gene→GO mapping
    gene2go: dict[int, set[str]] = {}
    with gene2go_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            tax, gene, go_id, *_ = line.strip().split("\t")
            if int(tax) != taxid:
                continue
            gene2go.setdefault(int(gene), set()).add(go_id)

    universe = {int(g) for g in gene_ids if str(g).isdigit()}

    return GOEnrichmentStudyNS(
        universe=universe,
        assoc=gene2go,
        godag=obodag,
        propagate_counts=False,
        methods=["fdr_bh"],
    )


# ────────────────────────────────────────────────────────────
# Public API
# ────────────────────────────────────────────────────────────

def go_assessment(results: Dict[str, List], gene_ids: Iterable[str], taxid: int = 9606, alpha: float = 0.05) -> pd.DataFrame:
    """Run GO enrichment for every bicluster and return a tidy DataFrame.

    Parameters
    ----------
    results
        Mapping of algorithm name → list of bicluster objects.
    gene_ids
        All gene identifiers present in the uploaded matrix (background set).
    taxid
        NCBI taxonomy ID (default: 9606 for human).
    alpha
        FDR threshold; only GO terms with ``p_fdr_bh <= alpha`` are kept.
    """
    # Convert to tuple so it is hashable for the cache key
    goea = _init_goea(taxid, tuple(gene_ids))

    records = []
    for algo, bicls in results.items():
        for idx, bc in enumerate(bicls):
            genes = getattr(bc, "genes", getattr(bc, "rows", []))
            genes_int = [int(g) for g in genes if str(g).isdigit()]
            if not genes_int:
                continue

            for res in goea.run_study(genes_int):
                if res.p_fdr_bh is None or res.p_fdr_bh > alpha:
                    continue

                records.append({
                    "algorithm": algo,
                    "bicluster_id": idx,
                    "GO_ID": res.GO,
                    "name": res.name,
                    "namespace": res.NS,
                    "p_FDR": res.p_fdr_bh,
                    "study_hits": f"{res.study_count}/{res.study_n}",
                    "pop_hits": f"{res.pop_count}/{res.pop_n}",
                })

    cols = ["algorithm", "bicluster_id", "GO_ID", "name", "namespace", "p_FDR", "study_hits", "pop_hits"]
    return pd.DataFrame(records, columns=cols).sort_values(["algorithm", "bicluster_id", "p_FDR"]).reset_index(drop=True)
