from __future__ import annotations
"""Unified wrappers that make every biclustering algorithm
return `List[Bicluster]` with the same signature ⇢ easier downstream code.

Usage
-----
>>> import numpy as np, pandas as pd
>>> from bicluster_adapters import ALL_ALGOS
>>> X = pd.read_csv("matrix.csv", index_col=0).values
>>> biclusters = {}
>>> for name, fn in ALL_ALGOS.items():
...     for i, bic in enumerate(fn(X), 1):
...         biclusters[f"{name}-{i}"] = bic
"""

from collections import namedtuple
from typing import List, Dict, Callable

import numpy as np

from LAS_algorithm import las_with_significance
from Chen_Church_algorithm import run_chen_church
from Bivisu_algorithm import bivisu
from ISA_algorithm import ISA_multi_seed
from OPSM_algorithm import run_opsm

# ---------------------------------------------------------------------------
# Common return type  ➜  rows, cols, score
# ---------------------------------------------------------------------------
Bicluster = namedtuple("Bicluster", ["rows", "cols", "score"])

# ---------------------------------------------------------------------------
# Helper to ensure np.ndarray[int]
# ---------------------------------------------------------------------------

def _as_idx_array(x) -> np.ndarray:
    """Convert list/set to 1‑D np.ndarray[int] (no copy if already array)."""
    if isinstance(x, np.ndarray):
        return x.astype(int, copy=False)
    return np.fromiter(x, dtype=int)

# ---------------------------------------------------------------------------
# LAS adapter
# ---------------------------------------------------------------------------

def wrap_las(X: np.ndarray, *, max_iter: int = 100, alpha: float = 0.05) -> List[Bicluster]:
    """Run LAS; always returns a **list** with ≤1 Bicluster."""
    h = las_with_significance(X, X.shape[0], X.shape[1], max_iter=int(max_iter), alpha=alpha)
    if h["rows"] is None or h["cols"] is None:
        return []
    rows = _as_idx_array(h["rows"])
    cols = _as_idx_array(h["cols"])
    score = float(h["z_score"])
    return [Bicluster(rows, cols, score)]

# ---------------------------------------------------------------------------
# Chen & Church adapter
# ---------------------------------------------------------------------------

def wrap_church(X: np.ndarray, *, delta: float = 0.6, alpha: float = 0.05, max_biclusters: int = 10) -> List[Bicluster]:
    bics = run_chen_church(X, msr_threshold=delta, alpha=alpha, max_biclusters=int(max_biclusters))
    out: List[Bicluster] = []
    for rows, cols, msr in bics:
        out.append(Bicluster(_as_idx_array(rows), _as_idx_array(cols), float(-msr)))  # negative MSR => higher better
    return out

# ---------------------------------------------------------------------------
# ISA adapter
# ---------------------------------------------------------------------------

def wrap_isa(
    X: np.ndarray,
    *,
    n_seeds: int = 500,
    seed_size: int = 5,
    t_g: float = 2.0,
    t_c: float = 2.0,
) -> List[Bicluster]:
    mods = ISA_multi_seed(X, n_seeds=int(n_seeds), seed_size=int(seed_size), t_g=t_g, t_c=t_c)
    out: List[Bicluster] = []
    for m in mods:
        rows = np.where(m["s_g"])[0]
        cols = np.where(m["s_c"])[0]
        out.append(Bicluster(rows, cols, float(m["score"])))
    return out

# ---------------------------------------------------------------------------
# OPSM adapter
# ---------------------------------------------------------------------------

def wrap_opsm(X: np.ndarray, *, k: int | None = None, restarts: int = 10) -> List[Bicluster]:
    rows, cols = run_opsm(X, k=None if k is None else int(k), restarts=int(restarts))
    return [Bicluster(_as_idx_array(rows), _as_idx_array(cols), float(rows.size))]

# ---------------------------------------------------------------------------
# BiVisu adapter
# ---------------------------------------------------------------------------

def wrap_bivisu(
    X: np.ndarray,
    *,
    model: str = "auto",
    eps: float = 0.3,
    thr: float = 0.05,
    min_rows: int = 5,
    min_cols: int = 2,
    max_iters: int = 5,
) -> List[Bicluster]:
    bics = bivisu(
        X,
        model=model,
        eps=eps,
        thr=thr,
        min_rows=int(min_rows),
        min_cols=int(min_cols),
        max_iters=int(max_iters),
    )
    out: List[Bicluster] = []
    for rows, cols in bics:
        out.append(Bicluster(_as_idx_array(list(rows)), _as_idx_array(list(cols)), float(len(rows))))
    return out

# ---------------------------------------------------------------------------
# Registry users can import in one line
# ---------------------------------------------------------------------------
ALL_ALGOS: Dict[str, Callable[[np.ndarray], List[Bicluster]]] = {
    "LAS": wrap_las,
    "Chen and Church": wrap_church,
    "ISA": wrap_isa,
    "OPSM": wrap_opsm,
    "Bivisu": wrap_bivisu,
}
