"""Chen & Church biclustering – classic MSR version
================================================
This single file contains **your original functions unchanged** plus a minimal
pipeline and a small CLI driver so you can run everything with
```
python Chen&Church_algorithm.py my_matrix.tsv --delta 0.6 --alpha 0.05
```
If you omit the matrix path, a random 50 × 40 demo matrix is generated.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
from typing import List, Tuple

import numpy as np

# ---------------------------------------------------------------------------
#  YOUR ORIGINAL FUNCTIONS (verbatim)
# ---------------------------------------------------------------------------

def seed_generator(X: np.ndarray):
    n, m = X.shape
    selected_genes = np.random.choice(n, size=2, replace=False)
    selected_conditions = np.random.choice(m, size=2, replace=False)
    return selected_genes, selected_conditions


def purging(
    X: np.ndarray,
    selected_genes,
    selected_conditions,
    msr_threshold: float,
    a: float,
):
    deleted_cells = 0
    original_cells = len(selected_genes) * len(selected_conditions)

    while True:
        sub = X[np.ix_(selected_genes, selected_conditions)]
        selected_gene_means = sub.mean(axis=1, keepdims=True)
        selected_conditions_means = sub.mean(axis=0, keepdims=True)
        sub_mean = sub.mean()
        residues = (
            sub - selected_gene_means - selected_conditions_means + sub_mean
        )
        msr = (residues ** 2).mean()
        if msr <= msr_threshold:
            return selected_genes, selected_conditions

        genes_msr = (residues ** 2).mean(axis=1)
        conditions_msr = (residues ** 2).mean(axis=0)

        if genes_msr.max() >= conditions_msr.max():
            drop_idx = int(np.argmax(genes_msr))
            deleted_cells += len(selected_conditions)
            selected_genes = np.delete(selected_genes, drop_idx)
        else:
            drop_idx = int(np.argmax(conditions_msr))
            deleted_cells += len(selected_genes)
            selected_conditions = np.delete(selected_conditions, drop_idx)

        if (
            deleted_cells > a * original_cells
            or min(len(selected_genes), len(selected_conditions)) < 2
        ):
            return np.array([]), np.array([])  # abandon seed


# ---------------------------------------------------------------------------
#  Helper – mean‑squared residue (unchanged math)
# ---------------------------------------------------------------------------

def _msr(block: np.ndarray) -> float:
    row_means = block.mean(axis=1, keepdims=True)
    col_means = block.mean(axis=0, keepdims=True)
    overall = block.mean()
    residues = block - row_means - col_means + overall
    return (residues ** 2).mean()


# ---------------------------------------------------------------------------
#  Expansion – greedy growth (adds rows/cols while MSR ≤ δ)
# ---------------------------------------------------------------------------

def expand_bicluster(
    X: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    msr_threshold: float,
):
    rows = rows.astype(int)
    cols = cols.astype(int)

    all_rows = np.arange(X.shape[0])
    all_cols = np.arange(X.shape[1])

    remaining_rows = np.setdiff1d(all_rows, rows, assume_unique=True)
    remaining_cols = np.setdiff1d(all_cols, cols, assume_unique=True)

    improved = True
    while improved:
        improved = False
        best_gain = -np.inf
        best_is_row = True
        best_idx = None

        current_block = X[np.ix_(rows, cols)]
        msr_current = _msr(current_block)

        # test rows
        for r in remaining_rows:
            test_rows = np.append(rows, r)
            msr_new = _msr(X[np.ix_(test_rows, cols)])
            if msr_new <= msr_threshold and msr_current - msr_new > best_gain:
                best_gain, best_is_row, best_idx = msr_current - msr_new, True, r

        # test columns
        for c in remaining_cols:
            test_cols = np.append(cols, c)
            msr_new = _msr(X[np.ix_(rows, test_cols)])
            if msr_new <= msr_threshold and msr_current - msr_new > best_gain:
                best_gain, best_is_row, best_idx = msr_current - msr_new, False, c

        if best_idx is not None:
            improved = True
            if best_is_row:
                rows = np.append(rows, best_idx)
                remaining_rows = remaining_rows[remaining_rows != best_idx]
            else:
                cols = np.append(cols, best_idx)
                remaining_cols = remaining_cols[remaining_cols != best_idx]

    return rows, cols


# ---------------------------------------------------------------------------
#  Full pipeline (seed → prune → expand → mask → repeat)
# ---------------------------------------------------------------------------

def run_chen_church(
    X: np.ndarray,
    msr_threshold: float = 0.6,
    alpha: float = 0.05,
    max_biclusters: int = 100,
    random_state: int | None = None,
):
    X_work = X.copy()
    biclusters: List[Tuple[np.ndarray, np.ndarray, float]] = []

    attempts_without_success = 0
    max_attempts = 50

    while len(biclusters) < max_biclusters and attempts_without_success < max_attempts:
        # seed (avoid NaNs)
        seed_rows, seed_cols = seed_generator(X_work)
        if np.isnan(X_work[np.ix_(seed_rows, seed_cols)]).any():
            attempts_without_success += 1
            continue

        # pruning
        core_rows, core_cols = purging(X_work, seed_rows, seed_cols, msr_threshold, alpha)
        if core_rows.size == 0:
            attempts_without_success += 1
            continue

        # expansion
        full_rows, full_cols = expand_bicluster(X_work, core_rows, core_cols, msr_threshold)
        msr_val = _msr(X_work[np.ix_(full_rows, full_cols)])
        biclusters.append((full_rows, full_cols, msr_val))

        # mask found bicluster
        X_work[np.ix_(full_rows, full_cols)] = np.nan
        attempts_without_success = 0

    return biclusters


# ---------------------------------------------------------------------------
#  CLI driver – no changes to your functions above
# ---------------------------------------------------------------------------

def _load_matrix(path: pathlib.Path) -> np.ndarray:
    if path.suffix.lower() == ".csv":
        return np.loadtxt(path, delimiter=",")
    return np.loadtxt(path, delimiter="\t")


def main(argv: List[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Classic Cheng & Church biclustering (your original code + driver)")
    p.add_argument("matrix", nargs="?", type=pathlib.Path,
                   help="TSV or CSV expression matrix; if omitted, a random demo matrix is used.")
    p.add_argument("--delta", type=float, default=0.6, help="MSR threshold δ")
    p.add_argument("--alpha", type=float, default=0.05, help="Deletion fraction α")
    p.add_argument("--max", dest="max_biclusters", type=int, default=10, help="Maximum biclusters to find")
    p.add_argument("--seed", type=int, default=0, help="Random seed")
    args = p.parse_args(argv)

    if args.matrix is None:
        print("[demo] generating random 50×40 matrix…", file=sys.stderr)
        X = np.random.default_rng(args.seed).standard_normal((50, 40))
    else:
        X = _load_matrix(args.matrix)
        print(f"Loaded {args.matrix} with shape {X.shape}", file=sys.stderr)

    bics = run_chen_church(X, args.delta, args.alpha, args.max_biclusters, args.seed)
    for i, (r, c, msr_val) in enumerate(bics, 1):
        print(f"bic{i}: {len(r)}×{len(c)}  MSR={msr_val:.5f}")
    if not bics:
        print("No biclusters found within the given parameters.")


if __name__ == "__main__":
    main()

