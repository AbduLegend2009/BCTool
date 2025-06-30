#!/usr/bin/env python3
"""
Greedy OPSM (Order-Preserving Sub-Matrix) search.

  • One self-contained file — no external helpers
  • Optional column cap  (--k)
  • Optional # restarts  (--restarts)
  • Optional RNG seed    (--seed)
  • Tiny “demo mode”     (python OPSM_algorithm.py --demo)

Author: you & ChatGPT, 2025-06-30
"""
from __future__ import annotations

import argparse
import random
from typing import Sequence, Tuple

import numpy as np


# ────────────────────────────────────────────────────────────
# basic helpers
# ────────────────────────────────────────────────────────────
def pick_seed_column(n_cols: int, rng: random.Random) -> int:
    """Return a random column index in 0…n_cols-1."""
    return rng.randint(0, n_cols - 1)  # inclusive


def score_submatrix(
    X: np.ndarray,
    row_mask: np.ndarray,
    cols: Sequence[int],
) -> Tuple[int, np.ndarray]:
    """
    Count how many `row_mask==True` rows are *strictly increasing*
    across the ordered columns in `cols`.

    Returns
    -------
    support : int
        Number of surviving rows.
    survives : np.ndarray[bool]
        New Boolean mask the same length as rows.
    """
    sub = X[row_mask][:, cols]  # shrink to alive rows, chosen columns
    survives = (np.diff(sub, axis=1) > 0).all(axis=1)
    return int(survives.sum()), survives


# ────────────────────────────────────────────────────────────
# greedy grower for **one** seed
# ────────────────────────────────────────────────────────────
def greedy_expand(
    X: np.ndarray,
    k: int | None = None,
    rng: random.Random | None = None,
) -> Tuple[np.ndarray, list[int]]:
    """
    Grow a column permutation one step at a time, keeping the column
    that leaves the largest survivor pool.

    Parameters
    ----------
    X : np.ndarray
        Data matrix (rows × columns).
    k : int | None
        Maximum columns to keep.  None = grow until no column helps.
    rng : random.Random | None
        Random-number generator (so caller can control seeding).

    Returns
    -------
    surviving_rows : np.ndarray[int]
        Indices of rows that remain in the final pattern.
    cols : list[int]
        Ordered column indices of the OPSM found in this run.
    """
    if rng is None:
        rng = random.Random()

    n_cols = X.shape[1]
    seed_col = pick_seed_column(n_cols, rng)

    cols: list[int] = [seed_col]                       # permutation so far
    row_mask = np.ones(X.shape[0], dtype=bool)         # everyone alive
    best_support = row_mask.sum()                      # m rows

    while k is None or len(cols) < k:
        best_gain = -1
        best_col = None
        best_row_mask = None

        for c in range(n_cols):
            if c in cols:
                continue

            support, survives = score_submatrix(X, row_mask, cols + [c])
            gain = support - best_support
            if gain > best_gain:
                best_gain = gain
                best_col = c
                best_row_mask = survives

        # no column improves (or even matches) → stop
        if best_gain <= 0 or best_col is None:
            break

        cols.append(best_col)
        row_mask[row_mask] = best_row_mask  # shrink alive set **in-place**
        best_support += best_gain

    surviving_rows = np.where(row_mask)[0]
    return surviving_rows, cols


# ────────────────────────────────────────────────────────────
# public API with multiple restarts
# ────────────────────────────────────────────────────────────
def run_opsm(
    X: np.ndarray,
    k: int | None = None,
    restarts: int = 10,
    seed: int | None = None,
) -> Tuple[np.ndarray, list[int]]:
    """
    Run greedy OPSM `restarts` times from different seeds; keep the best pattern.

    Returns
    -------
    best_rows : np.ndarray[int]
    best_cols : list[int]
    """
    global_rng = random.Random(seed)
    best_rows: np.ndarray | None = None
    best_cols: list[int] | None = None
    best_score = -1

    for _ in range(restarts):
        sub_seed = global_rng.randint(0, 2**31 - 1)
        rng = random.Random(sub_seed)
        rows, cols = greedy_expand(X, k=k, rng=rng)
        if rows.size > best_score:
            best_rows, best_cols, best_score = rows, cols, rows.size

    assert best_rows is not None and best_cols is not None
    return best_rows, best_cols


# ────────────────────────────────────────────────────────────
# tiny CLI / demo
# ────────────────────────────────────────────────────────────
def _demo():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(20, 10))  # 20 × 10 random matrix
    rows, cols = run_opsm(X, k=None, restarts=5, seed=0)
    print("Survivors :", rows.size)
    print("Row idxs  :", rows)
    print("Col order :", cols)


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Greedy OPSM search")
    p.add_argument("file", nargs="?", help="CSV file (rows × cols) to load")
    p.add_argument("-k", type=int, default=None, help="max columns to keep")
    p.add_argument("--restarts", type=int, default=10, help="# random restarts")
    p.add_argument("--seed", type=int, default=None, help="RNG seed")
    p.add_argument("--demo", action="store_true", help="run built-in demo")
    return p.parse_args()


def _load_csv(path: str) -> np.ndarray:
    return np.genfromtxt(path, delimiter=",", dtype=float)


def main() -> None:
    args = _parse_args()

    if args.demo:
        _demo()
        return

    if not args.file:
        raise SystemExit("Error: must pass a CSV file or use --demo")

    X = _load_csv(args.file)
    rows, cols = run_opsm(X, k=args.k, restarts=args.restarts, seed=args.seed)

    print("Rows surviving:", rows.size)
    print("Row indices   :", rows)
    print("Column order  :", cols)


if __name__ == "__main__":
    main()


