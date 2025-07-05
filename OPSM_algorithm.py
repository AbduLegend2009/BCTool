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
import argparse
import random
from typing import Sequence, Tuple, Optional

import numpy as np


# ────────────────────────────────────────────────────────────
# basic helpers
# ────────────────────────────────────────────────────────────
def pick_seed_column(n_cols: int, rng: random.Random) -> int:
    """Return a random column index in 0…n_cols-1."""
    return rng.randint(0, n_cols - 1)


def score_submatrix(
    X: np.ndarray,
    row_mask: np.ndarray,
    cols: Sequence[int],
) -> Tuple[int, np.ndarray]:
    """
    Count how many `row_mask==True` rows are *strictly increasing*
    across the ordered columns in `cols`.
    """
    sub = X[row_mask][:, cols]
    survives = (np.diff(sub, axis=1) > 0).all(axis=1)
    return int(survives.sum()), survives


# ────────────────────────────────────────────────────────────
# greedy grower for **one** seed
# ────────────────────────────────────────────────────────────
def greedy_expand(
    X: np.ndarray,
    k: Optional[int] = None,
    rng: Optional[random.Random] = None,
) -> Tuple[np.ndarray, list[int]]:
    """
    Grow a column permutation one step at a time, keeping the column
    that leaves the largest survivor pool.
    """
    if rng is None:
        rng = random.Random()

    n_cols = X.shape[1]
    seed_col = pick_seed_column(n_cols, rng)

    cols: list[int] = [seed_col]
    row_mask = np.ones(X.shape[0], dtype=bool)
    best_support = row_mask.sum()

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

        if best_gain <= 0 or best_col is None:
            break

        cols.append(best_col)
        row_indices = np.where(row_mask)[0]
        row_mask[:] = False
        row_mask[row_indices[best_row_mask]] = True
        best_support += best_gain

    surviving_rows = np.where(row_mask)[0]
    return surviving_rows, cols


# ────────────────────────────────────────────────────────────
# public API with multiple restarts
# ────────────────────────────────────────────────────────────
def run_opsm(
    X: np.ndarray,
    k: Optional[int] = None,
    restarts: int = 10,
    seed: Optional[int] = None,
) -> Tuple[np.ndarray, list[int]]:
    """
    Run greedy OPSM `restarts` times from different seeds; keep the best pattern.
    """
    global_rng = random.Random(seed)
    best_rows: Optional[np.ndarray] = None
    best_cols: Optional[list[int]] = None
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
    seed = 0
    rng = np.random.default_rng(seed)
    random.seed(seed)  # ensure consistency with random module
    X = rng.normal(size=(20, 10))
    rows, cols = run_opsm(X, k=None, restarts=5, seed=seed)
    print("Survivors :", rows.size)
    print("Row idxs  :", rows)
    print("Col order :", cols)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Greedy OPSM search")
    parser.add_argument("file", nargs="?",
                        help="CSV file (rows × cols) to load")
    parser.add_argument("-k", type=int, default=None,
                        help="max columns to keep")
    parser.add_argument("--restarts", type=int, default=10,
                        help="# random restarts")
    parser.add_argument("--seed", type=int, default=None, help="RNG seed")
    parser.add_argument("--demo", action="store_true",
                        help="run built-in demo")
    return parser.parse_args()


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
    rows, cols = run_opsm(X, k=args.k, restarts=args.restarts, seed=args.s_
