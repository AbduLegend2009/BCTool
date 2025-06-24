import numpy as np
from itertools import combinations
from collections import defaultdict

# ------------------------------------------------------------
# Core helper: fit additive or multiplicative model on a block
# ------------------------------------------------------------

def fit_block(X_sub, *, model="add", eps=1e-9):
    """Return MSR, reconstructed block, row factors a, column factors b."""
    # 1. optional log‑transform
    if model == "mult":
        Y = np.log(X_sub + eps)
    else:
        Y = X_sub

    # 2. additive fit in Y space
    mu = Y.mean()
    a = Y.mean(axis=1, keepdims=True) - mu  # (m,1)
    b = Y.mean(axis=0, keepdims=True) - mu  # (1,n)
    recon = mu + a + b                      # broadcast to (m,n)
    msr = np.mean((Y - recon) ** 2)         # mean‑squared residue

    # 3. map back if multiplicative
    if model == "mult":
        recon = np.exp(recon)
        a = np.exp(a).ravel()               # row scales 1‑D
        b = np.exp(b).ravel()               # col scales 1‑D
    else:
        a = a.ravel()                       # row shifts 1‑D
        b = b.ravel()                       # col shifts 1‑D

    return msr, recon, a, b

# ------------------------------------------------------------
# 1‑D histogram helper for signature binning
# ------------------------------------------------------------

def quantise(vec, *, eps):
    """Bin 1‑D vector with width eps → {bin_label: [row indices]}"""
    keys = np.round(vec / eps)
    bins = defaultdict(list)
    for idx, k in enumerate(keys):
        bins[k].append(idx)
    return bins

# ------------------------------------------------------------
# Generate candidate (row set, column set) seeds
# ------------------------------------------------------------

def generate_candidate_pairs(X, *, eps=0.3, min_rows=5, min_cols=2):
    m, n = X.shape
    row_signature = defaultdict(set)  # {frozenset(rows): set(cols)}

    for p, q in combinations(range(n), 2):
        diff = X[:, p] - X[:, q]                      # additive signature
        ratio = np.log(X[:, p] + 1e-9) - np.log(X[:, q] + 1e-9)  # multiplicative

        for vec in (diff, ratio):
            for _, rows in quantise(vec, eps=eps).items():
                if len(rows) >= min_rows:
                    row_signature[frozenset(rows)].update((p, q))

    seeds = []
    for rows, cols in row_signature.items():
        if len(cols) >= min_cols:
            seeds.append((set(rows), set(cols)))
    return seeds

# ------------------------------------------------------------
# Main BiVisu dispatcher
# ------------------------------------------------------------

def bivisu(
    X,
    *,
    model="auto",          # 'add' | 'mult' | 'auto'
    eps=0.3,                # bin width (z‑scored units)
    thr=0.05,               # MSR threshold for coherence
    min_rows=5,
    min_cols=2,
    max_iters=5,
):
    """Return list of biclusters as (row_set, col_set) tuples."""
    X = np.asarray(X, float)
    seeds = generate_candidate_pairs(
        X, eps=eps, min_rows=min_rows, min_cols=min_cols
    )

    biclusters = []
    for _ in range(max_iters):
        new_winners = []
        for R, C in seeds:
            X_sub = X[np.ix_(list(R), list(C))]

            # choose additive vs multiplicative
            if model == "auto":
                msr_add, *_ = fit_block(X_sub, model="add")
                msr_mul, *_ = fit_block(X_sub, model="mult")
                msr, chosen = (
                    (msr_add, "add") if msr_add < msr_mul else (msr_mul, "mult")
                )
            else:
                msr, *_ = fit_block(X_sub, model=model)
                chosen = model

            # global block passes?
            if msr <= thr:
                # prune noisy rows once
                good_rows = {
                    i
                    for i in R
                    if fit_block(X[i, list(C)][None, :], model=chosen)[0] <= thr
                }
                if len(good_rows) >= min_rows:
                    new_winners.append((good_rows, C))

        # convergence check
        if not new_winners or new_winners == biclusters:
            break
        biclusters = new_winners
        seeds = biclusters  # recycle winners for next iteration

    return biclusters

# ------------------------------------------------------------
# Example usage (comment out when importing as library)
# ------------------------------------------------------------
if __name__ == "__main__":
    # small synthetic example
    X_demo = np.array([
        [10, 12, 14, 18],
        [ 6,  8, 10, 12],
        [15, 18, 21, 24],
    ], dtype=float)

    bicls = bivisu(X_demo, model="auto", eps=0.25, thr=0.04, min_rows=2, min_cols=2)
    for k, (R, C) in enumerate(bicls, 1):
        print(f"Bicluster {k}: rows {sorted(R)} cols {sorted(C)}")
