
from __future__ import annotations
from typing import List, Dict
import numpy as np
__all__ = [
    "ISA_single_seed",
    "ISA_multi_seed",
]


def ISA_single_seed(X: np.ndarray, t_g, t_c, seed_size=5, max_iter=100):
    # Determine rows and columns
    n_genes, m_conditions = X.shape
    if seed_size > m_conditions:
        raise ValueError(
            f"seed_size ({seed_size}) exceeds the number of conditions "
            f"({m_conditions}); choose seed_size ≤ m_conditions.")
# update
    # Create new nparray with all zeros for unselected conditions
    s_c_old = np.zeros(m_conditions, dtype=int)
    # Choose random indices with same amount of selected conditions as the seed size
    seed_indices = np.random.choice(
        m_conditions, size=seed_size, replace=False)
    # Mark with 1 to show that conditions have been selected
    s_c_old[seed_indices] = 1
    # For genes to show that none have been selected yet
    s_g_old = np.zeros(n_genes, dtype=int)
    for i in range(max_iter):
        # Sum each gene's values over the currently selected conditions
        gene_scores = X.dot(s_c_old)
        k = s_c_old.sum()
        k = max(k, 1)
        gene_scores = gene_scores/np.sqrt(k)
        # We need to find the mean of the gene_scores to find the z_score
        gene_scores_mean = gene_scores.mean()
        # We need to find the standard deviation of the gene_scores to find the z_score
        gene_scores_standard_deviation = gene_scores.std()
        # If the standard deviation is 0, then all of the values are the same
        # (and the mean thus the same as all the values). Furthermore, the
        # z_scores would all be 0 as the formula for z_scores is
        # (gene_scores-gene_scores_mean)/gene_scores_standard_deviation
        if gene_scores_standard_deviation == 0:
            gene_scores_z_scores = np.zeros_like(gene_scores)
        # Otherwise, we calculate the z_scores normally
        else:
            gene_scores_z_scores = (
                gene_scores-gene_scores_mean)/gene_scores_standard_deviation
        # We now create an array that is full of 1s and 0s to see which genes
        # have been selected because they passed the gene z_score threshold
        s_g_new = (gene_scores_z_scores >= t_g).astype(int)
        # X.T changes the way the data is viewed: it was (n_genes,m_conditions),
        # but is now (m_condition,n_genes)
        condition_scores = X.T.dot(s_g_new)
        k = s_g_new.sum()
        k = max(k, 1)
        condition_scores = condition_scores/np.sqrt(k)
        # Repeating process as what we did for gene_scores
        condition_scores_mean = condition_scores.mean()

        condition_scores_standard_deviation = condition_scores.std()

        if condition_scores_standard_deviation == 0:

            condition_scores_z_scores = np.zeros_like(condition_scores)
        else:

            condition_scores_z_scores = (
                condition_scores-condition_scores_mean)/condition_scores_standard_deviation

        s_c_new = (condition_scores_z_scores >= t_c).astype(int)
        if s_g_new.sum() == 0 or s_c_new.sum() == 0:
            raise ValueError(
                "ISA: module collapsed — thresholds too high or seed too small")

        if np.array_equal(s_g_new, s_g_old) and np.array_equal(s_c_new, s_c_old):

            break
        s_g_old = s_g_new.copy()
        s_c_old = s_c_new.copy()
    return s_g_new, s_c_new


def _jaccard(mask_a: np.ndarray, mask_b: np.ndarray) -> float:
    intersect = np.logical_and(mask_a, mask_b).sum()
    union = np.logical_or(mask_a, mask_b).sum()
    return 0.0 if union == 0 else intersect/union


def ISA_multi_seed(
    X: np.ndarray,
    *,
    n_seeds: int = 500,
    seed_size: int = 5,
    t_g: float = 2.0,
    t_c: float = 2.0,
    max_iter: int = 100,
    jaccard_thresh: float = 0.90,
) -> List[Dict[str, np.ndarray]]:
    np.random.seed(0)
    modules: List[Dict[str, np.ndarray]] = []
    for _ in range(n_seeds):
        try:
            s_g_new, s_c_new = ISA_single_seed(
                X,
                seed_size=seed_size,
                t_g=t_g,
                t_c=t_c,
                max_iter=max_iter,
            )
        except ValueError:
            continue
        duplicate = False
        for m in modules:
            if _jaccard(s_g_new, m["s_g"]) >= jaccard_thresh and _jaccard(s_c_new, m["s_c"]) >= jaccard_thresh:
                duplicate = True
                break
        if duplicate:
            continue
        sub = X[np.ix_(s_g_new.astype(bool),   # rows  (genes in module)
                       s_c_new.astype(bool))]  # cols  (conds in module)
        score = np.abs(sub.mean(axis=1)).mean()
        modules.append({"s_g": s_g_new, "s_c": s_c_new, "score": score})
    modules.sort(key=lambda m: m["score"], reverse=True)
    return modules
