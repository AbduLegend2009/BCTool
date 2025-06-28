import numpy as np
import random
def random_condition_selection(X:np.ndarray):
    m,n = X.shape
    random_condition_index=random.randint(0,n-1)
    return random_condition_index, X[:,random_condition_index]
def score_submatrix(surviving_genes: np.ndarray, chosen_conditions: list[int], X:np.ndarray):
    sub_matrix=X[surviving_genes][:,chosen_conditions]
    surviving_genes= (np.diff(sub, axis=1) > 0).all(axis=1)
    return surviving_genes.sum(),surviving_genes

