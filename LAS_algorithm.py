import numpy as np
import math
from scipy.stats import norm


def las_with_significance(
    test_array,
    testarray_rows,
    testarray_cols,
    *,
    max_iter=100,
    alpha=0.05,
    k_rows=10,
    k_cols=20,
):
   
    #variables
    num_rows = testarray_rows
    standard_deviation = np.nanstd(test_array)
    array_mean = np.nanmean(test_array)
    alpha_level = alpha
    k_rows = int(k_rows)
    k_cols = int(k_cols)

    best_z_score = -np.inf
    best_p_score = 1.0
    best_submatrix = None
    best_rows = None
    best_cols = None

    # this is the loop tha creates random submatrices and calculates the means of them
    for _ in range(max_iter):
        rows = np.random.choice(num_rows, k_rows, replace=False)
        column_means = np.nanmean(test_array[rows, :], axis=0)
        # sorts the top k columns by means
        cols = np.argsort(column_means)[-k_cols:]
        row_means = np.nanmean(test_array[:, cols], axis=1)
        # sorts the top k columns by rows
        rows = np.argsort(row_means)[-k_rows:]
        new_submatrix = test_array[np.ix_(rows, cols)]  # new submatrix
        new_submatrix_mean = np.nanmean(new_submatrix)  # submatrix mean

        # the actual formula for calculating significance (calculating the required z scores)
        standard_error = math.sqrt(k_rows * k_cols)
        denominator = standard_deviation / standard_error
        numerator = new_submatrix_mean - array_mean
        z_score = numerator / denominator

        # coverting z score to p value
        p_value = 2 * (1 - norm.cdf(abs(z_score)))

        # testing
        if p_value < alpha_level and z_score > best_z_score:
            best_z_score = z_score
            best_p_score = p_value
            best_submatrix = new_submatrix.copy()
            best_rows = rows
            best_cols = cols

    # return results
    return {
        'z_score': best_z_score,
        'p_value': best_p_score,
        'significant': best_p_score < alpha_level,
        'submatrix': best_submatrix,
        'rows': best_rows,
        'cols': best_cols
    }

# -- new multi bicluster LAS --

def las_multi(
    X: np.ndarray,
    *,
    n_biclusters: int = 10,
    max_iter: int = 100,
    alpha: float = 0.05,
    k_rows: int = 10,
    k_cols: int = 20,
    random_state: int | None = None,
):
    """Find multiple biclusters using repeated LAS runs.

    Parameters
    ----------
    X : np.ndarray
        Input expression matrix.
    n_biclusters : int, optional
        Maximum number of biclusters to return.
    max_iter : int, optional
        Iterations per LAS run.
    alpha : float, optional
        Significance threshold for bicluster acceptance.
    k_rows : int, optional
        Number of rows in candidate submatrices.
    k_cols : int, optional
        Number of columns in candidate submatrices.
    random_state : int | None, optional
        Random seed for reproducibility.
    """
    if random_state is not None:
        np.random.seed(random_state)
    X_work = X.copy()
    bics = []
    for _ in range(n_biclusters):
        if np.isnan(X_work).all():
            break
        res = las_with_significance(
            X_work,
            X_work.shape[0],
            X_work.shape[1],
            max_iter=max_iter,
            alpha=alpha,
            k_rows=k_rows,
            k_cols=k_cols,
        )
        rows = res['rows']
        cols = res['cols']
        if rows is None or cols is None:
            break
        bics.append(res)
        X_work[np.ix_(rows, cols)] = np.nan
    return bics
