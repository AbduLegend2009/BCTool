import numpy as np
import math
from scipy.stats import norm


def las_with_significance(test_array, testarray_rows, testarray_cols, max_iter=100, alpha=0.05):
   
    #variables
    num_rows = testarray_rows
    num_columns = testarray_cols
    standard_deviation = test_array.std()
    array_mean = test_array.mean()
    alpha_level = alpha
    k_rows = 10
    k_cols = 20

    best_z_score = -np.inf
    best_p_score = 1.0
    best_submatrix = None
    best_rows = None
    best_cols = None

    # this is the loop tha creates random submatrices and calculates the means of them
    for _ in range(max_iter):
        rows = np.random.choice(num_rows, k_rows, replace=False)
        column_means = test_array[rows, :].mean(axis=0)
        # sorts the top k columns by means
        cols = np.argsort(column_means)[-k_cols:]
        row_means = test_array[:, cols].mean(axis=1)
        # sorts the top k columns by rows
        rows = np.argsort(row_means)[-k_rows:]
        new_submatrix = test_array[np.ix_(rows, cols)]  # new submatrix
        new_submatrix_mean = new_submatrix.mean()  # submatrix mean

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
