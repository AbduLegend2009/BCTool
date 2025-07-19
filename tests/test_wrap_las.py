import numpy as np
import adapter


def test_wrap_las_returns_bicluster():
    np.random.seed(0)
    X = np.zeros((15, 25))
    X[:10, :20] = 10
    result = adapter.wrap_las(X, max_iter=10, alpha=0.05)
    assert isinstance(result, list)
    assert len(result) == 1
    bic = result[0]
    assert bic.rows.size > 0 and bic.cols.size > 0
    assert bic.score > 0
