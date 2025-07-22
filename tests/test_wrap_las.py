from pathlib import Path
import pandas as pd
import adapter


def test_wrap_las_returns_bicluster():
    data_path = Path(__file__).with_name("data") / "sample_matrix.csv"
    X = pd.read_csv(data_path, index_col=0).values
    result = adapter.wrap_las(X, max_iter=10, alpha=0.05)
    assert isinstance(result, list)
    assert len(result) == 1
    bic = result[0]
    assert bic.rows.size > 0 and bic.cols.size > 0
    assert bic.score > 0
