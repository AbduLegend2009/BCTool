import streamlit as st
import pandas as pd
from GO_assessment import go_assessment
from adapter import ALL_ALGOS, Bicluster
def _matrix_uploader():
    uploaded = st.sidebar.file_uploader(
        "Upload expression matrix",
        type=["csv", "tsv", "txt", "xls", "xlsx"],
        help="Rows = genes, columns = conditions; first column will be treated as gene IDs",
    )
    if not uploaded:
        return None
    try:
        df = load_matrix(uploaded)
    except Exception as err:
        st.sidebar.error(f"❌ {err}")
        return None
    st.sidebar.success(f"✓ Loaded **{uploaded.name}**  →  shape {df.shape[0]}×{df.shape[1]}")
    return df

def load_matrix(file_obj):
    """Return a numeric *genes × conditions* DataFrame from a CSV/TSV/TXT/Excel.

    The parser auto‑detects tab vs comma delimiters for text files and
    coerces everything to *float* (non‑numeric cells → NaN → dropped if whole
    row/column is missing).  Unsupported formats or empty numeric tables raise
    *ValueError* so the GUI can surface a clear message.
    """
    name = file_obj.name.lower()

    # ── flat text (CSV / TSV / TXT) ──────────────────────────
    if name.endswith((".csv", ".tsv", ".txt")):
        sample = file_obj.read(4096).decode("utf‑8", errors="ignore")
        file_obj.seek(0)  # rewind for real read
        sep = "\t" if sample.count("\t") > sample.count(",") else ","
        df = pd.read_csv(file_obj, sep=sep, index_col=0)

    # ── Excel ────────────────────────────────────────────────
    elif name.endswith((".xls", ".xlsx")):
        df = pd.read_excel(file_obj, index_col=0)

    else:
        raise ValueError("Unsupported format – please upload CSV, TSV/TXT, or Excel.")

    # ensure numeric
    df = df.apply(pd.to_numeric, errors="coerce")
    df.dropna(axis=0, how="all", inplace=True)
    df.dropna(axis=1, how="all", inplace=True)
    if df.empty:
        raise ValueError("No numeric data found in the uploaded file.")

    df.index = df.index.astype(str)  # gene IDs as strings
    return df
def run_biclusters():
    st.title("BCTool – Streamlit prototype")

    df = _matrix_uploader()
    if df is None:
        st.info("Please upload a matrix to begin.")
        st.stop()

    st.header("Algorithm selection")
    chosen = st.multiselect(
        "Select biclustering algorithms to run",
        options=list(ALL_ALGOS.keys()),
        default=list(ALL_ALGOS.keys())[:1],
    )

    if st.button("Run selected algorithms"):
        with st.spinner("Running biclustering algorithms…"):
            results = {}
            for name in chosen:
                model = ALL_ALGOS[name]()
                biclusters = model.fit_transform(df.values)
                results[name] = biclusters
        st.success("Finished!")

        for name, bicls in results.items():
            st.subheader(f"{name}: {len(bicls)} biclusters")

        if st.checkbox("Assess GO enrichment of results"):
            with st.spinner("Calculating GO terms…"):
                go_summary = go_assessment(results, df.index)
            st.write(go_summary)


if __name__ == "__main__":
    run_biclusters()

