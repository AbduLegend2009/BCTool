import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path

from adapter import ALL_ALGOS  # {name: callable -> list[Bicluster]}
from GO_assessment import go_assessment

# ────────────────────────────────────────────────────────────
# Helper: robust file‑loader (CSV / TSV / Excel)
# ────────────────────────────────────────────────────────────

def load_matrix(file_obj) -> pd.DataFrame:
    """Return a numeric *genes × conditions* DataFrame regardless of
    whether the user uploads CSV, TSV/TXT, or XLS/XLSX.
    The first column is taken as the gene identifier / index.
    """
    # Determine file suffix
    suffix = Path(file_obj.name).suffix.lower()

    if suffix in {".xls", ".xlsx"}:
        df = pd.read_excel(file_obj, index_col=0, engine="openpyxl")
    else:
        # Let pandas sniff delimiter automatically (`sep=None`) with the
        # slower python engine – fine for uploads <= 200 MB.
        df = pd.read_csv(file_obj, index_col=0, sep=None, engine="python")

    # Coerce to numeric and drop non‑numeric cols if any
    df = df.apply(pd.to_numeric, errors="coerce")
    df.dropna(axis=1, how="all", inplace=True)

    # Sanity‑check
    if df.empty:
        st.error("Uploaded file contains no numeric expression values.")
        st.stop()

    # Normalise index to str so adapter / GO tools are happy
    df.index = df.index.map(str)
    return df.astype(np.float32, copy=False)  # memory‑friendly


# ────────────────────────────────────────────────────────────
# Streamlit app
# ────────────────────────────────────────────────────────────

st.set_page_config(page_title="BCTool", layout="wide")

# persistent storage across widget interactions
for key in ("df", "results", "go_df"):
    if key not in st.session_state:
        st.session_state[key] = None

st.title("🧬 BCTool – Biclustering & GO‑enrichment playground")

# 1️⃣ Upload matrix -------------------------------------------------------
uploaded = st.file_uploader("Upload gene‑expression matrix (CSV, TSV, TXT, XLS/XLSX)",
                            type=["csv", "tsv", "txt", "xls", "xlsx"],
                            help="First column = gene IDs. Remaining columns = samples / conditions.")

if uploaded:
    st.session_state.df = load_matrix(uploaded)
    st.success(f"Loaded matrix → {st.session_state.df.shape[0]} genes × {st.session_state.df.shape[1]} conditions")

# Stop early if no data yet
if st.session_state.df is None:
    st.stop()

df = st.session_state.df

# 2️⃣ Algorithm picker ----------------------------------------------------
chosen = st.multiselect("Choose biclustering algorithms to run",
                       options=list(ALL_ALGOS.keys()),
                       default=list(ALL_ALGOS.keys())[:2])

# 3️⃣ Run algorithms ------------------------------------------------------
if st.button("🚀 Run selected algorithms"):
    if not chosen:
        st.warning("Select at least one algorithm!")
        st.stop()

    with st.spinner("Running biclustering algorithms …"):
        results = {name: ALL_ALGOS[name](df.values) for name in chosen}

    st.session_state.results = results

    with st.spinner("Calculating GO‑term enrichment … (first time may download data)"):
        go_df = go_assessment(results, df.index)
    st.session_state.go_df = go_df

    st.success("Finished!")

# Stop until user clicks run
if st.session_state.results is None:
    st.stop()

results = st.session_state.results

# 4️⃣ Display biclusters --------------------------------------------------
for algo, bicls in results.items():
    st.subheader(f"{algo} — {len(bicls)} biclusters")
    for i, bc in enumerate(bicls, start=1):
        genes   = getattr(bc, "genes", getattr(bc, "rows", []))
        conds   = getattr(bc, "conditions", getattr(bc, "cols", []))
        header  = f"Bicluster {i}: {len(genes)} genes × {len(conds)} conditions"
        with st.expander(header):
            left, right = st.columns(2)
            left.markdown("**Genes (first 20)**")
            left.write(list(genes)[:20])
            right.markdown("**Conditions**")
            right.write(list(conds))

# 5️⃣ Show GO‑enrichment --------------------------------------------------
go_df = st.session_state.go_df
if go_df is not None:
    st.subheader("GO‑term enrichment across all biclusters")
    st.dataframe(go_df, use_container_width=True, hide_index=True)
    st.download_button("Download enrichment CSV",
                       go_df.to_csv(index=False).encode(),
                       "go_enrichment.csv",
                       mime="text/csv")
