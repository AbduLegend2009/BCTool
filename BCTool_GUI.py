import streamlit as st
import pandas as pd
import numpy as np

from adapter import ALL_ALGOS  # {name: callable -> list[Bicluster]}
from GO_assessment import go_assessment

# ────────────────────────────────────────────────────────────
# Helper: universal file‑loader
# ────────────────────────────────────────────────────────────

def load_matrix(file_obj) -> pd.DataFrame:
    """Return a *genes × conditions* numeric DataFrame.

    Accepts CSV, TSV/TXT (auto‑detected delimiter) or Excel (xls/xlsx).
    The first column is treated as gene IDs / row index.
    """
    fname = file_obj.name.lower()

    # Excel branch ─────────────────────
    if fname.endswith((".xls", ".xlsx")):
        df = pd.read_excel(file_obj, index_col=0)

    # Text branch ──────────────────────
    else:
        # Peek at first 4 KiB to guess delimiter
        head = file_obj.read(4096)
        file_obj.seek(0)
        comma = head.count(b",")
        tab   = head.count(b"\t")
        delim = "," if comma >= tab else "\t"
        df = pd.read_csv(file_obj, sep=delim, index_col=0)

    # Coerce to numeric, drop all‑NaN rows/cols
    df = df.apply(pd.to_numeric, errors="coerce")
    df.dropna(how="all", inplace=True)
    df.dropna(axis=1, how="all", inplace=True)

    if df.empty:
        raise ValueError("Parsed DataFrame is empty after numeric coercion.")

    # Sanitise index to string gene IDs
    df.index = df.index.map(str)
    return df

# ────────────────────────────────────────────────────────────
# Main Streamlit app
# ────────────────────────────────────────────────────────────

def run_biclusters():
    st.set_page_config(page_title="BCTool", layout="wide")
    st.title("🧬 BCTool – Biclustering & GO‑enrichment playground")

    # 1️⃣ Upload matrix ────────────────────────────────────
    upl = st.file_uploader(
        "Upload gene‑expression matrix (CSV, TSV, TXT, XLS/XLSX)",
        type=["csv", "tsv", "txt", "xls", "xlsx"],
        help="Rows = genes, columns = samples/conditions. First column must contain gene IDs.",
    )
    if upl is None:
        st.info("⬆️ Drag a file to begin.")
        return

    try:
        df = load_matrix(upl)
    except Exception as exc:
        st.error(f"❌ File‑read error: {exc}")
        return

    st.success(f"Loaded matrix → {df.shape[0]:,} genes × {df.shape[1]:,} conditions")

    # 2️⃣ Algorithm picker ─────────────────────────────────
    algo_choices = list(ALL_ALGOS.keys())
    default_pick = algo_choices[: min(len(algo_choices), 3)]
    selected = st.multiselect(
        "Choose biclustering algorithms to run",
        options=algo_choices,
        default=default_pick,
    )
    if not selected:
        st.warning("Select ≥1 algorithm to proceed.")
        return

    # 3️⃣ Execute button ───────────────────────────────────
    if st.button("🚀 Run selected algorithms"):
        with st.spinner("Running biclustering algorithms…"):
            results = {}
            for name in selected:
                try:
                    biclusters = ALL_ALGOS[name](df.values)
                    results[name] = biclusters
                except Exception as e:
                    st.error(f"{name} failed: {e}")
            st.success("Finished!")

        # 3a. Display per‑algo counts
        for name, bicls in results.items():
            st.subheader(f"**{name}** — {len(bicls)} biclusters")

        # 4️⃣ Optional GO enrichment ─────────────────────────
        if st.checkbox("Assess GO enrichment of results"):
            with st.spinner("Calculating GO terms…"):
                go_df = go_assessment(results, df.index)

            if go_df.empty:
                st.info("No significant GO terms at FDR ≤ 0.05.")
            else:
                st.dataframe(go_df, use_container_width=True)
                st.download_button(
                    "Download GO enrichment (CSV)",
                    go_df.to_csv(index=False).encode(),
                    "go_enrichment.csv",
                    mime="text/csv",
                )


if __name__ == "__main__":
    run_biclusters()
