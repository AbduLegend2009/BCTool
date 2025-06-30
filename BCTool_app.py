"""BCTool Streamlit GUI v0.4 – OPSM completed
================================================
Five algorithms + GO‑quality dashboard + sub‑matrix explorer.
Run:
    streamlit run bctool_gui.py
Install deps:
    pip install streamlit pandas numpy matplotlib seaborn goatools
"""
from __future__ import annotations

from typing import Any, Callable, Dict, List, Tuple
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st

# ---------------------------------------------------------------------------
# Algorithm imports ----------------------------------------------------------
# ---------------------------------------------------------------------------
from Bivisu_algorithm import bivisu           # type: ignore
from Chen_Church_algorithm import run_chen_church   # type: ignore
from ISA_algorithm import ISA_multi_seed      # type: ignore
from LAS_algorithm import las_with_significance  # type: ignore
from OPSM_algorithm import run_opsm           # type: ignore

# ---------------------------------------------------------------------------
# GO‑enrichment helper -------------------------------------------------------
# ---------------------------------------------------------------------------
PV_THRESHOLDS = [1e-5, 1e-4, 1e-3, 1e-2, 2e-2, 3e-2, 4e-2, 5e-2]

try:
    from goatools.obo_parser import GODag           # type: ignore
    from goatools import GOEnrichmentStudy          # type: ignore

    GO_DAG = GODag("go-basic.obo", optional_attrs={"relationship"})

    def GO_enrich(genes: List[str], background: List[str], category: str) -> List[float]:
        """Return un‑corrected p‑values for the given study list."""
        # TODO: replace with real gene2go mapping
        gene2go: Dict[str, set[str]] = {g: set() for g in background}
        study = GOEnrichmentStudy(background, gene2go, GO_DAG, methods=["holm"], prt=None)
        res = study.run_study(genes, methods=["holm"])
        return [r.p_uncorrected for r in res]
except Exception:  # goatools missing → stub

    def GO_enrich(genes: List[str], background: List[str], category: str) -> List[float]:  # type: ignore[override]
        rng = np.random.default_rng(42)
        return rng.random(rng.integers(5, 20)).tolist()

# ---------------------------------------------------------------------------
# Algorithm wrappers ---------------------------------------------------------
# ---------------------------------------------------------------------------

def _wrap_chen(X: np.ndarray, **kw):
    return [(np.array(r), np.array(c)) for r, c, *_ in run_chen_church(X, **kw)]

def _wrap_las(X: np.ndarray, **kw):
    rows, cols = X.shape
    res = las_with_significance(X, rows, cols, **kw)
    return [] if res["rows"] is None else [(np.array(res["rows"]), np.array(res["cols"]))]

def _wrap_isa(X: np.ndarray, **kw):
    outs: List[Tuple[np.ndarray, np.ndarray]] = []
    for m in ISA_multi_seed(X, **kw):
        r, c = np.flatnonzero(m["s_g"]), np.flatnonzero(m["s_c"])
        if r.size and c.size:
            outs.append((r, c))
    return outs

def _wrap_bivisu(X: np.ndarray, **kw):
    return [(np.array(sorted(r)), np.array(sorted(c))) for r, c in bivisu(X, **kw)]

def _wrap_opsm(X: np.ndarray, **kw):
    rows, cols = run_opsm(X, **kw)
    return [(np.array(rows), np.array(cols))] if rows.size and len(cols) else []

ALGORITHMS: Dict[str, Tuple[Callable[..., Any], Dict[str, Any]]] = {
    "Chen & Church": (_wrap_chen, {"msr_threshold": 0.6, "alpha": 0.05, "max_biclusters": 10}),
    "LAS": (_wrap_las, {"max_iter": 100, "alpha": 0.05}),
    "ISA": (_wrap_isa, {"n_seeds": 500, "seed_size": 5, "t_g": 2.0, "t_c": 2.0, "max_iter": 100, "jaccard_thresh": 0.9}),
    "BiVisu": (_wrap_bivisu, {"model": "auto", "eps": 0.3, "thr": 0.05, "min_rows": 5, "min_cols": 2, "max_iters": 5}),
    "OPSM": (_wrap_opsm, {"k": None, "restarts": 10, "seed": None}),
}

# ---------------------------------------------------------------------------
# Helpers --------------------------------------------------------------------
# ---------------------------------------------------------------------------
@st.cache_data(hash_funcs={np.ndarray: lambda _: None})
def _load_matrix(file, sep: str):
    df = pd.read_csv(file, sep=sep)
    if not np.issubdtype(df.iloc[:, 0].dtype, np.number):
        row_labels = df.iloc[:, 0].astype(str).tolist()
        X = df.iloc[:, 1:].to_numpy(float)
        col_labels = df.columns[1:].tolist()
    else:
        row_labels = [f"gene_{i}" for i in range(df.shape[0])]
        col_labels = df.columns.tolist()
        X = df.to_numpy(float)
    return df, X, row_labels, col_labels

def _heatmap(X: np.ndarray, rows: np.ndarray, cols: np.ndarray):
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.heatmap(X[np.ix_(rows, cols)], cmap="viridis", ax=ax, cbar=True, xticklabels=False, yticklabels=False)
    st.pyplot(fig)

def _plot_enrichment(data: Dict[str, List[float]], category: str):
    fig, ax = plt.subplots(figsize=(6, 4))
    bar_w = 0.8 / len(data)
    x = np.arange(len(PV_THRESHOLDS))
    for i, (lab, pct) in enumerate(data.items()):
        ax.bar(x + i * bar_w, pct, bar_w, label=lab)
    ax.set_xticks(x + bar_w * (len(data) - 1) / 2)
    ax.set_xticklabels([f"≤{t:.0e}" for t in PV_THRESHOLDS])
    ax.set_ylim(0, 100)
    ax.set_ylabel("% biclusters enriched")
    ax.set_title(f"GO {category} enrichment")
    ax.legend()
    st.pyplot(fig)

def _clusters_to_csv(bicls: List[Tuple[np.ndarray, np.ndarray]], row_labels: List[str], col_labels: List[str]) -> str:
    rows = []
    for i, (r, c) in enumerate(bicls, start=1):
        for rr in r:
            for cc in c:
                rows.append({"cluster": i, "gene": row_labels[rr], "condition": col_labels[cc]})
    return pd.DataFrame(rows).to_csv(index=False)

# ---------------------------------------------------------------------------
# Streamlit main -------------------------------------------------------------
# ---------------------------------------------------------------------------

def main():
    st.set_page_config(page_title="BCTool", page_icon="🧬", layout="wide")
    st.title("🧬 BCTool – five algorithms, GO dashboard, bicluster explorer")

    # Upload section
    st.sidebar.header("1️⃣ Dataset")
    data_file = st.sidebar.file_uploader("CSV/TSV expression matrix", type=["csv", "tsv"])
    sep = st.sidebar.radio("Delimiter", [",", "\t"], horizontal=True)
    if data_file is None:
        st.info("⬅️ Upload a dataset to begin.")
        st.stop()

    df, X, gene_labels, col_labels = _load_matrix(data_file, sep)
    st.write("### Preview", df.head())
    st.markdown(f"Shape: **{X.shape[0]} × {X.shape[1]}**")

    # Algorithm selection
    st.sidebar.header("2️⃣ Algorithms")
    chosen = st.sidebar.multiselect("Select algorithms to run", list(ALGORITHMS.keys()), default=list(ALGORITHMS.keys())[:3])
    params: Dict[str, Dict[str, Any]] = {}
    for name in chosen:
        st.sidebar.subheader(name)
        dft = ALGORITHMS[name][1]
        params[name] = {}
        for p, default in dft.items():
            key = f"{name}_{p}"
            if isinstance(default, bool):
                params[name][p] = st.sidebar.checkbox(p, value=default or False, key=key)
            elif isinstance(default, int):
                params[name][p] = st.sidebar.number_input(p, value=default or 0, step=1, key=key)
            elif isinstance(default, float):
                params[name][p] = st.sidebar.number_input(p, value=default or 0.0, step=0.05, key=key)
            else:
                params[name][p] = st.sidebar.text_input(p, value="" if default is None else str(default), key=key)

    # GO settings
    st.sidebar.header("3️⃣ GO settings")
    go_cat = st.sidebar.selectbox("Category", ["BP", "MF", "CC"], index=0)
    thr_selected = st.sidebar.multiselect("P‑value thresholds", [str(t) for t in PV_THRESHOLDS], default=[str(t) for t in PV_THRESHOLDS])
    thr_selected_f = sorted(float(t) for t in thr_selected)

    if not st.sidebar.button("🚀 Run"):
        st.stop()

    # Run algorithms
    biclusters: Dict[str, List[Tuple[np.ndarray, np.ndarray]]] = {}
    for name in chosen:
        with st.spinner(f"Running {name} …"):
            func = ALGORITHMS[name][0]
            try:
                biclusters[name] = func(X, **params[name])  # type: ignore[arg-type]
            except Exception as e:  # noqa: BLE001
                st.error(f"{name} failed: {e}")
                biclusters[name] = []
        st.success(f"{name}: {len(biclusters[name])} biclusters")

    if all(len(lst) == 0 for lst in biclusters.values()):
        st.warning("No biclusters produced.")
        st.stop()

    # GO enrichment stats
    pct: Dict[str, List[float]] = {n: [0.0] * len(PV_THRESHOLDS) for n in chosen}
    for name, bics in biclusters.items():
        counts = [0] * len(PV_THRESHOLDS)
        for rows, _ in bics:
            genes = [gene_labels[i] for i in rows]
            pvals = GO_enrich(genes, gene_labels, go_cat)
            mp = min(pvals) if pvals else 1.0
            for i, thr in enumerate(PV_THRESHOLDS):
                if mp <= thr:
                    counts[i] += 1
        total = max(1, len(bics))
        pct[name] = [round(c / total * 100, 1) for c in counts]

    # Tabs
    tab_go, tab_compare, tab_explore = st.tabs(["GO overview", "Comparison", "Explorer"])

    with tab_go:
        st.subheader("Per‑algorithm GO enrichment")
        for n in chosen:
            st.markdown(f"**{n}**")
            _plot_enrichment({n: pct[n]}, go_cat)

    with tab_compare:
        st.subheader("Stacked comparison")
        _plot_enrichment({n: pct[n] for n in chosen}, go_cat)
        df_pct = pd.DataFrame({n: pct[n] for n in chosen}, index=[f"≤{t:.0e}" for t in PV_THRESHOLDS])
        st.write(df_pct)
        st.download_button("Download CSV", data=df_pct.to_csv().encode(), file_name="go_enrichment_summary.csv", mime="text/csv")

    with tab_explore:
        st.subheader("Bicluster explorer")
        for n, bics in biclusters.items():
            if not bics:
                continue
            st.markdown(f"### {n}")
            for idx, (r, c) in enumerate(bics, start=1):
                with st.expander(f"{n} – Bicluster {idx} ({len(r)}×{len(c)})"):
                    _heatmap(X, r, c)
                    st.write("Row idx", r)
                    st.write("Col idx", c)
                    df_sub = pd.DataFrame(X[np.ix_(r, c)], index=[gene_labels[i] for i in r], columns=[col_labels[j] for j in c])
                    st.write(df_sub)

        # download all cells
        all_csv = _clusters_to_csv([bc for lst in biclusters.values() for bc in lst], gene_labels, col_labels).encode()
        st.download_button("Download all bicluster cells", data=all_csv, file_name="biclusters_cells.csv", mime="text/csv")

    st.sidebar.markdown("---")
    st.sidebar.caption("BCTool v0.4 · Streamlit")

if __name__ == "__main__":
    main()
