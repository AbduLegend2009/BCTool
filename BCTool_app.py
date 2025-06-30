"""BCTool Streamlit GUI – algorithms *plus* GO‑enrichment dashboard
===============================================================
This single‑file Streamlit app now handles **four biclustering algorithms**
(Chen‑&‑Church, LAS, ISA, BiVisu) *and* a light‑weight Gene Ontology
(GO) quality‑assessment pipeline:

1. Upload a gene‑expression matrix.
2. Choose any subset of algorithms and tweak their parameters.
3. Click **Run biclustering + GO** – the app:
   * executes each algorithm;
   * displays every bicluster (heat‑map, row/col indices);
   * performs GO enrichment for each bicluster (stubbed – see below);
   * shows bar‑charts comparing % enriched clusters across eight p‑value
     thresholds.
4. Optionally upload external bicluster sets (CSV / JSON) so the same
   GO dashboard can rank them side‑by‑side with in‑house results.

Dependencies (install once):
    pip install streamlit pandas numpy matplotlib seaborn goatools

*`goatools`* is used if available; otherwise a tiny **dummy_enrich()**
function returns random p‑values so the GUI works for layout tests.
Replace `GO_enrich()` with your real pipeline when ready.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple

import json
import math
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st

# ----------------------------------------------------------------------------
# Algorithm imports (same folder) --------------------------------------------
# ----------------------------------------------------------------------------
from Bivisu_algorithm import bivisu  # type: ignore
from Chen_Church_algorithm import run_chen_church  # type: ignore
from ISA_algorithm import ISA_multi_seed  # type: ignore
from LAS_algorithm import las_with_significance  # type: ignore

# ----------------------------------------------------------------------------
# GO‑enrichment helper --------------------------------------------------------
# ----------------------------------------------------------------------------

# Eight thresholds used in the PDF spec
PV_THRESHOLDS = [1e-5, 1e-4, 1e-3, 1e-2, 2e-2, 3e-2, 4e-2, 5e-2]

try:
    # Lazy import to keep goatools optional
    from goatools.obo_parser import GODag  # type: ignore
    from goatools import GOEnrichmentStudy  # type: ignore

    # Pre‑load minimal GO DAG once
    GO_DAG = GODag("go-basic.obo", optional_attrs={"relationship"})

    def GO_enrich(gene_symbols: List[str], background: List[str], category: str) -> List[float]:
        """Return a list of p‑values for enriched terms given *gene_symbols*.

        *category* must be one of {BP, MF, CC}. Uses goatools' *fisher* test
        with Holm multiple‑testing correction.
        """
        # Example background mapping: every gene → {GO terms}
        # In practice you’d load a GAF/GOA file aligned to your organism.
        # Here we fake a dictionary of empty lists so goatools runs.
        gene2go: Dict[str, set[str]] = {g: set() for g in background}
        study = GOEnrichmentStudy(
            background,
            gene2go,
            GO_DAG,
            methods=["holm"],
            prt=None,
            alpha=0.05,
            significance_levels=None,
        )
        res = study.run_study(gene_symbols, methods=["holm"])
        pvals = [r.p_uncorrected for r in res]
        return pvals

except Exception:  # noqa: BLE001

    def GO_enrich(gene_symbols: List[str], background: List[str], category: str) -> List[float]:  # type: ignore[override]
        """Fallback stub – returns random p‑values so charts render."""
        rng = np.random.default_rng(seed=42)
        return list(rng.random(size=rng.integers(5, 30)))

# ----------------------------------------------------------------------------
# Algorithm wrappers (return List[Tuple[np.ndarray, np.ndarray]]) -------------
# ----------------------------------------------------------------------------


def _wrap_chen(X: np.ndarray, **kw) -> List[Tuple[np.ndarray, np.ndarray]]:
    bicls = run_chen_church(X, **kw)
    return [(np.array(r), np.array(c)) for r, c, *_ in bicls]


def _wrap_las(X: np.ndarray, **kw) -> List[Tuple[np.ndarray, np.ndarray]]:
    rows, cols = X.shape
    res = las_with_significance(X, rows, cols, **kw)
    if res["rows"] is None or res["cols"] is None:
        return []
    return [(np.array(res["rows"]), np.array(res["cols"]))]


def _wrap_isa(X: np.ndarray, **kw) -> List[Tuple[np.ndarray, np.ndarray]]:
    modules = ISA_multi_seed(X, **kw)
    out: List[Tuple[np.ndarray, np.ndarray]] = []
    for m in modules:
        r = np.flatnonzero(m["s_g"])
        c = np.flatnonzero(m["s_c"])
        if r.size and c.size:
            out.append((r, c))
    return out


def _wrap_bivisu(X: np.ndarray, **kw) -> List[Tuple[np.ndarray, np.ndarray]]:
    bics = bivisu(X, **kw)
    return [(np.array(sorted(r)), np.array(sorted(c))) for r, c in bics]

# Registry: display name → (callable, default‑kwargs)
ALGORITHMS: Dict[str, Tuple[Callable[..., Any], Dict[str, Any]]] = {
    "Chen & Church": (_wrap_chen, {"msr_threshold": 0.6, "alpha": 0.05, "max_biclusters": 10}),
    "LAS": (_wrap_las, {"max_iter": 100, "alpha": 0.05}),
    "ISA": (_wrap_isa, {"n_seeds": 500, "seed_size": 5, "t_g": 2.0, "t_c": 2.0, "max_iter": 100, "jaccard_thresh": 0.9}),
    "BiVisu": (_wrap_bivisu, {"model": "auto", "eps": 0.3, "thr": 0.05, "min_rows": 5, "min_cols": 2, "max_iters": 5}),
}

# ----------------------------------------------------------------------------
# Utility functions ----------------------------------------------------------
# ----------------------------------------------------------------------------

@st.cache_data(hash_funcs={np.ndarray: lambda _: None})
def _load_matrix(file, sep: str) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    df = pd.read_csv(file, sep=sep)
    # If first column is non‑numeric, treat as row labels (gene names)
    if not np.issubdtype(df.iloc[:, 0].dtype, np.number):
        row_labels = df.iloc[:, 0].astype(str).tolist()
        X = df.iloc[:, 1:].to_numpy(float)
    else:
        row_labels = [f"gene_{i}" for i in range(df.shape[0])]
        X = df.to_numpy(float)
    return df, X, row_labels


def _heatmap(X: np.ndarray, rows: np.ndarray, cols: np.ndarray, cmap: str = "viridis") -> None:
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.heatmap(X[np.ix_(rows, cols)], cmap=cmap, cbar=True, ax=ax, xticklabels=False, yticklabels=False)
    st.pyplot(fig)


def _plot_enrichment_matrix(data: Dict[str, List[float]], category: str):
    """Stacked bar chart – % enriched biclusters vs p‑value threshold."""
    fig, ax = plt.subplots(figsize=(6, 4))

    width = 0.9 / len(data)  # bars-per‑group
    x = np.arange(len(PV_THRESHOLDS))

    for i, (algo, pct) in enumerate(data.items()):
        ax.bar(x + i * width, pct, width=width, label=algo)

    ax.set_xticks(x + width * (len(data) - 1) / 2)
    ax.set_xticklabels([f"≤{t:.0e}" for t in PV_THRESHOLDS])
    ax.set_ylim(0, 100)
    ax.set_ylabel("% biclusters enriched")
    ax.set_title(f"GO {category} enrichment across algorithms")
    ax.legend()
    st.pyplot(fig)

# ----------------------------------------------------------------------------
# Streamlit main -------------------------------------------------------------
# ----------------------------------------------------------------------------

def main() -> None:
    st.set_page_config(page_title="BCTool", layout="wide", page_icon="🧬")
    st.title("🧬 BCTool – Interactive Biclustering & GO Quality Dashboard")

    # ---- Sidebar: dataset upload -----------------------------------------
    st.sidebar.header("1️⃣ Dataset")
    data_file = st.sidebar.file_uploader("CSV/TSV expression matrix", type=["csv", "tsv"])
    sep = st.sidebar.radio("Delimiter", [",", "\t"], horizontal=True)

    if data_file is None:
        st.info("⬅️ Upload a dataset to begin.")
        st.stop()

    df, X, gene_labels = _load_matrix(data_file, sep)
    st.write("### Data preview", df.head())
    st.markdown(f"Shape: **{X.shape[0]} genes × {X.shape[1]} conditions**")

    # ---- Sidebar: algorithm multi-select ---------------------------------
    st.sidebar.header("2️⃣ Algorithms")
    algo_choices = st.sidebar.multiselect("Select algorithms to run", list(ALGORITHMS.keys()), default=list(ALGORITHMS.keys())[:2])

    algo_params: Dict[str, Dict[str, Any]] = {}
    for name in algo_choices:
        st.sidebar.subheader(name)
        defaults = ALGORITHMS[name][1]
        user_vals: Dict[str, Any] = {}
        for p, default in defaults.items():
            if isinstance(default, bool):
                user_vals[p] = st.sidebar.checkbox(p, value=default, key=f"{name}_{p}")
            elif isinstance(default, int):
                user_vals[p] = st.sidebar.number_input(p, value=default, step=1, key=f"{name}_{p}")
            elif isinstance(default, float):
                user_vals[p] = st.sidebar.number_input(p, value=default, step=0.05, key=f"{name}_{p}")
            else:
                user_vals[p] = st.sidebar.text_input(p, value=str(default), key=f"{name}_{p}")
        algo_params[name] = user_vals

    # ---- Sidebar: GO settings ---------------------------------------------
    st.sidebar.header("3️⃣ GO assessment")
    go_cat = st.sidebar.selectbox("GO category", ["BP", "MF", "CC"], index=0)
    thresholds_selected = st.sidebar.multiselect(
        "P‑value thresholds", options=[str(t) for t in PV_THRESHOLDS], default=[str(t) for t in PV_THRESHOLDS]
    )
    chosen_thresholds = [float(t) for t in thresholds_selected]
    chosen_thresholds.sort()

    run_button = st.sidebar.button("🚀 Run biclustering + GO")
    if not run_button:
        st.stop()

    # ---- Execute algorithms ----------------------------------------------
    bicluster_results: Dict[str, List[Tuple[np.ndarray, np.ndarray]]] = {}
    for name in algo_choices:
        with st.spinner(f"Running {name} …"):
            func = ALGORITHMS[name][0]
            try:
                bicluster_results[name] = func(X, **algo_params[name])  # type: ignore[arg-type]
            except Exception as exc:  # noqa: BLE001
                st.error(f"{name} failed: {exc}")
                bicluster_results[name] = []
        st.success(f"{name}: {len(bicluster_results[name])} bicluster(s) found.")

    if all(len(bcls) == 0 for bcls in bicluster_results.values()):
        st.warning("No biclusters produced by any algorithm with current settings.")
        st.stop()

    # ---- GO enrichment evaluation ----------------------------------------
    enrichment_pct: Dict[str, List[float]] = {name: [0.0 for _ in PV_THRESHOLDS] for name in algo_choices}

    background = gene_labels  # simplistic; ideally supply organism‑wide list
    for name, bicls in bicluster_results.items():
        enriched_counts = [0 for _ in PV_THRESHOLDS]
        for rows, _ in bicls:
            genes = [gene_labels[i] for i in rows]
            pvals = GO_enrich(genes, background, go_cat)
            min_p = min(pvals) if pvals else 1.0
            for i, thr in enumerate(PV_THRESHOLDS):
                if min_p <= thr:
                    enriched_counts[i] += 1
        total = max(1, len(bicls))
        enrichment_pct[name] = [round(c / total * 100, 1) for c in enriched_counts]

    # ---- Dashboard layout -------------------------------------------------
    tab_overview, tab_compare, tab_biclusters = st.tabs(["GO overview", "Algorithm comparison", "Bicluster explorer"])

    with tab_overview:
        st.subheader(f"GO {go_cat} enrichment – per algorithm")
        for name in algo_choices:
            st.markdown(f"**{name}**")
            _plot_enrichment_matrix({name: enrichment_pct[name]}, go_cat)

    with tab_compare:
        st.subheader("Combined view")
        _plot_enrichment_matrix({n: enrichment_pct[n] for n in algo_choices}, go_cat)
        df_metrics = pd.DataFrame(enrichment_pct, index=[f"≤{t:.0e}" for t in PV_THRESHOLDS])
        st.write("### Numeric table", df_metrics)
        st.download_button("Download CSV", data=df_metrics.to_csv().encode(), mime="text/csv", file_name="go_enrichment_summary.csv")

    with tab_biclusters:
        st.subheader("Interactive bicluster explorer")
        for name, bicls in bicluster_results.items():
            if not bicls:
                continue
            st.markdown(f"### {name}")
            expander_base = st.expander(f"{name}: {len(bicls)} biclusters")
            with expander_base:
                for idx, (rows, cols) in enumerate(bicls, start=1):
                    exp = st.expander(f"Bicluster {idx} – {len(rows)}×{len(cols)}")
                    with exp:
                        _heatmap(X, rows, cols)
                        st.write("Row indices", rows)
                        st.write("Column indices", cols)

    # ---- External bicluster upload ---------------------------------------
    st.sidebar.header("4️⃣ External bicluster set (optional)")
    ext_file = st.sidebar.file_uploader("CSV/JSON of biclusters", type=["csv", "json"], key="ext")
    ext_label = st.sidebar.text_input("Label for external set", value="External")
    if ext_file is not None:
        try:
            if ext_file.name.endswith(".csv"):
                df_ext = pd.read_csv(ext_file)
                bicls_ext: List[Tuple[np.ndarray, np.ndarray]] = []
                for cl_id, grp in df_ext.groupby("cluster"):
                    rows = grp["row"].astype(int).unique()
                    cols = grp["col"].astype(int).unique()
                    bicls_ext.append((np.array(rows), np.array(cols)))
            else:
                payload = json.load(ext_file)
                bicls_ext = [(np.array(b["rows"]), np.array(b["cols"])) for b in payload]
            bicluster_results[ext_label] = bicls_ext
            algo_choices.append(ext_label)
            enrichment_pct[ext_label] = [0.0 for _ in PV_THRESHOLDS]
            # Simple enrichment calc for external set
            enriched_counts = [0 for _ in PV_THRESHOLDS]
            for rows, _ in bicls_ext:
                genes = [gene_labels[i] for i in rows if i < len(gene_labels)]
                pvals = GO_enrich(genes, background, go_cat)
                if pvals:
                    for i, thr in enumerate(PV_THRESHOLDS):
                        if min(pvals) <= thr:
                            enriched_counts[i] += 1
            total_ext = max(1, len(bicls_ext))
            enrichment_pct[ext_label] = [round(c / total_ext * 100, 1) for c in enriched_counts]
            st.sidebar.success(f"Loaded {len(bicls_ext)} external biclusters.")
        except Exception as e:  # noqa: BLE001
            st.sidebar.error(f"Failed to parse external set: {e}")

    st.sidebar.markdown("---")
    st.sidebar.caption("BCTool v0.2 · Streamlit")


if __name__ == "__main__":
    main()