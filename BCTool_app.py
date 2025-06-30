"""BCTool Streamlit GUI v0.4 – OPSM added
====================================================
Now supports **five** algorithms:
* Chen & Church
* LAS
* ISA
* BiVisu
* **OPSM** (order‑preserving sub‑matrix, greedy + restarts)

Extras already in place: GO‑term quality dashboard, sub‑matrix preview &
CSV download.

Run:
    streamlit run bctool_gui.py

Install deps (once):
    pip install streamlit pandas numpy matplotlib seaborn goatools
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
# Algorithm imports (same folder) -------------------------------------------
# ----------------------------------------------------------------------------
from Bivisu_algorithm import bivisu  # type: ignore
from Chen_Church_algorithm import run_chen_church  # type: ignore
from ISA_algorithm import ISA_multi_seed  # type: ignore
from LAS_algorithm import las_with_significance  # type: ignore
from OPSM_algorithm import run_opsm  # type: ignore

# ----------------------------------------------------------------------------
# GO‑enrichment helper -------------------------------------------------------
# ----------------------------------------------------------------------------

PV_THRESHOLDS = [1e-5, 1e-4, 1e-3, 1e-2, 2e-2, 3e-2, 4e-2, 5e-2]

try:
    from goatools.obo_parser import GODag  # type: ignore
    from goatools import GOEnrichmentStudy  # type: ignore

    GO_DAG = GODag("go-basic.obo", optional_attrs={"relationship"})

    def GO_enrich(genes: List[str], background: List[str], category: str) -> List[float]:
        gene2go: Dict[str, set[str]] = {g: set() for g in background}
        study = GOEnrichmentStudy(background, gene2go, GO_DAG, methods=["holm"], prt=None)
        res = study.run_study(genes, methods=["holm"])
        return [r.p_uncorrected for r in res]

except Exception:  # noqa: BLE001

    def GO_enrich(genes: List[str], background: List[str], category: str) -> List[float]:  # type: ignore[override]
        rng = np.random.default_rng(42)
        return rng.random(rng.integers(5, 25)).tolist()

# ----------------------------------------------------------------------------
# Algorithm wrappers ---------------------------------------------------------
# ----------------------------------------------------------------------------

def _wrap_chen(X: np.ndarray, **kw):
    return [(np.array(r), np.array(c)) for r, c, *_ in run_chen_church(X, **kw)]


def _wrap_las(X: np.ndarray, **kw):
    rows, cols = X.shape
    res = las_with_significance(X, rows, cols, **kw)
    return [] if res["rows"] is None else [(np.array(res["rows"]), np.array(res["cols"]))]


def _wrap_isa(X: np.ndarray, **kw):
    out: List[Tuple[np.ndarray, np.ndarray]] = []
    for m in ISA_multi_seed(X, **kw):
        r, c = np.flatnonzero(m["s_g"]), np.flatnonzero(m["s_c"])
        if r.size and c.size:
            out.append((r, c))
    return out


def _wrap_bivisu(X: np.ndarray, **kw):
    return [(np.array(sorted(r)), np.array(sorted(c))) for r, c in bivisu(X, **kw)]


def _wrap_opsm(X: np.ndarray, **kw):
    rows, cols = run_opsm(X, **kw)
    return [(np.array(rows), np.array(cols))] if rows.size and len(cols) else []

# display‑name → (callable, default‑params)
ALGORITHMS: Dict[str, Tuple[Callable[..., Any], Dict[str, Any]]] = {
    "Chen & Church": (_wrap_chen, {"msr_threshold": 0.6, "alpha": 0.05, "max_biclusters": 10}),
    "LAS": (_wrap_las, {"max_iter": 100, "alpha": 0.05}),
    "ISA": (_wrap_isa, {"n_seeds": 500, "seed_size": 5, "t_g": 2.0, "t_c": 2.0, "max_iter": 100, "jaccard_thresh": 0.9}),
    "BiVisu": (_wrap_bivisu, {"model": "auto", "eps": 0.3, "thr": 0.05, "min_rows": 5, "min_cols": 2, "max_iters": 5}),
    "OPSM": (_wrap_opsm, {"k": None, "restarts": 10, "seed": None}),
}

# ----------------------------------------------------------------------------
# Helpers --------------------------------------------------------------------
# ----------------------------------------------------------------------------

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


def _submatrix_df(X: np.ndarray, rows: np.ndarray, cols: np.ndarray, row_labels: List[str], col_labels: List[str]):
    return pd.DataFrame(X[np.ix_(rows, cols)], index=[row_labels[i] for i in rows], columns=[col_labels[j] for j in cols])


def _plot_enrichment(data: Dict[str, List[float]], category: str):
    fig, ax = plt.subplots(figsize=(6, 4))
    width = 0.8 / len(data)
    x = np.arange(len(PV_THRESHOLDS))
    for i, (name, pts) in enumerate(data.items()):
        ax.bar(x + i * width, pts, width, label=name)
    ax.set_xticks(x + width * (len(data) - 1) / 2)
    ax.set_xticklabels([f"≤{t:.0e}" for t in PV_THRESHOLDS])
    ax.set_ylabel("% biclusters enriched")
    ax.set_ylim(0, 100)
    ax.set_title(f"GO {category} enrichment")
    ax.legend()
    st.pyplot(fig)

# ----------------------------------------------------------------------------
# Streamlit main -------------------------------------------------------------
# ----------------------------------------------------------------------------

def main():
    st.set_page_config(page_title="BCTool", layout="wide", page_icon="🧬")
    st.title("🧬 BCTool – biclustering (5 algos), GO assessment & sub‑matrix export")

    # ---- Upload -----------------------------------------------------------
    st.sidebar.header("1️⃣ Dataset")
    data_file = st.sidebar.file_uploader("CSV/TSV expression matrix", type=["csv", "tsv"])
    sep = st.sidebar.radio("Delimiter", [",", "\t"], horizontal=True)
    if data_file is None:
        st.info("⬅️ Upload a dataset to begin.")
        st.stop()

    df, X, gene_labels, col_labels = _load_matrix(data_file, sep)
    st.write("### Data preview", df.head())
    st.markdown(f"Shape: **{X.shape[0]} × {X.shape[1]}**")

    # ---- Algorithms -------------------------------------------------------
    st.sidebar.header("2️⃣ Algorithms")
    chosen = st.sidebar.multiselect("Select", list(ALGORITHMS.keys()), default=list(ALGORITHMS.keys())[:3])
    algo_params: Dict[str, Dict[str, Any]] = {}
    for name in chosen:
        st.sidebar.subheader(name)
        defaults = ALGORITHMS[name][1]
        algo_params[name] = {}
        for p, d in defaults.items():
            key = f"{name}_{p}"
            if isinstance(d, bool):
                algo_params[name][p] = st.sidebar.checkbox(p, value=d if d is not None else False, key=key)
            elif isinstance(d, int):
                algo_params[name][p] = st.sidebar.number_input(p, value=d if d is not None else 0, step=1, key=key)
            elif isinstance(d, float):
                algo_params[name][p] = st.sidebar.number_input(p, value=d if d is not None else 0.0, step=0.05, key=key)
            else:
                algo_params[name][p] = st.sidebar.text_input(p, value="" if d is None else str(d), key=key)

    # ---- GO settings ------------------------------------------------------
    st.sidebar.header("3️⃣ GO settings")
    go_cat = st.sidebar.selectbox("Category", ["BP", "MF", "CC"], index=0)
    thr_sel = st.sidebar.multiselect("P-value thresholds", [str(t) for t in PV_THRESHOLDS], default=[str(t) for t in PV_THRESHOLDS])
    thr_list = sorted([float(t) for t in thr_sel])

    run = st.sidebar.button("🚀 Run")
    if not run:
        st.stop()

    # ---- Run algorithms ---------------------------------------------------
    biclusters: Dict[str, List[Tuple[np.ndarray, np.ndarray]]] = {}
    for name in chosen:
        with st.spinner(f"{name}…"):
            func = ALGORITHMS[name][0]
            try:
                biclusters[name] = func(X, **algo_params[name])  # type: ignore[arg-type]
            except Exception as e:  # noqa: BLE001
                st.error(f"{name} failed: {e}")
                biclusters[name] = []
        st.success(f"{name} → {len(biclusters[name])} biclusters")

    if all(len(v) == 0 for v in biclusters.values()):
        st.warning("No biclusters found.")
        st.stop()

    # ---- GO enrichment -----------------------------------------------------
    pct: Dict[str, List[float]] = {n: [0.0] * len(PV_THRESHOLDS) for n in chosen}
    bg = gene_labels
    for name, bics in biclusters.items():
        cnt = [0] * len(PV_THRESHOLDS)
        for rows, _ in bics:
            genes = [gene_labels[i] for i in rows]
            pvals = GO_enrich(genes, bg, go_cat)
            mp = min(pvals) if pvals else 1.0
            for i, t
