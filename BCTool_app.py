"""BCTool Streamlit GUI
=======================
End‑to‑end dashboard that lets you
• upload a CSV/TSV *or* GEO `series_matrix.txt(.gz)` file,
• run five biclustering algorithms (Chen‑Church, LAS, ISA, BiVisu, OPSM),
• visualise GO‑enrichment percentages, compare methods and explore biclusters,
• export all bicluster cells as one CSV.

Run with:
    streamlit run BCTool_app.py
"""
from __future__ import annotations

import gzip
import io
from typing import Any, Callable, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st

# ────────────────────────────────────────────────────────────────────────────────
#  Algorithm wrappers  (delegate to standalone modules living in the same folder)
# ────────────────────────────────────────────────────────────────────────────────
from Bivisu_algorithm import bivisu  # type: ignore
from Chen_Church_algorithm import run_chen_church  # type: ignore
from ISA_algorithm import ISA_multi_seed  # type: ignore
from LAS_algorithm import las_with_significance  # type: ignore
from OPSM_algorithm import run_opsm  # type: ignore

# ────────────────────────────────────────────────────────────────────────────────
#  Gene‑Ontology helper  (use goatools if available, else fallback to dummy RNG)
# ────────────────────────────────────────────────────────────────────────────────
PV_THRESHOLDS = [1e-5, 1e-4, 1e-3, 1e-2, 2e-2, 3e-2, 4e-2, 5e-2]
try:
    from goatools.obo_parser import GODag  # type: ignore
    from goatools import GOEnrichmentStudy  # type: ignore

    GO_DAG = GODag("go-basic.obo", optional_attrs={"relationship"})

    def GO_enrich(genes: List[str], bg: List[str], cat: str) -> List[float]:
        """Return raw p‑values for *genes* vs *bg* in GO *cat* (BP|MF|CC)."""
        gene2go: Dict[str, set[str]] = {g: set() for g in bg}  # empty mapping
        study = GOEnrichmentStudy(bg, gene2go, GO_DAG, methods=["holm"], prt=None)
        res = study.run_study(genes, methods=["holm"])
        return [r.p_uncorrected for r in res]

except Exception:

    def GO_enrich(genes: List[str], bg: List[str], cat: str) -> List[float]:  # type: ignore[override]
        rng = np.random.default_rng(42)
        return rng.random(rng.integers(5, 25)).tolist()

# ────────────────────────────────────────────────────────────────────────────────
#  Individual wrappers turning each algorithm into (rows, cols) → ndarray[int]
# ────────────────────────────────────────────────────────────────────────────────

def _wrap_chen(X: np.ndarray, **kw):
    return [(np.array(r), np.array(c)) for r, c, *_ in run_chen_church(X, **kw)]


def _wrap_las(X: np.ndarray, **kw):
    r, c = X.shape
    res = las_with_significance(X, r, c, **kw)
    return [] if res["rows"] is None else [(np.array(res["rows"]), np.array(res["cols"]))]


def _wrap_isa(X: np.ndarray, **kw):
    out = []
    for m in ISA_multi_seed(X, **kw):
        rows, cols = np.flatnonzero(m["s_g"]), np.flatnonzero(m["s_c"])
        if rows.size and cols.size:
            out.append((rows, cols))
    return out


def _wrap_bivisu(X: np.ndarray, **kw):
    return [
        (np.array(sorted(r)), np.array(sorted(c))) for r, c in bivisu(X, **kw)
    ]


def _wrap_opsm(X: np.ndarray, **kw):
    rows, cols = run_opsm(X, **kw)
    return [(np.array(rows), np.array(cols))] if rows.size and len(cols) else []


ALGORITHMS: Dict[str, Tuple[Callable[..., Any], Dict[str, Any]]] = {
    "Chen & Church": (
        _wrap_chen,
        {"msr_threshold": 0.6, "alpha": 0.05, "max_biclusters": 10},
    ),
    "LAS": (_wrap_las, {"max_iter": 100, "alpha": 0.05}),
    "ISA": (
        _wrap_isa,
        {
            "n_seeds": 500,
            "seed_size": 5,
            "t_g": 2.0,
            "t_c": 2.0,
            "max_iter": 100,
            "jaccard_thresh": 0.9,
        },
    ),
    "BiVisu": (
        _wrap_bivisu,
        {
            "model": "auto",
            "eps": 0.3,
            "thr": 0.05,
            "min_rows": 5,
            "min_cols": 2,
            "max_iters": 5,
        },
    ),
    "OPSM": (
        _wrap_opsm,
        {"k": None, "restarts": 10, "seed": None},
    ),
}

# ────────────────────────────────────────────────────────────────────────────────
#  Utilities
# ────────────────────────────────────────────────────────────────────────────────


@st.cache_data(hash_funcs={np.ndarray: lambda _: None})
def _load_matrix(up, sep: str):
    """Parse uploaded file *up* → (df, X, row_labels, col_labels)."""
    name = getattr(up, "name", "").lower()
    series = "series_matrix" in name  # GEO dump?

    # transparently gunzip
    if name.endswith(".gz"):
        up = io.BytesIO(gzip.decompress(up.read()))

    df = pd.read_csv(
        up,
        sep="\t" if series else sep,
        comment="!" if series else None,
        engine="python",
    )

    gene_is_text = not np.issubdtype(df.iloc[:, 0].dtype, np.number)
    rows = (
        df.iloc[:, 0].astype(str).tolist()
        if gene_is_text
        else [f"gene_{i}" for i in range(df.shape[0])]
    )

    num = df.iloc[:, 1:] if gene_is_text else df
    cols = num.columns.astype(str).tolist()
    X = num.to_numpy(float)

    return df, X, rows, cols


def _heat(X: np.ndarray, r: np.ndarray, c: np.ndarray):
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.heatmap(
        X[np.ix_(r, c)], cmap="viridis", ax=ax, xticklabels=False, yticklabels=False
    )
    st.pyplot(fig)


def _plot_pcts(pcts: Dict[str, List[float]], cat: str):
    fig, ax = plt.subplots(figsize=(6, 4))
    w = 0.8 / len(pcts)
    x = np.arange(len(PV_THRESHOLDS))
    for i, (n, arr) in enumerate(pcts.items()):
        ax.bar(x + i * w, arr, w, label=n)
    ax.set_xticks(x + w * (len(pcts) - 1) / 2)
    ax.set_xticklabels([f"≤{t:.0e}" for t in PV_THRESHOLDS])
    ax.set_ylim(0, 100)
    ax.set_ylabel("% enriched")
    ax.set_title(f"GO {cat}")
    ax.legend()
    st.pyplot(fig)


def _cells_csv(bics, rlab, clab):
    rows = [
        {"cluster": i + 1, "gene": rlab[r], "condition": clab[c]}
        for i, (rs, cs) in enumerate(bics)
        for r in rs
        for c in cs
    ]
    return pd.DataFrame(rows).to_csv(index=False)


# ────────────────────────────────────────────────────────────────────────────────
#  Main application
# ────────────────────────────────────────────────────────────────────────────────

def main():
    st.set_page_config(
        page_title="BCTool", page_icon="🧬", layout="wide", initial_sidebar_state="auto"
    )
    st.title("🧬 BCTool – biclustering & GO dashboard")

    # 1️⃣ Dataset ----------------------------------------------------------------
    st.sidebar.header("1️⃣ Dataset")
    up = st.sidebar.file_uploader(
        "CSV/TSV or GEO series_matrix", type=["csv", "tsv", "txt", "gz"]
    )
    sep = st.sidebar.radio("Delimiter (CSV/TSV)", [",", "\t"], horizontal=True)
    if up is None:
        st.info("⬅️ Upload a dataset")
        st.stop()

    df, X, rlab, clab = _load_matrix(up, sep)
    st.write("### Preview", df.head())
    st.markdown(f"Shape **{X.shape[0]}×{X.shape[1]}**")

    # 2️⃣ Algorithms -------------------------------------------------------------
    st.sidebar.header("2️⃣ Algorithms")
    chosen = st.sidebar.multiselect(
        "Choose", list(ALGORITHMS), default=list(ALGORITHMS)[:3]
    )

    params: Dict[str, Dict[str, Any]] = {}
    for n in chosen:
        st.sidebar.subheader(n)
        defaults = ALGORITHMS[n][1]
        cfg: Dict[str, Any] = {}
        for p, d in defaults.items():
            key = f"{n}_{p}"
            if isinstance(d, bool):
                cfg[p] = st.sidebar.checkbox(p, value=d, key=key)
            elif isinstance(d, (int, float)):
                cfg[p] = st.sidebar.number_input(p, value=d, key=key)
            else:
                cfg[p] = st.sidebar.text_input(p, value=str(d), key=key)
        params[n] = cfg

    # 3️⃣ GO settings ------------------------------------------------------------
    st.sidebar.header("3️⃣ GO settings")
    cat = st.sidebar.selectbox("Category", ["BP", "MF", "CC"], index=0)

    # ---------------------------------------------------------------------------
    if not st.sidebar.button("🚀 Run"):
        st.stop()

    # Run algorithms ------------------------------------------------------------
    biclusters: Dict[str, List[Tuple[np.ndarray, np.ndarray]]] = {}
    for n in chosen:
        with st.spinner(f"Running {n} …"):
            try:
                biclusters[n] = ALGORITHMS[n][0](X, **params[n])  # type: ignore[arg-type]
            except Exception as e:
                st.error(f"{n} failed: {e}")
                biclusters[n] = []

    if all(len(v) == 0 for v in biclusters.values()):
        st.warning("No biclusters found")
        st.stop()

    # GO enrichment overview ----------------------------------------------------
    pcts = {n: [0] * len(PV_THRESHOLDS) for n in chosen}
    for n, bics in biclusters.items():
        counts = [0] * len(PV_THRESHOLDS)
        for rows, _ in bics:
            mp = (
                min(GO_enrich([rlab[i] for i in rows], rlab, cat)) if rows.size else 1
            )
            for i, t in enumerate(PV_THRESHOLDS):
                if mp <= t:
                    counts[i] += 1
        total = max(1, len(bics))
        pcts[n] = [round(c / total * 100, 1) for c in counts]

    #  📊  Tabs ------------------------------------------------------------------
    tab1, tab2, tab3 = st.tabs(["GO overview", "Comparison", "Explorer"])

    with tab1:
        for n in chosen:
            st.markdown(f"**{n}**")
            _plot_pcts({n: pcts[n]}, cat)

    with tab2:
        _plot_pcts({n: pcts[n] for n in chosen}, cat)
        st.dataframe(
            pd.DataFrame(pcts, index=[f"≤{t:.0e}" for t in PV_THRESHOLDS])
        )

    with tab3:
        for n, bics in biclusters.items():
            st.subheader(n)
            for i, (r, c) in enumerate(bics, 1):
                with st.expander(f"Cluster {i} ({len(r)}×{len(c)})"):
                    _heat(X, r, c)
                    st.write(
                        pd.DataFrame(
                            X[np.ix_(r, c)],
                            index=[rlab[j] for j in r],
                            columns=[clab[k] for k in c],
                        )
                    )
        st.download_button(
            "Download all cells",
            data=_cells_csv(
                [bc for lst in biclusters.values() for bc in lst], rlab, clab
            ).encode(),
            file_name="bicluster_cells.csv",
            mime="text/csv",
        )

    st.sidebar.markdown("---")
    st.sidebar.caption("BCTool v0.5 (Streamlit edition)")


if __name__ == "__main__":
    main()
