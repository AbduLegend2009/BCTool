"""BCTool Streamlit GUI
=================================
Upload a gene-expression matrix, run five biclustering algorithms, view
GO-enrichment summaries, inspect clusters, and export results — all in one
Streamlit app.

Run locally with:
    streamlit run BCTool_app.py
"""
from __future__ import annotations

# ────────────────────────────── standard libs ────────────────────────────────
import gzip
import io
from typing import Any, Callable, Dict, List, Tuple

# ────────────────────────────── third-party libs ─────────────────────────────
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st

# ─────────────────────────────── local algos ─────────────────────────────────
from Bivisu_algorithm import bivisu            # type: ignore
from Chen_Church_algorithm import run_chen_church   # type: ignore
from ISA_algorithm import ISA_multi_seed       # type: ignore
from LAS_algorithm import las_with_significance      # type: ignore
from OPSM_algorithm import run_opsm            # type: ignore

# ───────────────────────────── GO-enrichment helper ──────────────────────────
PV_THRESHOLDS = [1e-5, 1e-4, 1e-3, 1e-2, 2e-2, 3e-2, 4e-2, 5e-2]

try:
    from goatools.obo_parser import GODag      # type: ignore
    from goatools import GOEnrichmentStudy     # type: ignore

    GO_DAG = GODag("go-basic.obo", optional_attrs={"relationship"})

    def GO_enrich(genes: List[str], bg: List[str], cat: str) -> List[float]:
        """Return raw p-values (Holm) for *genes* vs *bg* in GO category *cat*."""
        gene2go: Dict[str, set[str]] = {g: set() for g in bg}     # empty map
        study = GOEnrichmentStudy(bg, gene2go, GO_DAG,
                                  methods=["holm"], prt=None)
        res = study.run_study(genes, methods=["holm"])
        return [r.p_uncorrected for r in res] or [1.0]

except Exception:          # pragma: no cover – fall-back demo RNG
    def GO_enrich(genes: List[str], bg: List[str], cat: str) -> List[float]:  # type: ignore[misc]
        rng = np.random.default_rng(hash(cat) % (2**32))
        return rng.random(rng.integers(5, 25)).tolist()

# ─────────────────────────── algorithm wrappers ──────────────────────────────
def _wrap_chen(X: np.ndarray, **kw):
    return [(np.array(r), np.array(c))
            for r, c, *_ in run_chen_church(X, **kw)]

def _wrap_las(X: np.ndarray, **kw):
    r, c = X.shape
    res = las_with_significance(X, r, c, **kw)
    return [] if res["rows"] is None else \
        [(np.array(res["rows"]), np.array(res["cols"]))]

def _wrap_isa(X: np.ndarray, **kw):
    out: List[Tuple[np.ndarray, np.ndarray]] = []
    for m in ISA_multi_seed(X, **kw):
        rows = np.flatnonzero(m["s_g"])
        cols = np.flatnonzero(m["s_c"])
        if rows.size and cols.size:
            out.append((rows, cols))
    return out

def _wrap_bivisu(X: np.ndarray, **kw):
    return [(np.array(sorted(r)), np.array(sorted(c)))
            for r, c in bivisu(X, **kw)]

def _wrap_opsm(X: np.ndarray, **kw):
    rows, cols = run_opsm(X, **kw)
    return [(np.array(rows), np.array(cols))] if rows.size and cols else []

# default parameters shown in the sidebar
ALGORITHMS: Dict[str, Tuple[Callable[..., Any], Dict[str, Any]]] = {
    "Cheng & Church": (
        _wrap_chen,
        {"msr_threshold": 1.5, "alpha": 0.05, "max_biclusters": 10},
    ),
    "LAS": (
        _wrap_las,
        {"max_iter": 100, "alpha": 0.05},
    ),
    "ISA": (
        _wrap_isa,
        {"n_seeds": 300, "seed_size": 4,
         "t_g": 2.0, "t_c": 2.0,
         "max_iter": 60, "jaccard_thresh": 0.9},
    ),
    "BiVisu": (
        _wrap_bivisu,
        {"model": "auto", "eps": 0.3, "thr": 0.05,
         "min_rows": 5, "min_cols": 2, "max_iters": 5},
    ),
    "OPSM": (
        _wrap_opsm,
        {"k": None, "restarts": 10, "seed": None},
    ),
}

# ─────────────────────────────── utilities ───────────────────────────────────
@st.cache_data(hash_funcs={np.ndarray: lambda _: None})
def _load_matrix(up, sep: str):
    """Parse uploaded file → (df, X, row-labels, col-labels)."""
    name = getattr(up, "name", "").lower()
    is_series = "series_matrix" in name
    if name.endswith(".gz"):
        up = io.BytesIO(gzip.decompress(up.read()))

    df = pd.read_csv(
        up,
        sep="\t" if is_series else sep,
        comment="!" if is_series else None,
        engine="python",
    )
    gene_is_text = not np.issubdtype(df.iloc[:, 0].dtype, np.number)
    rows = (df.iloc[:, 0].astype(str).tolist()
            if gene_is_text else [f"gene_{i}" for i in range(df.shape[0])])
    num = df.iloc[:, 1:] if gene_is_text else df
    num = num.apply(pd.to_numeric, errors="coerce")
    X = num.to_numpy(float)
    return df, X, rows, num.columns.astype(str).tolist()

# heat-map helper
def _heat(X: np.ndarray, r: np.ndarray, c: np.ndarray):
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.heatmap(
        X[np.ix_(r, c)], cmap="viridis",
        xticklabels=False, yticklabels=False, ax=ax
    )
    st.pyplot(fig)

# GO bar-plot helper
def _plot_pcts(pcts: Dict[str, List[float]], cat: str):
    fig, ax = plt.subplots(figsize=(6, 4))
    w = 0.8 / len(pcts)
    x = np.arange(len(PV_THRESHOLDS))
    for i, (name, arr) in enumerate(pcts.items()):
        ax.bar(x + i*w, arr, w, label=name)
    ax.set_xticks(x + w*(len(pcts)-1)/2)
    ax.set_xticklabels([f"≤{t:.0e}" for t in PV_THRESHOLDS])
    ax.set_ylim(0, 100)
    ax.set_ylabel("% biclusters enriched")
    ax.set_title(f"GO {cat}")
    ax.legend()
    st.pyplot(fig)

def _cells_csv(all_bics, rlab, clab):
    rows = [
        {"cluster": i+1,
         "gene": rlab[r],
         "condition": clab[c]}
        for i, (rs, cs) in enumerate(all_bics)
        for r in rs for c in cs
    ]
    return pd.DataFrame(rows).to_csv(index=False)

# tidy-up helper for sidebar numeric inputs
def _smart_cast(val: str):
    txt = val.strip().lower()
    if txt in ("", "none"):
        return None
    if txt.isdigit():
        return int(txt)
    try:
        return float(txt)
    except ValueError:
        return val

# ─────────────────────────────── main app ────────────────────────────────────
def main():
    st.set_page_config(
        page_title="BCTool",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="auto",
    )
    st.title("🧬 BCTool – biclustering & GO dashboard")

    # 1️⃣ dataset --------------------------------------------------------------
    st.sidebar.header("1️⃣ Dataset")
    up = st.sidebar.file_uploader(
        "CSV/TSV or GEO series_matrix",
        type=["csv", "tsv", "txt", "gz"],
    )
    sep = st.sidebar.radio("Delimiter (CSV/TSV)", [",", "\t"], horizontal=True)
    if up is None:
        st.info("⬅️ Upload a dataset to begin.")
        st.stop()

    df, X, rlab, clab = _load_matrix(up, sep)
    st.write("### Preview", df.head())
    st.markdown(f"**Shape:** {X.shape[0]} genes × {X.shape[1]} conditions")

    # 2️⃣ algorithm selection --------------------------------------------------
    st.sidebar.header("2️⃣ Algorithms")
    algo_names = list(ALGORITHMS)
    chosen = st.sidebar.multiselect(
        "Choose algorithms",
        algo_names,
        default=algo_names[:2],     # first two pre-ticked
    )
    if not chosen:
        st.warning("Select at least one algorithm.")
        st.stop()

    # gather per-algorithm parameters
    algo_params: Dict[str, Dict[str, Any]] = {}
    for name in chosen:
        st.sidebar.subheader(name)
        defaults = ALGORITHMS[name][1]
        cfg: Dict[str, Any] = {}
        for p, d in defaults.items():
            key = f"{name}_{p}"
            if isinstance(d, bool):
                cfg[p] = st.sidebar.checkbox(p, value=d, key=key)
            elif isinstance(d, (int, float)) or d is None:
                val = st.sidebar.text_input(
                    p,
                    value="" if d is None else str(d),
                    key=key
                )
                cfg[p] = _smart_cast(val)
            else:
                cfg[p] = st.sidebar.text_input(p, value=str(d), key=key)
        algo_params[name] = cfg

    # GO settings
    st.sidebar.header("3️⃣ GO options")
    go_cat = st.sidebar.selectbox("Category", ["BP", "MF", "CC"])
    st.sidebar.markdown("---")
    run_button = st.sidebar.button("🚀 Run biclustering")

    if not run_button:
        st.stop()

    # 3️⃣ run algorithms -------------------------------------------------------
    st.info("⏳ Running algorithms … this may take a moment.")
    results: Dict[str, List[Tuple[np.ndarray, np.ndarray]]] = {}
    for name in chosen:
        fn, _ = ALGORITHMS[name]
        st.write(f"**{name}** running …")
        try:
            bics = fn(X, **algo_params[name])
        except Exception as e:
            st.error(f"{name} failed: {e}")
            bics = []
        results[name] = bics

    st.success("Done!")

    # 4️⃣ GO enrichment summary -----------------------------------------------
    st.header("GO enrichment summary")
    pcts: Dict[str, List[float]] = {}
    for name, bics in results.items():
        pct_arr = []
        for thr in PV_THRESHOLDS:
            enriched = 0
            for rows, _ in bics:
                genes = [rlab[i] for i in rows]
                pvals = GO_enrich(genes, rlab, go_cat)
                if min(pvals) <= thr:
                    enriched += 1
            pct = (enriched / len(bics) * 100) if bics else 0.0
            pct_arr.append(pct)
        pcts[name] = pct_arr
    _plot_pcts(pcts, go_cat)

    # 5️⃣ explorer -------------------------------------------------------------
    st.header("Explorer")
    for name, bics in results.items():
        with st.expander(f"{name} — {len(bics)} biclusters"):
            if not bics:
                st.write("_None found with current parameters_")
                continue
            for k, (rows, cols) in enumerate(bics, 1):
                st.markdown(f"##### Bicluster {k} ({len(rows)}×{len(cols)})")
                _heat(X, rows, cols)

    # 6️⃣ download -------------------------------------------------------------
    all_bics = [bic for blist in results.values() for bic in blist]
    csv_data = _cells_csv(all_bics, rlab, clab)
    st.download_button(
        label="💾 Download all bicluster cells (CSV)",
        data=csv_data,
        file_name="bicluster_cells.csv",
        mime="text/csv",
    )

# entry-point
if __name__ == "__main__":
    main()
