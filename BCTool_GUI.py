import streamlit as st
import pandas as pd
import numpy as np
from GO_assessment import go_assessment
import adapter

def uploaded_data(data):
    if data is None:
        return None, None
    suffix = data.name.split(".")[-1].lower()
    if suffix in {"csv", "tsv"}:
        df = pd.read_csv(data, sep="\t" if suffix == "tsv" else ",", index_col=0)
    elif suffix in {"xls", "xlsx"}:
        df = pd.read_excel(data, engine="openpyxl", index_col=0)
    else:
        st.error(f"unsupported data type: {suffix}")
        return None, None
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.isna().any().any():
        st.warning("Non-numeric values found; coerced to NaN and dropped.")
        df = df.dropna(axis=0, how="any")
    gene_ids = df.index.tolist()        # for GO assessment
    matrix   = df.to_numpy(dtype=float) # for the biclustering algos
    return matrix, gene_ids

algorithms = ["LAS", "Chen and Church", "ISA", "OPSM", "Bivisu"]


def summarize_biclusters(bics, gene_ids, alg_name=None):
    """Create a summary ``DataFrame`` from biclustering results.

    Parameters
    ----------
    bics : list | dict
        Either a list of ``Bicluster`` objects from a single algorithm or a
        dictionary mapping algorithm names to such lists.
    gene_ids : list[str]
        Mapping of row indices to gene identifiers.
    alg_name : str | None
        Optional algorithm name when ``bics`` is a plain list.
    """

    rows = []

    def _append_rows(bic_list, alg):
        for i, bic in enumerate(bic_list):
            rows.append({
                "ID": f"{alg}_{i}",
                "Algorithm": alg,
                "Rows": len(bic.rows),
                "Cols": len(bic.cols),
                "Genes": ", ".join(str(gene_ids[j]) for j in bic.rows[:3]) + " …",
                "Score": getattr(bic, "score", np.nan),
            })

    if isinstance(bics, dict):
        for alg, bic_list in bics.items():
            _append_rows(bic_list, alg)
    else:
        _append_rows(bics, alg_name or "")

    return pd.DataFrame(rows)



def main():
    st.title("BCTool 🧬")
    data = st.file_uploader(
        "Please upload expression matrix (CSV, TSV, or Excel)",
        type=["csv", "tsv", "xls", "xlsx"],
    )
    with st.sidebar:
        st.header("Customization")
        matrix, gene_ids = uploaded_data(data)
        taxid = st.number_input("What is the taxonomy identifier?")
        sel_alg = st.multiselect("Please select algorithms", algorithms)
        if "LAS" in sel_alg:
            st.header("LAS")
            alpha_LAS = st.number_input("What is the p-value cutoff for biclusters?")
            max_iter_LAS = st.number_input("What is the maximum amount of iterations?")
            k_rows_LAS = st.number_input(
                "How many rows in candidate submatrices?", value=10, step=1, min_value=1
            )
            k_cols_LAS = st.number_input(
                "How many columns in candidate submatrices?", value=20, step=1, min_value=1
            )
        if "Chen and Church" in sel_alg:
            st.header("Chen and Church")
            delta_C_C = st.number_input("What is the MSR threshold for biclusters?")
            alpha_C_C = st.number_input("What fraction of the cells are we allowed to prune?")
            max_bi_C_C = st.number_input("What is the maximum number of biclusters?")
        if "ISA" in sel_alg:
            st.header("ISA")
            n_seeds_ISA = st.number_input("How many starting seeds?")
            seed_size_ISA = st.number_input("How many starting conditions in each seed?")
            t_g_ISA = st.number_input("What is the gene z-score threshold?")
            t_c_ISA = st.number_input("What is the condition z-score threshold?")
        if "OPSM" in sel_alg:
            st.header("OPSM")
            k_OPSM = st.number_input("What is the maximum amount of conditions in each bicluster?")
            restarts_OPSM = st.number_input("How many restarts?")
        if "Bivisu" in sel_alg:
            st.header("Bivisu")
            model_Bivisu = st.selectbox("Which model do you want to be used?", ["add", "mult", "auto"], index = 2)
            eps_Bivisu = st.number_input("What is the bin width when quantasising the signature algorithms?")
            msr_Bivisu = st.number_input("What is the MSR cutoff?")
            min_genes_Bivisu = st.number_input("What is the lowest number of genes per bicluster?")
            min_cond_Bivisu = st.number_input("What is the lowest number of conditions per bicluster?")
        if st.button("Run algorithms 🚀"):
            if not sel_alg:
                st.error("No algorithms found. Please ensure that algorithms are available.")
                st.stop()
            s = st.session_state.get("Biclusters", {})

            if isinstance(s, list):
                # backward compatibility if previous state was stored as a list
                s = {name: bic_list for name, bic_list in s}

            for _ in sel_alg:
                if _ == "LAS":
                    LAS = adapter.wrap_las(
                        matrix,
                        max_iter=int(max_iter_LAS),
                        alpha=alpha_LAS,
                        k_rows=int(k_rows_LAS),
                        k_cols=int(k_cols_LAS),
                    )
                    s["LAS"] = LAS
                elif _ == "Chen and Church":
                    Chen_and_Church = adapter.wrap_church(
                        matrix,
                        delta=delta_C_C,
                        alpha=alpha_C_C,
                        max_biclusters=int(max_bi_C_C),
                    )
                    s["Chen and Church"] = Chen_and_Church
                elif _ == "ISA":
                    ISA = adapter.wrap_isa(
                        matrix,
                        n_seeds=int(n_seeds_ISA),
                        seed_size=int(seed_size_ISA),
                        t_g=t_g_ISA,
                        t_c=t_c_ISA,
                    )
                    s["ISA"] = ISA
                elif _ == "OPSM":
                    opsm = adapter.wrap_opsm(
                        matrix,
                        k=int(k_OPSM) if k_OPSM else None,
                        restarts=int(restarts_OPSM),
                    )
                    s["OPSM"] = opsm
                elif _ == "Bivisu":
                    Bi = adapter.wrap_bivisu(
                        matrix,
                        model=model_Bivisu,
                        eps=eps_Bivisu,
                        thr=msr_Bivisu,
                        min_rows=int(min_genes_Bivisu),
                        min_cols=int(min_cond_Bivisu),
                    )

                    s["Bivisu"] = Bi

            # Persist biclusters for later visualization
            st.session_state["Biclusters"] = s

    if "Biclusters" in st.session_state and matrix is not None:
        gene_universe = set(gene_ids)
        p_vals = (0.05, 0.01, 0.001)
        all_enrich = []
        alg_names = list(st.session_state["Biclusters"].keys())
        tabs = st.tabs(alg_names)

        for (alg, bic_list), tab in zip(st.session_state["Biclusters"].items(), tabs):
            with tab:
                st.header(alg)
                st.write(summarize_biclusters(bic_list, gene_ids, alg))

                bic_gene_lists = [[gene_ids[i] for i in bic.rows] for bic in bic_list]
                enrich = go_assessment(int(taxid), bic_gene_lists, gene_universe,
                                       p_vals=p_vals)

                row = {"Algorithm": alg}
                for pv in p_vals:
                    row[str(pv)] = 100 * enrich[pv]
                all_enrich.append(row)

        if all_enrich:
            enrich_df = pd.DataFrame(all_enrich).set_index("Algorithm")
            enrich_df = enrich_df[[str(p) for p in p_vals]]
            st.bar_chart(enrich_df)

if __name__=="__main__":
    main()
