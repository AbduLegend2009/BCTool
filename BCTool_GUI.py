import streamlit as st
import pandas as pd
import io, numpy as np
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


def summarize_biclusters(bic_dict, gene_ids):
    """
    bic_dict : {"LAS": [Bicluster,…], "ISA": …}
    gene_ids : list[str]  – index → gene symbol
    returns   : pd.DataFrame
    """
    rows = []
    for alg, bics in bic_dict.items():
        for i, bic in enumerate(bics):
            rows.append({
                "ID"       : f"{alg}_{i}",
                "Algorithm": alg,
                "Rows"     : len(bic.rows),
                "Cols"     : len(bic.cols),
                "Genes"    : ", ".join(gene_ids[j] for j in bic.rows[:3]) + " …",
                "Score"    : getattr(bic, "score", np.nan),   # not all wrappers set score
            })
    return pd.DataFrame(rows)



def main():
    st.title("BCTool 🧬")
    data = st.file_uploader("Please upload expression matrix (CSV, TSV, or Excel)", type = ["csv", "tsv", "xls", "xlsx"])
    with st.sidebar:
        st.header("Customization")
        matrix, gene_ids = uploaded_data(data)
        taxid = st.number_input("What is the taxonomy identifier?")
        sel_alg = st.multiselect("Please select algorithms", algorithms)
        if "LAS" in sel_alg:
            st.header("LAS")
            alpha_LAS = st.number_input("What is the p-value cutoff for biclusters?")
            max_iter_LAS = st.number_input("What is the maximum amount of iterations?")
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
            k_OPSM = st.number_input("What is the maximum amoung of conditions in each bicluster?")
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
            s = []
            for _ in sel_alg:
                if _ == "LAS":
                    LAS = adapter.wrap_las(matrix, max_iter_LAS, alpha_LAS)
                    s.append(["LAS",LAS])
                elif _ == "Chen and Church":
                    Chen_and_Church = adapter.wrap_church(matrix, delta_C_C, delta_C_C, max_bi_C_C)
                    s.append(["Chen_and_Church", Chen_and_Church])
                elif _ == "ISA":
                    ISA = adapter.wrap_isa(matrix, n_seeds_ISA, seed_size_ISA, t_g_ISA, t_c_ISA)
                    s.append(["ISA", ISA])
                elif _ == "OPSM":
                    opsm = adapter.wrap_opsm(matrix, k_OPSM, restarts_OPSM)
                    s.append(["OPSM", opsm])
                elif _ == "Bivisu":
                    Bi = adapter.wrap_bivisu(matrix, model_Bivisu, eps_Bivisu, msr_Bivisu, min_genes_Bivisu, min_cond_Bivisu)
                    s.append(["Bivisu", Bi])
            st.session_state["Biclusters"] = s
    for sub in s:
        st.write(f"{sub[0]}\n{summarise_biclusters(sub[1], gene_ids)}")
    






 

    



if __name__=="__main__":
    main()