import streamlit as st
import pandas as pd
import io, numpy as np
import GO_assessment
from adapter import ALL_ALGOS

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

def main():
    st.title("BCTool 🧬")
    with st.sidebar:
        st.header("settings")
        data = st.file_uploader("Please upload expression matrix (CSV, TSV, or Excel)", type = ["csv", "tsv", "xls", "xlsx"])
        matrix, gene_ids = uploaded_data(data)
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

    if st.button("Run Algorithms 🚀"):
