import streamlit as st
import pandas as pd
import io, numpy as np
import GO_assessment
from adapter import ALL_ALGOS

def uploaded_data():
    data = st.file_uploader("Please upload expression matrix (CSV, TSV, or Excel)", type = ["csv", "tsv", "xls", "xlsx"])
    suffix = data.name.split(".")[-1].lower # type: ignore
    if suffix in {"csv", "tsv"}:
        df = pd.read_csv(data, sep="\t" if suffix == "tsv" else ",", index_col=0) # type: ignore
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


def main():
    st.title("BCTool 🧬")
    with st.sidebar:
        st.header("settings")
        matrix, gene_ids = uploaded_data()
    

