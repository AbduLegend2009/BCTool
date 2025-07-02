import streamlit as st
import pandas as pd
from LAS_algorithm import las_with_significance
from Chen_Church_algorithm import run_chen_church
from Bivisu_algorithm import bivisu
from ISA_algorithm import ISA_multi_seed
from OPSM_algorithm import run_opsm
from GO_assessment import go_assessment
def upload_matrix():
    """Handles the CSV uploader and populates session_state."""
    file = st.sidebar.file_uploader("Upload gene-expression matrix (CSV)")
    if file:
        df = pd.read_csv(file, index_col=0)       # genes as rows
        st.session_state["expr_df"] = df
        st.session_state["universe"] = set(df.index)   # ← THE universe
        st.success(f"Loaded {df.shape[0]:,} genes × {df.shape[1]:,} samples")
def run_biclusters():
    if "expr_df" not in st.session_state:
        st.error("Please upload data first")
        return {}
    df = st.session_state["expr_df"]
    biclusters = {}
    alg_options = ["LAS", "ISA", "OPSM", "Bivisu", "Chen and Church"]
    sel_algs = st.sidebar.multiselect("Choose algorithms", alg_options, default = [])
    st.session_state["sel_algs"]=sel_algs
    algo_map = {"LAS":las_with_significance, "Chen_and_Church":run_chen_church, "Bivisu": bivisu, "ISA": ISA_multi_seed, "OPSM": run_opsm}
    
