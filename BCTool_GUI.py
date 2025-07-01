import streamlit as st
import pandas as pd
from GO_assessment import go_assessment
def upload_matrix():
    """Handles the CSV uploader and populates session_state."""
    file = st.sidebar.file_uploader("Upload gene-expression matrix (CSV)")
    if file:
        df = pd.read_csv(file, index_col=0)       # genes as rows
        st.session_state["expr_df"] = df
        st.session_state["universe"] = set(df.index)   # ← THE universe
        st.success(f"Loaded {df.shape[0]:,} genes × {df.shape[1]:,} samples")