BCTool — Biclustering Tool

BCTool is a comprehensive graphical interface for analyzing large-scale gene expression datasets. It takes in CSV, TSV, XLSV, and XLS datasets and find bicluster submatrices using various methods.

🔧 Features Overview

📥 Data Handling & Conversion

Multi-Format Compatibility: Supports CSV, TSV, Excel (.xls, .xlsx)

Automatic Validation: Performs file integrity and format checks during loading

🔬 Analytical Tools

Integrated Biclustering Algorithms:

  LAS (Large Average Submatrices)

  Chen & Church Method

  ISA (Iterative Signature Algorithm)

  OPSM (Order-Preserving Submatrix)

  BiVisu with visualization support



Output:

  All biclusters found for the selected algorithms




⚙️ Installation Instructions

Requirements

Python 3.8 or newer

pip package manager

Quick Start

pip install -r requirements.txt
streamlit run BCTool_GUI.py



Manual Setup

pip install -r requirements.txt
streamlit run BCTool_GUI.py

Required Python Packages

pandas >= 1.3.0

numpy >= 1.21.0

matplotlib >= 3.5.0

seaborn >= 0.11.0

scikit-learn >= 1.0.0

scipy >= 1.7.0

openpyxl >= 3.0.0







🖥️ Using the Application

Launching

streamlit run BCTool_GUI.py


Loading Your Data

Method 1: Expression Matrix Files

Select the Data Loading tab

Choose Load CSV/TSV or Load Excel

The matrix should have gene identifiers in the first column and conditions/samples in the remaining columns



Method 2: HCL File Processing

Click Load HCL File to select a compatible HCL file

Press Convert to Matrix to transform HCL content into matrix form

Click Use HCL Matrix to proceed to analysis

Running Analyses

Navigate to the Analysis tab

Select the biclustering algorithms to apply

Start the process with Run Analysis

Results appear in the output panel

Visual Exploration

Open the Visualization tab

Available plot types:

Heatmaps for expression overview

Visual comparison of biclusters

Exporting Outputs

Head to the Results tab

Select the desired file format (CSV, JSON, Excel)

Export biclusters or matrices






📄 HCL File Format Example

Supported HCL files should resemble the following:

version = "1.0"
dataset_name = "gene_study"

gene_data "MYC" {
  expression_values = [4.2, 3.8, 2.1, 5.6]
  chromosome = "8"
}

gene_data "EGFR" {
  expression_values = [3.1, 2.9, 4.0, 3.7]
  chromosome = "7"
}

experimental_conditions {
  condition_1 = "Control"
  condition_2 = "Treatment"
  condition_3 = "Recovery"
  condition_4 = "Extended"
}

metadata {
  organism = "Homo sapiens"
  platform = "RNA-seq"
}

📁 Project Layout

BCTool/
├── BCTool_GUI.py          # Streamlit interface
├── adapter.py             # Algorithm adapters
├── LAS_algorithm.py
├── Chen_Church_algorithm.py
├── ISA_algorithm.py
├── OPSM_algorithm.py
├── Bivisu_algorithm.py
├── requirements.txt
└── README.md

🛠️ Algorithms Included

Biclustering

LAS: Locates large, high-mean submatrices with significance testing

Chen & Church: Classic biclustering via mean squared residue minimization

ISA: Iterative method to reveal co-expressed gene modules

OPSM: Discovers order-preserving submatrices in expression data

BiVisu: Combines biclustering with visual output capabilities


🛠️ Troubleshooting Guide

Issue

Recommended Solution


Import Errors

Verify dependencies with pip install -r requirements.txt


Memory Errors with Large Data

Subsample the dataset or upgrade system RAM (8GB+ recommended)


HCL Parsing Problems

Review HCL file structure against provided example




Additional Help

streamlit run BCTool_GUI.py

Check terminal messages for error details.



⚡ Optimization Tips

For datasets with more than 5,000 genes:

Start with a reduced dataset

Run algorithms separately to manage memory usage

Consider systems with higher RAM capacity


🧪 Running Tests
----------------

The repository includes a small test suite powered by `pytest`.
After installing the dependencies you can run all tests with:

```bash
pip install -r requirements.txt pytest
pytest
```



👍 Contributing to BCTool

To incorporate new biclustering algorithms:

Implement the algorithm in the algorithms/ directory following the adapter format

Register it within adapter.py under the ALL_ALGOS dictionary

Update GUI options accordingly



📜 License & Citation

BCTool integrates multiple published algorithms under their respective licenses. The overall project is released under the [MIT License](LICENSE). Review individual files for licensing details.

If you publish work utilizing BCTool, please acknowledge the source algorithms accordingly.


✅ Final Notes

BCTool simplifies complex gene expression analysis with a visual, user-friendly environment. Nevertheless, always perform appropriate validation to ensure research integrity.

