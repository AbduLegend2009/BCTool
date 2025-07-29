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




🖥️ Using the Application

Launching

Visit the following website:

bctoolgui.streamlit.app






Running Analyses

Navigate to the sidebar

Select the biclustering algorithms to apply

Start the process with Run Analysis

Results appear below


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

LAS: Locates large, high-mean submatrices with significance testing

Chen & Church: Classic biclustering via mean squared residue minimization

ISA: Iterative method to reveal co-expressed gene modules

OPSM: Discovers order-preserving submatrices in expression data

BiVisu: Combines biclustering with visual output capabilities



⚡ Optimization Tips

For datasets with more than 5,000 genes:

Start with a reduced dataset

Run algorithms separately to manage memory usage


📜 License & Citation

BCTool integrates multiple published algorithms under their respective licenses. The overall project is released under the [MIT License](LICENSE). Review individual files for licensing details.

If you publish work utilizing BCTool, please acknowledge the source algorithms accordingly.


✅ Final Notes

BCTool simplifies complex gene expression analysis with a visual, user-friendly environment. Nevertheless, always perform appropriate validation to ensure research integrity.


📣 Wrap Up

Thank you for exploring **BCTool**. If you encounter issues or wish to contribute enhancements, open a GitHub issue or pull request. Your feedback helps improve the project.
