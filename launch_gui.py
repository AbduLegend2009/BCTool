#!/usr/bin/env python3
"""
Gene Expression Analysis GUI Launcher
=====================================
This script launches the Gene Expression Matrix Analyzer GUI application.
It handles dependency checking and provides helpful error messages.
Usage:
    python launch_gui.py
Features:
- Load and analyze gene expression matrices
- Convert HCL text files to expression matrices
- Run multiple biclustering algorithms
- Perform GO enrichment analysis
- Visualize results with interactive plots
"""

import sys
import subprocess
import importlib
import tkinter as tk
from tkinter import messagebox
from pathlib import Path


def check_dependencies():
    """Check if required dependencies are installed"""
    required_packages = [
        'pandas', 'numpy', 'matplotlib', 'seaborn',
        'sklearn', 'scipy', 'openpyxl'
    ]

    missing_packages = []

    for package in required_packages:
        try:
            if package == 'sklearn':
                importlib.import_module('sklearn')
            else:
                importlib.import_module(package)
        except ImportError:
            missing_packages.append(package)

    return missing_packages


def install_packages(packages):
    """Install missing packages using pip"""
    for package in packages:
        try:
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", package])
            print(f"Successfully installed {package}")
        except subprocess.CalledProcessError:
            print(f"Failed to install {package}")
            return False
    return True


def check_files():
    """Check if required files exist"""
    required_files = [
        'gene_expression_gui.py',
        'adapter.py',
        'GO_assessment.py',
        'hcl_parser.py'
    ]

    missing_files = []
    for file in required_files:
        if not Path(file).exists():
            missing_files.append(file)

    return missing_files


def main():
    """Main launcher function"""
    print("Gene Expression Analysis GUI Launcher")
    print("=" * 40)

    # Check if required files exist
    missing_files = check_files()
    if missing_files:
        print("Error: Missing required files:")
        for file in missing_files:
            print(f"  - {file}")
        print("\nPlease ensure all required files are in the current directory.")
        return

    # Check dependencies
    print("Checking dependencies...")
    missing_packages = check_dependencies()

    if missing_packages:
        print(f"Missing packages: {', '.join(missing_packages)}")

        # Ask user if they want to install missing packages
        if len(sys.argv) > 1 and '--auto-install' in sys.argv:
            install_missing = True
        else:
            response = input(
                "Would you like to install missing packages? (y/n): ").lower()
            install_missing = response in ['y', 'yes']

        if install_missing:
            print("Installing missing packages...")
            if not install_packages(missing_packages):
                print("Failed to install some packages. Please install manually:")
                print(f"pip install {' '.join(missing_packages)}")
                return
        else:
            print("Please install missing packages manually:")
            print(f"pip install {' '.join(missing_packages)}")
            return

    # Launch the GUI
    print("Starting Gene Expression Analysis GUI...")
    try:
        from gene_expression_gui import GeneExpressionGUI

        root = tk.Tk()
        app = GeneExpressionGUI(root)

        # Set up proper window closing
        def on_closing():
            root.quit()
            root.destroy()

        root.protocol("WM_DELETE_WINDOW", on_closing)

        print("GUI launched successfully!")
        print("Close this terminal window or press Ctrl+C to exit.")

        root.mainloop()

    except ImportError as e:
        print(f"Import error: {e}")
        print(
            "Please ensure all required files are present and dependencies are installed.")
    except Exception as e:
        print(f"Error launching GUI: {e}")


def show_help():
    """Show help information"""
    help_text = """

This application provides a comprehensive interface for analyzing gene expression data.
Features:
---------
1. Data Loading:
   - Load CSV/TSV files with gene expression data
   - Load Excel files (.xlsx, .xls)
   - Convert HCL text files to expression matrices
2. Analysis:
   - Multiple biclustering algorithms (LAS, Chen & Church, ISA, OPSM, BiVisu)
   - GO enrichment analysis
   - Statistical summaries
3. Visualization:
   - Expression heatmaps
   - Bicluster analysis plots
   - GO enrichment charts
4. Export:
   - Export biclusters to JSON/CSV
   - Export GO results
   - Export expression matrices
File Formats:
------------
- CSV/TSV: First column should contain gene IDs, remaining columns are conditions
- Excel: Same format as CSV/TSV
- HCL: Hierarchical configuration language format for gene expression data
Sample HCL Format:
-----------------
gene_data "GENE_NAME" {
  expression_values = [1.2, 2.3, 1.8, ...]
  chromosome = "1"
  start_position = 12345
}
Usage:
------
python launch_gui.py [options]
Options:
  --auto-install    Automatically install missing packages
  --help           Show this help message
Requirements:
------------
- Python 3.7+
- pandas, numpy, matplotlib, seaborn
- scikit-learn, scipy, openpyxl
- tkinter (usually included with Python)
"""
    print(help_text)


if __name__ == "__main__":
    if '--help' in sys.argv or '-h' in sys.argv:
        show_help()
    else:
        main()
