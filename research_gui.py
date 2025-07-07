import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import seaborn as sns
from pathlib import Path
import threading
import json
import re
from typing import Dict, List, Optional, Tuple

from adapter import ALL_ALGOS
from GO_assessment import go_assessment, OBO, G2G
from hcl_parser import HCLParser


class GeneExpressionGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Gene Expression Matrix Analyzer")
        self.root.geometry("1400x900")

        # Data storage
        self.expression_data = None
        self.hcl_data = None
        self.biclusters = None
        self.go_results = None

        # Setup UI
        self.setup_ui()

    def setup_ui(self):
        # Create main notebook for tabs
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Create tabs
        self.create_data_tab()
        self.create_analysis_tab()
        self.create_visualization_tab()
        self.create_results_tab()

    def create_data_tab(self):
        # Data Loading and Conversion Tab
        data_frame = ttk.Frame(self.notebook)
        self.notebook.add(data_frame, text="Data Loading")

        # HCL File Conversion Section
        hcl_frame = ttk.LabelFrame(
            data_frame, text="HCL Text File Conversion", padding=10)
        hcl_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(hcl_frame, text="Load HCL File",
                   command=self.load_hcl_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(hcl_frame, text="Convert to Matrix",
                   command=self.convert_hcl_to_matrix).pack(side=tk.LEFT, padx=5)

        self.hcl_info_label = ttk.Label(hcl_frame, text="No HCL file loaded")
        self.hcl_info_label.pack(side=tk.LEFT, padx=10)

        # Expression Data Section
        expr_frame = ttk.LabelFrame(
            data_frame, text="Gene Expression Data", padding=10)
        expr_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # File loading buttons
        button_frame = ttk.Frame(expr_frame)
        button_frame.pack(fill=tk.X, pady=5)

        ttk.Button(button_frame, text="Load CSV/TSV",
                   command=self.load_expression_csv).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Load Excel",
                   command=self.load_expression_excel).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Use HCL Matrix",
                   command=self.use_hcl_matrix).pack(side=tk.LEFT, padx=5)

        # Data info
        self.data_info_text = scrolledtext.ScrolledText(expr_frame, height=15)
        self.data_info_text.pack(fill=tk.BOTH, expand=True, pady=5)

    def create_analysis_tab(self):
        # Analysis Configuration Tab
        analysis_frame = ttk.Frame(self.notebook)
        self.notebook.add(analysis_frame, text="Analysis")

        # Algorithm Selection
        algo_frame = ttk.LabelFrame(
            analysis_frame, text="Biclustering Algorithms", padding=10)
        algo_frame.pack(fill=tk.X, padx=5, pady=5)

        self.algo_vars = {}
        row = 0
        for i, algo_name in enumerate(ALL_ALGOS.keys()):
            var = tk.BooleanVar(value=True)
            self.algo_vars[algo_name] = var
            ttk.Checkbutton(algo_frame, text=algo_name, variable=var).grid(
                row=row//2, column=row % 2, sticky=tk.W, padx=10, pady=2)
            row += 1

        # Analysis Controls
        control_frame = ttk.LabelFrame(
            analysis_frame, text="Analysis Controls", padding=10)
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(control_frame, text="Run Analysis", command=self.run_analysis,
                   style="Accent.TButton").pack(side=tk.LEFT, padx=5)
        self.progress_var = tk.StringVar(value="Ready")
        ttk.Label(control_frame, textvariable=self.progress_var).pack(
            side=tk.LEFT, padx=10)

        # Progress bar
        self.progress_bar = ttk.Progressbar(
            control_frame, mode='indeterminate')
        self.progress_bar.pack(side=tk.RIGHT, padx=5)

        # Analysis Results Summary
        summary_frame = ttk.LabelFrame(
            analysis_frame, text="Analysis Summary", padding=10)
        summary_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.summary_text = scrolledtext.ScrolledText(summary_frame, height=15)
        self.summary_text.pack(fill=tk.BOTH, expand=True)

    def create_visualization_tab(self):
        # Visualization Tab
        viz_frame = ttk.Frame(self.notebook)
        self.notebook.add(viz_frame, text="Visualization")

        # Plot controls
        control_frame = ttk.Frame(viz_frame)
        control_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(control_frame, text="Heatmap",
                   command=self.plot_heatmap).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="Bicluster Visualization",
                   command=self.plot_biclusters).pack(side=tk.LEFT, padx=5)
        ttk.Button(control_frame, text="GO Enrichment",
                   command=self.plot_go_enrichment).pack(side=tk.LEFT, padx=5)

        # Matplotlib canvas
        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        self.canvas = FigureCanvasTkAgg(self.fig, viz_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def create_results_tab(self):
        # Results Export Tab
        results_frame = ttk.Frame(self.notebook)
        self.notebook.add(results_frame, text="Results")

        # Export controls
        export_frame = ttk.LabelFrame(
            results_frame, text="Export Results", padding=10)
        export_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(export_frame, text="Export Biclusters",
                   command=self.export_biclusters).pack(side=tk.LEFT, padx=5)
        ttk.Button(export_frame, text="Export GO Terms",
                   command=self.export_go_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(export_frame, text="Export Matrices",
                   command=self.export_matrices).pack(side=tk.LEFT, padx=5)

        # Results display
        self.results_tree = ttk.Treeview(results_frame, columns=(
            'Algorithm', 'Biclusters', 'Avg_Size'), show='tree headings')
        self.results_tree.heading('#0', text='ID')
        self.results_tree.heading('Algorithm', text='Algorithm')
        self.results_tree.heading('Biclusters', text='# Biclusters')
        self.results_tree.heading('Avg_Size', text='Avg Size')
        self.results_tree.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    def load_hcl_file(self):
        """Load HCL text file for conversion"""
        file_path = filedialog.askopenfilename(
            title="Select HCL Text File",
            filetypes=[("Text files", "*.txt"),
                       ("HCL files", "*.hcl"), ("All files", "*.*")]
        )

        if file_path:
            try:
                self.hcl_parser = HCLParser()
                self.hcl_data = self.hcl_parser.parse_file(file_path)
                self.hcl_info_label.config(
                    text=f"HCL file loaded: {Path(file_path).name}")
                messagebox.showinfo("Success", "HCL file loaded successfully!")
            except Exception as e:
                messagebox.showerror(
                    "Error", f"Failed to load HCL file: {str(e)}")

    def convert_hcl_to_matrix(self):
        """Convert HCL data to gene expression matrix format"""
        if not hasattr(self, 'hcl_parser') or not self.hcl_data:
            messagebox.showwarning("Warning", "Please load an HCL file first")
            return

        try:
            # Use the enhanced parser to convert to matrix
            self.hcl_matrix = self.hcl_parser.to_expression_matrix()
            self.hcl_metadata = self.hcl_parser.get_metadata()
            self.hcl_annotations = self.hcl_parser.get_gene_annotations()

            info_msg = f"HCL converted to matrix: {self.hcl_matrix.shape[0]} genes × {self.hcl_matrix.shape[1]} conditions"
            if self.hcl_metadata:
                info_msg += f"\nDataset: {self.hcl_metadata.get('dataset_name', 'Unknown')}"

            messagebox.showinfo("Success", info_msg)

        except Exception as e:
            messagebox.showerror("Error", f"Failed to convert HCL: {str(e)}")

    def load_expression_csv(self):
        """Load gene expression data from CSV/TSV"""
        file_path = filedialog.askopenfilename(
            title="Select Expression Data File",
            filetypes=[("CSV files", "*.csv"), ("TSV files",
                                                "*.tsv"), ("Text files", "*.txt")]
        )

        if file_path:
            try:
                # Try to detect delimiter
                with open(file_path, 'r') as f:
                    first_line = f.readline()
                    delimiter = '\t' if '\t' in first_line else ','

                self.expression_data = pd.read_csv(
                    file_path, index_col=0, sep=delimiter)
                self.expression_data = self.expression_data.apply(
                    pd.to_numeric, errors='coerce')
                self.expression_data.dropna(axis=1, how='all', inplace=True)

                self.update_data_info()
                messagebox.showinfo(
                    "Success", "Expression data loaded successfully!")

            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {str(e)}")

    def load_expression_excel(self):
        """Load gene expression data from Excel"""
        file_path = filedialog.askopenfilename(
            title="Select Excel File",
            filetypes=[("Excel files", "*.xlsx"), ("Excel files", "*.xls")]
        )

        if file_path:
            try:
                self.expression_data = pd.read_excel(file_path, index_col=0)
                self.expression_data = self.expression_data.apply(
                    pd.to_numeric, errors='coerce')
                self.expression_data.dropna(axis=1, how='all', inplace=True)

                self.update_data_info()
                messagebox.showinfo(
                    "Success", "Expression data loaded successfully!")

            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {str(e)}")

    def use_hcl_matrix(self):
        """Use converted HCL matrix as expression data"""
        if not hasattr(self, 'hcl_matrix'):
            messagebox.showwarning(
                "Warning", "Please convert HCL file to matrix first")
            return

        self.expression_data = self.hcl_matrix.copy()
        self.update_data_info()
        messagebox.showinfo("Success", "Using HCL matrix as expression data!")

    def update_data_info(self):
        """Update data information display"""
        if self.expression_data is not None:
            info = f"Data Shape: {self.expression_data.shape[0]} genes × {self.expression_data.shape[1]} conditions\n\n"
            info += f"Gene IDs (first 10):\n{list(self.expression_data.index[:10])}\n\n"
            info += f"Conditions (first 10):\n{list(self.expression_data.columns[:10])}\n\n"
            info += f"Data Statistics:\n{self.expression_data.describe()}\n\n"
            info += f"Missing Values: {self.expression_data.isnull().sum().sum()}"

            self.data_info_text.delete(1.0, tk.END)
            self.data_info_text.insert(1.0, info)

    def run_analysis(self):
        """Run biclustering analysis on loaded data"""
        if self.expression_data is None:
            messagebox.showwarning(
                "Warning", "Please load expression data first")
            return

        selected_algos = [name for name,
                          var in self.algo_vars.items() if var.get()]
        if not selected_algos:
            messagebox.showwarning(
                "Warning", "Please select at least one algorithm")
            return

        # Start analysis in separate thread
        threading.Thread(target=self._run_analysis_thread,
                         args=(selected_algos,), daemon=True).start()

    def _run_analysis_thread(self, selected_algos):
        """Run analysis in background thread"""
        try:
            self.progress_var.set("Running biclustering algorithms...")
            self.progress_bar.start()

            # Convert to numpy array
            data_matrix = self.expression_data.values.astype(np.float32)

            # Run selected algorithms
            self.biclusters = {}
            for algo_name in selected_algos:
                self.progress_var.set(f"Running {algo_name}...")
                try:
                    algo_func = ALL_ALGOS[algo_name]
                    biclusters = algo_func(data_matrix)
                    self.biclusters[algo_name] = biclusters
                except Exception as e:
                    print(f"Error running {algo_name}: {e}")
                    self.biclusters[algo_name] = []

            # Prompt for taxonomy ID
            def ask_taxid():
                taxid = tk.simpledialog.askstring(
                    "Taxonomy ID",
                    "Enter NCBI Taxonomy ID for your organism (e.g., 9606 for human):",
                    initialvalue="9606"
                )
                if taxid is None or not taxid.strip().isdigit():
                    return 9606  # Default to human
                return int(taxid.strip())

            # Ask for taxid in main thread
            taxid = [9606]
            event = threading.Event()
            def get_taxid():
                taxid[0] = ask_taxid()
                event.set()
            self.root.after(0, get_taxid)
            event.wait()

            # Prepare biclusters and universe for GO_assessment
            # Flatten all biclusters from all algorithms into a list of lists of gene indices
            all_biclusters = []
            for algo_bics in self.biclusters.values():
                for bic in algo_bics:
                    # Convert .rows to list of gene indices (as integers or as gene IDs)
                    # If expression_data.index is gene names, map indices to names
                    if hasattr(self.expression_data, 'index'):
                        genes = [self.expression_data.index[i] for i in bic.rows]
                    else:
                        genes = list(bic.rows)
                    all_biclusters.append(genes)
            # Universe is all gene IDs in the data
            universe = set(self.expression_data.index)

            # Run GO enrichment analysis
            self.progress_var.set("Running GO enrichment analysis...")
            try:
                self.go_results = go_assessment(
                    taxid[0], all_biclusters, universe, OBO=OBO, G2G=G2G)
            except Exception as e:
                print(f"Error in GO analysis: {e}")
                self.go_results = None

            # Update UI
            self.root.after(0, self._analysis_complete)

        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror(
                "Error", f"Analysis failed: {str(e)}"))
        finally:
            self.root.after(0, lambda: self.progress_bar.stop())

    def _analysis_complete(self):
        """Handle analysis completion"""
        self.progress_var.set("Analysis complete!")
        self.update_analysis_summary()
        self.update_results_tree()
        messagebox.showinfo("Success", "Analysis completed successfully!")

    def update_analysis_summary(self):
        """Update analysis summary display"""
        if not self.biclusters:
            return

        summary = "BICLUSTERING ANALYSIS RESULTS\n"
        summary += "=" * 50 + "\n\n"

        total_biclusters = 0
        for algo_name, biclusters in self.biclusters.items():
            summary += f"{algo_name}:\n"
            summary += f"  • Found {len(biclusters)} biclusters\n"

            if biclusters:
                avg_genes = np.mean([len(bc.rows) for bc in biclusters])
                avg_conditions = np.mean([len(bc.cols) for bc in biclusters])
                avg_score = np.mean([bc.score for bc in biclusters])

                summary += f"  • Average genes per bicluster: {avg_genes:.1f}\n"
                summary += f"  • Average conditions per bicluster: {avg_conditions:.1f}\n"
                summary += f"  • Average score: {avg_score:.3f}\n"

            summary += "\n"
            total_biclusters += len(biclusters)

        summary += f"Total biclusters found: {total_biclusters}\n\n"

        if self.go_results is not None:
            summary += f"GO ENRICHMENT ANALYSIS\n"
            summary += "=" * 50 + "\n"
            summary += f"Found {len(self.go_results)} enriched GO terms\n"
            if len(self.go_results) > 0:
                summary += f"Top enriched terms:\n"
                for _, row in self.go_results.head().iterrows():
                    summary += f"  • {row.get('go_term', 'Unknown')}: p={row.get('p_value', 'N/A')}\n"

        self.summary_text.delete(1.0, tk.END)
        self.summary_text.insert(1.0, summary)

    def update_results_tree(self):
        """Update results tree view"""
        # Clear existing items
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)

        if not self.biclusters:
            return

        for algo_name, biclusters in self.biclusters.items():
            avg_size = 0
            if biclusters:
                avg_size = np.mean([len(bc.rows) * len(bc.cols)
                                   for bc in biclusters])

            self.results_tree.insert('', 'end', text=algo_name,
                                     values=(algo_name, len(biclusters), f"{avg_size:.1f}"))

    def plot_heatmap(self):
        """Plot expression data heatmap"""
        if self.expression_data is None:
            messagebox.showwarning("Warning", "No data to plot")
            return

        self.ax.clear()

        # Sample data if too large
        data_to_plot = self.expression_data
        if data_to_plot.shape[0] > 100:
            data_to_plot = data_to_plot.sample(n=100)

        if data_to_plot.shape[1] > 50:
            data_to_plot = data_to_plot.iloc[:, :50]

        sns.heatmap(data_to_plot, ax=self.ax, cmap='viridis', cbar=True)
        self.ax.set_title("Gene Expression Heatmap")
        self.ax.set_xlabel("Conditions")
        self.ax.set_ylabel("Genes")

        self.canvas.draw()

    def plot_biclusters(self):
        """Plot bicluster visualization"""
        if not self.biclusters:
            messagebox.showwarning("Warning", "No biclusters to plot")
            return

        self.ax.clear()

        # Plot bicluster sizes
        algo_names = []
        bicluster_counts = []
        avg_sizes = []

        for algo_name, biclusters in self.biclusters.items():
            algo_names.append(algo_name)
            bicluster_counts.append(len(biclusters))
            if biclusters:
                avg_sizes.append(
                    np.mean([len(bc.rows) * len(bc.cols) for bc in biclusters]))
            else:
                avg_sizes.append(0)

        x = np.arange(len(algo_names))
        width = 0.35

        self.ax.bar(x - width/2, bicluster_counts, width,
                    label='# Biclusters', alpha=0.8)
        self.ax.bar(x + width/2, avg_sizes, width, label='Avg Size', alpha=0.8)

        self.ax.set_xlabel('Algorithm')
        self.ax.set_ylabel('Count / Size')
        self.ax.set_title('Bicluster Analysis Results')
        self.ax.set_xticks(x)
        self.ax.set_xticklabels(algo_names, rotation=45)
        self.ax.legend()

        self.canvas.draw()

    def plot_go_enrichment(self):
        """Plot GO enrichment results"""
        if self.go_results is None or len(self.go_results) == 0:
            messagebox.showwarning("Warning", "No GO enrichment data to plot")
            return

        self.ax.clear()

        # Plot top 10 GO terms
        top_terms = self.go_results.head(10)

        if 'p_value' in top_terms.columns and 'go_term' in top_terms.columns:
            y_pos = np.arange(len(top_terms))
            p_values = -np.log10(top_terms['p_value'].astype(float))

            self.ax.barh(y_pos, p_values)
            self.ax.set_yticks(y_pos)
            self.ax.set_yticklabels(top_terms['go_term'].str[:50])
            self.ax.set_xlabel('-log10(p-value)')
            self.ax.set_title('Top GO Terms Enrichment')

        self.canvas.draw()

    def export_biclusters(self):
        """Export bicluster results to file"""
        if not self.biclusters:
            messagebox.showwarning("Warning", "No biclusters to export")
            return

        file_path = filedialog.asksaveasfilename(
            title="Save Biclusters",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("CSV files", "*.csv")]
        )

        if file_path:
            try:
                if file_path.endswith('.json'):
                    # Export as JSON
                    export_data = {}
                    for algo_name, biclusters in self.biclusters.items():
                        export_data[algo_name] = []
                        for i, bc in enumerate(biclusters):
                            export_data[algo_name].append({
                                'id': i,
                                'genes': bc.rows.tolist() if hasattr(bc.rows, 'tolist') else list(bc.rows),
                                'conditions': bc.cols.tolist() if hasattr(bc.cols, 'tolist') else list(bc.cols),
                                'score': float(bc.score)
                            })

                    with open(file_path, 'w') as f:
                        json.dump(export_data, f, indent=2)

                else:
                    # Export as CSV
                    rows = []
                    for algo_name, biclusters in self.biclusters.items():
                        for i, bc in enumerate(biclusters):
                            rows.append({
                                'Algorithm': algo_name,
                                'Bicluster_ID': i,
                                'Num_Genes': len(bc.rows),
                                'Num_Conditions': len(bc.cols),
                                'Score': bc.score,
                                'Genes': ','.join(map(str, bc.rows)),
                                'Conditions': ','.join(map(str, bc.cols))
                            })

                    pd.DataFrame(rows).to_csv(file_path, index=False)

                messagebox.showinfo(
                    "Success", f"Biclusters exported to {file_path}")

            except Exception as e:
                messagebox.showerror("Error", f"Failed to export: {str(e)}")

    def export_go_results(self):
        """Export GO enrichment results"""
        if self.go_results is None:
            messagebox.showwarning("Warning", "No GO results to export")
            return

        file_path = filedialog.asksaveasfilename(
            title="Save GO Results",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")]
        )

        if file_path:
            try:
                self.go_results.to_csv(file_path, index=False)
                messagebox.showinfo(
                    "Success", f"GO results exported to {file_path}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to export: {str(e)}")

    def export_matrices(self):
        """Export expression matrices"""
        if self.expression_data is None:
            messagebox.showwarning("Warning", "No matrix data to export")
            return

        file_path = filedialog.asksaveasfilename(
            title="Save Expression Matrix",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx")]
        )

        if file_path:
            try:
                if file_path.endswith('.xlsx'):
                    self.expression_data.to_excel(file_path)
                else:
                    self.expression_data.to_csv(file_path)

                messagebox.showinfo(
                    "Success", f"Matrix exported to {file_path}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to export: {str(e)}")


if __name__ == "__main__":
    root = tk.Tk()
    app = GeneExpressionGUI(root)
    root.mainloop()
