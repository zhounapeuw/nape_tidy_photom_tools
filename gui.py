import tkinter as tk
from tkinter import filedialog, simpledialog, messagebox, ttk
from tkinter.scrolledtext import ScrolledText
import pandas as pd
import threading
import os
import plotly.graph_objects as go
import plotly.offline as pyo
import sys
import glob

from functions.tdt_analysis import fp_preprocess  # Ensure this is in the same directory


# --- Redirect stdout to GUI ---
class StdoutRedirector:
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, s):
        self.text_widget.insert(tk.END, s)
        self.text_widget.see(tk.END)
        self.text_widget.update()

    def flush(self):
        pass


class FiberPhotometryGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Fiber Photometry Preprocessing GUI")

        self.log_fp_path = None
        self.overwrite = 0
        self.df = None  # Cached DataFrame
        self.fiber_id = tk.StringVar(value="1")
        self.processed_data_path = None
        self.processed_files = []

        self.traces = [
            "signal", "control", "control_fitted", "delta_signal_fitted_control", "poly_signal", "poly_control",
            "delta_signal_poly_zscore", "delta_control_poly_zscore"
        ]
        self.selected_traces = {trace: tk.BooleanVar(value=True) for trace in self.traces}

        self.create_widgets()

    def create_widgets(self):
        # Create notebook for tabs
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        # Preprocessing tab
        self.preprocessing_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.preprocessing_frame, text="Preprocessing")

        # Visualization tab
        self.visualization_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.visualization_frame, text="Data Visualization")

        self.create_preprocessing_widgets()
        self.create_visualization_widgets()

    def create_preprocessing_widgets(self):
        # File picker
        frame_top = tk.Frame(self.preprocessing_frame)
        frame_top.pack(fill=tk.X, padx=10, pady=5)

        self.log_fp_label = tk.Label(frame_top, text="No log_fp file selected")
        self.log_fp_label.pack(side=tk.LEFT)

        tk.Button(frame_top, text="Browse log_fp", command=self.browse_log_fp).pack(side=tk.RIGHT)

        # Overwrite parameter dialog
        param_frame = tk.Frame(self.preprocessing_frame)
        param_frame.pack(pady=5)
        tk.Label(param_frame, text="Overwrite (0 or 1):").pack(side=tk.LEFT)
        self.overwrite_entry = tk.Entry(param_frame, width=5)
        self.overwrite_entry.insert(0, "0")
        self.overwrite_entry.pack(side=tk.LEFT)

        # Fiber selection dropdown
        fiber_frame = tk.Frame(self.preprocessing_frame)
        fiber_frame.pack(padx=10, pady=5, fill=tk.X)
        tk.Label(fiber_frame, text="Select Fiber:").pack(side=tk.LEFT)
        self.fiber_dropdown = ttk.Combobox(fiber_frame, textvariable=self.fiber_id, values=["1"], state="readonly", width=10)
        self.fiber_dropdown.pack(side=tk.LEFT)

        # Trace selection checkboxes
        self.trace_frame = tk.LabelFrame(self.preprocessing_frame, text="Select traces to plot")
        self.trace_frame.pack(padx=10, pady=5, fill=tk.X)

        for i, trace in enumerate(self.traces):
            cb = tk.Checkbutton(self.trace_frame, text=trace, variable=self.selected_traces[trace])
            cb.grid(row=i // 3, column=i % 3, sticky="w")

        # Run and refresh buttons
        button_frame = tk.Frame(self.preprocessing_frame)
        button_frame.pack(pady=10)

        tk.Button(button_frame, text="Run Preprocessing", command=self.run_pipeline).pack(side=tk.LEFT, padx=10)
        tk.Button(button_frame, text="Refresh Plot", command=self.refresh_plot).pack(side=tk.LEFT, padx=10)

        # Console output
        self.console = ScrolledText(self.preprocessing_frame, height=10)
        self.console.pack(padx=10, pady=5, fill=tk.BOTH, expand=True)
        sys.stdout = StdoutRedirector(self.console)

    def create_visualization_widgets(self):
        # Processed data selection frame
        data_frame = tk.LabelFrame(self.visualization_frame, text="Select Processed Data")
        data_frame.pack(fill=tk.X, padx=10, pady=5)

        # Browse button for processed data folder
        browse_frame = tk.Frame(data_frame)
        browse_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.processed_data_label = tk.Label(browse_frame, text="No processed data folder selected")
        self.processed_data_label.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        tk.Button(browse_frame, text="Browse Processed Data", command=self.browse_processed_data).pack(side=tk.RIGHT)

        # File selection dropdown
        file_frame = tk.Frame(data_frame)
        file_frame.pack(fill=tk.X, padx=5, pady=5)
        
        tk.Label(file_frame, text="Select File:").pack(side=tk.LEFT)
        self.file_dropdown = ttk.Combobox(file_frame, values=[], state="readonly", width=50)
        self.file_dropdown.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(5, 0))
        self.file_dropdown.bind('<<ComboboxSelected>>', self.on_file_selected)

        # Load and plot buttons
        button_frame = tk.Frame(data_frame)
        button_frame.pack(pady=5)
        
        tk.Button(button_frame, text="Load Selected File", command=self.load_selected_file).pack(side=tk.LEFT, padx=5)
        tk.Button(button_frame, text="Refresh File List", command=self.refresh_file_list).pack(side=tk.LEFT, padx=5)

        # Fiber selection for visualization
        viz_fiber_frame = tk.Frame(self.visualization_frame)
        viz_fiber_frame.pack(padx=10, pady=5, fill=tk.X)
        tk.Label(viz_fiber_frame, text="Select Fiber:").pack(side=tk.LEFT)
        self.viz_fiber_dropdown = ttk.Combobox(viz_fiber_frame, textvariable=self.fiber_id, values=["1"], state="readonly", width=10)
        self.viz_fiber_dropdown.pack(side=tk.LEFT)

        # Trace selection checkboxes for visualization
        self.viz_trace_frame = tk.LabelFrame(self.visualization_frame, text="Select traces to plot")
        self.viz_trace_frame.pack(padx=10, pady=5, fill=tk.X)

        for i, trace in enumerate(self.traces):
            cb = tk.Checkbutton(self.viz_trace_frame, text=trace, variable=self.selected_traces[trace])
            cb.grid(row=i // 3, column=i % 3, sticky="w")

        # Plot button for visualization
        plot_button_frame = tk.Frame(self.visualization_frame)
        plot_button_frame.pack(pady=10)
        
        tk.Button(plot_button_frame, text="Plot Selected Data", command=self.plot_visualization_data).pack(padx=10)

    def browse_log_fp(self):
        path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
        if path:
            self.log_fp_path = path
            self.log_fp_label.config(text=os.path.basename(path))

    def run_pipeline(self):
        if not self.log_fp_path:
            messagebox.showerror("Missing Input", "Please select a log_fp file.")
            return
        try:
            self.overwrite = int(self.overwrite_entry.get())
        except ValueError:
            messagebox.showerror("Invalid Input", "Overwrite must be 0 or 1.")
            return

        threading.Thread(target=self.run_processing, daemon=True).start()

    def run_processing(self):
        try:
            dir_extracted = os.path.join(os.path.dirname(self.log_fp_path), "extracted")
            dir_processed = os.path.join(os.path.dirname(self.log_fp_path), "processed")
            log_fp = pd.read_csv(self.log_fp_path)

            print("Running preprocessing...")
            fp_preprocess(dir_extracted, dir_processed, log_fp, overwrite=self.overwrite)

            output_path = os.path.join(dir_processed, log_fp['blockname'][0] + "_streams_session.csv")
            if os.path.exists(output_path):
                df = pd.read_csv(output_path)
                self.df = df  # Cache DataFrame

                # Populate fiber ID dropdown
                fiber_ids = sorted(df['fiber_id'].dropna().astype(str).unique())
                self.fiber_dropdown['values'] = fiber_ids
                self.fiber_id.set(fiber_ids[0])

                print(f"Available fibers: {fiber_ids}")
                self.plot_selected_traces(df)
            else:
                print(f"No output file found at {output_path}")
        except Exception as e:
            print(f"Error: {e}")

    def refresh_plot(self):
        if self.df is not None:
            self.plot_selected_traces(self.df)
        else:
            print("No data to refresh. Please run preprocessing first.")

    def plot_selected_traces(self, df):
        # Filter by selected fiber
        fiber = self.fiber_id.get()
        df = df[df['fiber_id'].astype(str) == fiber]

        fig = go.Figure()
        for trace, var in self.selected_traces.items():
            if var.get() and trace in df.columns:
                y = df[trace].dropna()
                if not y.empty:
                    fig.add_trace(go.Scatter(x=df["time"], y=y, mode='lines', name=trace))

        fig.update_layout(
            title=f"Fiber Photometry Traces (Fiber {fiber})",
            xaxis_title="Time (s)",
            yaxis_title="AU or z-score"
        )

        html_path = "temp_plot.html"
        pyo.plot(fig, filename=html_path, auto_open=True)

    def browse_processed_data(self):
        path = filedialog.askdirectory(title="Select Processed Data Folder")
        if path:
            self.processed_data_path = path
            self.processed_data_label.config(text=os.path.basename(path))
            self.refresh_file_list()

    def refresh_file_list(self):
        if not self.processed_data_path:
            messagebox.showwarning("No Folder", "Please select a processed data folder first.")
            return
        
        try:
            # Find all CSV files in the processed data folder and subdirectories
            csv_files = []
            for root, dirs, files in os.walk(self.processed_data_path):
                for file in files:
                    if file.endswith('.csv'):
                        # Get relative path from the processed data folder
                        rel_path = os.path.relpath(os.path.join(root, file), self.processed_data_path)
                        csv_files.append(rel_path)
            
            self.processed_files = csv_files
            self.file_dropdown['values'] = csv_files
            
            if csv_files:
                self.file_dropdown.set(csv_files[0])
                print(f"Found {len(csv_files)} CSV files in processed data folder")
            else:
                print("No CSV files found in the selected folder")
                
        except Exception as e:
            print(f"Error refreshing file list: {e}")

    def on_file_selected(self, event):
        selected_file = self.file_dropdown.get()
        if selected_file:
            print(f"Selected file: {selected_file}")

    def load_selected_file(self):
        selected_file = self.file_dropdown.get()
        if not selected_file:
            messagebox.showwarning("No File", "Please select a file from the dropdown.")
            return
        
        try:
            file_path = os.path.join(self.processed_data_path, selected_file)
            df = pd.read_csv(file_path)
            self.df = df  # Cache DataFrame
            
            # Populate fiber ID dropdown if fiber_id column exists
            if 'fiber_id' in df.columns:
                fiber_ids = sorted(df['fiber_id'].dropna().astype(str).unique())
                self.viz_fiber_dropdown['values'] = fiber_ids
                if fiber_ids:
                    self.fiber_id.set(fiber_ids[0])
                print(f"Available fibers: {fiber_ids}")
            else:
                print("No 'fiber_id' column found in the selected file")
            
            print(f"Successfully loaded {selected_file}")
            print(f"Data shape: {df.shape}")
            
        except Exception as e:
            print(f"Error loading file: {e}")
            messagebox.showerror("Error", f"Failed to load file: {e}")

    def plot_visualization_data(self):
        if self.df is not None:
            self.plot_selected_traces(self.df)
        else:
            print("No data loaded. Please select and load a file first.")


# --- Run the App ---
if __name__ == "__main__":
    root = tk.Tk()
    app = FiberPhotometryGUI(root)
    root.geometry("1000x800")
    root.mainloop()