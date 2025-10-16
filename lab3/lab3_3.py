import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math

# --- Core Bioinformatics Functions (Unchanged) ---

def calculate_tm_basic(sequence: str) -> float:
    """Calculates Tm using the basic formula: 4(G+C) + 2(A+T)."""
    seq = sequence.upper()
    g_count = seq.count('G')
    c_count = seq.count('C')
    a_count = seq.count('A')
    t_count = seq.count('T')
    return float(4 * (g_count + c_count) + 2 * (a_count + t_count))

def calculate_tm_advanced(sequence: str, na_concentration: float) -> float:
    """Calculates Tm using the salt-adjusted formula."""
    seq_len = len(sequence)
    if seq_len == 0:
        return 0.0
    seq = sequence.upper()
    gc_percent = ((seq.count('G') + seq.count('C')) / seq_len) * 100
    tm = (81.5 +
          (16.6 * math.log10(na_concentration)) +
          (0.41 * gc_percent) -
          (600 / seq_len))
    return tm

def parse_fasta(filepath: str) -> str:
    """Parses a FASTA file and returns the concatenated DNA sequence."""
    sequence = ""
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence += line.strip()
    except Exception as e:
        messagebox.showerror("File Error", f"Could not read the file:\n{e}")
        return ""
    return sequence.upper()

# --- Main Application Class ---

class TmScannerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("🧬 Advanced DNA Melting Temperature Scanner")
        self.geometry("900x750")

        self.dna_sequence = ""
        self.window_size = 9

        # --- Create UI Frames ---
        top_frame = tk.Frame(self)
        top_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)
        
        control_frame = tk.LabelFrame(top_frame, text="Controls", padx=10, pady=10)
        control_frame.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        results_frame = tk.LabelFrame(top_frame, text="Statistics", padx=10, pady=10)
        results_frame.pack(side=tk.RIGHT, fill=tk.X, expand=True)

        chart_frame = tk.Frame(self)
        chart_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        # --- Control Widgets ---
        self.load_button = tk.Button(control_frame, text="Load FASTA File", command=self.load_file)
        self.load_button.grid(row=0, column=0, padx=5, pady=2, sticky='w')

        self.file_label = tk.Label(control_frame, text="No file loaded.", fg="gray", width=40, anchor='w')
        self.file_label.grid(row=0, column=1, columnspan=3, padx=5, sticky='w')

        tk.Label(control_frame, text="Na+ Conc (M):").grid(row=1, column=0, padx=5, pady=2, sticky='w')
        self.na_entry = tk.Entry(control_frame, width=10)
        self.na_entry.insert(0, "0.05")
        self.na_entry.grid(row=1, column=1, padx=5, sticky='w')

        tk.Label(control_frame, text="Tm Threshold (°C):").grid(row=1, column=2, padx=5, pady=2, sticky='w')
        self.threshold_entry = tk.Entry(control_frame, width=10)
        self.threshold_entry.insert(0, "40")
        self.threshold_entry.grid(row=1, column=3, padx=5, sticky='w')

        self.analyze_button = tk.Button(control_frame, text="Analyze Sequence", font=('Helvetica', 10, 'bold'), command=self.run_analysis)
        self.analyze_button.grid(row=0, column=4, rowspan=2, padx=10, ipady=10, sticky='e')
        
        # --- Results Display Widget ---
        self.results_label = tk.Label(results_frame, text="Min/Max values will appear here.", justify=tk.LEFT, anchor='w')
        self.results_label.pack(anchor='w')

        # --- Matplotlib Chart (with two subplots) ---
        self.fig, (self.ax1, self.ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
        self.canvas = FigureCanvasTkAgg(self.fig, master=chart_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.setup_charts()
        
    def setup_charts(self):
        """Initial setup for the charts' appearance."""
        # Top chart (Signals)
        self.ax1.set_title("Melting Temperature Profile")
        self.ax1.set_ylabel("Tm (°C)")
        self.ax1.grid(True, linestyle='--', alpha=0.6)
        
        # Bottom chart (Threshold Regions)
        self.ax2.set_xlabel(f"Window Start Position (Window Size = {self.window_size} bp)")
        self.ax2.set_ylabel("Regions\nAbove Thr.")
        self.ax2.set_yticks([0.5, 1.5])
        self.ax2.set_yticklabels(['Advanced', 'Basic'])
        self.ax2.set_ylim(0, 2)
        
        self.fig.tight_layout()
        self.canvas.draw()
        
    def load_file(self):
        filepath = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")])
        if not filepath: return
        self.dna_sequence = parse_fasta(filepath)
        if self.dna_sequence:
            self.file_label.config(text=f"{filepath.split('/')[-1]} ({len(self.dna_sequence)} bp)", fg="black")

    def run_analysis(self):
        if not self.dna_sequence:
            messagebox.showwarning("No Data", "Please load a FASTA file first.")
            return

        try:
            na_concentration = float(self.na_entry.get())
            threshold = float(self.threshold_entry.get())
            if na_concentration <= 0: raise ValueError("Concentration must be positive.")
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter valid positive numbers for concentration and threshold.")
            return

        # Perform the sliding window analysis
        positions, basic_tms, advanced_tms = self.sliding_window_analysis(na_concentration)
        
        if not positions:
            messagebox.showwarning("Sequence Too Short", f"Sequence is shorter than the window size ({self.window_size} bp).")
            return
            
        # Update statistics
        self.update_statistics(basic_tms, advanced_tms)
        
        # Find regions above the threshold
        basic_regions = self.find_regions_above_threshold(basic_tms, threshold)
        advanced_regions = self.find_regions_above_threshold(advanced_tms, threshold)

        # Plot the results
        self.plot_data(positions, basic_tms, advanced_tms, threshold, basic_regions, advanced_regions)

    def sliding_window_analysis(self, na_concentration: float):
        positions, basic_tms, advanced_tms = [], [], []
        if len(self.dna_sequence) < self.window_size: return [], [], []
            
        for i in range(len(self.dna_sequence) - self.window_size + 1):
            sub_sequence = self.dna_sequence[i : i + self.window_size]
            positions.append(i)
            basic_tms.append(calculate_tm_basic(sub_sequence))
            advanced_tms.append(calculate_tm_advanced(sub_sequence, na_concentration))
        return positions, basic_tms, advanced_tms

    def find_regions_above_threshold(self, data: list, threshold: float) -> list:
        """Identifies contiguous segments in the data that are above a threshold."""
        regions, in_region, start_idx = [], False, 0
        for i, value in enumerate(data):
            is_above = value >= threshold
            if is_above and not in_region:
                in_region = True
                start_idx = i
            elif not is_above and in_region:
                in_region = False
                regions.append((start_idx, i)) # region is from start_idx to i-1
        if in_region: # If the sequence ends while in a region
            regions.append((start_idx, len(data)))
        return regions

    def update_statistics(self, basic_tms, advanced_tms):
        stats_text = (
            f"Basic Tm:\tMin: {min(basic_tms):.2f}°C\tMax: {max(basic_tms):.2f}°C\n"
            f"Advanced Tm:\tMin: {min(advanced_tms):.2f}°C\tMax: {max(advanced_tms):.2f}°C"
        )
        self.results_label.config(text=stats_text)

    def plot_data(self, positions, basic_tms, advanced_tms, threshold, basic_regions, advanced_regions):
        # Clear previous plots
        self.ax1.cla()
        self.ax2.cla()

        # Plot main signals on the top chart
        self.ax1.plot(positions, basic_tms, label="Basic Formula (Tm)", color="blue", alpha=0.8)
        self.ax1.plot(positions, advanced_tms, label="Advanced Formula (Tm)", color="red", alpha=0.8)
        self.ax1.axhline(y=threshold, color='green', linestyle='--', label=f'Threshold ({threshold}°C)')
        self.ax1.legend()

        # Plot horizontal bars on the bottom chart
        for start, end in basic_regions:
            self.ax2.broken_barh([(start, end - start)], (1.25, 0.5), facecolors='blue', alpha=0.7)
        for start, end in advanced_regions:
            self.ax2.broken_barh([(start, end - start)], (0.25, 0.5), facecolors='red', alpha=0.7)
        
        # Redraw charts with proper labels
        self.setup_charts()

if __name__ == "__main__":
    app = TmScannerApp()
    app.mainloop()