import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import math

# --- Core Bioinformatics Functions ---

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
    g_count = seq.count('G')
    c_count = seq.count('C')
    gc_percent = ((g_count + c_count) / seq_len) * 100

    # Formula: Tm = 81.5 + 16.6(log10([Na+])) + 0.41*(%GC) – 600/length
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
        self.title("🧬 DNA Melting Temperature Scanner")
        self.geometry("800x600")

        self.dna_sequence = ""
        self.window_size = 9

        # --- Create UI Frames ---
        control_frame = tk.Frame(self, padx=10, pady=10)
        control_frame.pack(side=tk.TOP, fill=tk.X)

        chart_frame = tk.Frame(self)
        chart_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        # --- Control Widgets ---
        self.load_button = tk.Button(control_frame, text="Load FASTA File", command=self.load_file)
        self.load_button.pack(side=tk.LEFT, padx=5)

        self.file_label = tk.Label(control_frame, text="No file loaded.", fg="gray")
        self.file_label.pack(side=tk.LEFT, padx=5)

        tk.Label(control_frame, text="Na+ Conc (M):").pack(side=tk.LEFT, padx=(20, 0))
        self.na_entry = tk.Entry(control_frame, width=10)
        self.na_entry.insert(0, "0.05") # Default value
        self.na_entry.pack(side=tk.LEFT, padx=5)

        self.analyze_button = tk.Button(control_frame, text="Analyze", command=self.run_analysis)
        self.analyze_button.pack(side=tk.LEFT, padx=5)

        # --- Matplotlib Chart ---
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=chart_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.setup_chart()
        
    def setup_chart(self):
        """Initial setup for the chart appearance."""
        self.ax.set_title("Melting Temperature Profile")
        self.ax.set_xlabel(f"Window Start Position (Window Size = {self.window_size} bp)")
        self.ax.set_ylabel("Melting Temperature (Tm) °C")
        self.ax.grid(True, linestyle='--', alpha=0.6)
        self.fig.tight_layout()
        self.canvas.draw()
        
    def load_file(self):
        """Opens a file dialog to select a FASTA file and loads the sequence."""
        filepath = filedialog.askopenfilename(
            title="Select a FASTA file",
            filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")]
        )
        if not filepath:
            return
        
        self.dna_sequence = parse_fasta(filepath)
        if self.dna_sequence:
            filename = filepath.split('/')[-1]
            self.file_label.config(text=f"Loaded: {filename} ({len(self.dna_sequence)} bp)", fg="black")

    def run_analysis(self):
        """Performs the sliding window analysis and plots the results."""
        if not self.dna_sequence:
            messagebox.showwarning("No Data", "Please load a FASTA file first.")
            return

        try:
            na_concentration = float(self.na_entry.get())
            if na_concentration <= 0:
                raise ValueError("Concentration must be positive.")
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter a valid positive number for Na+ concentration.")
            return

        if len(self.dna_sequence) < self.window_size:
            messagebox.showwarning("Sequence Too Short", f"The sequence is shorter than the window size of {self.window_size} bp.")
            return
            
        # Perform the analysis
        positions, basic_tms, advanced_tms = self.sliding_window_analysis(na_concentration)
        
        # Plot the results
        self.plot_data(positions, basic_tms, advanced_tms)

    def sliding_window_analysis(self, na_concentration: float):
        """Slides a window across the sequence and calculates Tm at each step."""
        positions = []
        basic_tms = []
        advanced_tms = []

        for i in range(len(self.dna_sequence) - self.window_size + 1):
            sub_sequence = self.dna_sequence[i : i + self.window_size]
            
            positions.append(i)
            basic_tms.append(calculate_tm_basic(sub_sequence))
            advanced_tms.append(calculate_tm_advanced(sub_sequence, na_concentration))
            
        return positions, basic_tms, advanced_tms

    def plot_data(self, positions, basic_tms, advanced_tms):
        """Clears the old chart and plots the new data."""
        self.ax.cla() # Clear the previous plot
        self.ax.plot(positions, basic_tms, label="Basic Formula", color="blue", alpha=0.8)
        self.ax.plot(positions, advanced_tms, label="Advanced Formula", color="red", alpha=0.8)
        self.ax.legend()
        self.setup_chart() # Re-apply labels and grid

if __name__ == "__main__":
    app = TmScannerApp()
    app.mainloop()