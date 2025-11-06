# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 14:40:36 2025

@author: Antonio
"""

import matplotlib.pyplot as plt
import re

def parse_multi_fasta(filename):
    """
    Parses a multi-FASTA file and returns a dictionary of sequences.
    Keys are sequence headers (e.g., 'NC_026431.1')
    Values are the complete sequence strings.
    """
    sequences = {}
    current_sequence_name = None
    current_sequence = []

    print(f"Reading file: {filename}")
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue  # Skip empty lines
                
                if line.startswith('>'):
                    # If we are already building a sequence, save it first
                    if current_sequence_name:
                        sequences[current_sequence_name] = "".join(current_sequence)
                    
                    # Start the new sequence
                    # Clean up the header to get a short, usable name
                    try:
                        # Assumes format like >ID | description
                        current_sequence_name = line.split('|')[0].lstrip('>').strip()
                    except Exception:
                        # Fallback for simple headers like >ID
                        current_sequence_name = line.lstrip('>').strip()
                    
                    current_sequence = []
                else:
                    # This is a sequence line
                    if current_sequence_name:  # Ensure we're inside a sequence block
                        current_sequence.append(line)
        
        # Don't forget to save the very last sequence in the file
        if current_sequence_name:
            sequences[current_sequence_name] = "".join(current_sequence)

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return {}
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return {}
    
    if not sequences:
        print("Warning: No sequences were parsed from the file.")
    
    return sequences

def restriction_digest(sequence, enzyme_site="GAATTC"):
    """
    Performs an in-silico digest on a linear DNA sequence.
    The default site is EcoRI (GAATTC), which cuts after the 'G'.
    """
    # finditer finds all non-overlapping matches
    # We add +1 because EcoRI cuts after the 'G' at position 1
    cut_sites = [match.start() + 1 for match in re.finditer(enzyme_site, sequence.upper())]
    
    # Add the start and end of the sequence to the list of cut locations
    all_cuts = [0] + cut_sites + [len(sequence)]
    
    # Calculate fragment lengths
    fragments = []
    for i in range(1, len(all_cuts)):
        fragment_len = all_cuts[i] - all_cuts[i-1]
        if fragment_len > 0:
            fragments.append(fragment_len)
            
    return fragments

def plot_multi_gel(lanes_data, lane_names, ladder_bands, ladder_labels):
    """
    Plots a multi-lane gel simulation.
    lanes_data: A list of lists, where each inner list contains fragment sizes.
    lane_names: A list of names for each lane.
    """
    num_lanes = len(lanes_data)
    
    # Create a figure wide enough to accommodate all lanes
    # Adjust width dynamically based on the number of lanes
    fig_width = max(8, num_lanes * 1.0 + 2) # +2 for ladder and margins
    fig, ax = plt.subplots(figsize=(fig_width, 8))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')

    # Set up the plot area
    ax.set_xlim(0, num_lanes + 1)
    ax.set_ylim(50, 10000)  # Standard range, log scale
    ax.set_yscale('log')
    ax.invert_yaxis()  # Large fragments at top

    # --- Draw Lane 1 (Ladder) ---
    lane_x = 0.5
    lane_width = 0.8
    # Draw well
    ax.add_patch(plt.Rectangle((lane_x - lane_width/2, 9500), lane_width, 1000, 
                               facecolor='#555', edgecolor='#888'))
    ax.text(lane_x, 11000, "Ladder", color='white', ha='center', fontsize=9)
    
    # Draw ladder bands
    for band_size in ladder_bands:
        lw = 2
        if band_size == 3000 or band_size == 1000:
            lw = 4  # Brighter band
        ax.plot([lane_x - lane_width/2, lane_x + lane_width/2], 
                [band_size, band_size], color='white', linewidth=lw, solid_capstyle='butt')

    # --- Draw Sample Lanes ---
    for i, (fragments, name) in enumerate(zip(lanes_data, lane_names)):
        lane_x = i + 1.5  # Position each sample lane (starting from 1.5)
        # Draw well
        ax.add_patch(plt.Rectangle((lane_x - lane_width/2, 9500), lane_width, 1000, 
                                   facecolor='#555', edgecolor='#888'))
        # Draw lane name (rotated for clarity)
        ax.text(lane_x, 11000, name, color='white', ha='center', fontsize=8, rotation=45, va='bottom')
        
        # Draw fragment bands
        for frag_len in fragments:
            ax.plot([lane_x - lane_width/2, lane_x + lane_width/2], 
                    [frag_len, frag_len], color='white', linewidth=3, solid_capstyle='butt')

    # --- Add Labels ---
    for band_size, label_text in ladder_labels.items():
        # Position labels to the left of the ladder
        ax.text(0.1, band_size, label_text, color='white', fontsize=10, ha='right', va='center')

    # --- Clean up axes ---
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    ax.set_title("Simulated EcoRI Restriction Digest", color='white', fontsize=16)

    plt.tight_layout()
    plt.savefig('multi_gel_simulation.png', dpi=150)
    print(f"\nSuccessfully generated image: multi_gel_simulation.png")

# --- Main script execution ---
def main():
    # 1. DEFINE FILES AND PARAMETERS
    fasta_file = "sequences.fasta"
    enzyme_site = "GAATTC"  # EcoRI
    
    # Define the DNA ladder (based on a common 1kb ladder)
    ladder_bands = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250, 100]
    ladder_labels = {
        3000: '3000 bp -',
        1500: '1500 bp -',
        500: '500 bp -'
    }

    # 2. PROCESS FILE AND PERFORM DIGEST
    sequences = parse_multi_fasta(fasta_file)
    if not sequences:
        print("No sequences found in file. Exiting.")
        return

    all_fragments = []
    lane_names = []
    digest_results = {}
    
    print(f"Found {len(sequences)} sequences. Starting in-silico restriction digest...")
    
    for name, seq in sequences.items():
        print(f"  Digesting {name} (Length: {len(seq)} bp)...")
        fragments = restriction_digest(seq, enzyme_site)
        
        # Store results for plotting and analysis
        all_fragments.append(fragments)
        lane_names.append(name)
        digest_results[name] = fragments

    # 3. ANALYZE AND COMPARE
    print("\n--- Analysis Results ---")
    most_fragments = -1
    genome_with_most_fragments = "None"
    
    for name, frags in digest_results.items():
        num_fragments = len(frags)
        print(f"  {name}: {num_fragments} segments")
        
        if num_fragments > most_fragments:
            most_fragments = num_fragments
            genome_with_most_fragments = name
            
    if most_fragments == -1:
        most_fragments = 0 # Handle case where no fragments were found
    
    print("------------------------")
    print(f"**Sequence with the most DNA segments:** {genome_with_most_fragments} ({most_fragments} segments)")
    print("------------------------")

    # 4. PLOT THE GEL
    plot_multi_gel(all_fragments, lane_names, ladder_bands, ladder_labels)

if __name__ == "__main__":
    main()