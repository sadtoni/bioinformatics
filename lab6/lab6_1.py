# -*- coding: utf-8 -*-
"""
Created on Thu Nov  6 09:56:54 2025

@author: Antonio
"""

import random
import matplotlib.pyplot as plt
import numpy as np
import re

# --- Step 1: Parse the FASTA file ---
# Assumes the file 'covid.fasta' is in the same directory
file_path = 'covid.fasta'

sequence_lines = []
full_genome = ""
genome_length = 0

try:
    with open(file_path, 'r') as f:
        for line in f:
            # Skip the header line
            if not line.startswith('>'):
                sequence_lines.append(line.strip())
    
    # Join all sequence lines into one large string
    full_genome = "".join(sequence_lines)
    genome_length = len(full_genome)
    
    if genome_length == 0:
        print("Warning: Genome sequence is empty. Check FASTA file.")
        # Use a fallback length if parsing failed
        genome_length = 30000 

except Exception as e:
    print(f"Error reading file: {e}. Using a dummy length.")
    # Fallback length in case of file read error
    genome_length = 30000

print(f"Parsed genome. Total length: {genome_length} bp")

# --- Step 2 & 3: Get 10 random samples (lengths) ---
# The prompt asks for 10 samples, 100-3000 bases.
# For this simulation, we only need the *lengths* of the fragments.

fragment_lengths = []
for _ in range(10):
    # Generate a random length between 100 and 3000 (inclusive)
    length = random.randint(100, 3000)
    fragment_lengths.append(length)

# Sort them from largest to smallest to mimic gel loading order (optional)
fragment_lengths.sort(reverse=True)
print(f"Generated 10 random fragment lengths (bp): {fragment_lengths}")

# --- Step 4: Simulate and visualize the gel ---

# Define the DNA ladder (based on a common 1kb ladder)
# These are the sizes (in bp) of the bands in the ladder lane
ladder_bands = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250, 100]

# Define labels to match the user's example image
ladder_labels = {
    3000: '3000 bp -',
    1500: '1500 bp -',
    500: '500 bp -'
}

# --- Create the plot ---
fig, ax = plt.subplots(figsize=(4, 7))
# Set the background color to black, like a gel box
fig.patch.set_facecolor('black')
ax.set_facecolor('black')

# Set up the plot area
ax.set_xlim(0, 1)
# Set Y-axis limits (base pair sizes) with extra space for wells
ax.set_ylim(50, 12000) 
# Use a logarithmic scale for the Y-axis, as DNA migrates logarithmically
ax.set_yscale('log') 
# Invert the Y-axis so large fragments (slow) are at the top
ax.invert_yaxis() 

# --- Draw the lanes (wells) ---
lane_width = 0.3
ladder_x = 0.3 # X-coordinate for the ladder lane
sample_x = 0.7 # X-coordinate for the sample lane
well_y_start = 11000 # Y-coordinate for the top of the well
well_height = 1000  # Height of the well

# Ladder well (a gray box)
well_patch_ladder = plt.Rectangle((ladder_x - lane_width/2, well_y_start), 
                                  lane_width, well_height, 
                                  facecolor='#555555', 
                                  edgecolor='#888888')
ax.add_patch(well_patch_ladder)

# Sample well (a gray box)
well_patch_sample = plt.Rectangle((sample_x - lane_width/2, well_y_start), 
                                  lane_width, well_height, 
                                  facecolor='#555555', 
                                  edgecolor='#888888')
ax.add_patch(well_patch_sample)


# --- Draw Lane 1 (Ladder Bands) ---
for band_size in ladder_bands:
    # Make some bands brighter (e.g., 1000 and 3000)
    lw = 2 # default linewidth
    if band_size == 3000 or band_size == 1000:
        lw = 4 # Brighter band (thicker line)
        
    # Draw a white line for the band
    ax.plot([ladder_x - lane_width/2, ladder_x + lane_width/2], 
            [band_size, band_size], 
            color='white', 
            linewidth=lw, 
            solid_capstyle='butt') # Use square line ends

# --- Draw Lane 2 (Sample Bands) ---
for frag_len in fragment_lengths:
    # Draw a white line for each sample fragment
    ax.plot([sample_x - lane_width/2, sample_x + lane_width/2], 
            [frag_len, frag_len], 
            color='white', 
            linewidth=3, 
            solid_capstyle='butt')

# --- Add Labels ---
for band_size, label_text in ladder_labels.items():
    # Add text to the left of the ladder lane
    ax.text(ladder_x - lane_width/2 - 0.05, # X-position
            band_size,                     # Y-position (matches band)
            label_text,                    # The text itself
            color='white', 
            fontsize=12, 
            fontfamily='sans-serif',
            ha='right',                    # Horizontal alignment
            va='center')                   # Vertical alignment

# --- Clean up the plot (remove axes lines and ticks) ---
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])

# --- Save the final image ---
plt.tight_layout()
plt.savefig('gel_electrophoresis_simulation.png', dpi=150)

print("Generated image: gel_electrophoresis_simulation.png")