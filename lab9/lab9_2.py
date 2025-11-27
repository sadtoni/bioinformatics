# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 15:16:53 2025

@author: Antonio
"""

import sys

def main():
    sequences = {}
    current_header = None
    try:
        with open("sequences.fasta", "r") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    current_header = line[1:]
                    sequences[current_header] = ""
                else:
                    if current_header:
                        sequences[current_header] += line
    except FileNotFoundError:
        print("Error: sequences.fasta not found.")
        return

    enzymes = [
        ("EcoRI", "GAATTC", 1),
        ("BamHI", "GGATCC", 1),
        ("HindIII", "AAGCTT", 1),
        ("TaqI", "TCGA", 1),
        ("HaeIII", "GGCC", 2)
    ]

    all_strains_data = {}
    
    for header, sequence in sequences.items():
        strain_data = {}
        seq_len = len(sequence)
        
        for name, site, offset in enzymes:
            cut_positions = []
            search_start = 0
            while True:
                idx = sequence.find(site, search_start)
                if idx == -1:
                    break
                cut_positions.append(idx + offset)
                search_start = idx + 1
            
            fragments = []
            last_cut = 0
            sorted_cuts = sorted(cut_positions)
            for cut in sorted_cuts:
                fragments.append(cut - last_cut)
                last_cut = cut
            fragments.append(seq_len - last_cut)
            
            strain_data[name] = sorted(fragments, reverse=True)
        all_strains_data[header] = strain_data

    common_bands = {name: set() for name, _, _ in enzymes}
    
    for name, _, _ in enzymes:
        first_strain = True
        for header in sequences:
            current_frags = set(all_strains_data[header][name])
            if first_strain:
                common_bands[name] = current_frags
                first_strain = False
            else:
                common_bands[name] = common_bands[name].intersection(current_frags)

    merged_bands = {name: [] for name, _, _ in enzymes}
    max_frag_len = 0

    for header in sequences:
        for name, _, _ in enzymes:
            for frag in all_strains_data[header][name]:
                if frag > max_frag_len:
                    max_frag_len = frag
                if frag not in common_bands[name]:
                    merged_bands[name].append(frag)

    print(f"Analyzed {len(sequences)} sequences.")
    print("Generating Differential Electrophoresis Gel (Conserved bands removed)...")
    print("")

    print("BP       | " + " | ".join(f"{name:^9}" for name, _, _ in enzymes))
    print("-" * (11 + 12 * len(enzymes)))

    step = 50
    if max_frag_len > 2000: step = 100
    if max_frag_len > 5000: step = 200

    for marker in range(max_frag_len + step, 0, -step):
        line = f"{marker:<8} | "
        for name, _, _ in enzymes:
            frags = merged_bands[name]
            has_band = False
            for frag in frags:
                if marker - (step/2) <= frag < marker + (step/2):
                    has_band = True
                    break
            symbol = "=========" if has_band else "         "
            line += f"{symbol} | "
        print(line)

if __name__ == "__main__":
    main()