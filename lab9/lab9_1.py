# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 11:50:33 2025

@author: Antonio
"""

import sys

def main():
    try:
        sequence = ""
        with open("bacteria.fasta", "r") as f:
            for line in f:
                if not line.startswith(">"):
                    sequence += line.strip()
    except FileNotFoundError:
        print("Error: bacteria.fasta file not found.")
        return

    seq_len = len(sequence)
    print(f"DNA Sequence Length: {seq_len} bp\n")

    enzymes = [
        ("EcoRI", "GAATTC", 1),
        ("BamHI", "GGATCC", 1),
        ("HindIII", "AAGCTT", 1),
        ("TaqI", "TCGA", 1),
        ("HaeIII", "GGCC", 2)
    ]

    enzyme_results = {}

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
        
        enzyme_results[name] = {
            "cuts": sorted_cuts,
            "fragments": sorted(fragments, reverse=True)
        }

        print(f"--- {name} ---")
        print(f"Number of cleavages: {len(sorted_cuts)}")
        print(f"Cleavage positions: {sorted_cuts}")
        print(f"Fragment lengths: {enzyme_results[name]['fragments']}")
        print("")

    print("--- Simulated Electrophoresis Gel ---")
    print("Base Pairs (bp) | " + " | ".join(f"{name:^9}" for name, _, _ in enzymes))
    print("-" * (18 + 12 * len(enzymes)))

    step = seq_len // 40
    if step < 50: step = 50

    for marker in range(seq_len, 0, -step):
        line = f"{marker:<15} | "
        for name, _, _ in enzymes:
            frags = enzyme_results[name]["fragments"]
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