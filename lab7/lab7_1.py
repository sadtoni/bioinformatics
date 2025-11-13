# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 14:01:28 2025

@author: Antonio
"""

import sys
import os
from collections import defaultdict

def extract_sequence(fasta_file):
    """
    Reads a FASTA file and returns the concatenated sequence, 
    stripping newline characters and ignoring the header.
    """
    try:
        with open(fasta_file, 'r') as f:
            lines = [line.strip() for line in f if not line.startswith('>')]
        
        # Check if the sequence is empty
        if not lines:
            print(f"Error: No sequence found in {fasta_file}. The file might be empty or improperly formatted.")
            return None
            
        # Concatenate and convert to uppercase for consistency
        sequence = "".join(lines).upper()
        
        # Simple validation for a DNA sequence (A, T, C, G)
        if any(base not in 'ATCG' for base in sequence):
            print("Warning: The sequence contains non-standard DNA bases. Proceeding, but results might be unexpected.")

        return sequence

    except FileNotFoundError:
        print(f"Error: The file '{fasta_file}' was not found.")
        return None
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None

def find_tandem_repeats(sequence, min_len=3, max_len=6):
    """
    Detects tandem repetitions (consecutive identical motifs) 
    in the DNA sequence for motif lengths between min_len (3) and max_len (6).
    The minimum number of repetitions is 2.
    """
    # A dictionary to store results: {motif_length: {motif: [list_of_locations]}}
    # Each list_of_locations holds tuples: (start_index, repeat_count)
    all_repeats = defaultdict(lambda: defaultdict(list))
    
    if not sequence:
        return all_repeats

    seq_len = len(sequence)
    print(f"\nSequence length: {seq_len} nucleotides.")
    
    # Iterate through all required motif lengths (3, 4, 5, 6)
    for k in range(min_len, max_len + 1):
        # Slide a window of size k across the sequence
        for i in range(seq_len - k + 1):
            # The current k-mer (motif)
            motif = sequence[i:i + k]
            
            # The expected location of the second repeat immediately following the first
            next_start_index = i + k
            
            # Check if there's enough sequence left for at least one more repetition (min 2 copies)
            if next_start_index + k <= seq_len:
                next_motif = sequence[next_start_index:next_start_index + k]
                
                # Check for a repetition (a tandem repeat of length k)
                if motif == next_motif:
                    # Found a repeat! Check for a longer run of this same motif
                    repeat_count = 2 # Start with 2 copies (motif + next_motif)
                    current_end = next_start_index + k
                    
                    # Continue checking for additional consecutive repeats
                    while current_end + k <= seq_len and sequence[current_end:current_end + k] == motif:
                        repeat_count += 1
                        current_end += k
                    
                    # Store the repeat if it hasn't been stored yet
                    # We only store the *first* occurrence of a contiguous block
                    start_index = i
                    
                    # Accesses the start_index (the int) of the last stored repeat run
                    if not all_repeats[k][motif] or start_index > all_repeats[k][motif][-1][0]:
                        
                        # Add the start index and the count of repeats
                        all_repeats[k][motif].append((start_index, repeat_count))
    
    return all_repeats

def print_results(results):
    """
    Prints the detected tandem repeats in a readable format.
    """
    has_repeats = False
    for k, motifs in sorted(results.items()):
        if motifs:
            has_repeats = True
            print(f"\n--- 📏 Motifs of Length {k} ---")
            
            # Sort by motif sequence for cleaner output
            for motif, locations in sorted(motifs.items()):
                print(f"**Motif: {motif}**")
                for start, count in locations:
                    end = start + (k * count) - 1
                    print(f"  - Location: {start} to {end} (0-indexed). Repeats: {count} copies. Total length: {k*count}b")
                print("-" * 20) # Separator for motifs within the same k-length

    if not has_repeats:
        print("\n✅ No tandem repeats of lengths 3 to 6 (min 2 copies) were found in the sequence.")
    else:
        print("\n--- Detection Complete ---")


# --- Main Execution Block ---
if __name__ == "__main__":
    FASTA_FILE_NAME = "covid.fasta"
    
    # 1. Read the DNA sequence
    dna_sequence = extract_sequence(FASTA_FILE_NAME)
    
    if dna_sequence:
        # 2. Find the repeats, now explicitly targeting 3 to 6 bases
        # The default arguments in the function definition now handle the 3-6 range.
        repeat_data = find_tandem_repeats(dna_sequence)
        
        # 3. Print the results
        print_results(repeat_data)