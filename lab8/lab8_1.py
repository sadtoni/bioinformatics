# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 14:14:11 2025

@author: Antonio
"""

import random

# ==========================================
# PART 1: SIMULATION & BIOLOGICAL MODELING
# ==========================================

def generate_random_dna(length):
    return "".join(random.choice("ACGT") for _ in range(length))

class Transposon:
    def __init__(self, name, internal_seq, ir_seq):
        self.name = name
        self.ir_seq = ir_seq  # Inverted Repeat
        # In real biology, the right IR is usually the reverse complement of the left
        # For this simulation, we will simplify and use the string reverse or identical for detection ease
        self.ir_seq_rev = ir_seq[::-1] 
        self.internal = internal_seq
        # The full transposon is IR + Internal + Reverse IR
        self.full_sequence = self.ir_seq + self.internal + self.ir_seq_rev
        self.length = len(self.full_sequence)

def insert_transposon(host_dna, transposon, position):
    """
    Simulates insertion as shown in the image.
    1. Cuts host at position.
    2. Generates Target Site Duplication (TSD) (The 'Direct Repeats' in the image).
    3. Inserts Transposon.
    """
    # Biology logic: TSDs are usually 2-8bp of the host duplicated at both ends
    tsd_len = 4 
    
    # If position is too close to end/start, adjust
    if position < tsd_len or position > len(host_dna) - tsd_len:
        position = len(host_dna) // 2
        
    # The TSD is the sequence at the insertion site that gets duplicated
    tsd_seq = host_dna[position:position+tsd_len]
    
    print(f"--- Inserting {transposon.name} at index {position} ---")
    print(f"    Target Site Duplication (TSD) created: {tsd_seq}")
    
    # Construct new DNA: Up to pos + TSD + Transposon + TSD + Rest of DNA
    # Note: In the image, the TSD appears on both sides of the transposon.
    
    left_flank = host_dna[:position+tsd_len] # Includes the sequence to be duplicated
    right_flank = host_dna[position+tsd_len:]
    
    # The simulation result: 
    # ...Host... [TSD] [TRANSPOSON] [TSD] ...Host...
    new_dna = host_dna[:position] + tsd_seq + transposon.full_sequence + tsd_seq + host_dna[position+tsd_len:]
    
    return new_dna

# ==========================================
# PART 2: DETECTION ALGORITHM
# ==========================================

def detect_transposons(genome, library):
    """
    Scans the genome for the transposons in the library.
    Handles contiguous TEs and Intersected (Split) TEs.
    """
    print("\n" + "="*40)
    print("RUNNING DETECTION ALGORITHM")
    print("="*40)
    
    results = []

    for te in library:
        # Strategy 1: Exact Match (The TE is intact)
        # We search for the full sequence in the genome
        start_idx = genome.find(te.full_sequence)
        
        if start_idx != -1:
            end_idx = start_idx + len(te.full_sequence)
            results.append({
                "name": te.name,
                "status": "Intact",
                "start": start_idx,
                "end": end_idx,
                "details": "Found contiguous sequence."
            })
        else:
            # Strategy 2: Split Match (Intersection/Nesting detection)
            # If the full sequence isn't found, it might have been disrupted by another TE.
            # We look for the "Heads" (Left IR) and "Tails" (Right IR) separately.
            
            head_seq = te.ir_seq + te.internal[:5] # Use IR + bit of internal as a signature
            tail_seq = te.internal[-5:] + te.ir_seq_rev
            
            head_pos = genome.find(head_seq)
            tail_pos = genome.find(tail_seq)
            
            if head_pos != -1 and tail_pos != -1:
                # Calculate where the tail ends
                tail_end_pos = tail_pos + len(tail_seq)
                
                # Check the distance. If distance > original length, something was inserted inside.
                observed_len = tail_end_pos - head_pos
                gap_size = observed_len - te.length
                
                if gap_size > 0:
                    results.append({
                        "name": te.name,
                        "status": "Intersected/Nested",
                        "start": head_pos,
                        "end": tail_end_pos,
                        "details": f"TE was split! Gap of {gap_size}bp detected inside."
                    })
                else:
                     results.append({
                        "name": te.name,
                        "status": "Fragmented",
                        "start": head_pos,
                        "end": tail_end_pos,
                        "details": "Found ends, but length is inconsistent."
                    })
            else:
                 results.append({
                        "name": te.name,
                        "status": "Not Found",
                        "details": "Significant fragments missing."
                    })

    return results

# ==========================================
# MAIN EXECUTION
# ==========================================

# 1. Define Library of Transposons (Distinct signatures for detection)
te_library = [
    Transposon("Tn_Alpha", "CCCCCCCCCC", "TTTT"),  # Poly-C internal, Poly-T repeats
    Transposon("Tn_Beta",  "GGGGGGGGGG", "AAAA"),  # Poly-G internal, Poly-A repeats
    Transposon("Tn_Gamma", "ATATATATAT", "CGCG")   # AT internal, GC repeats
]

# 2. Create Host Genome (300bp)
genome = generate_random_dna(300)
print(f"Original Host DNA Length: {len(genome)}bp")

# 3. Perform Insertions
# Insertion 1: Tn_Alpha
genome = insert_transposon(genome, te_library[0], 50)

# Insertion 2: Tn_Beta (Placed far away)
genome = insert_transposon(genome, te_library[1], 200)

# Insertion 3: Tn_Gamma (INTERSECTION CASE)
# We deliberately insert Gamma *inside* Alpha.
# Alpha is at 50. Alpha length is 4 (IR) + 10 (Int) + 4 (IR) = 18bp.
# We insert Gamma at index 60 (inside Alpha).
genome = insert_transposon(genome, te_library[2], 60)

print(f"Final Mutated Genome Length: {len(genome)}bp")

# 4. Detect
detections = detect_transposons(genome, te_library)

# 5. Report
print("\nDetection Report:")
print(f"{'TE Name':<10} | {'Status':<20} | {'Position (Start-End)':<20} | {'Notes'}")
print("-" * 80)
for res in detections:
    pos_str = f"{res.get('start', 'N/A')}-{res.get('end', 'N/A')}"
    print(f"{res['name']:<10} | {res['status']:<20} | {pos_str:<20} | {res['details']}")

print("-" * 80)