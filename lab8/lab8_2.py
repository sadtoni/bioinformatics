import sys
import time

def get_reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def read_fasta(file_path):
    """
    Reads a FASTA file and yields (header, sequence) tuples.
    Handles multi-line FASTA formatting.
    """
    with open(file_path, 'r') as f:
        header = None
        sequence = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(sequence)
                header = line[1:] # Remove >
                sequence = []
            else:
                sequence.append(line)
        if header:
            yield header, "".join(sequence)

def find_inverted_repeats(sequence, min_k=4, max_k=6, min_gap=100, max_gap=5000):
    """
    Detects potential Transposable Elements (TEs) based on Inverted Repeats (IRs).
    
    Args:
        sequence (str): The genomic sequence.
        min_k (int): Minimum size of IR (4 bases).
        max_k (int): Maximum size of IR (6 bases).
        min_gap (int): Minimum distance between IRs (spacer length).
        max_gap (int): Maximum distance between IRs.
        
    Returns:
        list: A list of dictionaries containing TE details.
    """
    candidates = []
    seq_len = len(sequence)
    sequence = sequence.upper()

    # Iterate through allowed IR lengths (4, 5, 6)
    for k in range(min_k, max_k + 1):
        # Optimization: Build a hash map (dictionary) of k-mer positions
        # kmer_positions[sequence_string] = [list of start indices]
        kmer_positions = {}
        
        for i in range(seq_len - k + 1):
            kmer = sequence[i : i+k]
            # Skip sequences with N (unknown bases)
            if 'N' in kmer:
                continue
                
            if kmer not in kmer_positions:
                kmer_positions[kmer] = []
            kmer_positions[kmer].append(i)

        # Search for IR pairs using the hash map
        # logic: for every kmer, calculate its RC. If RC exists in map, check positions.
        
        # We use a set to avoid processing the same pair twice (e.g. Palindromes or self-matches)
        checked_kmers = set()

        for kmer, start_indices in kmer_positions.items():
            if kmer in checked_kmers:
                continue
            
            rc_kmer = get_reverse_complement(kmer)
            
            # If the Reverse Complement exists in the genome
            if rc_kmer in kmer_positions:
                end_indices = kmer_positions[rc_kmer]
                
                # Check every combination of start (Left IR) and end (Right IR)
                # This O(N*M) loop is constrained by max_gap, so it's efficient in practice.
                for s_idx in start_indices:
                    for e_idx in end_indices:
                        # Basic validation: 
                        # 1. The Right IR (e_idx) must be after the Left IR (s_idx)
                        # 2. They must not be the exact same physical atoms (though start < end handles this)
                        
                        if e_idx > s_idx:
                            # Calculate distance between the END of Left IR and START of Right IR
                            # Structure: [IR Left] --- spacer --- [IR Right]
                            # s_idx is start of Left. e_idx is start of Right.
                            
                            # Distance calculation:
                            # Left IR ends at: s_idx + k
                            # Right IR starts at: e_idx
                            dist = e_idx - (s_idx + k)
                            
                            if min_gap <= dist <= max_gap:
                                # Calculate total TE length (start of Left to end of Right)
                                total_len = (e_idx + k) - s_idx
                                
                                candidates.append({
                                    'start': s_idx + 1, # 1-based indexing for biology standard
                                    'end': e_idx + k,   # 1-based inclusive
                                    'length': total_len,
                                    'ir_seq': kmer,
                                    'rc_seq': rc_kmer,
                                    'ir_len': k,
                                    'gap_size': dist
                                })
            
            checked_kmers.add(kmer)

    # Sort candidates by start position
    candidates.sort(key=lambda x: x['start'])
    return candidates

def analyze_overlaps(candidates):
    """
    Analyzes the raw candidates to specifically highlight nesting and overlaps.
    (Note: The raw detection already includes them, this just tags them for display).
    """
    results = []
    for i, current in enumerate(candidates):
        status = "Distinct"
        
        # Check against other candidates
        # This is a simplistic O(N^2) check suitable for report generation
        for j, other in enumerate(candidates):
            if i == j: continue
            
            # Check for Nesting (Involving)
            if current['start'] > other['start'] and current['end'] < other['end']:
                status = f"Nested inside TE at {other['start']}"
                break
            
            # Check for Overlap (partial)
            elif (current['start'] < other['end'] and current['start'] > other['start']) or \
                 (current['end'] > other['start'] and current['end'] < other['end']):
                status = f"Overlaps with TE at {other['start']}"
                break
        
        current['status'] = status
        results.append(current)
    return results

def main():
    input_file = "bacteria.fasta"
    
    print(f"--- Processing {input_file} ---")
    print("Constraints: IR 4-6bp | Spacing 100-5000bp")
    print("----------------------------------------------------")

    try:
        for header, sequence in read_fasta(input_file):
            genome_id = header.split()[0] # Take first word of header
            print(f"\nAnalyzing Genome: {genome_id} (Length: {len(sequence)} bp)")
            
            start_time = time.time()
            
            # Detect IRs
            raw_candidates = find_inverted_repeats(sequence)
            
            # Tag overlaps/nesting
            final_candidates = analyze_overlaps(raw_candidates)
            
            elapsed = time.time() - start_time
            
            # Output results
            if not final_candidates:
                print("No transposons detected matching criteria.")
            else:
                print(f"Found {len(final_candidates)} potential elements in {elapsed:.2f}s")
                print(f"{'Start':<10} {'End':<10} {'Length':<8} {'IR Seq':<10} {'Status'}")
                print("-" * 60)
                
                # Limit output to first 50 to prevent flooding console if too many hits
                # (4bp IRs can generate many random hits)
                limit = 50 
                for i, te in enumerate(final_candidates):
                    if i >= limit:
                        print(f"... and {len(final_candidates) - limit} more.")
                        break
                    print(f"{te['start']:<10} {te['end']:<10} {te['length']:<8} {te['ir_seq']:<10} {te['status']}")

    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found. Please upload bacteria.fasta.")

if __name__ == "__main__":
    main()