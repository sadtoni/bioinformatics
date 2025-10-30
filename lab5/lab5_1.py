import random

def parse_fasta(filename="covid.fasta"):
    """
    Parses a FASTA file and returns the first sequence as a single string.
    Removes newline characters.
    """
    print(f"--- 1. Parsing {filename} ---")
    sequence = ""
    try:
        with open(filename, 'r') as f:
            # Skip the header line
            header = f.readline()
            if not header.startswith(">"):
                print("Error: This does not appear to be a valid FASTA file.")
                return None
            print(f"Found sequence: {header.strip()}")
            
            # Read the rest of the file and build the sequence
            for line in f:
                sequence += line.strip()
                
    except FileNotFoundError:
        print(f"Error: '{filename}' not found.")
        print("Please make sure 'covid.fasta' is in the same directory as this script.")
        return None
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None
        
    if not sequence:
        print("Error: No sequence data found in the file.")
        return None
        
    print(f"Successfully parsed {len(sequence)} total bases.\n")
    return sequence

def take_samples(sequence, num_samples, min_len, max_len):
    """
    Takes random samples (reads) from a given sequence.
    """
    print(f"--- 2. Taking {num_samples} random samples ---")
    samples = []
    seq_len = len(sequence)
    
    for _ in range(num_samples):
        # 1. Choose a random length for this sample
        sample_len = random.randint(min_len, max_len)
        
        # 2. Choose a random starting position
        # Ensure the sample doesn't go past the end of the sequence
        max_start_index = seq_len - sample_len
        if max_start_index < 0:
            print("Warning: Sequence is shorter than sample length.")
            continue
            
        start_index = random.randint(0, max_start_index)
        
        # 3. Extract the sample and add to our list
        sample = sequence[start_index : start_index + sample_len]
        samples.append(sample)
        
    print(f"Generated {len(samples)} samples (e.g., '{samples[0]}...').\n")
    return samples

def rebuild_sequence_greedy(samples, min_overlap):
    """
    Attempts to rebuild the sequence using a simple greedy overlap algorithm.
    This demonstrates the flaws of the approach.
    """
    print(f"--- 3. Attempting to rebuild sequence (min_overlap={min_overlap}) ---")
    
    if not samples:
        print("No samples to rebuild from.")
        return ""

    # 1. Start our assembly with the first sample and remove it from the pool
    # We use a set for efficient removal of used indices
    rebuilt_sequence = samples[0]
    used_indices = {0}
    
    # We'll try to extend the right side of our sequence
    while True:
        best_overlap_len = -1
        best_sample_index = -1
        
        # This is the "tail" we are trying to find an overlap for
        # We only check the end of our current sequence for efficiency
        search_tail = rebuilt_sequence[-150:] 
        
        # 2. Iterate through all *unused* samples
        for i in range(len(samples)):
            if i in used_indices:
                continue
                
            current_sample = samples[i]
            
            # 3. Check for the best possible overlap with this sample
            # We check from the largest possible overlap down to the minimum
            for overlap_len in range(min(len(search_tail), len(current_sample)), min_overlap - 1, -1):
                
                if search_tail.endswith(current_sample[:overlap_len]):
                    # Found a valid overlap!
                    # Is it the best one we've found *in this pass*?
                    if overlap_len > best_overlap_len:
                        best_overlap_len = overlap_len
                        best_sample_index = i
                    
                    # We found the longest overlap for *this sample*,
                    # so we can stop checking shorter overlaps for it.
                    break 
        
        # 4. Did we find any valid overlap in this entire pass?
        if best_sample_index != -1:
            # Yes! Add the new, non-overlapping part to our sequence
            best_sample = samples[best_sample_index]
            rebuilt_sequence += best_sample[best_overlap_len:]
            
            # Mark this sample as used
            used_indices.add(best_sample_index)
        else:
            # No. We are "stuck". No more samples overlap with the end.
            # This is a "coverage gap" or an unresolvable repeat.
            print("Assembly stuck. No more forward overlaps found.")
            break
            
    print(f"Assembly finished. Found {len(used_indices)} total overlapping samples.\n")
    return rebuilt_sequence

def main():
    # --- Parameters ---
    FASTA_FILE = "covid.fasta"
    SEQUENCE_LENGTH = 3000
    NUM_SAMPLES = 2000
    MIN_SAMPLE_LEN = 100
    MAX_SAMPLE_LEN = 150
    MIN_OVERLAP = 10
    
    # --- Step 1: Parse and Shorten ---
    original_full_sequence = parse_fasta(FASTA_FILE)
    if not original_full_sequence:
        return
        
    original_sequence = original_full_sequence[:SEQUENCE_LENGTH]
    
    # --- Step 2: Take Samples ---
    samples = take_samples(original_sequence, NUM_SAMPLES, MIN_SAMPLE_LEN, MAX_SAMPLE_LEN)
    
    # --- Step 3: Rebuild ---
    # We shuffle the samples to better simulate the "random"
    # nature of which sample we pick first.
    random.shuffle(samples) 
    rebuilt_sequence = rebuild_sequence_greedy(samples, MIN_OVERLAP)
    
    # --- Step 4: Compare Results ---
    print("--- 4. Final Comparison ---")
    print(f"Original Sequence Length: {len(original_sequence)}")
    print(f"Rebuilt Sequence Length:  {len(rebuilt_sequence)}")
    
    print("\nOriginal Sequence (first 100 bases):")
    print(original_sequence[:100])
    
    print("\nRebuilt Sequence (first 100 bases):")
    print(rebuilt_sequence[:100])
    
    if len(original_sequence) != len(rebuilt_sequence) or original_sequence != rebuilt_sequence:
        print("\n*** RESULT: The rebuilt sequence DOES NOT match the original. ***")
        print("This demonstrates the algorithm's failure due to coverage gaps and/or repeat ambiguity.")
    else:
        print("\n*** RESULT: The rebuilt sequence MATCHES the original. ***")
        print("This is highly unlikely, but possible with a lucky random seed.")


if __name__ == "__main__":
    main()
