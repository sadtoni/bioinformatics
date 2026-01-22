import math

def get_index(base):
    if base == 'A': return 0
    if base == 'C': return 1
    if base == 'G': return 2
    if base == 'T': return 3
    return -1

def count_transitions(seq):
    counts = [[0.0 for _ in range(4)] for _ in range(4)]
    for i in range(len(seq) - 1):
        idx_from = get_index(seq[i])
        idx_to = get_index(seq[i+1])
        if idx_from != -1 and idx_to != -1:
            counts[idx_from][idx_to] += 1
    return counts

def to_probabilities(counts):
    probs = [[0.0 for _ in range(4)] for _ in range(4)]
    for i in range(4):
        row_sum = sum(counts[i])
        if row_sum > 0:
            for j in range(4):
                probs[i][j] = counts[i][j] / row_sum
    return probs

s1 = "ATCGATTCGATATCATACACGTAT"
s2 = "CTCGACTAGTATGAAGTCCACGCTTG"
s_test = "CAGGTTGGAAACGTAA"
bases = ['A', 'C', 'G', 'T']

c1 = count_transitions(s1)
c2 = count_transitions(s2)

p1 = to_probabilities(c1)
p2 = to_probabilities(c2)

log_matrix = [[0.0 for _ in range(4)] for _ in range(4)]
epsilon = 1e-10

print("Log-Likelihood Matrix:")
print("      A      C      G      T")
for i in range(4):
    row_str = f"{bases[i]}  "
    for j in range(4):
        p_plus = p1[i][j]
        p_minus = p2[i][j]
        
        if p_plus == 0: p_plus = epsilon
        if p_minus == 0: p_minus = epsilon
        
        val = math.log(p_plus / p_minus) / math.log(2)
        log_matrix[i][j] = val
        row_str += f"{val:6.2f} "
    print(row_str)

total_score = 0
print("\nTesting Sequence:", s_test)
for i in range(len(s_test) - 1):
    idx_from = get_index(s_test[i])
    idx_to = get_index(s_test[i+1])
    score = log_matrix[idx_from][idx_to]
    total_score += score

print(f"Total Log-Likelihood Score: {total_score:.4f}")

if total_score > 0:
    print("Result: Sequence belongs to a CpG island.")
else:
    print("Result: Sequence does NOT belong to a CpG island.")