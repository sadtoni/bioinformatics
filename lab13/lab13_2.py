# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 14:38:39 2026

@author: Antonio
"""

import random
import json

nucleotides = ['A', 'C', 'G', 'T']
sequence = ''.join(random.choices(nucleotides, k=50))

counts = {n1: {n2: 0 for n2 in nucleotides} for n1 in nucleotides}

for i in range(len(sequence) - 1):
    current_char = sequence[i]
    next_char = sequence[i + 1]
    counts[current_char][next_char] += 1

transition_matrix = {}

for n1 in nucleotides:
    total_transitions = sum(counts[n1].values())
    transition_matrix[n1] = {}
    for n2 in nucleotides:
        if total_transitions > 0:
            transition_matrix[n1][n2] = counts[n1][n2] / total_transitions
        else:
            transition_matrix[n1][n2] = 0.0

with open('transition_matrix.json', 'w') as f:
    json.dump(transition_matrix, f, indent=4)