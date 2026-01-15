# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 15:17:29 2026

@author: Antonio
"""

import json
import random
import numpy as np

with open('word_transition_matrix.json', 'r') as f:
    matrix_data = json.load(f)

symbols = sorted(matrix_data.keys())
current_symbol = random.choice(symbols)
generated_sequence = [current_symbol]

for _ in range(20):
    probs = [matrix_data[current_symbol][s] for s in symbols]
    
    if sum(probs) == 0:
        current_symbol = random.choice(symbols)
    else:
        current_symbol = np.random.choice(symbols, p=probs)
    
    generated_sequence.append(current_symbol)

print("Generated Sequence (Symbols):", " ".join(generated_sequence))