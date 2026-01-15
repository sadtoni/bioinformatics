# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 14:38:49 2026

@author: Antonio
"""

import random
import json

word_pool = [
    "the", "quick", "brown", "fox", "jumps", "over", "lazy", "dog",
    "alice", "was", "beginning", "to", "get", "very", "tired",
    "of", "sitting", "by", "her", "sister", "on", "bank", "and"
]

text_list = []
current_length = 0
while current_length < 300:
    word = random.choice(word_pool)
    text_list.append(word)
    current_length += len(word) + 1

unique_words = sorted(list(set(text_list)))
word_to_char = {word: chr(33 + i) for i, word in enumerate(unique_words)}

counts = {s: {s2: 0 for s2 in word_to_char.values()} for s in word_to_char.values()}

for i in range(len(text_list) - 1):
    curr_sym = word_to_char[text_list[i]]
    next_sym = word_to_char[text_list[i+1]]
    counts[curr_sym][next_sym] += 1

transition_matrix = {}
for row_sym, transitions in counts.items():
    total = sum(transitions.values())
    transition_matrix[row_sym] = {}
    for col_sym, count in transitions.items():
        if total > 0:
            transition_matrix[row_sym][col_sym] = count / total
        else:
            transition_matrix[row_sym][col_sym] = 0.0

with open('word_transition_matrix.json', 'w') as f:
    json.dump(transition_matrix, f, indent=4)

print("Word to Symbol Mapping:")
print(word_to_char)