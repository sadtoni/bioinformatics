# -*- coding: utf-8 -*-
"""
Created on Thu Jan 22 14:43:22 2026

@author: Antonio
"""

import math
import re
from collections import defaultdict

def get_tokens(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        text = f.read().lower()
    text = re.sub(r'[^\w\s]', '', text)
    return text.split()

def train_model(tokens):
    transitions = defaultdict(lambda: defaultdict(int))
    counts = defaultdict(int)
    
    for i in range(len(tokens) - 1):
        prev, curr = tokens[i], tokens[i+1]
        transitions[prev][curr] += 1
        counts[prev] += 1
        
    probs = defaultdict(dict)
    for prev, followers in transitions.items():
        for curr, count in followers.items():
            probs[prev][curr] = count / counts[prev]
    return probs

def get_log_likelihood(w1, w2, m_eminescu, m_stanescu):
    p_e = m_eminescu.get(w1, {}).get(w2, 0)
    p_s = m_stanescu.get(w1, {}).get(w2, 0)
    
    if p_e == 0 and p_s == 0:
        return 0
    if p_e == 0:
        return -10 
    if p_s == 0:
        return 10
        
    return math.log2(p_e / p_s)

tokens_eminescu = get_tokens('eminescu.txt')
tokens_stanescu = get_tokens('stanescu.txt')
tokens_suspicious = get_tokens('combinatie.txt')

model_eminescu = train_model(tokens_eminescu)
model_stanescu = train_model(tokens_stanescu)

window_size = 5
print(f"{'Index':<10} {'Score':<10} {'Attribution':<15} {'Text Segment'}")
print("-" * 80)

for i in range(len(tokens_suspicious) - window_size):
    window = tokens_suspicious[i : i + window_size]
    score = 0
    
    for j in range(len(window) - 1):
        w1, w2 = window[j], window[j+1]
        score += get_log_likelihood(w1, w2, model_eminescu, model_stanescu)
    
    attribution = "NEUTRAL"
    if score > 0:
        attribution = "EMINESCU"
    elif score < 0:
        attribution = "STANESCU"
        
    segment = " ".join(window)
    print(f"{i:<10} {score:<10.2f} {attribution:<15} {segment}")