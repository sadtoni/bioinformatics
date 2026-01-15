# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 14:01:58 2026

@author: Antonio
"""

import numpy as np

def predict_steps(matrix, vector):
    A = np.array(matrix)
    v = np.array(vector)

    for i in range(5):
        v = A @ v
        print(f"Step {i + 1}: {v}")

initial_matrix = [
    [0.8, 0.3, 0.1],
    [0.1, 0.6, 0.1],
    [0.1, 0.1, 0.8]
]

initial_vector = [10.0, 5.0, 20.0]

predict_steps(initial_matrix, initial_vector)