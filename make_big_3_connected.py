#!/usr/bin/python3
import igraph as ig
import math
import re
import random
import matplotlib.pyplot as plt
import os
import pdb
from sympy import geometry as geo
import cs695util

def main():
    n = 10

    input_filename  = 'triangle_graph.txt'
    output_filename = 'triangle_graph-{0}.txt'.format(6 * (n+1))
    small_g = cs695util.load_tsv_edges("data/" + input_filename, False)

    new_edges = [[e.source, e.target] for e in small_g.es]

    to_add = list(new_edges)

    for i in range(n):
        next_new = max(set([v[0] for v in new_edges] +
                           [v[1] for v in new_edges])) + 1
        new_edges.extend([[x[0] + next_new, x[1] + next_new] for x in to_add])
        next_old = next_new - 3
        for j in range(3):
            new_edges.append([next_old, next_new])
            next_old += 1
            next_new += 1

    with open(os.path.join('data', output_filename), 'w') as out:
        for e in new_edges:
            out.write('{0} {1}\n'.format(*e))

if __name__=="__main__":
    main()
    
