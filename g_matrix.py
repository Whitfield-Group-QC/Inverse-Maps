"""
User level functions for dealing with encodings

"""

import numpy as np
from graph_generators.sierpinski import unfilled_sierpinski_with_children
from utils.encoding_graph_utils import G_matrix_from_graph, greedy_remove_edges

def pruned_sierpinski_G_matrix(N):
    return G_matrix_from_graph(unfilled_sierpinski_with_children(N), N)


def optimized_sierpinski_G_matrix(N):
    # This is equivalent to Ternary
    Graph = unfilled_sierpinski_with_children(N)
    Graph_opt, _, _ = greedy_remove_edges(Graph, N)

    return G_matrix_from_graph(Graph_opt, N)


def fenwick_G_matrix(N):
    # TODO
    pass

def random_enc_G_matrix(N):
    '''
    Makes random matrix in GLNF2: Corresponds to random encoding via 
    conjugating Jordan Wigner by CNOTS
    '''
    mat = np.random.randint(0,2,(N,N))
    while int(np.linalg.det(mat)) % 2 != 1:
        mat = np.random.randint(0,2,(N,N))
    return mat