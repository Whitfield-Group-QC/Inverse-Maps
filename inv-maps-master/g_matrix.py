"""
User level functions for dealing with encodings

"""

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
    # TODO
    pass