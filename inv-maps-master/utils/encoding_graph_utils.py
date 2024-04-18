import copy
import networkx as nx

from scipy.sparse import identity

from .majoranas_utils import avg_pauli_weight_gmatrix


def G_matrix_from_graph(Graph, N=None, cascaded=True):
    if not cascaded:
        Graph = cascade_children(Graph)
    G = nx.adjacency_matrix(Graph)[:N,:N] + identity(N)
    return G


def greedy_remove_edges(Graph, N):
    # NOTE: Requires Graph to be uncascaded
     
    best_weight = avg_pauli_weight_digraph(Graph, N, cascaded=False)
    current_weight = best_weight
    edges_removed = []

    while True:
        for edge in Graph.edges:
            temp = copy.deepcopy(Graph)
            temp.remove_edge(*edge)
            temp = cascade_children(temp)
            
            temp_weight = avg_pauli_weight_digraph(temp, N)
            if temp_weight < best_weight:
                best_weight = temp_weight
                edge_cut = edge
                
        # All edges checked and there was no edge 
        # we could delete to improve avg weight
        if current_weight == best_weight:
            break
        
        edges_removed.append(edge_cut)
        Graph.remove_edge(*edge_cut)
        current_weight = best_weight

    return Graph, current_weight, edges_removed


def cascade_children(Graph):
    # This returns the "transitive closure" of the
    # input DiGraph
    while True:
        cascaded = copy.deepcopy(Graph)
        for u, v in Graph.edges:
            for children in Graph.successors(v):
                cascaded.add_edge(u, children)
            
        # A full iteration and no edges were added
        if cascaded.edges == Graph.edges:
            break
        else:
            Graph = cascaded
    return cascaded


def avg_pauli_weight_digraph(Graph, N, cascaded=True):

    G = G_matrix_from_graph(Graph, N, cascaded)

    return avg_pauli_weight_gmatrix(G, N)