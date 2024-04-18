import networkx as nx


def full_sierpinski(depth):
    G = nx.DiGraph()
    N = pow(3, depth)
    G.add_nodes_from(range(N))

    sierpinski_recursive(G, 0, N-1)

    return G


def unfilled_sierpinski(N):
    depth = 1
    while pow(3, depth) < N:
        depth += 1
    
    G = full_sierpinski(depth)
    for l in range(pow(3, depth)-1, N-1, -1):
        G.remove_node(l)

    return G


def prune_right_branches(depth):
    G = nx.DiGraph()
    N = pow(3, depth)
    G.add_nodes_from(range(N))
    

    root = int((N - 1)/2)
    hard_code = []
    for k in range(depth-1, -1, -1):
        root += pow(3, k)
        hard_code.append(root)
    
    sierpinski_recursive(G, 0, N - 1, hard_coded_roots=hard_code)

    return G


def unfilled_sierpinski_pruned(N):
    depth = 1
    while pow(3, depth) < N:
        depth += 1
    
    G = prune_right_branches(depth)
    for l in range(pow(3, depth)-1, N-1, -1):
        G.remove_node(l)

    return G


def full_sierpinski_with_children(depth):
    G = nx.DiGraph()
    N = pow(3, depth)
    G.add_nodes_from(range(N))

    sierpinski_with_children(G, 0, N-1)

    return G


def unfilled_sierpinski_with_children(N):
    depth = 1
    while pow(3, depth) < N:
        depth += 1
    
    G = prune_right_branches_with_children(depth)
    for l in range(pow(3, depth)-1, N-1, -1):
        G.remove_node(l)

    return G
    

def prune_right_branches_with_children(depth):
    G = nx.DiGraph()
    N = pow(3, depth)
    G.add_nodes_from(range(N))
    

    root = int((N - 1)/2)
    hard_code = []
    for k in range(depth-1, -1, -1):
        root += pow(3, k)
        hard_code.append(root)
    
    sierpinski_with_children(G, 0, N - 1, hard_coded_roots=hard_code)

    return G


def sierpinski_recursive(G, S, E, hard_coded_roots=[]):
    """
    args: 
        G (nx.DiGraph): The Graph object we'll be updating
        S (int):        The index of the starting node
        E (int):        The index of the ending node   

    NOTE:
        This adds the condition that if there exists an edge
        i -> j,
        there must also be edges from i to every child of j
    """
    if S != E:
        L = int( S + ((E-S-2)/6 ) ) # "left" point, midpoint of first third of [S,E]
        C = int( (S+E)/2 )          # "center" point, midpoint of second third
        R = int( E - (L-S) )        # "right" point, midpoint of last third
        G.add_edge(C, L)
        if R not in hard_coded_roots:
            G.add_edge(C, R)
        if L not in hard_coded_roots:
            G.add_edge(C, L)


        # Divide interval into thirds, apply function to each third
        T = 2*(L-S)
        sierpinski_recursive(G, S, S+T, hard_coded_roots=hard_coded_roots)
        sierpinski_recursive(G, S+T+1, S+2*T+1, hard_coded_roots=hard_coded_roots)
        sierpinski_recursive(G, S+2*T+2, E, hard_coded_roots=hard_coded_roots)
    else:
        pass


def sierpinski_with_children(G, S, E, hard_coded_roots=[]):
    """
    args: 
        G (nx.DiGraph): The Graph object we'll be updating
        S (int):        The index of the starting node
        E (int):        The index of the ending node   

    NOTE:
        This adds the condition that if there exists an edge
        i -> j,
        there must also be edges from i to every child of j
    """
    if S != E:
        L = int( S + ((E-S-2)/6 ) ) # "left" point, midpoint of first third of [S,E]
        C = int( (S+E)/2 )          # "center" point, midpoint of second third
        R = int( E - (L-S) )        # "right" point, midpoint of last third
        G.add_edge(C, L)
        if R not in hard_coded_roots:
            G.add_edge(C, R)
            for j in G.predecessors(C):
                G.add_edge(j, R)
        if L not in hard_coded_roots:
            G.add_edge(C, L)
            for j in G.predecessors(C):
                G.add_edge(j, L)


        # Divide interval into thirds, apply function to each third
        T = 2*(L-S)
        sierpinski_with_children(G, S, S+T, hard_coded_roots=hard_coded_roots)
        sierpinski_with_children(G, S+T+1, S+2*T+1, hard_coded_roots=hard_coded_roots)
        sierpinski_with_children(G, S+2*T+2, E, hard_coded_roots=hard_coded_roots)
    else:
        pass


def sierpinski_positions(depth=None, N=None):
    if N is not None:
        depth = 0
        while pow(3, depth) < N:
            depth += 1

    t = 1
    pos = {0: (0,0), 1: (1,2), 2:(2, 0)}

    while t < depth:
        addition = pow(3, t)
        second_pos, third_pos = {}, {}
        for node, position in pos.items():
            second_pos[node + addition] = (position[0]+pow(2,t), position[1]+pow(2,t+1) - 1)
            third_pos[node+2*addition] = (position[0] + pow(2, t+1), position[1])
        
        pos.update(second_pos)
        pos.update(third_pos)
        t += 1

    return pos


def sierpinski_positions_2Np1(depth=None, N=None):
    if N is not None:
        depth = 0
        while pow(3, depth) < N:
            depth += 1

    t = 1
    pos = {0: (0,0), 1: (1,2), 2:(2, 0)}

    layers = [pow(3, depth-1) - 1]
    shifts = [0]
    for m in range(depth-2, -1, -1):
        layers.append(layers[-1] + pow(3, m))
        shifts.append(shifts[-1] + pow(2, m))
    layers.append(pow(3, 2*N+1))
    shifts.append(shifts[-1] + 1)

    s = 0
    while t < depth:
        addition = pow(3, t)
        second_pos, third_pos = {}, {}
        for node, position in pos.items():
            if node + addition > layers[s]:
                s += 1
            second_pos[node + addition] = (position[0]+pow(2,t) - shifts[s], position[1]+pow(2,t+1) )
        
        for node, position in pos.items():
            third_pos[node+2*addition] = (position[0] + pow(2, t+1) - shifts[s], position[1])
        
        pos.update(second_pos)
        pos.update(third_pos)
        t += 1

    return pos