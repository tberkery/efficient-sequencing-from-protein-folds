from collections import defaultdict
import bisect
seq = 'abcdefabcdefefefefefefefefefabcdefgzzzzz'
len_seq = len(seq)
print(len_seq)
sep = '_'
num_proposals = 3
max_path_len = 10

def add_edge(graph, idx_layer, u, v, new_flow=1):
    layer = graph[idx_layer]
    if u not in layer:
        layer[u] = [0, defaultdict(int)]    
    layer[u][0] += new_flow
    layer[u][1][v] += new_flow

"""
Max edge from an adjacency dict with capacities
"""
def get_max_edge(dict_edges):
    c_max = float('-inf')
    e_max = None
    for e, c in dict_edges.items():
        if c > c_max:
            c_max = c
            e_max = e
    return e, c

"""
Use the "incoming" graph to backtrace a heavy path.
    graph_start is inclusive
    graph_end is exclusive
"""
def backtrace_heavy_path(seed, idx_curr, graph_start, graph_end, graph, graph_incoming, m=int(num_proposals*0.1)):
    path = str(seed)
    capacities = []
    if idx_curr > graph_start and idx_curr < graph_end-1:
        #print(idx_curr, seed)
        capacities.append(graph[idx_curr][seed][0])
    
    # walk left
    seed_left = seed
    for skip_left in range(0, idx_curr-graph_start):
        idx_left = idx_curr-skip_left
        capacity, data = graph_incoming[idx_left][seed_left]
        seed_next, c_edge = get_max_edge(data)
        #print(graph_incoming[idx_left][seed_left])
        bisect.insort(capacities, capacity)
        path = str(seed_next) + path
        #capacities.insert(0, graph_incoming[idx_left][seed_left][0])
        seed_left = seed_next
    
    # walk right
    seed_right = seed
    for idx_right in range(idx_curr, graph_end):
        capacity, data = graph[idx_right][seed_right]
        seed_next, c_edge = get_max_edge(data)
        bisect.insort(capacities, capacity)
        if seed_next is None:
            break
        path += str(seed_next)
        seed_right = seed_next

    if len(capacities) > m:
        bottleneck = capacities[m]
    elif len(capacities) > 0:
        bottleneck = capacities[-1]
    else:
        bottleneck = 0
    return bottleneck, path

def propose_n(proposals, graph_sorted, graph, graph_incoming, n:int, min_passes:int):
    next_best = 0
    len_graph = len(graph_sorted)
    while(len(proposals) < n and next_best < n):
        passes = 0
        for idx in range(0, len_graph):
            layer = graph_sorted[idx]
            if next_best >= len(layer):
                continue
            else:
                passes += 1
            bottleneck, heavy_path = backtrace_heavy_path(layer[next_best][1], idx, graph_start=0, graph_end=len_graph, graph=graph, graph_incoming=graph_incoming)
            proposals.add(heavy_path)
        if passes == 0 or passes < min_passes:
            return proposals
        next_best += 1
    return proposals

def mlse_viterbi(threshold=0.01):
    len_graph = max_path_len
    graph_keys = range(len_graph)
    # flow graph
    graph = list()
    # incoming edge graph (may act as residual graph under certain considerations)
    graph_incoming = list()
    # all edges per layer, sorted from highest to lowest capacity
    graph_sorted = list()
    graph_incoming_sorted = list()
    
    for idx in graph_keys:
        graph.append(dict())
        graph_incoming.append(dict())
        graph_sorted.append(list())
        graph_incoming_sorted.append(list())
    total_flow = 0
    for idx_seq in range(0, len_seq - max_path_len+1, max_path_len):
        add_edge(graph_incoming, 0, seq[idx_seq], None)
        for path_offset in range(0, len_graph-1):
            u = seq[idx_seq+path_offset]
            v = seq[idx_seq+path_offset+1]
            add_edge(graph, path_offset, u, v)
            add_edge(graph_incoming, path_offset+1, v, u)
        add_edge(graph, path_offset+1, v, None)
        total_flow += 1
    # some triming to make them line up
    graph_incoming = graph_incoming

    # REMOVE LOW-CAPACITY EDGES
    required_flow = total_flow * threshold
    noise = set()
    #noise_rev = set() # sanity check
    for idx_layer in graph_keys:
        layer = graph[idx_layer]
        layer_rev = graph_incoming[len_graph-idx_layer-1]
        incomings = list(layer.keys())
        for u in incomings:
            outgoings = list(layer[u][1].keys())
            for v in outgoings:
                if layer[u][1][v] < required_flow:
                    noise.add((idx_layer, u, v, layer[u][1].pop(v)))
                    layer_rev[v][1].pop(u)
                    #noise_rev.add((idx_layer, u, v, test))



    # SORT NODES IN EACH LAYER BY HIGHEST CAPACITY
    for idx_layer in graph_keys:
        layer = graph[idx_layer]
        for key, data in layer.items():
            capacity, edges = data
            bisect.insort(graph_sorted[idx_layer], [capacity, key, edges], key=lambda x: -x[0])
        layer_rev = graph_incoming[idx_layer]
        for key, data in layer_rev.items():
            capacity, edges = data
            bisect.insort(graph_incoming_sorted[idx_layer], [capacity, key, edges], key=lambda x: -x[0])

    # MAKE N PROPOSALS AND REMOVE PROBABLE DUPLICATES
    proposals = set()
    propose_n(proposals, graph_sorted, graph, graph_incoming, n=num_proposals, min_passes=0.5*len_graph)
    print(proposals)
    
    '''
    for idx_layer in graph_keys:
        print(idx_layer)
        layer = graph_sorted[idx_layer]
        for data in layer:
            print(' ', data)
    print('__')
    for idx_layer in graph_keys:
        print(idx_layer)
        layer = graph_incoming_sorted[idx_layer]
        for data in layer:
            print(' ', data)
    layer = graph_sorted[0]
    for data in layer:
        print(' ', data)
    print('__')
    layer_rev = graph_incoming_sorted[-1]
    for data in layer_rev:
        print(' ', data)
    '''


mlse_viterbi()