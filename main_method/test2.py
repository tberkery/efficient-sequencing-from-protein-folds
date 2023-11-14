from collections import defaultdict
import bisect
seq = 'abcdefabcdefefefefefefefefefabcdefgzzzzz'
len_seq = len(seq)
sep = '_'
num_proposals = 100
max_path_len = 8

def add_edge(graph, idx_layer, u, v, new_flow=1):
    layer = graph[idx_layer]
    if u not in layer:
        layer[u] = [0, defaultdict(int)]    
    layer[u][0] += new_flow
    layer[u][1][v] += new_flow

def layer_dist_max(layer):
    u_max = ''
    v_max = ''
    flow_max = float('-inf')
    for v in layer.keys():
        flow, u = layer[v]
        if flow > flow_max:
            u_max = u
            v_max = v
            flow_max = flow
    return flow_max, u_max, v_max

def update_dist(path, new_flow, graph, terminator=None):
    for idx in range(0, len(graph)-1):
        u = path[idx]
        v = path[idx+1]
        graph[idx][u][1][v] += new_flow
    u = path[-1]
    graph[-1][u][1][terminator] += new_flow
    return

"""
Use the distance graph to deduce the heaviest path (dynamic programming).
    graph_start, graph_end are INCLUSIVE
"""
def backtrace_heavy_path(len_graph, graph, dist_from_source, m=0):
    if len_graph == 0:
        return 0, ''
    
    graph_start = 0
    graph_end = len_graph - 1
    path = ''
    flows = []
    
    # walk left and complete the heavy path
    flow, u, v = layer_dist_max(dist_from_source[graph_end])
    bisect.insort(flows, flow)
    if v is not None:
        path = str(v) + path
    for idx in range(graph_end, 0, -1):
        v = u
        flow, u = dist_from_source[idx-1][v]
        bisect.insort(flows, flow)
        if v is not None:
            path = str(v) + path
    if u is not None:
        path = str(u) + path
    if len(flows) > m:
        bottleneck = flows[m] # a less guaranteed version where you use the mth lowest flow as the bottleneck (but you still use the lowest flow for path deduction)
    elif len(flows) > 0:
        bottleneck = flows[0]
    else:
        bottleneck = 0

    # update the distances
    update_dist(path, -1*bottleneck, graph)

    return bottleneck, path


def propose_n_candidates(proposals, seen, graph, n:int, min_walks:int):

    len_graph = len(graph)
    len_proposals = len(proposals)
    prev_len = None
    consec_repeats = 0
    while(len_proposals < n):
        # find longest path to each node in this topologically sorted DAG
        dist_from_source = build_dist_graph(graph)
        bottleneck, proposal = backtrace_heavy_path(len_graph, graph, dist_from_source)
        if proposal not in seen:
            consec_repeats = 0
            seen.add(proposal)
            bisect.insort(proposals, (bottleneck, proposal), key=lambda x: -x[0])
            #idx_add += 1
            len_proposals += 1
        else:
            # too many re-walks of the same path
            #   (this case is usually IMPOSSIBLE, unless the algorithm to deduct heaviest path changes)
            consec_repeats += 1
            if consec_repeats > min_walks:
                return
            
        if (prev_len and prev_len == len_proposals):
            return
        prev_len = len_proposals
    return

def build_dist_graph(graph):
    len_graph = len(graph)
    dist_from_source = []
    for idx in range(0, len_graph):
        dist_from_source.append(dict())
        for u, data in graph[idx].items():
            _, edges_outgoing = data
            for v, c in edges_outgoing.items():
                if idx == 0:
                    d = c
                else:
                    d = c + dist_from_source[idx-1][u][0]
                if v not in dist_from_source[idx] or d > dist_from_source[idx][v][0]:
                    dist_from_source[idx][v] = (d, u)
    return dist_from_source

def propose_n(graph, n:int, min_walks:int):
    # construct proposals and track seen ("visted"? oops) paths
    proposals = list()
    seen = set()
    prev_len = len(proposals)
    propose_n_candidates(proposals, seen, graph, n, min_walks)
    new_len = len(proposals)
    while (new_len < n and new_len != prev_len):
        prev_len = new_len
        propose_n_candidates(proposals, seen, graph, n, min_walks)
        new_len = len(proposals)
    return proposals

def mlse_viterbi(threshold=0.01):
    len_graph = max_path_len
    graph_keys = range(len_graph)
    # flow graph
    graph = list()
    
    for idx in graph_keys:
        graph.append(dict())
    total_flow = 0
    for idx_seq in range(0, len_seq - max_path_len + 1, 1): # max_path_len):
        for path_offset in range(0, len_graph-1):
            u = seq[idx_seq+path_offset]
            v = seq[idx_seq+path_offset+1]
            add_edge(graph, path_offset, u, v)
        add_edge(graph, path_offset+1, v, None)
        total_flow += 1

    # REMOVE LOW-CAPACITY EDGES
    required_flow = total_flow * threshold
    noise = set()
    for idx_layer in graph_keys:
        layer = graph[idx_layer]
        incomings = list(layer.keys())
        for u in incomings:
            outgoings = list(layer[u][1].keys())
            for v in outgoings:
                if layer[u][1][v] < required_flow:
                    noise.add((idx_layer, u, v, layer[u][1].pop(v)))
    
    # MAKE N PROPOSALS AND REMOVE PROBABLE DUPLICATES
    #proposals = propose_n(graph_sorted, graph, graph_incoming, n=num_proposals, min_passes=0.5*len_graph)
    proposals = propose_n(graph, n=num_proposals, min_walks=0.5*len_graph)
    
    print("proposal", "exact_matches", sep='\t')
    for num_exact_matches, proposal in proposals:
        print(proposal, num_exact_matches, sep='\t')


mlse_viterbi()