"""
Compress an input sequence an propose the n most frequent sequence paths within 
length k, using a Viterbi algorithm inspired by Maximum Likelihood Sequence 
Estimators. 
"""
from collections import defaultdict
import bisect
import networkx as nx
from matplotlib import pyplot as plt

"""
Compress single-character runs. 

Return:
    shortSeq - the compressed sequence; "AAABBC" --> "ABC"
    runLengths - the runs indexed to shortSeq; "ABC" --> [3,2,1]
"""
def compress_sequence(fileName):
    seq = ""
    with open(fileName, 'r') as ifile:
        for line in ifile:
            seq += line.strip()
    
    shortSeq = seq[0]
    runLengths = [1]
    j = 0 # index in shortSeq
    
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            runLengths[j] += 1
        else:
            shortSeq = shortSeq + seq[i]
            runLengths.append(1)
            j += 1
    
    return shortSeq, runLengths

"""
Propose the n most frequent sequence paths of <= k length. 
"""
class MLSE_Propose:
    def __init__(self, seq, num_proposals, max_path_len, sep='_', threshold=0.01, verbose=True):
        self.seq = seq
        self.len_seq = len(seq)
        self.num_proposals = num_proposals
        self.max_path_len = max_path_len
        self.sep = sep # char separator for node/edge labels
        self.verbose=verbose

        # flow graph
        self.len_graph = 0
        self.graph = list()
        self.threshold = threshold # ignore a paths with less than a certain percentage of flow

        self.mlse_viterbi()
    
    def get_proposals(self):
        return self.proposals

    def report(self):
        print("proposal", "exact_matches", sep='\t')
        for num_exact_matches, proposal in self.proposals:
            print(proposal, num_exact_matches, sep='\t')

    def add_edge(self, idx_layer, u, v, new_flow=1):
        layer = self.graph[idx_layer]
        if u not in layer:
            layer[u] = [0, defaultdict(int)]    
        layer[u][0] += new_flow
        layer[u][1][v] += new_flow

    '''
    layer is a soft copy of a graph layer (equal topology)
    '''
    def layer_dist_max(self, layer):
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

    def update_dist(self, path, new_flow, terminator=None):
        for idx in range(0, self.len_graph-1):
            u = path[idx]
            v = path[idx+1]
            self.graph[idx][u][1][v] += new_flow
        u = path[-1]
        self.graph[-1][u][1][terminator] += new_flow
        return

    """
    Use the distance graph to deduce the heaviest path (dynamic programming).
    """
    def backtrace_heavy_path(self, dist_from_source, m=0):
        if self.len_graph == 0:
            return 0, ''
        
        graph_end = self.len_graph - 1 # inclusive end
        path = ''
        flows = []
        
        # walk left and complete the heavy path
        flow, u, v = self.layer_dist_max(dist_from_source[graph_end])
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
        self.update_dist(path, -1*bottleneck)

        return bottleneck, path
    
    def build_dist_graph(self):
        dist_from_source = []
        for idx in range(0, self.len_graph):
            dist_from_source.append(dict())
            for u, data in self.graph[idx].items():
                _, edges_outgoing = data
                for v, c in edges_outgoing.items():
                    if idx == 0:
                        d = c
                    else:
                        d = c + dist_from_source[idx-1][u][0]
                    if v not in dist_from_source[idx] or d > dist_from_source[idx][v][0]:
                        dist_from_source[idx][v] = (d, u)
        return dist_from_source

    def propose_n_candidates(self, proposals, seen, n:int, min_walks:int):
        len_proposals = len(proposals)
        prev_len = None
        consec_repeats = 0
        while(len_proposals < n):
            # find longest path to each node in this topologically sorted DAG
            dist_from_source = self.build_dist_graph()
            bottleneck, proposal = self.backtrace_heavy_path(dist_from_source)
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

    def propose_n(self, n:int, min_walks:int):
        # construct proposals and track seen ("visted"? oops) paths
        proposals = list()
        seen = set()
        prev_len = len(proposals)
        self.propose_n_candidates(proposals, seen, n, min_walks)
        new_len = len(proposals)
        while (new_len < n and new_len != prev_len):
            prev_len = new_len
            self.propose_n_candidates(proposals, seen, n, min_walks)
            new_len = len(proposals)
        return proposals

    def mlse_viterbi(self):
        self.len_graph = self.max_path_len
        graph_keys = range(self.len_graph)
        
        for idx in graph_keys:
            self.graph.append(dict())
        total_flow = 0
        for idx_seq in range(0, self.len_seq - self.len_graph + 1, 1): # max_path_len):
            for path_offset in range(0, self.len_graph-1):
                u = self.seq[idx_seq+path_offset]
                v = self.seq[idx_seq+path_offset+1]
                self.add_edge(path_offset, u, v)
            self.add_edge(path_offset+1, v, None)
            total_flow += 1

        # REMOVE LOW-CAPACITY EDGES
        required_flow = total_flow * self.threshold
        noise = set()
        for idx_layer in graph_keys:
            layer = self.graph[idx_layer]
            incomings = list(layer.keys())
            for u in incomings:
                outgoings = list(layer[u][1].keys())
                for v in outgoings:
                    if layer[u][1][v] < required_flow:
                        noise.add((idx_layer, u, v, layer[u][1].pop(v)))
        
        # MAKE N PROPOSALS AND REMOVE PROBABLE DUPLICATES
        #proposals = propose_n(graph_sorted, graph, graph_incoming, n=num_proposals, min_passes=0.5*len_graph)
        self.proposals = self.propose_n(n=self.num_proposals, min_walks=0.5*self.len_graph)
        
        if self.verbose:
            report()

"""
Use networkx functionality to make a simple visualization of what our graphs look like
"""
class MLSE_Plot:
    def __init__(self, seq, num_proposals, max_path_len):
        self.seq = seq
        self.len_seq = len(seq)
        self.num_proposals = num_proposals
        self.max_path_len = max_path_len
        self.sep = '_' # char separator for node/edge labels
        
        self.G = None
        self.visualize_viterbi()
    
    def add_edge_flow(self, u, v, new_flow=1):
        if self.G.has_edge(u, v):
            self.G[u][v]['capacity'] += new_flow
        else:
            self.G.add_edge(u, v, capacity=new_flow)
        return

    def decode_path(self, path, sep='_', keep_ends=True):
        len_path = len(path)
        if len_path == 0:
            return ''
        seq = ''
        for idx in range(len_path):
            end = path[idx].find(sep)
            if end > 0:
                seq += path[idx][ 0:end]
            else:
                seq += path[idx]
        if keep_ends:
            return seq
        return seq[1:len_path-1-1]

    """
    Simple visualization of our Viterbi MLSE graph (sanity check and 
    for seeing the size of the graph, number of edges, etc.)
    """
    def visualize_viterbi(self, source_node='s', sink_node='t'):
        self.G = nx.DiGraph()

        for idx in range(0, self.len_seq - self.max_path_len + 1): #self.max_path_len):
            # add the sliding window as a path
            self.add_edge_flow(source_node, self.seq[idx] + self.sep + str(0)) # source to first node
        for path_offset in range(0, self.max_path_len-1): # second to penultimate
            u = self.seq[idx + path_offset] + self.sep + str(path_offset)
            v = self.seq[idx + path_offset + 1] + self.sep + str(path_offset+1)
            self.add_edge_flow(u, v)
        
        self.add_edge_flow(self.seq[idx + self.max_path_len-1] + self.sep + str(self.max_path_len-1), sink_node) # source to first node
    
        for layer, nodes in enumerate(nx.topological_generations(self.G)):
            # `multipartite_layout` expects the layer as a node attribute, so add the
            # numeric layer value as a node attribute
            for node in nodes:
                self.G.nodes[node]["layer"] = layer
        pos = nx.multipartite_layout(self.G, subset_key="layer")
        fig, ax = plt.subplots()
        nx.draw_networkx(self.G, pos=pos, ax=ax)
        ax.set_title("DAG layout in topological order")
        fig.tight_layout()
        plt.show()
        return

def main():
    seq = 'abcdefabcdefefefefefefefefefabcdefgzzzzz'
    len_seq = len(seq)
    sep = '_'
    num_proposals = 3
    max_path_len = 8

    plotter = MLSE_Plot(seq, num_proposals, max_path_len)
    
    return

if __name__ == '__main__':
    main()