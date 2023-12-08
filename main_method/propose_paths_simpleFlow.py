"""
Compress an input sequence an propose the n most frequent sequence paths within 
length k, using a Viterbi algorithm inspired by Maximum Likelihood Sequence 
Estimators. 
"""
from collections import defaultdict
import bisect
import networkx as nx
from matplotlib import pyplot as plt
from numpy import log

"""
Compress single-character runs. 

Return:
    shortSeq - the compressed sequence; "AAABBC" --> "ABC"
    runLengths - the runs indexed to shortSeq; "ABC" --> [3,2,1]
"""
def compress_runs(fileName):
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
    len_seq = j
    return shortSeq, runLengths, len_seq

"""
Propose the n most frequent sequence paths of <= k length. 
"""
class MLSE_Propose:
    def __init__(self, seq, num_proposals, max_path_len, sep='_', threshold=0.01, verbose=True):
        #self.filename = filename
        self.seq = seq
        self.len_seq = len(seq)
        self.num_proposals = num_proposals
        self.max_path_len = max_path_len
        self.sep = sep # char separator for node/edge labels
        self.verbose=verbose

        # flow graph
        self.len_graph = 0
        self.graph = list()
        self.graph_keys = None
        self.threshold = threshold # ignore a paths with less than a certain percentage of flow
        self.total_flow = 0

        self.mlse_viterbi()
    
    def get_proposals(self):
        return self.proposals

    def report(self):
        print("proposal", "exact_matches", sep='\t')
        for num_exact_matches, proposal in self.proposals:
            print(proposal, num_exact_matches, sep='\t')
    def print_graph(self):
        for i in range(0, self.len_graph):
            print(i)
            for node, data in self.graph[i].items():
                flow_in, dict_edges = data
                print(f'\t{node} {flow_in}')
                for out, flow_out in dict_edges.items():
                    print(f'\t\t{out}, {flow_out}')

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

    def update_dist(self, path, terminator=None):
        bottleneck = float('inf')
        for idx in range(0, self.len_graph-1):
            u = path[idx]
            v = path[idx+1]
            bottleneck = min(bottleneck, self.graph[idx][u][1][v])

        for idx in range(0, self.len_graph-1):
            u = path[idx]
            v = path[idx+1]
            self.graph[idx][u][1][v] = self.graph[idx][u][1][v] - bottleneck
        u = path[-1]
        #self.graph[-1][u][1][terminator] += new_flow
        return bottleneck

    """
    Use the distance graph to deduce the heaviest path (dynamic programming).
    """
    def backtrace_heavy_path(self, dist_from_source):
        if self.len_graph == 0:
            return 0, ''
        
        graph_end = self.len_graph - 1 # inclusive end
        path = ''
        
        # walk left and complete the heavy path
        dist_layer, u, v = self.layer_dist_max(dist_from_source[graph_end])
        if v is not None:
            path = str(v) + path
        for idx in range(graph_end, 0, -1):
            v = u
            dist_layer, u = dist_from_source[idx-1][v]
            if v is not None:
                path = str(v) + path

        if u is not None:
            path = str(u) + path

        # update the distances
        bottleneck = self.update_dist(path)
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
            bottleneck, proposed_path = self.backtrace_heavy_path(dist_from_source)
            if proposed_path not in seen:
                consec_repeats = 0
                seen.add(proposed_path)
                bisect.insort(proposals, (bottleneck, proposed_path), key=lambda x: -x[0])
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
    
    def update_graph(self, subseq):
        for path_offset in range(0, self.len_graph-1):
            u = subseq[path_offset]
            v = subseq[path_offset+1]
            self.add_edge(path_offset, u, v)
        self.add_edge(path_offset+1, v, None)
        self.total_flow += 1
    
    def seq_to_graph(self):
        # (re)-initialize graph
        self.len_graph = self.max_path_len
        self.graph_keys = range(self.len_graph)
        self.total_flow = 0
        for idx in self.graph_keys:
            self.graph.append(dict())
        # parse seq to graph
        seq_idx = 0
        len_window = 0
        window = ''
        while seq_idx < self.len_seq and len_window < self.len_graph:
            next = str(self.seq[seq_idx]).strip()
            seq_idx += 1
            if len_window == 0 or next != window[-1]: # compress duplicate runs
                window += next
                len_window += 1

        if len_window == self.len_graph:
            self.update_graph(window)
            while seq_idx < self.len_seq:
                next = str(self.seq[seq_idx]).strip()
                seq_idx += 1
                if next and len(next) > 0:
                    if next != window[-1]:
                        window = window[1:]
                        window += next
                        self.update_graph(window)
                else:
                    break

    def mlse_viterbi(self):
        self.seq_to_graph()
        
        """
        # parse seq directly from file
        with open(self.filename, 'r') as fh:
            len_window = 0
            window = ''
            while len_window < self.len_graph:
                next = str(fh.read(1)).strip()
                if len_window == 0 or next != window[-1]: # compress duplicate runs
                    window += next
                    len_window += 1

            if len_window == self.len_graph:
                self.update_graph(window)
                while True:
                    next = str(fh.read(1)).strip()
                    if next and len(next) > 0:
                        if next != window[-1]:
                            window = window[1:]
                            window += next
                            self.update_graph(window)
                    else:
                        break
        """
        
        # REMOVE LOW-CAPACITY EDGES (one signal processing paper just does beam search instead)
        required_flow = self.total_flow * self.threshold
        noise = set()
        for idx_layer in self.graph_keys:
            layer = self.graph[idx_layer]
            incomings = list(layer.keys())
            for u in incomings:
                outgoings = list(layer[u][1].keys())
                for v in outgoings:
                    if layer[u][1][v] < required_flow:
                        noise.add((idx_layer, u, v, layer[u][1].pop(v)))
        
        """
        for idx_layer in self.graph_keys:
            layer = self.graph[idx_layer]
            incomings = list(layer.keys())
            for u in incomings:
                outgoings = list(layer[u][1].keys())
                for v in outgoings:
                    layer[u][1][v] = log(layer[u][1][v])
        """

        # MAKE N PROPOSALS AND REMOVE PROBABLE DUPLICATES
        #proposals = propose_n(graph_sorted, graph, graph_incoming, n=num_proposals, min_passes=0.5*len_graph)
        self.proposals = self.propose_n(n=self.num_proposals, min_walks=0.5*self.len_graph)
        
        if self.verbose:
            self.report()

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

    def add_path_flow(self, path, new_flow, create:bool):
        len_path = len(path)
        for idx in range(len_path-1):
            u = path[idx]
            v = path[idx+1]

            if not create and not self.G.has_edge(u, v):
                raise ValueError("Cannot add path flow, as an edge does not exist and create=False is specified")
            print('prev capacity:', G[u][v]['capacity'], end=';')
            self.add_edge_flow(u, v, new_flow=new_flow)
            print('new capacity:', self.G[u][v]['capacity'])

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

    def get_heavy_path(self, source, target):
        flow_value, flow_dict = nx.maximum_flow(self.G, source, target)
        # Find the path with maximum flow
        path = nx.shortest_path(self.G, source=source, target=target, weight='capacity')
        return flow_value, path

    def visualize_viterbi(self, source_node='s', sink_node='t', visualize=True):
        # CREATE FLOW GRAPH
        self.G = nx.DiGraph()

        for idx in range(0, self.len_seq - self.max_path_len + 1): #max_path_len):
            # add the sliding window as a path
            self.add_edge_flow(source_node, self.seq[idx] + self.sep + str(0)) # source to first node
            for path_offset in range(0, self.max_path_len-1): # second to penultimate
                u = self.seq[idx + path_offset] + self.sep + str(path_offset)
                v = self.seq[idx + path_offset + 1] + self.sep + str(path_offset+1)
                self.add_edge_flow(u, v)
            
            self.add_edge_flow(self.seq[idx + self.max_path_len-1] + self.sep + str(self.max_path_len-1), sink_node) # source to first node
        
        # IF SPECIFIED, VISUALIZE GRAPH
        if visualize:
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


def query_seq_file(filepath_seq, num_proposals, max_path_len, sep='_', threshold=0.01, verbose=True):
    seq, runs, len_seq = compress_runs(filepath_seq)
    return MLSE_Propose(seq, len_seq, num_proposals, max_path_len, sep=sep, threshold=threshold, verbose=verbose)

def query_seq(seq, num_proposals, max_path_len, sep='_', threshold=0.01, verbose=True):
    return MLSE_Propose(seq, num_proposals, max_path_len, sep=sep, threshold=threshold, verbose=verbose)

def main():
    #seqfile = '../generate_sequence/simulation.txt'
    seqfile = '../Langevin_clustering/sequence.txt'
    import os, psutil
    process = psutil.Process()
    m0 = process.memory_info().rss
    print("Starting memory:", m0)
    seq, runs, len_seq = compress_runs(seqfile) 
    print("length of seq:", len_seq)
    #tests = ['17383463846233278957234957324878926784302958903246803246478', 'abcdefghijklmnop','abcdefghghghghghabcdefedcbaghghghghabcdefghghgh', 'abcdefgh', 'abcdefabcdefefefefefefefefefabcdefgzzzzz']
    #seq = tests[0]
    #len_seq = len(seq)
    m1 = process.memory_info().rss
    print("Added memory:", m1 - m0)
    #len_seq = len(seq)
    num_proposals = 50
    max_path_len = 10

    #plotter = MLSE_Plot(seq, num_proposals, max_path_len)

    #proposer = MLSE_Propose(seqfile, num_proposals, max_path_len, sep='_', threshold=0.01, verbose=True)
    #proposer = MLSE_Propose(seq, num_proposals, max_path_len, sep='_', threshold=0.01, verbose=True)
    proposer = query_seq(seq, num_proposals, max_path_len)
    #proposer.print_graph()
    m2 = process.memory_info().rss
    print("Added memory:", m2-m1)
    return

if __name__ == '__main__':
    main()