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
    def __init__(self, filename, num_proposals, max_path_len, sep='_', threshold=0.01, verbose=True):
        self.filename = filename
        #self.seq = seq
        #self.len_seq = len_seq
        self.num_proposals = num_proposals
        self.max_path_len = max_path_len
        self.transition_key = sep # char separator for node/edge labels
        self.verbose=verbose

        # flow graph
        self.len_graph = 0
        self.graph = dict()
        self.graph_range = range(0, 0)
        self.threshold = threshold # ignore a paths with less than a certain percentage of flow
        self.total_flow = 0

        # special short fluctuations between conformations that were previously observed in a Viterbi parse
        #   key: int: the number of states in the fluctuation
        #   values: Set[str]: k-length fluctuation as a single string, e.g.: 1<->2 fluctuation adds '->12->' and '->21->'
        self.fluctuations = defaultdict(set)
        self.transition_key = '->' # delimiter to encode multi-chararcter fluctuation as its own "code" (e.g.: '->12->' is a single code, just like '1' or '2')
        self.len_key = len(self.transition_key)

        self.mlse_viterbi()
        self.mlse_viterbi(parse=2)

    def read_next(self, subseq, offset):
        end = offset + 1
        if subseq[offset:].find(self.transition_key) == 0:
            end = offset + 2*self.len_key + subseq[offset+self.len_key:].find(self.transition_key)
        return subseq[offset:end], end
    
    def is_transition(self, candidate):
        len_candidate = len(candidate)
        return len_candidate > 2*self.len_key and candidate.find(self.transition_key) == 0 and \
            candidate[self.len_key:].find(self.transition_key) == len_candidate - 2*self.len_key
    
    def strip(self, seq):
        if self.is_transition(seq):
            return seq[self.len_key:-1*self.len_key]
        return seq

    def get_proposals(self):
        return self.proposals
    
    def compress_transitions(self, seq):
        seq_out = ''
        offset = 0
        len_path = len(seq)
        while offset < len_path:
            next, offset = self.read_next(seq, offset)
            t = self.is_transition(next)
            next = self.strip(next)
            if seq_out and t:
                seq_out += next[1:]
            else:
                seq_out += next
        return seq_out

    def report(self, compact=True):
        print("proposal", "exact_matches", sep='\t')
        for num_exact_matches, proposal in self.proposals:
            if compact:
                proposal = self.compress_transitions(proposal)
            print(proposal, num_exact_matches, sep='\t')
    
    def add_edge(self, seed, idx_layer, u, v, new_flow=1):
        layer = self.graph[seed][idx_layer]
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
    
    def update_dist(self, seed, path, terminator=None):
        bottleneck = float('inf')

        u, u_end = self.read_next(path, 0)
        idx = 0
        while idx < self.len_graph-1:
            v, v_end = self.read_next(path, u_end)
            bottleneck = min(bottleneck, self.graph[seed][idx][u][1][v])
            
            u, u_end = v, v_end
            idx += 1
        
        u, u_end = self.read_next(path, 0)
        idx = 0
        while idx < self.len_graph-1:
            v, v_end = self.read_next(path, u_end)
            self.graph[seed][idx][u][1][v] = self.graph[seed][idx][u][1][v] - bottleneck
            u, u_end = v, v_end
            idx += 1
        
        #self.graph[-1][u][1][terminator] += new_flow
        return bottleneck

    """
    Use the distance graph to deduce the heaviest path (dynamic programming).
    """
    def backtrace_heavy_path(self, seed, dist_from_source):
        if self.len_graph == 0:
            return 0, ''
        
        graph_end = self.len_graph - 1 # inclusive end
        path = str(seed)
        
        # walk left and complete the heavy path
        dist_layer, u, v = self.layer_dist_max(dist_from_source[seed][graph_end])
        if v is not None:
            path = str(v) + path
        for idx in range(graph_end, 0, -1):
            v = u
            dist_layer, u = dist_from_source[seed][idx-1][v]
            if v is not None:
                path = str(v) + path

        if u is not None:
            path = str(u) + path

        # update the distances
        bottleneck = self.update_dist(seed, path)
        return bottleneck, path
    
    def build_dist_graph(self):
        dist_from_source = dict()
        for seed in self.graph:
            dist_from_source[seed] = list()
            for idx in range(0, self.len_graph):
                dist_from_source[seed].append(dict())
                for u, data in self.graph[seed][idx].items():
                    _, edges_outgoing = data
                    for v, c in edges_outgoing.items():
                        if idx == 0:
                            d = c
                        else:
                            d = c + dist_from_source[seed][idx-1][u][0]
                        if v not in dist_from_source[seed][idx] or d > dist_from_source[seed][idx][v][0]:
                            dist_from_source[seed][idx][v] = (d, u)
        return dist_from_source

    def propose_n_candidates(self, seed, proposals, seen, n:int, min_walks:int):
        len_proposals = len(proposals)
        prev_len = None
        consec_repeats = 0
        while(len_proposals < n):
            # find longest path to each node in this topologically sorted DAG
            dist_from_source = self.build_dist_graph()
            #print(self.graph[seed][0:2])
            #print(dist_from_source[seed][0:2])
            #assert(False)
            bottleneck, proposed_path = self.backtrace_heavy_path(seed, dist_from_source)
            if proposed_path not in seen:
                consec_repeats = 0
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

    def propose_n(self, n:int, num_seeds:int, min_walks:int):
        # construct proposals and track seen ("visted"? oops) paths
        proposals = list()
        seen = set()
        prev_len = len(proposals)

        n_seed = max(1, int(n/num_seeds))
        for seed in self.graph:
            seed_proposals = list()
            self.propose_n_candidates(seed, seed_proposals, seen, n_seed, min_walks)
            new_len = len(proposals)
            while (new_len < n and new_len != prev_len):
                prev_len = new_len
                self.propose_n_candidates(seed, seed_proposals, seen, n_seed, min_walks)
                new_len = len(proposals)

            for proposal in seed_proposals:
                bisect.insort(proposals, proposal, key=lambda x: -x[0])
        return proposals[:n]
    
    def update_graph(self, subseq):
        seed, seed_end = self.read_next(subseq, 0)
        if seed not in self.graph:
            # make subgraph starting with the letter
            self.graph[seed] = list()
            for idx in self.graph_range:
                self.graph[seed].append(dict())

        # add flow to path
        idx = 0
        u, u_end = seed, seed_end
        while idx < self.len_graph-1:
            v, v_end = self.read_next(subseq, u_end)
            #print(path_offset)
            self.add_edge(seed, idx, u, v)
            u, u_end = v, v_end
            idx += 1
        #print(path_offset)
        self.add_edge(seed, idx, v, None)
        self.total_flow += 1
    
    def record_fluctuations(self, proposed_path, len_fluctuation):
        counts = defaultdict(int)
        total = 0
        len_path = len(proposed_path)
        u, u_end = self.read_next(proposed_path, 0)
        idx = 0
        while idx < self.len_graph-1:
            v, v_end = self.read_next(proposed_path, u_end)
            counts[(u, v)] += 1
            total += 1
            u, u_end = v, v_end
            idx += 1

        for transition_pair, count in counts.items():
            if count >= total/2 * 0.5:
                u, v = transition_pair
                if len(u) > 0 and len(v) > 0 and u != v:
                    self.fluctuations[len_fluctuation].add(
                        self.transition_key + u + v + self.transition_key
                    )
                    self.fluctuations[len_fluctuation].add(
                        self.transition_key + v + u + self.transition_key
                    )

    """
    parse: the iteration of mlse_viterbi, either: {1 = first iteration, 2 = second iteration}
        on the kth iteration, the kth-order "codebook of fluctuations" will be used
        to test for folding paths that allow/minimize some long run of "ABABABAB" to '[AB]'
    """
    def mlse_viterbi(self, parse=1):
        self.graph = dict()
        self.len_graph = self.max_path_len
        self.graph_range = range(self.len_graph)
        self.total_flow = 0

        with open(self.filename, 'r') as fh:
            len_window = 0
            window = ''
            next = str(fh.read(1)).strip()
            prev = None
            while len_window < self.len_graph:
                if len_window > 0 and next == window[-1]:
                    next = str(fh.read(1)).strip()
                else:
                    if len_window == 0:
                        window += next
                        prev = next
                        next = str(fh.read(1)).strip()
                    elif next != prev: # compress duplicate runs
                        transition_chars = [prev, next]
                        transition = self.transition_key + window[-1] + next + self.transition_key
                        if parse > 1 and transition in self.fluctuations[parse]:
                            # compress this possible fluctuation into one code
                            window += transition
                            while next and next in transition_chars:
                                # skip long fluctuations of these particular states
                                prev = next
                                next = str(fh.read(1)).strip()
                        else:
                            window += next
                            prev = next
                            next = str(fh.read(1)).strip()
                    len_window += 1

            if len_window == self.len_graph:
                self.update_graph(window)
                prev = next
                next = str(fh.read(1)).strip()
                while True:
                    if not (next and len(next) > 0):
                        break
                    if next == prev:
                        next = str(fh.read(1)).strip()
                        continue
                    else:
                        cut, len_cut = self.read_next(window, 0) # find oldest "letter" to cut from window (may be represented with >1 character)

                        transition_chars = [prev, next]
                        transition = self.transition_key + prev + next + self.transition_key
                        if parse > 1 and transition in self.fluctuations[parse]:
                            # compress this possible fluctuation into one code
                            if len_cut > 6 or len_cut not in [1, 5, 6]:
                                assert(False)
                            window = window[len_cut:]
                            window += transition
                            while next and next in transition_chars:
                                # skip long fluctuations of these particular states
                                prev = next
                                next = str(fh.read(1)).strip()
                        else:
                            window += next
                            window = window[len_cut:]
                            prev = next
                            next = str(fh.read(1)).strip()
                        self.update_graph(window)
        
        # MAKE N PROPOSALS AND REMOVE PROBABLE DUPLICATES
        #proposals = propose_n(graph_sorted, graph, graph_incoming, n=num_proposals, min_passes=0.5*len_graph)
        num_seeds = len(self.graph)
        self.proposals = self.propose_n(n=self.num_proposals, num_seeds=num_seeds, min_walks=0.5*self.len_graph)

        len_fluctuation = parse + 1
        for bottleneck, proposed_path in self.proposals:
            if bottleneck >= 0.01 * self.total_flow: # significant flow
                self.record_fluctuations(proposed_path, len_fluctuation)
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

def query_seq(seq, len_seq, num_proposals, max_path_len, sep='_', threshold=0.01, verbose=True):
    return MLSE_Propose(seq, len_seq, num_proposals, max_path_len, sep=sep, threshold=threshold, verbose=verbose)

def main():
    seqfile = '../Langevin_clustering/sequence.txt'
    #seqfile = '../generate_sequence/simulation.txt'
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
    num_proposals = 30
    max_path_len = 10

    #plotter = MLSE_Plot(seq, num_proposals, max_path_len)

    proposer = MLSE_Propose(seqfile, num_proposals, max_path_len, sep='_', threshold=0.01, verbose=True)
    print(proposer.fluctuations)

    m2 = process.memory_info().rss
    print("Added memory:", m2-m1)
    return

if __name__ == '__main__':
    main()