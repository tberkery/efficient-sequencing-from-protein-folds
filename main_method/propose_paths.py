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

"""
Use networkx functionality to make a simple visualization of what our graphs look like
"""
class MLSE_plot:
    def __init__(self, seq, num_proposals, max_path_len):
        self.seq = seq
        self.len_seq = len(seq)
        self.num_proposals = num_proposals
        self.max_path_len = max_path_len
        self.sep = '_' # char separator for node/edge labels
    
    def add_edge_flow(G, u, v, new_flow=1):
        if G.has_edge(u, v):
            G[u][v]['capacity'] += new_flow
        else:
            G.add_edge(u, v, capacity=new_flow)
        return

    def decode_path(path, sep='_', keep_ends=True):
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
        G = nx.DiGraph()

        for idx in range(0, self.len_seq - self.max_path_len + 1, self.max_path_len):
            # add the sliding window as a path
            add_edge_flow(G, source_node, self.seq[idx] + self.sep + str(0)) # source to first node
        for path_offset in range(0, self.max_path_len-1): # second to penultimate
            u = self.seq[idx + path_offset] + self.sep + str(path_offset)
            v = self.seq[idx + path_offset + 1] + self.sep + str(path_offset+1)
            add_edge_flow(G, u, v)
        
        add_edge_flow(G, self.seq[idx + self.max_path_len-1] + sep + str(self.max_path_len-1), sink_node) # source to first node
    
        for layer, nodes in enumerate(nx.topological_generations(G)):
            # `multipartite_layout` expects the layer as a node attribute, so add the
            # numeric layer value as a node attribute
            for node in nodes:
                G.nodes[node]["layer"] = layer
        pos = nx.multipartite_layout(G, subset_key="layer")
        fig, ax = plt.subplots()
        nx.draw_networkx(G, pos=pos, ax=ax)
        ax.set_title("DAG layout in topological order")
        fig.tight_layout()
        plt.show()
        return
