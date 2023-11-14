import networkx as nx
from matplotlib import pyplot as plt

def add_edge_flow(G, u, v, new_flow=1):
    if G.has_edge(u, v):
        G[u][v]['capacity'] += new_flow
    else:
        G.add_edge(u, v, capacity=new_flow)
    return

def add_path_flow(G, path, new_flow, create:bool):
    len_path = len(path)
    for idx in range(len_path-1):
        u = path[idx]
        v = path[idx+1]

        if not create and not G.has_edge(u, v):
            raise ValueError("Cannot add path flow, as an edge does not exist and create=False is specified")
        print('prev capacity:', G[u][v]['capacity'], end=';')
        add_edge_flow(G, u, v, new_flow=new_flow)
        print('new capacity:', G[u][v]['capacity'])

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

def get_heavy_path(G, source, target):
    flow_value, flow_dict = nx.maximum_flow(G, source, target)
    # Find the path with maximum flow
    path = nx.shortest_path(G, source=source, target=target, weight='capacity')
    return flow_value, path

def max_flow(source_node='s', sink_node='t', visualize=False):
    # CREATE FLOW GRAPH
    G = nx.DiGraph()
    seq = 'abcdefabcdefefefefefefefefefabcdefgzzzzzzzz'
    len_seq = len(seq)
    sep = '_'
    num_proposals = 3
    max_path_len = 8

    for idx in range(0, len_seq - max_path_len + 1, max_path_len):
        # add the sliding window as a path
        add_edge_flow(G, source_node, seq[idx] + sep + str(0)) # source to first node
        for path_offset in range(0, max_path_len-1): # second to penultimate
            u = seq[idx + path_offset] + sep + str(path_offset)
            v = seq[idx + path_offset + 1] + sep + str(path_offset+1)
            add_edge_flow(G, u, v)
        
        add_edge_flow(G, seq[idx + max_path_len-1] + sep + str(max_path_len-1), sink_node) # source to first node
    
    # IF SPECIFIED, VISUALIZE GRAPH
    if visualize:
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
    
    # MAX FLOW, REMOVING FLOW FROM THE HEAVIEST PATH
    proposals = dict()
    for idx_proposal in range(num_proposals):
        max_flow, heaviest_path = get_heavy_path(G, source_node, sink_node)
        proposal = decode_path(heaviest_path, keep_ends=False)
        print('proposed path:', proposal, '; flow:', max_flow)
        # deduct flow from heaviest path
        add_path_flow(G, heaviest_path, new_flow=-1*max_flow, create=False)

def main():
    max_flow(visualize=True)

if __name__ == '__main__':
    main()