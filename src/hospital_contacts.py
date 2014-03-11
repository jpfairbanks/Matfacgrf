import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import argparse
import sys


def load_row(G, row, graphframe):
    src = row[1]
    dst = row[2]
    if src not in G:
        G.add_node(src, role=graphframe.ix[src][3])
    if dst not in G:
        G.add_node(dst, role=graphframe.ix[dst][4])
    G.add_edge(src, dst)


def mapcolors(nodelist):
    """Map each role to a color pass the values to nx.draw_networkx"""
    colors = dict()

    def color(role):
        """return the color string"""
        if 'PAT' in role:
            col = 'y'
        elif 'NUR' in role:
            col = 'r'
        elif 'MED' in role:
            col = 'g'
        elif 'ADM' in role:
            col = 'b'
        else:
            print(role)
            col = 'k'
        return col

    for vertex, properties in nodelist:
        colors[vertex] = color(properties["role"])
    return colors


def draw_graph(G,
               color_func=mapcolors,
               layout_func=nx.drawing.layout.spring_layout):
    """Draw the graph using the colors and layout that we picked"""
    colors = color_func(G.nodes(data=True))
    layout = layout_func(G)
    nx.draw_networkx(G, pos=layout, node_color=colors.values())


def get_args():
    """Produce the args namespace based on command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('--draw-graph', action='store_true')
    parser.add_argument('--drawfile', type=str)
    parser.add_argument('--matfile', type=str,
                        help='Store the graph matrices into this file in matlab format')
    arg = parser.parse_args()
    return arg


def write_adj_mat(G, fileobj=sys.stdout):
    """Write G to a sparse matrix format that Julia and Matlab can read."""
    lapmatrix = nx.laplacian_matrix(G)
    norm_lapl = nx.normalized_laplacian_matrix(G)
    adjmatrix = nx.adjacency_matrix(G)
    mdict = {'laplacian': lapmatrix,
             'norm_lapl': norm_lapl,
             'adjacency': adjmatrix}
    sio.savemat(fileobj, mdict)
    return mdict


if __name__ == '__main__':
    arg = get_args()
    graphframe = pd.read_table(arg.filename, header=None)
    G = nx.Graph()
    insertion_count = 0
    batchsize = 100
    for row_tuple in graphframe.iterrows():
        load_row(G, row_tuple[1], graphframe)
        insertion_count += 1
        if insertion_count % batchsize == 0:
            #do something for each batch
            pass
    if arg.draw_graph:
        draw_graph(G)
        if arg.drawfile:
            plt.savefig(arg.drawfile)
        else:
            plt.show()
    if arg.matfile:
        import scipy.io as sio
        write_adj_mat(G, arg.matfile)
