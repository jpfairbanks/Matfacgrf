import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import argparse


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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('--draw-graph', action='store_true')
    arg = parser.parse_args()
    graphframe = pd.read_table(arg.filename, header=None)
    G = nx.Graph()
    insertion_count = 0
    batchsize = 100
    for row_tuple in graphframe.iterrows():
        load_row(G, row_tuple[1], graphframe)
        insertion_count += 1
        if insertion_count % batchsize == 0:
            pass
    if arg.draw_graph:
        draw_graph(G)
        plt.show()
