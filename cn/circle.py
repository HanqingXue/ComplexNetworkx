import matplotlib.pyplot as plt
import networkx as nx
from gr_io import *

try:
    import pygraphviz
    from networkx.drawing.nx_agraph import graphviz_layout
except ImportError:
    try:
        import pydot
        from networkx.drawing.nx_pydot import graphviz_layout
    except ImportError:
        raise ImportError("This example needs Graphviz and either "
                          "PyGraphviz or pydot")

#G = nx.balanced_tree(3, 5)
#G=read_edgelist('testdata/humanNet.net')

#G = nx.balanced_tree(3, 5)
G=read_edgelist('testdata/humanNet.net')
seed =read_terminals('pathway/{0}.txt'.format(1))
sub = G.subgraph(seed)

degree_sequence = sorted([d for n, d in sub.degree()], reverse=True)
# print "Degree sequence", degree_sequence
dmax = max(degree_sequence)

plt.loglog(degree_sequence, 'b-', marker='o')
plt.title("Degree rank plot")
plt.ylabel("degree")
plt.xlabel("rank")

plt.axes([0.45, 0.45, 0.45, 0.45])
Gcc = sorted(nx.connected_component_subgraphs(sub), key=len, reverse=True)[0]
pos = nx.spring_layout(Gcc)
plt.axis('off')
nx.draw_networkx_nodes(Gcc, pos, node_size=20)
nx.draw_networkx_edges(Gcc, pos, alpha=0.4)

plt.show()



'''
seed =read_terminals('pathway/{0}.txt'.format(1))
sub = G.subgraph(seed)
mapper = {}

node = list(sub.node())
print type(node)
for i in range(0, len(node)):
    mapper[node[i]] = i

g = nx.Graph()

for u, v in sub.edges():
    g.add_edge(mapper[u], mapper[v])

pos = graphviz_layout(g, prog='neato')
plt.figure(figsize=(8, 8))
nx.draw(g, pos, node_size=20, alpha=0.5, node_color="blue", with_labels=False)

nx.draw_circular(g)
plt.axis('equal')
plt.show()


'''
# Author: Aric Hagberg (hagberg@lanl.gov)

#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
'''

import random

try:
    import pygraphviz
    from networkx.drawing.nx_agraph import graphviz_layout
except ImportError:
    try:
        import pydot
        from networkx.drawing.nx_pydot import graphviz_layout
    except ImportError:
        raise ImportError("This example needs Graphviz and either "
                          "PyGraphviz or pydot.")

import matplotlib.pyplot as plt

import networkx as nx
from networkx.algorithms.isomorphism.isomorph import graph_could_be_isomorphic as isomorphic
from networkx.generators.atlas import graph_atlas_g


def atlas6():
    """ Return the atlas of all connected graphs of 6 nodes or less.
        Attempt to check for isomorphisms and remove.
    """

    Atlas = graph_atlas_g()[0:208]  # 208
    # remove isolated nodes, only connected graphs are left
    U = nx.Graph()  # graph for union of all graphs in atlas
    for G in Atlas:
        zerodegree = [n for n in G if G.degree(n) == 0]
        for n in zerodegree:
            G.remove_node(n)
        U = nx.disjoint_union(U, G)

    # list of graphs of all connected components
    C = nx.connected_component_subgraphs(U)

    UU = nx.Graph()
    # do quick isomorphic-like check, not a true isomorphism checker
    nlist = []  # list of nonisomorphic graphs
    for G in C:
        # check against all nonisomorphic graphs so far
        if not iso(G, nlist):
            nlist.append(G)
            UU = nx.disjoint_union(UU, G)  # union the nonisomorphic graphs
    return UU


def iso(G1, glist):
    """Quick and dirty nonisomorphism checker used to check isomorphisms."""
    for G2 in glist:
        if isomorphic(G1, G2):
            return True
    return False


if __name__ == '__main__':
    G = atlas6()

    print("graph has %d nodes with %d edges"
          % (nx.number_of_nodes(G), nx.number_of_edges(G)))
    print(nx.number_connected_components(G), "connected components")

    plt.figure(1, figsize=(8, 8))
    # layout graphs with positions using graphviz neato
    pos = graphviz_layout(G, prog="neato")
    # color nodes the same in each connected subgraph
    C = nx.connected_component_subgraphs(G)
    for g in C:
        c = [random.random()] * nx.number_of_nodes(g)  # random color...
        nx.draw(g,
                pos,
                node_size=40,
                node_color=c,
                vmin=0.0,
                vmax=1.0,
                with_labels=False
               )
    plt.show()
    '''