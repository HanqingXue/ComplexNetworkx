import networkx as nx  
from collections import deque  
from collections import defaultdict

def bfs_edges(G,source):
    """Produce edges in a breadth-first-search starting at source."""
    # Based on http://www.ics.uci.edu/~eppstein/PADS/BFS.py
    # by D. Eppstein, July 2004.
    print visited
    visited=set([source])
    stack = [(source,iter(G[source]))]
    while stack:
        parent,children = stack[0]
        try:
            child = next(children)
            if child not in visited:
                print '{0}'.format(parent)
                print '{0}'.format(child)
                yield parent,child
                visited.add(child)
                stack.append((child,iter(G[child])))
        except StopIteration:
            stack.pop(0)


def bfs_tree(G, source):
    """Return directed tree of breadth-first-search from source."""
    return nx.DiGraph(bfs_edges(G,source))


def bfs_predecessors(G, source):
    """Return dictionary of predecessors in breadth-first-search from source."""
    return dict((t,s) for s,t in bfs_edges(G,source))


def bfs_successors(G, source):
    """Return dictionary of successors in breadth-first-search from source."""
    d=defaultdict(list)
    for s,t in bfs_edges(G,source):
        d[s].append(t)
    return dict(d)

def  main():
    time=0  
    G = nx.Graph()  
    G.add_edges_from([(1,2),(1,3),(2,4),(2,5),(3,6),(4,8),(5,8),(3,7)])  
    bfs_edges(G, 1)
    bfs_tree(G, 1)
    pass


if __name__ == '__main__':
    main()

