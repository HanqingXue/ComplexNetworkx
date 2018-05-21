'''
This program implement the approximate algorithm for node weighted
Steiner tree problem proposed by Klein P. and Ravi R. in 1995.
Since it aims to find a minimum-weight subnetwork, gene scores probably
need transformations to find a reasonable network.
Author: Siyuan Zheng, siyuan.zheng@vanderbilt.edu
Module for GenRev v1.0.
'''

import networkx as nx
import time
import operator
import numpy

def node2tree_link(G, node, T, asp):
    '''
    Calculate the cost of connecting node to a tree.
    Note that the cost of end points are not included.
    Return cost and merged tree
    G: global network
    node: the query node
    T: tree, graph object
    asp: all-pairs shortest path length of G
    '''
    tree=T.copy()
    vertices=T.nodes()
    linker=vertices[0]    #initial linker to query node.
    cost=None
    for x in vertices:
        w = 0
        pathVs = asp[node][x]     #'node is query'

        for vv in pathVs[1:-1]:   #exclude the end points
            w = w+G.node[vv]['weight']
        if cost != None:
            if w < cost:
                cost=w
                linker=x
        else:
            cost=w
            linker=x

    path=asp[node][linker]
    tree.add_path(path)
    return (cost,tree)

    
#############################################################

def NWConnSteiner(G,terminals,shortestPath=None):
    '''
    Find node weighted steiner tree solution from a
    connected graph. Vertex in the graph must have "weight" attribute.
    G: a node weighted graph
    terminals: the input seeds
    shortestPath: all pairs shortest paths
    '''
    gG=G.copy()
    terminals=set(terminals) & set(gG.nodes())  ##internal terminals
    net_vertices=set(gG.nodes())
    
    if not nx.is_connected(gG):
        raise Exception('Disjoint graph not allowed')
    if len(terminals) == 0:
        raise Exception('Input terminals not present in the graph')
    
    if shortestPath:
        asp=shortestPath
    else:
        asp=nx.shortest_path(gG)

    ##initialize the tree
    trees=nx.components.connected_component_subgraphs(nx.subgraph(gG,terminals))
    collector=set(terminals)    #collect terminals and all steiner nodes, for final subgraph extraction

    while len(trees) > 1:
        tree_nodes=[]
        [tree_nodes.extend(x.nodes()) for x in trees]    #tree nodes in current cycle
        non_tree_nodes = net_vertices - set(tree_nodes)  #non-tree nodes in current cycle
        ntvs = list(non_tree_nodes)
        steiner_node = False
        quotient_cost = False
        select_tree_index = False

        for ntv in ntvs:    #iterate the non-tree nodes, and find the candidate
            #print ntv
            cost_pool=[]
            #get the cost from this node to each tree
            for i in range(len(trees)):
                cst=node2tree_link(gG,ntv,trees[i],asp)[0]
                cost_pool.append((i,cst))

            #sort the node-tree cost
            cost_pool.sort(key=operator.itemgetter(1))
            cst_seq = [x[1] for x in cost_pool]
            #get the index of minimum cost from ntv to the trees
            #count the top n trees with minimal cost
            index=sum([1 for x in cst_seq if x==min(cst_seq)])
            #there must be at least two trees to calculate quotient cost
            if index==1: index=2   
            #calculate the quotient cost
            #total cost is the path cost plus node cost
            cost_sum = sum(cst_seq[:index]) + gG.node[ntv]['weight']
            ntv_quotient_cost=float(cost_sum)/index
            tree_index = [cost_pool[x][0] for x in range(index)]
            if not steiner_node:
                #initialize the metrics 
                steiner_node=ntv
                quotient_cost=ntv_quotient_cost
                select_tree_index=tree_index
            else:
                #update the metrics
                if ntv_quotient_cost < quotient_cost:
                    quotient_cost = ntv_quotient_cost
                    steiner_node=ntv
                    select_tree_index=tree_index
                if ntv_quotient_cost == quotient_cost:
                    #if two nodes have same quotient cost, prefer the light weight one.
                    if gG.node[ntv]['weight'] < gG.node[steiner_node]['weight']:
                        quotient_cost = ntv_quotient_cost
                        steiner_node=ntv
                        select_tree_index=tree_index
            
            #print ntv,ntv_quotient_cost,cost_pool, steiner_node
        #print steiner_node,cost_pool
        #merge the trees, given steiner node and tree index
        pretrees=[]    #used to store the updated treelist
        newtree=nx.Graph()  #the new formed tree
        for i in range(len(trees)):
            if i in select_tree_index:
                ct=trees[i]
                ct_update=node2tree_link(gG,steiner_node,ct,asp)[1]
                newtree=nx.compose(newtree,ct_update)
            else:
                pretrees.append(trees[i])
        pretrees.append(newtree)
        trees=pretrees    ##update the tree list
        collector=collector.union(newtree.nodes())
       
    subG=nx.subgraph(gG,collector)    #the returned subnetwork retains all attrs of global network. 
    return subG

def NWSteiner(G,terminals,shortestPath=None):
    '''
    Find a subnetwork from a set of terminals.
    This function implement the approximate algorithm for node weighted
    Steiner tree problem proposed by Klein P and Ravi R in 1995.
    G: the weighted network
    terminals: a list of input seed nodes
    shortestPath: all pairs shortest path
    '''
    if shortestPath:
        asp=shortestPath
    else:
        asp=nx.shortest_path(G)
    subG=nx.Graph()
    for gG in nx.components.connected_component_subgraphs(G):
        local_terminals=set(terminals) & set(gG.nodes())
        if len(local_terminals) > 0:
            subgraph=NWConnSteiner(gG,local_terminals,asp)
            subG=nx.compose(subG,subgraph)
    return subG
            

if __name__=='__main__':
    print 'Node weighted Steiner algorithm testing...'
    from gr_io import *
    import tools
    G=read_edgelist('testdata/ppi_id.txt')
    nm=read_nodes('testdata/case/gene.score.txt')
    term=read_terminals('testdata/case/top50.terminal')
    G=layNode2Graph(G,nm)
    term=set(term) & set(G.nodes())
    mapAttr(G,f=lambda x: 1.0/numpy.sqrt(x))
    t1=time.time()
    gG=tools.netShrink(G,term,bs=1)
    g=NWSteiner(gG,term)
    t2=time.time()
    print t2-t1
 

    

