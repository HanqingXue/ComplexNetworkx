'''
This program implements the heuristic searching algorithm proposed 
by Trey Ideker et al. published in Molecular Systems Biology in 2009.
Author: Siyuan Zheng, syzhenger@gmail.com
For GenRev v1.0.
'''

import networkx as nx
import time
import numpy

#use graph node 'weight' attr to calculate overall score.
def avg_score(G):
    s=0
    N=G.order()
    for item in G:
        w=G.node[item]['weight']
        s=s+w
    score=float(s)/N
    return score

def sum_score(G):
    s=0
    for item in G:
        w=G.node[item]['weight']
        s=s+w
    return s

def neighbor(G,nodes,d=1):
    '''
    Find neighbors for nodes in the order of d.
    '''
    nodeset=set()
    init_nodes=set(nodes)
    while d:
        for node in init_nodes:
            nn=G.neighbors(node)
            nodeset=nodeset.union(set(nn))
        init_nodes=nodeset
        d=d-1
    ngb = nodeset - set(nodes)
    return ngb

def node_tree_merge(G,node,T):
    T=T.copy()
    ntl=None
    path=None
    for n in T.nodes():
        P=nx.shortest_path(G,node,n)
        L=len(P)
        if ntl:
            if L < ntl:
                ntl = L
                path = P
        else:
            ntl=L
            path=P
    T=nx.subgraph(G,(T.nodes()+path))
    return T

################################################################

def seedQuery(G,seedG,scoreFun=sum_score,d=2,r=0.1):
    '''
    Return a graph object.
    This function is to exapnd the seed graph seedG.
    Seed graph must be connected.
    scoreFun: function to score the graph
    d: search radius
    r: the expansion factor
    '''
    
    if not nx.components.is_connected(seedG):
        raise Exception('Disjoint terminals not allowed')
    if r==1.0:
        r=0.99    #to prevent ZeroDivisionError
        
    while True:
        ##initial state of each iteration
        N=seedG.order()
        seed=seedG.nodes()
        #iterate for each radius value
        for rad in range(1,d+1):
            subsum=scoreFun(seedG)    #The seed network score 
            pot_nodes=neighbor(G,seed,rad)     #The nodes potentially added to seed 
            if len(pot_nodes)==0:    #No neighbors available
                break
            add_node=[]    #nodes potentially can be added to the network
            w=max([G.node[x]['weight'] for x in pot_nodes])
            for item in pot_nodes:    #find nodes with local maximum weight
                if G.node[item]['weight']==w:
                    add_node.append(item)
            
            #newset=seed[:]     #newset is a collection of nodes after expansion
            #newset.append(add_node[0])
            #tmp_subg=nx.subgraph(G,newset)    #candidate expanded network
            #subsum_u=scoreFun(tmp_subg)

            if w >= subsum*(r/(1-r)):    ##fulfill the increment rate
                for node in add_node:
                    seedG=node_tree_merge(G,node,seedG)
                break
            
            #if subsum_u > subsum*(1+r):    #addition of the max node is valid under policy
            #    for node in add_node:      #add all valid nodes to the network
            #        seedG=node_tree_merg(G,node,seedG)
            #    break
            
        if seedG.order()==N:    #after iterations, no valid nodes found. 
            break
    return seedG

def listQuery(G,terminals,scoreFun=sum_score,d=1,r=0.2):
    '''
    Find subnetworks from a set of terminals.
    Implement the algorithm proposed by Chuang, et al, Trey Ideker. 2007.
    G: global network
    terminals: input seed nodes
    scoreFun: function to score the graph
    d: search radius
    r: the score increment factor
    '''
    subG=nx.subgraph(G,terminals)
    resultG=nx.Graph()
    for sg in nx.components.connected_component_subgraphs(subG):
        g=seedQuery(G,sg,scoreFun,d,r)
        resultG=nx.compose(resultG,g)
    return resultG
    

if __name__=='__main__':
    print 'Now testing the heuristic search algorithm...'    
       
            
    
    
        
        
        
