'''
Analysis tool kit for GenRev 1.0.
'''
import networkx as nx
import numpy as np
from gr_io import *

def geneRank(G,by=set(['degree','betweenness','score']),top=10):
    '''
    This function ranks genes in the network by measures specified. 
    All genes are ranked decreasingly. Top N genes are returned as 
    well as their measures. 
    '''
    import operator    

    index=by     #measurements used to rank genes
    N=0          #how many genes will be returned
    
    if top >= G.order():
        N=G.order()
    else:
        N=top
    
    if index not in set(['degree','score','betweenness']):
        raise Exception('Rank measurement is not valid.')

    reservior=[]    
    measure={}    #a dictionary, key is node, value is measurement.
    if index=='degree':   #for weighted network, weight is omitted
        measure=nx.degree(G)
    if index=='score':    #node score
        for item in G.node:
            measure[item]=G.node[item]['score']
    if index=='betweenness':    #normalized betweenness
        measure=nx.algorithms.centrality.betweenness_centrality(G)
    
    for gene in measure:
        measurement=measure[gene]
        reservior.append((gene,measurement))
    
    reservior=sorted(reservior, key=operator.itemgetter(1),reverse=True)  #the sorted gene list

    return (reservior[0:N])


def MCL(G,r=2,order_limit=1000):
    '''
    Implement the Markov Cluster Algorithm (MCL,http://www.micans.org/mcl/).
    For big networks, this process will be omitted in GenRev.
    '''
    partition=[]

    if G.order() > order_limit:    #for big networks, GenRev doesn't provide MCL. 
        return None
     
    allnodes=G.nodes()
    adjm=nx.linalg.spectrum.adj_matrix(G) #use edge weight is applicable.
    for i in range(adjm.shape[0]):        ##M+I
        adjm[i,i]=1
    dgr=np.apply_along_axis(sum,0,adjm)
    Pm=(adjm/dgr)

    M1=Pm.copy()
    M2=np.zeros(M1.shape)
    i=0
    while ((M2-M1).any()):
        M2=np.dot(M1,M1)                                 #expansion
        M1=np.apply_along_axis(lambda x: x**r, 0, M2)    #inflation
        M1=M1/np.apply_along_axis(sum,0,M1)              #
        #M2=M2.round(3)
        #M1=M1.round(3)

    #print M1
    g=nx.convert.from_numpy_matrix(M1)
    g_parts=nx.components.connected_components(g)
    
    for part in g_parts:
        collect=[G.nodes()[x] for x in part]
        partition.append(collect)

    return partition


def modularity_ref(G,partition):
    '''
    Straightforward implementation of  modularity. For ref only. 
    ''' 
    m=float(G.size())
    dg=nx.degree(G)
    
    Q=0
    for part in partition:
        for gA in part:
            for gB in part:
                if gA != gB:
                    if G.has_edge(gA,gB):
                        A=1
                    else:
                        A=0
                    ki,kj=dg[gA],dg[gB]
                    v=A-(ki*kj)/(2*m)
                    Q=Q+v
                    #print gA, gB, v
    Q=Q/(2*m)
    return Q
      

def modularity(G,partition):
    '''
    Compute modularity for MCL partitions. 
    The modularity, Q, lies in the range [-1,1].
    Q > 0 means the number of edges within groups exceeds the number 
    expected by chance.
    Partition is a list of list, with each element is a node name. 
    '''
    #def fill_diag(M,v):
    #    for i in range(M.shape[0]):
    #        M[i,i]=v
         
    m=float(G.size())
    dg=nx.degree(G)
    Q=0
   
    for part in partition:
        n=len(part)
        M1=np.zeros((n,n))
        M2=np.zeros((n,n))
        for i in range(n):
            M1[i,]=dg[part[i]]
            M1[i,i]=0            #this will make sure the diagonal is 0. 
            M2[:,i]=dg[part[i]]
       
        MM=np.multiply(M1,M2)
        MM=MM/(2*m)
        expE=np.sum(MM)
        realE=2*(nx.subgraph(G,part).size())
        v=realE - expE
        Q=Q+v
    Q=Q/(2*m)
    return Q

def listNeighbors(G,nodes,d=1):
    '''
    Find neighbors for nodes in the order of d.
    If d=0, return nodes themselves. 
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

def netShrink(G,nodes,bs=1):
    '''
    Return a shrinked network
    bs: block size
    '''
    ngbs=listNeighbors(G,nodes,bs)
    totalnodes=ngbs.union(nodes)
    subg=nx.subgraph(G,totalnodes)
    return subg

def simplify(G):
    '''
    Remove self loops from networks.
    Inplace change.
    '''
    for node in G.edge:
        cont=G.edge[node]
        if cont.has_key(node):
            cont.pop(node)
    


if __name__=='__main__':
    print 'Testing functions...'
    G=read_edgelist('testdata/ppi_id.txt')
    print G.degree('983')
    print len(G.edge['983'])
    


