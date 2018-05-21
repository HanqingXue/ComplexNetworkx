'''
This program implements the limited k-walk algorithm proposed by Dupont P. in 2006. 
It also includes the k-walk algorithm, which is very fast for small networks. 
Author: Siyuan Zheng, siyuan.zheng@vanderbilt.edu
Module for GenRev 1.0.
'''
import warnings
warnings.filterwarnings("ignore")

import networkx as nx
import time
import numpy as np
import math
import os

#########################################################################

def matrix_power(M,n):
    m=M.copy()
    if n < 0:
        raise Exception('n should positive integer')
    if n==0:
        m=np.identity(m.shape[0])
        return m
    if n==1:
        return m
    while n > 1:
        m=np.dot(m,M)
        n=n-1
    return m

def match(seq,pool):
    index=[None]*len(seq)
    for i in xrange(len(pool)):
        if pool[i] in seq:
            for pos in range(len(seq)):
                if seq[pos]==pool[i]:
                    index[pos]=i
                    break
    return index
                

def kWalk(K,G,r=0.5):
    '''
    Straightfoward implementation of the k walks algorithm, proposed
    by P. Dupont et al in 2006.
    K, the set of seed nodes;
    G, the global network. Order should no bigger than 1000;
    r, a inclusion parameter, the proportion of the expected passage
    time of a node Vt to that of the node of interest x. The rationale is,
    if the expected passage number of x is much bigger than other candidate
    nodes, it is reasonable to keep the walk to x. Otherwise, if the possibility
    to other node is bigger than M*r, this node will be included. M is
    the expected passage time of x. 
    '''
    if r==0:
        r=0.001
    if not nx.components.connected.is_connected(G):
        raise Exception('G must be connected!')
    if (len(K)<2):
        raise Exception('There must be more than two seed nodes!')
    allnodes=G.nodes()
    adjm=nx.linalg.spectrum.adj_matrix(G)
    dgr=np.apply_along_axis(sum,0,adjm)
    Pm=(adjm/dgr).transpose()    #initial transition matrix
    collect=set()
    
    for NOI in K:
        select=set()
        Kprime=list(set(K)-set([NOI]))
        CandN=list(set(allnodes) - set(K))
        Qlabels=[NOI]+CandN    ##the first element in Qlabels would be NOI

        index=match(Qlabels, allnodes)
        Qx=Pm[index,:][:,index]    #the first row is NOI, then the candidate genes
        index2=match(Kprime,allnodes)
        Rx=Pm[index,:][:,index2]
        #print Qx
        #print Qlabels

        I=np.identity(Qx.shape[0])

        N=np.linalg.inv(I-Qx)
        #print NOI
        #print Qlabels
        #print N.round(2)

        tmpN=N.copy()[0,:]
        scoreVec=tmpN.tolist()[0]
        scores=list(set(scoreVec))
        scores.sort(reverse=True)
        w1=scores[0]
        w2=scores[1]

        flag=0
        if scoreVec[0]==scores[0]:
            flag=1
        for i in range(1,tmpN.shape[1]):
            node=Qlabels[i]
            if scoreVec[i]==w1:
                select.add(node)
            if scoreVec[i]==w2:
                if flag and scoreVec[i]*(1.0/r) >= w1:
                    select.add(node)

        collect=collect.union(select)
        #print I-Qx
        #print Qlabels
        #print N
        #print select
        #print '\n\n'
    finalset=collect.union(K)
    subG=nx.subgraph(G,finalset)
    return subG


#G=nx.components.strongly_connected_component_subgraphs(G)[0]
#K=['a','f']
#a=kWalk(K,G)

#################################################################################
#limited k walks
#G=nx.components.connected_component_subgraphs(G)[0]

def limkWalks(K,G,L=50,iteration=1):
    '''
    The limited k-walk algorithm, proposed
    by P. Dupont et al in 2006. This function is for connected network only.
    Otherwise, call limkSearch(). 
    '''
    if not nx.components.connected.is_connected(G):
        raise Exception ('G has to be connected!')
    #iteration=1
    allnodes=G.nodes()
    #L=10
    n=len(allnodes)
    adjm = nx.linalg.spectrum.adj_matrix(G)  


      #this will use 'edge weight' attribute.
    dgr=np.apply_along_axis(sum,0,adjm)
    Pm=(adjm/dgr).transpose()    #initial transition matrix

    collect=set(K)

    for iter_i in range(iteration):

        K=set(K)
        k=len(K)
        #print 'iter',iter_i,K
        for NOI in K:

            selected=NOI
            Kprime=list(K-set([NOI]))
            CandN=list(set(allnodes) - K)
            Qlabels=[NOI]+CandN    ##the first element in Qlabels would be NOI
                                   #Qlabels represents the transient states

            index=match(Qlabels, allnodes)
            Qx=Pm[index,:][:,index]    #the first row is NOI, then the candidate genes
            index2=match(Kprime,allnodes)
            Rx=Pm[index,:][:,index2]

            selected=NOI
            
            ##build two lattices, Lalpha and Lbeta

            ###############################################################################
            ##the row names of Lalpha part1 and part2 are Qlabels and Kprime respectively.
            #this matrix is the probability (i,l) of starting the walk in x and
            #reaching state i in l steps.
            #the basis probability
            Lalpha_part1=np.zeros((n-k+1,L))   #for the transient state
            Lalpha_part2=np.zeros((k-1,L))     #for the absorbing state
            Lalpha_part1[0,0]=1

            for l in range(1,L):
                vec_part1=np.dot(Lalpha_part1[:,l-1],Qx)
                Lalpha_part1[:,l]=vec_part1
                vec_part2=np.dot(Lalpha_part1[:,l-1],Rx)
                Lalpha_part2[:,l]=vec_part2

            Lalpha=np.vstack((Lalpha_part1,Lalpha_part2))

            #the following code is a straight implementation, just for confirmation purpose
            #for l in range(1,L):
            #    Lalpha_part1[:,l]=matrix_power(Qx,l)[0,:]
            #    Lalpha_part2[:,l]=np.dot(matrix_power(Qx,l-1),Rx)[0,:]

            ################################################################################
            ##the beta lattice. Row names are NOI + Candinate nodes.
            #I compared the results, they are identical with the straightforward
            #implementation.
            #
            Lbeta=np.zeros((n-k+1,L))
            Lbeta[:,0]=0    #when sl=L-1, the basis of Beta lattice

            for sl in range((L-1)-1,-1,-1):
                if sl==(L-1)-1:
                    #print np.apply_along_axis(sum, 1, np.nan_to_num(Rx))
                    cur_sum = np.nan_to_num(Rx).sum(axis=1)
                    cur_sum = np.array(cur_sum)
                    Lbeta[:,(L-1)-sl] = cur_sum.T[0]
                    continue
                else:
                    vec_b=np.dot(Qx,Lbeta[:,(L-1)-sl-1])
                    Lbeta[:,(L-1)-sl]=vec_b
            
            #the following code is a straight implementation, just for confirmation purpose
            #for sl in range((L-1)-1,-1,-1):
            #    print (L-1)-sl,sl
            #    tmp=np.dot(matrix_power(Qx,(L-1)-sl-1),Rx)
            #    vec=np.apply_along_axis(sum,1,tmp)
            #    Lbeta[:,(L-1)-sl]=vec
            
            ################################################################################
            ##calculate the expected passage times for each edge.
            #actually, we don't have to calculate all the edges

            edgeM=np.zeros((n-k+1,n-k+1))    ###expected passage times for edges between transient states
            for i in range(n-k+1):
                stack=[0]*(n-k+1)
                for j in range(n-k+1):
                    B=Lbeta[0,L-1]
                    if B==0:
                        stack[j]=0
                        continue
                    #if j < n-k+1:
                    A=sum(Lalpha[i,]*Qx[i,j]*Lbeta[j,])
                    stack[j]=A/B
                        #print 'transient',i,j
                    #else:
                    #    A=Lalpha[i,L-1]*Rx[i,j-(n-k+1)]
                    #    stack[j]=A/B
                        #print 'absorbing',i,j
                edgeM[i,]=stack

            ##automatic edge weight threshold selection
            #maximal theta to induce a connected subgraph from edges
            thset=np.apply_along_axis(max,1,edgeM)
            thset=list(set(thset))
            thset.sort(reverse=True)
            for theta in thset:
                if theta==0:
                    break
                else:
                    tmpM=edgeM.copy()
                    tmpM[tmpM < theta]=0
                    tmpM[tmpM >= theta]=1
                    tmpg=nx.convert.from_numpy_matrix(tmpM)
                    tmpg2=nx.convert.from_edgelist(tmpg.edges())
                    if nx.components.connected.is_connected(tmpg2):
                        break
                    else:
                        continue

            collect=collect.union(set([Qlabels[x] for x in tmpg2.nodes()]))
        K=collect    #for iteration purpose.
    subg=nx.subgraph(G,collect)        
    return subg


def limkSearch(K,G,L=50,iteration=1):
    '''
    Find subnetwork from a set of terminals using limited k-walk algorithm.
    K: terminals
    G: the edge weighted network
    '''
    subG=nx.Graph()
    for gG in nx.components.connected_component_subgraphs(G):
        local_terminals=set(K) & set(gG.nodes())
        if len(local_terminals) >= 2:
            subgraph=limkWalks(local_terminals,gG,L,iteration)
            subG=nx.compose(subG,subgraph)
    return subG


if __name__=='__main__':
    print 'Now testing limited k-walk algorithm now...'
    from gr_io import *
    G=read_edgelist('testdata/lesmis.net')

    term=read_terminals('testdata/lesmis.terminal')
    t1=time.time()
    g=limkSearch(term,G,L=10,iteration=1)
    #print g.nodes()
    #print g.edges()
    t2=time.time()
    print t2 - t1
