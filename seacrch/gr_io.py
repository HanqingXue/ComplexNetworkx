'''
IO module for GenRev v1.0.
'''

import networkx as nx
import numpy as np
import os,sys
import re
import time


def read_edgelist(path, comment='#',default_score=1):
    '''
    Return a graph object.
    Read the edgelist network file. 
    Input file is tab or space delimited. 
    First two columns are genes, and the third column is positive edge scores. 
    Larger score indicates more relevant relationship.
    If the third column is omitted, default value is 1. 
    '''
    infilename=path
    infile=open(infilename)
    basename=os.path.basename(infilename)
    pref=basename.split('.')[0]
    G=nx.Graph(name=pref)
    neg_warning=0
    for line in infile:
        line=line.strip()
        if line.startswith(comment):
            continue
        elif line:
            info=line.split()
            if len(info) >=3:
                try:
                    w=float(info[2])
                    if w <= 0:
                        neg_warning=1
                except ValueError:
                    raise Exception('The third columne of network file should be edge score!')
                    
            else:
                w=default_score
            G.add_edge(info[0],info[1],score=w)
    infile.close()
    if neg_warning:
        raise Exception('Edge scores should be positive!')
    else:
        return G


def read_sif(path, comment='#'):
    '''
    Return a graph object.
    For test only, since SIF format doesn't have edge information.
    '''
    infilename=path
    infile=open(infilename)
    basename=os.path.basename(infilename)
    pref=basename.split('.')[0]
    G=nx.Graph(name=pref)
    for line in infile:
        line=line.strip()
        if line.startswith(comment):
            continue
        elif line:
            info=line.split()
            gA, gB=info[0],info[2]
            G.add_edge(gA,gB)
    infile.close()
    return G
  

def read_nodes(path, comment='#',default_score=1):
    '''
    Return a dictionary, key is node, value is score.
    Read node information from space or tab delimited files.
    Second column is positive scores. If omitted, default value is 1. 
    '''
    infile=open(path)
    nodes={}
    neg_warning=0
    for line in infile:
        line=line.strip()
        if line.startswith(comment):
            continue
        else:
            info=line.split()
            if len(info)>=2:
                try:
                    w=float(info[1])
                    if w <= 0:
                        neg_warning=1
                except ValueError:
                    raise Exception('The second columne of node file shoud be node score!')
            else:
                w=default_score
            nodes[info[0]]=w
    infile.close()
    if neg_warning:
        raise Exception('Node scores should be positive!')
    else:
        return nodes
    

def read_terminals(path):
    '''
    Return a set of nodes.
    Read terminals from terminal file.
    Terminal file is single column, each line is a terminal. 
    '''
    terminals=set()
    infile=open(path)
    for line in infile:
        info=line.strip()
        if info:
            terminals.add(info)
    infile.close()
    return terminals


def layNode2Graph(G,node_info,attr='score'):
    '''
    Return the overlaid graph.
    Overlay node weights to a graph.
    If the node in G is not present in node_info, it will be omitted. 
    '''
    total=set(G.nodes())
    wN=set(node_info.keys())
    overlap=total & wN
    subg=nx.subgraph(G,overlap)
    for node in subg:
        subg.node[node][attr]=node_info[node]
    return subg


def mapAttr(G,entry='node',attr_from='score',attr_to='weight',f=lambda x:x):
    '''
    In place change the graph.
    Since User input edge and node score, I have to define a function to map 
    these scores to weights. 'weight' attribute is used for calculation, so it's better
    to seperate with 'score' attribute.
    '''
    if entry=='node':
        for node in G:
            av1=G.node[node][attr_from]
            av2=f(av1)
            G.node[node][attr_to]=av2
    elif entry=='edge':
        for gA in G:
            for gB in G.edge[gA]:
                av1=G.edge[gA][gB][attr_from]
                av2=f(av1)
                G.edge[gA][gB][attr_to]=av2   
    #return G    


def setNodeScore(G, default=1):
    '''
    In place change G.
    Set node score to default.
    '''
    for n in G:
        G.node[n]['score']=default
    #return G


def setEdgeScore(G, default=1):
    '''
    In place change G.
    Set edge weight to default.
    '''
    for n in G.edge:
        for neighbor in G.edge[n]:
            G.edge[n][neighbor]['score']=default
    #return G    
    

def autoDetectDir(path,prefix='GenRev_analysis'):
    '''
    Detect directories with name similar with GenRev_analysis12192010_1.
    The pattern is "prefix" + "date" + "_" + "num". For analyses in one day, 
    number will increase by 1 for each run. 
    This function will return the directory name for current analysis.
    '''
    T=time.localtime()
    mon=T.tm_mon
    mday=T.tm_mday
    year=T.tm_year
    if mon < 10:
        mon='0'+str(mon)
    else:
        mon=str(mon)
    if mday < 10:
        mday='0'+str(mday)
    else:
        mday=str(mday)
    year=str(year)
    presence=mon+mday+year
    existDirPat=re.compile('^'+prefix+presence+'_'+'(\d*)')
    runCol=[]
    for f in os.listdir(path):
        if existDirPat.search(f):
            d=int(existDirPat.search(f).group(1))
            runCol.append(d)
    else:
        runCol.append(0)
    runCol.sort()
    cur_run=str(runCol[-1]+1)
    directName=prefix+presence+'_'+cur_run
    return(directName)
    
    
def writeSif(G,path,edgeType='ee'):   ##
    '''Write graph to Cytoscape SIF format.'''
    outfile=open(path,'w')
    for edge in G.edges():
        gA,gB=edge
        outfile.write('%s\t%s\t%s\n'%(gA,edgeType,gB))
    outfile.close()


def writeEdgeAttr(G,path,edgeType='ee',attrName='EdgeCategory',attrValues=[],Class='String'):   ##eda format
    '''
    Write the edge attribute to output file. 
    Conventional file format is ".eda".
    attrsValues is a list of vlaues with equal length of network size.
    '''
    if len(attrValues)!=G.size():
        raise Exception('Number of attribute values is not equal to network size!')
    outfile=open(path,'w')
    outfile.write('%s class=%s\n'%(attrName,Class))
    es=G.edges()
    for i in xrange(G.size()):
        E=es[i]
        gA,gB=E
        outfile.write('%s (%s) %s = %s\n'%(gA,edgeType,gB,attrValues[i]))        
    outfile.close()

    
def writeNodeAttr(G,path,attrName='NodeCategory',attrValues=[], Class='String'):  ##noa
    '''
    Write the node attribute to output file. 
    Conventional file format is ".noa".
    attrValues is a list of values with equal length of network order.
    '''
    if len(attrValues)!=G.order():
        raise Exception('Number of attribute values is not equal to network order!')
    outfile=open(path,'w')
    outfile.write('%s class=%s\n'%(attrName,Class))
    ns=G.nodes()
    for i in xrange(G.order()):
        gen=ns[i]
        w=attrValues[i]
        outfile.write('%s = %s\n'%(gen,w))
    outfile.close()
   

def extNodeAttr(G,attr='weight'):
    '''
    Return the node attr list. The node attr value sequence is 
    corresponding to G.nodes()
    ''' 
    wvec=[]
    for N in G.nodes():
        w=G.node[N][attr]
        wvec.append(w)
    return wvec


def extEdgeAttr(G,attr='weight'):
    '''
    Return the edge attr list. The edge attr value sequence is 
    corresponding to G.edges()
    '''
    wvec=[]
    for E in G.edges():
        gA,gB=E
        w=G.edge[gA][gB][attr]
        wvec.append(w)
    return wvec
 

def determCat(G, subG=None, terminals=[]):
    '''
    Determine the categories of nodes and edges. 
    Node category includes "terminal","linker" and "other", and 
    accordingly edge categories include "terminal_terminal","terminal_linker", 
    "linker_linker" and "other". 
    G is global network, subG is extracted subnetwork,
    Terminals are input seeds.
    if subG is none, G is the subnetwork.
    Return a tuple of length 2, the first is node category list, the second is 
    edge category list. 
    '''
    nodeCategory=[]
    edgeCategory=[]
    #if subG is not None, then we will collect node and edge categorization info 
    #for the global network.
    if subG:  #
        for node in G:
            if node in terminals:
                cat='terminal'
            elif subG.has_node(node):
                cat='linker'
            else:
                cat='other'
            nodeCategory.append(cat)
        for edge in G.edges():
            gA,gB=edge
            if subG.has_edge(gA,gB):
                cat='subnetwork'
            else:
                cat='other'
            edgeCategory.append(cat)
    #if subG is None, then G is probably the extracted subnetwork itself. Then we 
    #will collect node and edge categorization info for it. 
    else:   #
        for node in G:
            if node in terminals:
                cat='terminal'
            else:
                cat='linker'
            nodeCategory.append(cat)
        for edge in G.edges():
            gA,gB=edge
            tmp=set([gA,gB]) & terminals
            if len(tmp)==2:
                cat='terminal_terminal'
            elif len(tmp)==1:
                cat='terminal_linker'
            else:
                cat='linker_linker'
            edgeCategory.append(cat)

    return (nodeCategory, edgeCategory)
            





