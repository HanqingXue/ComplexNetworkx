import warnings
warnings.filterwarnings("ignore")

import networkx as nx
import time
import numpy as np
import math
import os
from gr_io import *
import matplotlib.pyplot as plt

def get_cn_matrix():
	f = open('adj.txt', 'r')
	matrix = f.readlines()
	print len(matrix)

	out = open('cn.txt', 'w')
	for i in range(0, len(matrix)):
		for j in range(0, len(matrix)):
			cn = (matrix[i] or matrix[j]).count('1')
			out.write(str(cn)+',')
		print 'i:{}'.format(i)
		out.write('\n')

	out.close()
			

def get_string_matrix(matrix):
	row_num ,col_num = matrix.shape
	slist = []
	count = 0
	out = open('adj.txt', 'w')
	for i in range(0, row_num):
		s=''
		print count
		for j in range(0, col_num):
			if matrix[i, j] == 1:
				s += '1'
			else:
				s += '0'
		out.write(s+'\n')
		count += 1

	out.close()

def get_k_neighbor(G, L, node):
    seed = []
    pre_node =  nx.predecessor(G, node, None, L)
    for index in pre_node:
        seed += pre_node[index]

    return seed

def get_k_neighbor_term(G, term, L):
    seed = []
    for item in term:
        pre_node = get_k_neighbor(G, L, item)
        seed += pre_node

    return set(seed)


def main():
	
	G=read_edgelist('testdata/demo.net')
	M =  nx.to_numpy_matrix(G)
	print M.dot(M.T)
	print nx.jaccard_coefficient(G, ['a'])
	for u, v, p in nx.jaccard_coefficient(G, ['a']):
		print u, v, p


def link_pred(G, seeds, fname):
	print str(fname) + '.txt'
	expansion = []
	for node in G.node:
		G.node[node]['community'] = 1

	for seed in seeds:
		if seed in G.node.keys(): 
			G.node[seed]['community'] = 0

	count = 0
	cutoff = 10
	subgraph = G.subgraph(seeds)

	#seeds = list(set(subgraph.node) - set(seeds))


	preds = nx.cn_soundarajan_hopcroft(G, subgraph.edges())

	for u, v, p in preds:
		if p > cutoff:
			cutoff = p 

	print 'cutoff:{0}'.format(cutoff)
	nodes = []

	for seed in seeds:
		nodes += G.neighbors(seed)

	node = set()
	candiate_nodes = list(set(G.subgraph(nodes).edges()) - set(subgraph.edges()))
	preds = nx.cn_soundarajan_hopcroft(G, candiate_nodes)
	for u, v, p in preds:
		if p <  cutoff*8:
			continue
		expansion.append((u, v))
		node.add(u)
		node.add(v)

	#nx.write_edgelist(subgraph, "result/origingraph/{0}.edgelist".format(fname))
	origing_graph = open("result/origingraph/{0}.edgelist".format(fname), 'w')
	for u, v in subgraph.edges():
		origing_graph.write("{0}\t{1}\n".format(u, v))
	origing_graph.close()

	node_file = open('result/node/{0}.txt'.format(fname), 'w')
	for item in node:
		node_file.write('{0}\n'.format(item))
	node_file.close()

	edge_file = open('result/expanedge/{0}.txt'.format(fname), 'w')
	for u, v in expansion:
		edge_file.write('{0}\t{1}\n'.format(u, v))
	edge_file.close()
	
	print len(node)

if __name__ == '__main__':
	G=read_edgelist('../testdata/humanNet.net')

	for i in range(1, 299):
		seed =read_terminals('../pathway/{0}.txt'.format(i))
		link_pred(G, seed, str(i))