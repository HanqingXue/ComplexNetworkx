from gr_io import *
import networkx as nx  
from collections import deque  
import json
import matplotlib.pyplot as plt
import time

class SearchBaseExpansion(object):
	"""docstring for ClassName"""
	def __init__(self, cut_off, fname):
		file = open('{0}.json'.format(fname), 'r')
		self.fname = fname
		self.clustering = json.load(file)
		self.G = read_edgelist('../testdata/{0}'.format(fname))
		self.term = []
		self.result = []
		self.paths = []
		self.cut_off= cut_off
		file.close()

	def update_terminal(self, fname):
		return read_terminals('../pathway2/{0}.txt'.format(fname))

	def mul(self, x, y):
		return x*y 
	
	def get_clustering(self, key):
		return self.clustering[key]
	
	def bfs(self, source):		
		gen =  nx.dfs_edges(self.G, source, 2)
		expansion = []
		res = []
		paths = []

		while True:
			try:
				u, v = gen.next()
				path_from_source = nx.shortest_path(self.G, source, u)
				path_from_target = nx.shortest_path(self.G, source, v)
				likehood_source =  reduce(self.mul, map(self.get_clustering, path_from_source))
				likehood_target =  reduce(self.mul, map(self.get_clustering, path_from_target))
				if likehood_source > self.cut_off:
					if u not in self.term:
						res.append(u)
						paths.append(path_from_source)

				if likehood_source > self.cut_off:
					if v not in self.term:
						res.append(v)
						paths.append(path_from_target)

			except StopIteration as e:
				break
		return res, paths
		

	def run(self, fname):
		print fname
		start = time.time()
		self.term = self.update_terminal(fname)
		for gene in self.term:
			self.result.extend(self.bfs(gene)[0])
			self.paths.extend(self.bfs(gene)[1])

		node = open('./result/node/{0}.txt'.format(fname), 'w')
		for item in set(self.result):
			node.write('{0}\n'.format(item))
		node.close()

		edge = open('./result/expanedge/{0}.txt'.format(fname), 'w')
		for item in self.paths:
			edge.write(str(item)+'\n')
		edge.close()
		end = time.time()
		print 'time:{0}'.format(end - start)

	def plot_distribute(self):
		weigh = []
		count = 11 * [0]
		for u, v in self.clustering.items():
			index = int(v * 100) / 10
			count[index] += 1

		plt.bar(range(len(count)), count, ec='k', lw=1, hatch='o')
		plt.show()


	def calcuate(self):
		self.G.clustering
		pass

if __name__ == '__main__':
	algorithms = SearchBaseExpansion(0.7, 'Tumor_Net_Basic_p.net')
	'''
	f = open('{0}.json'.format(algorithms.fname), 'w')
	f.write(json.dumps(nx.clustering(algorithms.G)))
	f.close()
	'''
	
	
	for i in range(1, 299):
		algorithms.run(str(i))

