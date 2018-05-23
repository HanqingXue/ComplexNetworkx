import os

def load_id2symbol():
	gene_info = open('geneinfo_ncbihgnc.csv', 'r')
	mappers = {}

	for record in gene_info.readlines():
		record = record.split(',')
		keys = record[0]
		val = record[1][:-1]
		mappers[keys] = val

	return mappers

def net_converter(net, mappers={}):
	net = open(net, 'r')
	result = open('humanNet.net', 'w')
	count = 0 
	for line in net.readlines():
		line = line.split('\t')
		if line[0] in mappers.keys() and line[1] in mappers.keys():
			result.write('{0}\t{1}\n'.format( mappers[line[0]], mappers[line[1]]))
			count += 1

	result.close()
	return True

def main():
	
	net = open('humanNet.net', 'r')
	all_nodes = []
	for item in net:
		item = item[:-1].split('\t')
		all_nodes.append(item[0])
		all_nodes.append(item[-1])


	for i in range(1, 299):
		out = open('../pathway/{0}.txt'.format(i), 'w')
		f = open('../new_pathway/{0}.txt'.format(i), 'r')
	
		for line in f.readlines():
			if line[:-1] not in all_nodes:
				print line
				pass
			else:
				out.write(line)

		
		f.close()
		out.close()
	#net_converter("HumanNet.txt", mappers)

if __name__ == '__main__':
	main()