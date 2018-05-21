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
			print count 
			count += 1

	result.close()
	return True

def main():
	mappers = load_id2symbol()
	net_converter("HumanNet.txt", mappers)

if __name__ == '__main__':
	main()