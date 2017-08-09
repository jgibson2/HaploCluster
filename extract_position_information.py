import argparse

def __main__():
	
	parser = argparse.ArgumentParser('Extract position information from clustered haplotypes')
	parser.add_argument('-p', '--positions', nargs='+',  help='Positions to extract', required=True)
	parser.add_argument('-i', '--input', help='Input filename', required=True)
	args = parser.parse_args()

	with open(args.input) as f:
		header = None
		indexes = list()
		for line in f:
			line = line.strip()
			if line:
				if line[0] == '#':
					header = line[1:].split('\t')
					for position in args.positions:
						try:
							indexes.append(header.index(position))
						except ValueError:
							pass
					print('#' + '\t'.join([pos + '::{}'.format(indexes[i]+1) for i, pos in enumerate(args.positions)]))
				elif header is None:
					raise ValueError('Did not detect header with positions!')
				else:
					data = line.split('\t')
					hap = data[0].split(' ')
					pos_haps = list()
					#print(indexes, hap)
					for loc in indexes:
						pos_haps.append(hap[loc])
					print(' '.join(pos_haps), data[1], sep='\t')
			else:
				print('\n')
	
if __name__ == '__main__':
	__main__()

