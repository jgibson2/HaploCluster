import sys

with open(sys.argv[1], 'r') as f:
	for line in f:
		line = line.strip()
		if line[0] == '#':
			print(line)
		else:
			data = line.split('\t')
			data[1] = str(int(sys.argv[2]) + int(data[1]))
			print('\t'.join(data))
