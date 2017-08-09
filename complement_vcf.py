import sys
import vcf
"""
replacement_dict = {'A':'T', 'C':'G', 'T':'A', 'G':'C', '.':'.'}

with open(sys.argv[1]) as f1:
		for line in f1:
			if line[0] == '#':
				print(line.strip())
				continue
			data = line.strip().split('\t')
			data[3] = replacement_dict[data[3]]
			data[4] = replacement_dict[data[4]]
			print('\t'.join(data))

"""

"""
f1_dict = dict()

with open(sys.argv[1]) as f1:
	with open(sys.argv[2]) as f2:
		for line in f1:
			if line[0] == '#':
				continue
			data = line.strip().split('\t')
			f1_dict[(data[0],data[1])] = [data[3], data[4]]
		for line in f2:
			if line[0] == '#':
				print(line.strip())
				continue
			data = line.strip().split('\t')
			if (data[0],data[1]) in f1_dict and f1_dict[(data[0],data[1])][0] != data[3]:
				print('\t'.join(data[0:3] + f1_dict[(data[0],data[1])][0:1] + data[4:]))
			else:
				print(line.strip())

"""
reader = vcf.Reader(fsock=None, filename=sys.argv[1], prepend_chr=False, strict_whitespace=True, encoding='ascii')
writer = vcf.Writer(sys.stdout, reader)
for vcf_rec in reader:
	#if vcf_rec.is_snp:
	for call in vcf_rec.samples:
		if call.data.GT == '0/0':
			call.data= call.data._replace(GT='0')
		elif call.data.GT == '1/1':
			call.data= call.data._replace(GT='1')
		else:
			call.data= call.data._replace(GT='.')
	writer.write_record(vcf_rec)
writer.close()

