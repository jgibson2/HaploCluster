import vcf
import argparse

def __main__():
	parser = argparse.ArgumentParser()
	parser.add_argument('files', nargs='+', help='Files to combine')
	parser.add_argument('--classes', help='File containing positions and classes in tab delineated format (chr, pos, class)')
	args = parser.parse_args()

	position_sample_dict = dict()
	samples_list = list()
	class_dict = dict()

	if args.classes:
		with open(args.classes, 'r') as class_file:
			for line in class_file:
				line = line.strip()
				if line[0] == '#':
					continue
				data = line.split('\t')
				class_dict[(data[0], int(data[1]))] = data[2]

	for i, f in enumerate(args.files):
		vcf_reader = vcf.Reader(open(f, 'r'))
		samples_list.append(vcf_reader.samples)
		for record in vcf_reader:
			if (record.CHROM, record.POS) not in position_sample_dict:	
				position_sample_dict[(record.CHROM, record.POS)] = list()
			position_sample_dict[(record.CHROM, record.POS)].append((f,[str(sample['GT']) for sample in record.samples]))

	if args.classes:
		print('#CHROM\tPOS\tCLASS\t' + '\t'.join(['\t'.join(samples) for samples in samples_list]))
	else:
		print('#CHROM\tPOS\t' + '\t'.join(['\t'.join(samples) for samples in samples_list]))
	for chrom_pos_tuple, file_genotype_list in sorted(position_sample_dict.items()):
		print('\t'.join(map(lambda x: str(x), chrom_pos_tuple)), end='\t')
		if args.classes:
			print(class_dict.setdefault(chrom_pos_tuple, '.'), end='\t')
			for index, file in enumerate(args.files):
				index = min(index, len(file_genotype_list)-1)
				if not file == file_genotype_list[index][0]:
					print('\t'.join(['.' for i in samples_list[index]]), end='\t')
				else:
					print('\t'.join(file_genotype_list[index][1]), end='\t')
			print()
		else:	
			for index, file in enumerate(args.files):
				index = min(index, len(file_genotype_list)-1)
				if not file == file_genotype_list[index][0]:
					print('\t'.join(['.' for i in samples_list[index]]), end='\t')
				else:
					print('\t'.join(file_genotype_list[index][1]), end='\t')
			print()

__main__()			
		
			
