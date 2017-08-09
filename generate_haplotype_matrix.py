import vcf
import re
import argparse
from collections import Counter

def __main__():
	parser = argparse.ArgumentParser('Generate Haplotype Matrix from VCF File')
	parser.add_argument('vcf', help='VCF file to generate matrix from')
	parser.add_argument('--nucleotides', action='store_true', help='Use nucleotides instead of numbers')
	parser.add_argument('--haplotypes', action='store', help='Filename for storing  haplotypes and samples that contain them')
	parser.add_argument('-c', '--compressed', action='store_true', help='VCF file is compressed')
	parser.add_argument('-o', '--output', help='Output file (do not specify for stdout)')
	parser.add_argument('--assume-wild-type', action='store_true', help='Assume all unknown positions are wild type in final haplotype file', dest='assume_wt')
	
	args = parser.parse_args()
	
	#split only on tabs, specify filename instead of socket
	reader = vcf.Reader(fsock=None, filename=args.vcf, compressed=args.compressed, prepend_chr=False, strict_whitespace=True, encoding='ascii')

	header = []
	samples = []
	sample_dict = {}

	for vcf_rec in reader:
		#if vcf_rec.is_snp:
		if not samples:
			for call in vcf_rec.samples:
				samples.append(call.sample)
		header.append(vcf_rec.CHROM + ':' + str(vcf_rec.POS)) # 1-based position
		for call in vcf_rec.samples:
			if not call.sample in sample_dict:
				sample_dict[call.sample] = [] #if we haven't gotten to this sample yet, add it to the dict
			if args.nucleotides:
				s = call.gt_bases
				if s is None:
					s = '.'
				sample_dict[call.sample].append(tuple(re.split('[/\|]', s)))
			else:
				sample_dict[call.sample].append(tuple(re.split('[/\|]', call.data.GT)))


	if args.haplotypes:

		hap_dict = __generate_hap_dict__(sample_dict, assume_wt=args.assume_wt)
		
		if args.haplotypes:
			with open(args.haplotypes, 'w') as f:
				f.write('#' + '\t'.join(header) + '\n')
				for haplotype, sample_list in hap_dict.items():
					f.write(' '.join(['/'.join([str(x) for x in t]) for t in haplotype]) + '\t' + ' '.join(sample_list) + '\n')
		else:
			print('#' + '\t'.join(header))
			for haplotype, sample_list in hap_dict.items():
				print(' '.join(['/'.join([str(x) for x in t]) for t in haplotype]) + '\t' + ' '.join(sample_list))

	if args.output:
		with open(args.output, 'w') as f:
			f.write('#' + 'Sample\t' + '\t'.join(header) + '\n')
			for sample in samples:
				f.write(sample + '\t' + '\t'.join(['/'.join([str(x) for x in t]) for t in sample_dict[sample]]) + '\n')
	else:
		print('#' + 'Sample\t' + '\t'.join(header))
		for sample in samples:
			print(sample + '\t' + '\t'.join(['/'.join([str(x) for x in t]) for t in sample_dict[sample]]))

"""
	returns dict of haplotypes mapped to the samples that have them
"""
def __generate_hap_dict__(sample_dict, assume_wt=False):
	haps = set()
	haplotype_sample_map = dict()
	for sample, haplotype in sample_dict.items():
		if assume_wt:
			haplotype = [__replace_unknowns__(t) for t in haplotype]
		haplotype = tuple(haplotype)
		haps.add(haplotype)
		if not haplotype in haplotype_sample_map:
			haplotype_sample_map[haplotype] = []
		haplotype_sample_map[haplotype].append(sample)

	iter_haps = haps.copy()
	
	for h in iter_haps:
		to_remove = False
		new_hap = list(h)
		unknown_locs = [i for i,x in enumerate(h) if '.' in x]
		new_hap = [x for i,x in enumerate(new_hap) if i not in unknown_locs]
		
		#print(new_hap)

		for hap in haps:
			temp_hap = list(hap)
			temp_hap = [x for i,x in enumerate(temp_hap) if i not in unknown_locs]
			if temp_hap == new_hap and not hap == h:
				to_remove = True
				#print(hap, h)
				tmp = haplotype_sample_map[h]
				haplotype_sample_map[hap].extend(tmp)
		if to_remove:
			haps.remove(h)
			del haplotype_sample_map[h]
	
	return haplotype_sample_map


def __replace_unknowns__(allelic_calls : tuple):
	result = []
	for call in allelic_calls:
		if '.' in call:
			if any(filter(lambda x : '.' not in x, allelic_calls)):
				c = Counter(allelic_calls)
				c.subtract(Counter({'.':10000000}))
				c += Counter()
				#print(c, allelic_calls)
				call = c.most_common()[0][0]
			else:
				call = '0'
		result.append(call)
	#print(result)
	return tuple(result)

	
if  __name__ == '__main__':
	__main__()
