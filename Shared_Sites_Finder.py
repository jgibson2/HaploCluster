"""
This script detects all of the shared sites given a file listing all of the unique genotypes.
Coded By: Caleb Palagyi (2017)
A.M.D.G.
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", help="Declares the file from which the genotype data will be pulled", type=str, required=True)
parser.add_argument("-p", "--positions_file", help="Declares the file from which the positions will be pulled", type=str, required=True)
parser.add_argument("-o", "--output_file", help="Declares the file to which the data will be added", type=str, required=True)
args = parser.parse_args()
input_file_name = args.input_file
positions_file_name = args.positions_file
output_file_name = args.output_file
i = open(input_file_name, 'r')
p = open(positions_file_name, 'r')
o = open(output_file_name, 'w')
genotype_list = []
line_counter = 0
for line in i:
    genotype_as_list = line.strip().split(' ')
    genotype = str(genotype_as_list)[2:-2].replace('\', \'', '')
    genotype_list.append(genotype)
    line_counter = line_counter + 1
positions_list = []
for line in p:
    positions_list.append(line)
character_dict = {}
for position in positions_list:
    site_list = []
    for x in genotype_list:
        char = x[positions_list.index(position)]
        site_list.append(int(char))
    character_dict[position] = site_list
results_dict = {}
for key in character_dict.keys():
    check = sum(character_dict[key])
    if check/line_counter == 0:
        results_dict[key] = 0
    if check / line_counter == 1:
        results_dict[key] = 1
for key in results_dict.keys():
    o.write(str("chrX"))
    o.write('\t')
    o.write(str(key[5:]).strip())
    o.write('\t')
    int_key = int(key[5:])
    int_key = int_key+1
    o.write(str(int_key))
    o.write('\n')
