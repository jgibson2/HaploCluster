"""
This script takes a tab delimited input file containing genotype and sample numbers and another file containing
information about the samples and outputs the results as a sample number with the desired information.
Script by: Caleb Palagyi (2017)
A.M.D.G.
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-is", "--input_file_cluster", help="Declares the file from which the sample data will be pulled", type=str, required=True)
parser.add_argument("-ii", "--input_file_information", help="Declares the file from which the information data will be pulled", type=str, required=True)
parser.add_argument("-o", "--output_file", help="Declares the file to which the data will be added", type=str, required=True)
args = parser.parse_args()
input_file_name_cluster = args.input_file_cluster
input_file_name_information = args.input_file_information
output_file_name = args.output_file
cluster = open(input_file_name_cluster, 'r')
output = open(output_file_name, 'w')
information = open(input_file_name_information, 'r')
sample_info_dict = {}
for lines in information:
    datum = lines.strip().split('\t')
    sample_info_dict[datum[0]] = [datum[1], datum[2], datum[3]]
for line in cluster:
    sample_list = line.strip().split(' ')
    for sample in sample_list:
        output.write(str(sample))
        output.write('\t')
        output.write(str(sample_info_dict[sample]).replace('\'', '').replace('[', '').replace(']', '').replace(',', ''))
        output.write('\n')
