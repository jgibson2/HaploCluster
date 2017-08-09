"""
This script is designed to take a tab delimited input file of female haplotypes relating to the G6PD allele on the X
Chromosome and output a list of possible Male genotypes that led to this haplotype in a Tab Delimited File.

Coded By: Caleb Palagyi (2017)
A.M.D.G.
"""
import re
import math
import argparse
import matplotlib.pyplot
parser = argparse.ArgumentParser()
parser.add_argument("-im", "--input_file_male", help="Declares the file from which the male data will be pulled", type=str, required=True)
parser.add_argument("-if", "--input_file_female", help="Declares the file from which the female data will be pulled", type=str, required=True)
parser.add_argument("-o", "--output_file", help="Declares the file to which the data will be added", type=str, required=True)
args = parser.parse_args()
input_file_name_male = args.input_file_male
input_file_name_female = args.input_file_female
output_file_name = args.output_file
f = open(input_file_name_female, 'r')
m = open(input_file_name_male, 'r')
male_new = open('g6pd_males_new.hap', 'w')
female_new = open('g6pd_females_new.temp', 'w')
output = open(output_file_name, 'w')
test_file = open('zebra','w')
input_original_males = open('g6pd.merge.combined_positions.positive_strand.haps.males', 'r')
xsum = 0
ysum = 0
counter = {}
genotype_updated = ""


# This function converts genotype/genotype to either a 0, 1, or 2
def double_to_single(genotype_f_func, genotype_updated_func):
    for i in genotype_f_func:
        if i == '0/0':
            genotype_updated_func += '0'
        if i == '1/1':
            genotype_updated_func += '2'
        if i == '0/1':
            genotype_updated_func += '1'
        if i == '1/0':
            genotype_updated_func += '1'
    return genotype_updated_func


# This turns the male haplotypes into a string to make it easier for them to be compared to the female
# genotypes.  This is less readable than the other output file, that's why this temp file is created.
def male_rankings_temp():
    male_geno = open('g6pd_males_new.hap', 'r')
    male_geno_rankings_temp = open('g6pd_males_rankings.temp', 'w')
    for lines in male_geno:
        data_male_geno = lines.strip().split('\t')
        male_geno_list = data_male_geno[0].strip().split(' ')
        male_geno_string = ''.join(str(i) for i in male_geno_list)
        male_geno_rankings_temp.write(male_geno_string)
        male_geno_rankings_temp.write('\t')
        male_geno_rankings_temp.write(data_male_geno[1])
        male_geno_rankings_temp.write('\n')
    male_geno_rankings_temp.close()
    return


# This function takes in female data converted to 0/1/2, converts it to 0/*/1, generates the M1 Template, searches for
# a match in the Male Genotype Dict.  If no match exists, it spits out Error 1.  Then it subtracts match from the 0/1/2
# list, creating M2.  It searches the Male Genotype Dict again and tries to find a match.  If there is no match
def female_to_two_males(key_func):
    all_genotype_scores_males = {}
    results_tup_o = None
    error_list = []
    temp_string = ""
    male_geno_rankings_temp = open('g6pd_males_rankings.temp', 'r')
    zero_loc = [pos for pos, char in enumerate(key_func) if char == '0']
    zero_loc = list(zip(zero_loc, [0 for c in zero_loc]))
    one_loc = ([pos for pos, char in enumerate(key_func) if char == '1'])
    one_loc = list(zip(one_loc, [1 for c in one_loc]))
    two_loc = ([pos for pos, char in enumerate(key_func) if char == '2'])
    two_loc = list(zip(two_loc, [2 for c in two_loc]))
    loc_list = sorted(zero_loc + one_loc + two_loc, key=lambda q: q[0])
    for i in loc_list:
        temp_string += str(i[1])
    male_1_string = ""
    male_2_string = ""
    for x in temp_string:
        if x == '2':
            male_1_string += '1'
        if x == '1':
            male_1_string += '*'
        if x == '0':
            male_1_string += '0'
    matches_score_dict = {}
    for lines in male_geno_rankings_temp:
        data_mr = lines.strip().split('\t')
        all_genotype_scores_males[data_mr[0]] = data_mr[1]
        asterisk_positions = [match.start(0) for match in re.finditer('\*', male_1_string)]
        x = [True if i == j or pos in asterisk_positions else False for pos, (i, j) in enumerate(zip(list(data_mr[0]), list(male_1_string)))]
        if all(x):
            matches_score_dict[data_mr[0]] = data_mr[1]
    if matches_score_dict == {}:
        error_list.append('Error 1')
        print('Error 1')
        return error_list
    else:
        double_wammy = {}
        double_wammy_scores = {}
        for possible_matches in matches_score_dict.keys():
            possible_matches_list = [int(s) for s in possible_matches]
            key_list = [int(v) for v in key_func]
            possible_matches_mates_list = [a - b for a, b in zip(key_list, possible_matches_list)]
            possible_matches_mates_string = ''.join(str(i) for i in possible_matches_mates_list)
            male_geno_rankings_temp_2 = open('g6pd_males_rankings.temp', 'r')
            for lines in male_geno_rankings_temp_2:
                data_mr_2 = lines.strip().split('\t')
                x_two = [True if i == o else False for (i, o) in zip(list(data_mr_2[0]), list(possible_matches_mates_string))]
                if all(x_two):
                    double_wammy[data_mr_2[0]] = possible_matches
                    double_wammy_scores[data_mr_2[0]] = data_mr_2[1]
            if double_wammy == {}:
                v = list(matches_score_dict.values())
                k = list(matches_score_dict.keys())
                top_pick = matches_score_dict[max(matches_score_dict, key=matches_score_dict.get)]
                answer_one = [k for k, v in matches_score_dict.items() if v == top_pick]
                if len(answer_one) > 1:
                    error_list.append('Error 2')
                    return error_list
                else:
                    error_list.append('Error 3')
                    return error_list
            else:
                x = double_wammy.keys()
                y = double_wammy.values()
                rating = {}
                for i in x:
                    k = all_genotype_scores_males[i]
                    l = all_genotype_scores_males[double_wammy[i]]
                    rating[i] = (int(k)+int(l)) * (math.sqrt(int(k)) * math.sqrt(int(l))) / (math.sqrt(1 + int(k) + int(l)))
                max_val = rating[max(rating, key=rating.get)]
                answer = [k for k, v in rating.items() if v == max_val]
                if len(answer) > 1:
                    error_list.append('Error 2')
                    return error_list
                else:
                    g = str(answer)
                    g = g[2:-2]
                    answer_tuple = (g, double_wammy[g])
                    return answer_tuple

# This processes all of the data on males and returns the genotype as a key with the value as the number of samples with
# that genotype.
males_samples_dict = {}
for line in m:
    data_m = line.strip().split('\t')
    genotype_m = str(data_m[0])
    samples_count_m = len(data_m[1].split(' '))
    xsum = xsum + samples_count_m                 # Counts the number of occurences of a sequence
    male_new.write(genotype_m)
    male_new.write('\t')
    male_new.write(str(samples_count_m))
    male_new.write('\n')
    males_samples_dict[''.join(str(i) for i in genotype_m.replace(' ', ''))] = data_m[1]
male_new.close()

# This processes all of the data on females and returns the genotype as a key with the value as the number of samples
# with that genotype.
female_temp_file = open('female_geno_samples.temp', 'w')
females_samples_dict_two = {}
for line in f:
    samples_list = []
    if line[0] == '#':
        continue
    genotype_updated = " "
    data_f = line.strip().split('\t')
    sample = data_f[1]
    samples_count_f = len(data_f[1].split(' '))
    samples_list = sample.strip().split(' ')
    genotype_f = data_f[0].split(' ')
    genotype_updated = double_to_single(genotype_f, genotype_updated)
    test_file.write(str(genotype_updated))
    test_file.write('\n')
    female_new.write(genotype_updated)
    female_new.write('\n')
    females_samples_dict_two[genotype_updated.strip()] = data_f[1]
    for i in range(samples_count_f):
        female_temp_file.write(samples_list[i])
        female_temp_file.write('\t')
        female_temp_file.write(genotype_updated.strip())
        female_temp_file.write('\n')
female_temp_file.close()
female_new.close()
male_rankings_temp()
checker = []
female_new_dos = open('g6pd_females_new.temp', 'r')
for line in female_new_dos:
    data_f_two = line.strip().split('\t')
    key = str(data_f_two[0])
    results_tup = female_to_two_males(key)
    if results_tup is not None:
        if key in checker:
            continue
        else:
            checker.append(key)
            if type(results_tup) == list:
                output.write(key)
                output.write('\t')
                output.write(str(results_tup))
                output.write('\n')
            else:
                output.write(key)
                output.write('\t')
                output.write(results_tup[0])
                output.write('\t')
                output.write(results_tup[1])
                output.write('\n')
output.close()
output = open(output_file_name, 'r')
output_males_new_counts = open('new_g6pd_male_counts.txt', 'w')
haps_output = []
counter_output = {}
output_dict = {}
for line in output:
    data = line.strip().split('\t')
    if data[1] == '[\'Error 1\']' or data[1] == '[\'Error 2\']' or data[1] == '[\'Error 3\']':
        continue
    else:
        haps_output.append(data[1])
        haps_output.append(data[2])
        key_1 = data[1]
        key_2 = data[2]
        output_dict.setdefault(key_1, [])
        output_dict.setdefault(key_2, [])
        output_dict[key_1].append(data[0])
        output_dict[key_2].append(data[0])
for hap in haps_output:
    if hap in counter_output:
        counter_output[hap] += 1
    else:
        counter_output[hap] = 1
for h in counter_output.keys():
    output_males_new_counts.write(h)
    output_males_new_counts.write('\t')
    output_males_new_counts.write(str(counter_output[h]))
    output_males_new_counts.write('\n')
output_males_new_counts.close()
output.close()
line_counter = 0
haps_output_two = []
output = open('new_g6pd_male_counts.txt', 'r')
output_two = open(output_file_name, 'r')
new_total = 0
for line in output:
    data = line.strip().split('\t')
    haps_output_two.append(data[0])
    line_counter = line_counter + 1
for h in haps_output_two:
    if h in males_samples_dict:
        continue
    else:
        if h[0] == '!':
            h = h[1:]
final_dict = {}
counter_final = {}
output_three = open(output_file_name, 'r')
female_temp_file_read = open('female_geno_samples.temp', 'r')
female_final_dict = {}
for line in female_temp_file_read:
    data_fem = line.strip().split('\t')
    sample_id = data_fem[0]
    female_genotype = data_fem[1]
    female_final_dict[sample_id] = female_genotype
haps_output_final = open('g6pd_genotypes_results.res', 'w')
for line in output_three:
    final_data = line.strip().split('\t')
    if final_data[1] == '[\'Error 1\']' or final_data[1] == '[\'Error 2\']' or final_data[1] == '[\'Error 3\']':
        new_total = new_total + 1
        continue
    else:
        existing_value_1 = []
        existing_value_2 = []
        if final_data[1] in final_dict.keys():
            existing_value_1 = final_dict[final_data[1]]
            existing_value_1 = existing_value_1 + ' ' + final_data[0]
            final_dict[final_data[1]] = existing_value_1
        else:
            final_dict[final_data[1]] = final_data[0]
        if final_data[2] in final_dict.keys():
            existing_value_2 = final_dict[final_data[2]]
            existing_value_2 = existing_value_2 + ' ' + final_data[0]
            final_dict[final_data[2]] = existing_value_2
        else:
            final_dict[final_data[2]] = final_data[0]
haps_output_checker = []
for keyz in final_dict.keys():
    female_genotype_list = final_dict[keyz].split(' ')
    haps_output_final.write(str(" ".join(keyz)))
    haps_output_final.write('\t')
    sample_id_list = []
    for female_genotype_x in female_genotype_list:
        for linex in open("female_geno_samples.temp", 'r'):
            if female_genotype_x in linex:
                sample_data = linex.strip().split('\t')
                sample_id_list.append(sample_data[0])
    haps_output_final.write(str(sample_id_list)[2:-2].replace('\', \'', ' '))
    haps_output_final.write('\n')
haps_output_final.close()
haps_output_final = open('g6pd_genotypes_results.res', 'r')
final_output = open('Final_Results', 'w')
last_dict = {}
for linez in haps_output_final:
    data_female = linez.strip().split('\t')
    sample_id_females = data_female[1]
    female_sequence = data_female[0]
    last_dict[female_sequence] = sample_id_females
for keyy in last_dict.keys():
    input_original_males = open('g6pd.merge.combined_positions.positive_strand.haps.males', 'r')
    for line_input in input_original_males:
        male_sep = line_input.strip().split('\t')
        male = male_sep[0]
        male_id = male_sep[1]
        if keyy == male:
            value = last_dict[keyy]
            value = value + ' ' + male_id
            last_dict[keyy] = value
    final_output.write(keyy)
    final_output.write('\t')
    final_output.write(last_dict[keyy])
    final_output.write('\n')
final_output.close()
final_output = open('Final_Results', 'r')
input_original_males = open('g6pd.merge.combined_positions.positive_strand.haps.males', 'r')
last_checker = []
for liney in final_output:
    datum = liney.strip().split('\t')
    last_checker.append(datum[0])
final_output = open('Final_Results', 'a')
for line in input_original_males:
    data = line.strip().split('\t')
    if data[0] in last_checker:
        continue
    else:
         final_output.write(line)

import numpy as np
import matplotlib.pyplot as plt
import os.path
noobs = []
num_o_runs = []
if os.path.exists('data_log.txt'):
    data_log = open('data_log.txt', 'r')
    for line in data_log:
        data = line.strip().split('\t')
        num_o_runs.append(data[0])
        noobs.append(data[1])
    y = num_o_runs[-1:]
    y = int(y[0]) + 1
    num_o_runs.append(str(y))
    new_total = int(new_total)
    noobs.append(str(new_total))
    heights = [int(noobs[i]) for i in range(y)]
    plt.bar([int(x) for x in num_o_runs], height=heights)
    data_log.close()
    data_log_append = open('data_log.txt', 'a')
    data_log_append.write(str(num_o_runs[-1:][0]))
    data_log_append.write('\t')
    data_log_append.write(str(new_total) + '\n')
    data_log_append.close()
else:
    data_log_write = open('data_log.txt', 'w')
    x = np.arange(1)
    plt.xticks(x + .5, [1]);
    plt.bar(x, height=[new_total])
    data_log_write.write('1')
    data_log_write.write('\t')
    data_log_write.write(str(new_total))
    data_log_write.write('\n')
    data_log_write.close()
plt.xlabel('Number of Runs')
plt.ylabel('New Haplotypes')
plt.title(r'Histogram of Haplotypes')

plt.show()
