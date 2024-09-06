'''
This script was used to analyze the data show in Figures 2, 3b, 4b-g, and S3.
Given the results from "Fig2, Fig3b, Fig4b-g, FigS2a, FigS3 - analyze_insertion_compositions_at_Cas9_target_sites.py", this script compiles the results for every sample into a single CSV, for convenience.
'''

import pandas as pd

##### USER INPUTS #####

sample_numbers = [1, 3, 4, 5, 8, 9, 21, 22] #this corresponds with the way the output files from "Fig2, Fig3b, Fig4b-g, FigS2a, FigS3 - analyze_insertion_compositions_at_Cas9_target_sites.py" are named. Numbers in this list were adjusted to capture all samples.
range_of_indel_sizes = [-25, 30]

##### END USER INPUTS #####


single_sample_data = []
all_sample_data = []

for sample in sample_numbers:
    input_file_name = 'TdT4-quality_filtered-{}_indel_histogram'.format(sample) 
    sample_file = open('{}.csv'.format(input_file_name),'r').read().splitlines() 

    single_sample_data.append(sample)

    counter = 0
    for line in sample_file: 
        if counter == 0:
            counter +=1
            continue
        else:
            interpretable_line = line.split(',')
            single_sample_data.append(interpretable_line[1])


    input_file_name_2 = 'TdT4-quality_filtered-{}_insertion_compositions'.format(sample)
    sample_file_2 = open('{}.csv'.format(input_file_name_2),'r').read().splitlines() 

    counter = 0
    for line in sample_file_2:
        if counter == 0:
            counter += 1
            continue
        if counter == 1:
            interpretable_line_2 = line.split(',')
            single_sample_data.append(interpretable_line_2[7])
            single_sample_data.append(interpretable_line_2[8])
            single_sample_data.append(interpretable_line_2[9])
            single_sample_data.append(interpretable_line_2[10])
            single_sample_data.append(interpretable_line_2[11])
            counter += 1
        if counter > 1:
            break

    all_sample_data.append(single_sample_data)
    single_sample_data = []

columns = ['sample']
for n in range(range_of_indel_sizes[0], range_of_indel_sizes[1]+1):
    columns.append(str(n))

columns.append('avg ins length')
columns.append('ins / (ins + del)')
columns.append('editing efficiency')
columns.append('total n ins')
columns.append('overall A/T bias')

df_compiled_data = pd.DataFrame(all_sample_data, columns = columns)
df_compiled_data.to_csv('COMPILED_all_samples_results.csv')







        