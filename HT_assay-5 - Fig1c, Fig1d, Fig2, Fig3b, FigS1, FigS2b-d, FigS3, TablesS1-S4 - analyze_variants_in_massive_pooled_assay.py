'''
This script was used to generate data shown in Figures 1c, 1d, 2, 3b, S1, S2b-d, and S3. It also produced the data provided in Tables S1-S4.
The input for this script is NGS reads after merging F and R reads for each cluster, removing short reads, demultiplexing, and aligning to a reference sequence. 
To generate the input for this script, run HT_assay-1... > HT_assay-2... > HT_assay-3... > HT_assay-4...
The output of this script is the A/T bias, insertion fraction, insertion length data, and number of insertions for each TdT variant in the pooled high-throughput assay.
'''

import sys
import os
import numpy as np
import multiprocessing as mp
import datetime
import os
import time
import shutil
import re
import pandas as pd
import statistics
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
from Bio.Seq import translate
from statistics import mean
import seaborn as sns
sns.set(style="ticks", color_codes=True)
import matplotlib.pyplot as plt


##### USER-DEFINED INPUTS #####
max_mismatches = 50  #total number of allowable X's (mismatches) in a given read (this is the concatenated F + R read for each cluster, so it's a 500 bp sequence)
max_deletions = 251  #total number of allowable D's (deletions) in a given read (this is the concatenated F + R read for each cluster, so it's a 500 bp sequence)
max_insertions = 15  #total number of allowable contiguous I's (insertions) (remember 6 I's on either side due to barcodes) in a given read (this is the concatenated F + R read for each cluster, so it's a 500 bp sequence)
NNK_start = 47  #what's the position of the 1st N in the NNKs?  47 for library 1.  37 for library 2.  90 for library 3.
NNK_stop = 55  #what's the position of the last K in the NNKs?  55 for library 1.  48 for library 2.  113 for library 3.
req_NNK_Ms = 2  #how many "M" DXMI characters do you require on either side of the NNKs? This is for an alignment before extracting NNK sequences.
NNK_DXMI = 'XXXXXXXXX'  #what DXMI characters are expected around the NNKs? 'XXXXXXXXX' for library 1. 'XXXMMMXXXXXX' for library 2. 'XXXMMMMMMXXXMMMMMMMMMXXX' for library 3.
cutsite = 38  #position of Cas9 cut site (distance from the 3' end of NGS amplicon).  38 for library 1.  37 for library 2.  125 for library 3.
req_ins_Ms = 2  #how many "M" DXMI characters do you require on either side of the ins? This is for an alignment before extracting edit outcomes at the Cas9 target site.
ins_window = 2  #how many nt on either side of the Cas9 cutsite do you want to check for insertions and deletions? This is for an alignment before extracting edit outcomes at the Cas9 target site.
wt_window = 2  #how many nt on either side of the cutsite need to be 'M' for it to be counted as wt? This is for an alignment before extracting edit outcomes at the Cas9 target site.
activity_cutoff = 0  #what minimum insertion fraction should be the cutoff for TdT variants that are reported at the end of this analysis?
ins_cutoff = 0 #how many insertion sequences associated with a given TdT variant are required for the variant to be reported at the end of this analysis?
##### END OF USER-DEFINED INPUTS #####


#FIRST QUALITY CHECK; THROWS OUT READS WITH TOO MANY X's, D's, or I's
init_seqs = []  #this will hold the DNA sequence for each read
mapped_seqs = []  #this will hold the DXMI characters for each read
low_quality_seqs = []  #this will hold the DNA sequences that get filtered out due to too many sequencing errors / bad alignment
low_quality_DXMI = []  #this will hold the DXMI characters for seqs that get filtered out

#GRABBING ALIGNED READS WITH DXMI NOTATION
input_file = 'aligned_all_TdT_reads_PassedFilter_Gs_all021.txt'  #put the input file name here. My file names during the original analysis used "021", "022", and "023" to denote libraries 1, 2, and 3, respectively
with open(input_file,'r') as aligned_seqs:

    counter = 0
    for line in aligned_seqs:
        counter += 1  
        if counter%2 == 1:
            try:
                seq = line.split("\t")[1][:-1]
            except:
                print('could not read the following line: ' + line)
            if seq.count("X") <= max_mismatches and seq.count("D") <= max_deletions:
                init_seqs.append(line.split('\t')[0])
                mapped_seqs.append(seq)
            else:
                low_quality_seqs.append(line.split('\t')[0])
                low_quality_DXMI.append(seq)

print(str(counter/2) + " total aligned reads in this sample\n")
print(str(len(mapped_seqs)) + " good reads are in this sample, based on DXMI\n")
print(str(len(low_quality_DXMI)) + " bad reads are in this sample, based on DXMI\n")


#BINNING THE READS ACCORDING TO NNK SEQUENCES AND GRABBING INSERTION SEQUENCES, DELETIONS, AND WT READS
low_quality_NNKs = 0
low_quality_ins = 0
no_Ms = 0
binned_by_NNKs = {}									#creating a dictionary that will have NNK sequences as the keys, and insertions as the values
binned_deletions = {}                               #creating a dictionary that will have NNK sequences as the keys, and reads with Cas9-induced deletions as the values
binned_wt = {}                                      #creating a dictionary that will have NNK sequences as the keys, and reads with wt Cas9 targets as the values
binned_other = {}                                   #creating a dictionary that will keep track of Cas9 cut sites that couldn't be classifed as ins, del, or wt
all_NNKs = []
long_insertions = []                                #creating a list that will keep track of insertions that were filtered out because they were too long
bad_NNKs = []                                     

for c1,c2 in zip(init_seqs,mapped_seqs):			                #adjusting the position from which NNKs need to be grabbed due to sequencing errors
    sequencing_insertions = c2[:NNK_start].count('I')
    sequencing_deletions = c2[:NNK_start].count('D')
    
    if sequencing_insertions == 6 and sequencing_deletions == 0: #indicates no sequencing indels upstream of the NNKs; therefore, it can grab the NNKs from expected positions. 6 "insertions" are from the multiplexing barcode.
        potential_NNKs = c1[NNK_start-1:NNK_stop]						#expected NNKs locus
        potential_DXMIs = c2[NNK_start-1:NNK_stop]
        NNK_window = c1[NNK_start-req_NNK_Ms-1:NNK_stop+req_NNK_Ms]		#expected NNKs locus + user-defined # nts on either side
        DXMI_window = c2[NNK_start-req_NNK_Ms-1:NNK_stop+req_NNK_Ms]	#DXMIs corresponding to NNKs window
    
    elif sequencing_insertions == 7 and sequencing_deletions == 0: #indicates one insertion from sequencing error upstream of the NNKs; grab NNKs from modified position
        potential_NNKs = c1[NNK_start:NNK_stop+1]
        potential_DXMIs = c2[NNK_start:NNK_stop+1]
        NNK_window = c1[NNK_start-req_NNK_Ms:NNK_stop+req_NNK_Ms+1]		
        DXMI_window = c2[NNK_start-req_NNK_Ms:NNK_stop+req_NNK_Ms+1]	
    
    elif sequencing_insertions == 6 and sequencing_deletions == 1: #indicates one deletion from sequencing error upstream of the NNKs; grab NNKs from modified position
        potential_NNKs = c1[NNK_start-2:NNK_stop-1]						
        potential_DXMIs = c2[NNK_start-1:NNK_stop]
        NNK_window = c1[NNK_start-req_NNK_Ms-2:NNK_stop+req_NNK_Ms-1]		
        DXMI_window = c2[NNK_start-req_NNK_Ms-1:NNK_stop+req_NNK_Ms]

    else: #reads with more than one insertion or deletion upstream of the NNKs aren't counted. They fail the above alignment steps.
        continue

    upstream_DXMI = DXMI_window[0:req_NNK_Ms]						#DXMIs for nts immediately upstream of NNKs
    downstream_DXMI = DXMI_window[-req_NNK_Ms:]						#DXMIs for nts immediately downstream of NNKs
    matches = 'M'*req_NNK_Ms

    
    if upstream_DXMI == matches and downstream_DXMI == matches:	#ensures there is a user-defined number of "M's" upstream and downstream of the NNKs
        if potential_DXMIs == NNK_DXMI:                 #makes sure NNK locus does not have mismatches inside the region that should be constant (most important for libraries 2 and 3)
            
            if NNK_DXMI == 'XXXMMMXXXXXX': #for library 2, filter out TdT variants that have a mutation between NNKs
                if potential_NNKs[3:6] != 'GCT':
                    continue
            if NNK_DXMI == 'XXXMMMMMMXXXMMMMMMMMMXXX':
                if potential_NNKs[3:9] != 'CAGTTT' or potential_NNKs[12:21] != 'AGAGACCTC': #for library 3, filter out TdT variants that have a mutation between NNKs
                    continue

            if potential_NNKs not in all_NNKs:          #keeping track of all NNKs evaluated (with proper alignment of the NNKs)
                all_NNKs.append(potential_NNKs)

            ins_region = c2[-cutsite-ins_window-1:-cutsite+ins_window]   #now that the NNKs have been verified, check if there is an insertion

            if 'I' in ins_region:                 
                no_Ms = 0                                      #keeps track of how many reads had no 'M' in the ins_window, meaning the alignment failed
                scanner = 0                                         #finds where the insertion starts relative to the 3' end of the read
                if c2[-cutsite+ins_window+1-scanner] != 'I':         
                    while c2[-cutsite+ins_window-scanner] != 'I':
                        scanner += 1
                    first_i = -cutsite+ins_window-scanner
                    
                    ins = 0                                         #counts how long the insertion is
                    while c2[first_i-ins] == 'I':          
                        ins +=1

                    ins_seq = c1[first_i-ins+1:first_i+1]                         #grabs the DNA sequence of the insertion
                    ins_DXMI = c2[first_i-ins-req_ins_Ms+1:first_i+req_ins_Ms+1]       #grabs the DXMI characters of the insertion + user-defined window around it
                    ins_upstream_DXMI = ins_DXMI[0:req_ins_Ms]                  #grabs the DXMI characters upstream of the ins
                    ins_downstream_DXMI = ins_DXMI[-req_ins_Ms:]                #grabs the DXMI characters downstream of the ins
                    ins_matches = 'M'*req_ins_Ms                                #calculates how many M's are required on either side of the ins

                    if ins >= max_insertions:                                   #filtering out insertions that are too long (according to user-defined length cutoff)
                        long_insertions.append(potential_NNKs + '/t' + ins_seq + '/n')
                        continue

                    if ins_upstream_DXMI == ins_matches and ins_downstream_DXMI == ins_matches:     #if the alignment is good around the ins, adds NNKs and ins to dictionary
                        if potential_NNKs not in binned_by_NNKs:	
                            binned_by_NNKs[potential_NNKs] = [ins_seq]
                        else:
                            binned_by_NNKs[potential_NNKs].append(ins_seq)
                    else:
                        low_quality_ins += 1
                        
                else:
                    no_Ms += 1

            elif 'D' in ins_region:         #keeping track of reads that had a deletion at the Cas9 cut site
                if potential_NNKs not in binned_deletions:
                    binned_deletions[potential_NNKs] = [c2]
                else:
                    binned_deletions[potential_NNKs].append(c2) 

            else:
                wt_Ms = wt_window*'M'*2                                             #keeping track of reads that had wt Cas9 cut site
                wt_region = c2[-cutsite-wt_window:-cutsite+wt_window]               #looking at a user-defined window around the cut site to determine if it is truly wt
                if wt_region == wt_Ms:                                      
                    if potential_NNKs not in binned_wt:
                        binned_wt[potential_NNKs] = [c2]
                    else:
                        binned_wt[potential_NNKs].append(c2)

                else:                                                               #keeping track of reads that couldn't be classified as ins, del, or wt
                    if potential_NNKs not in binned_other:
                        binned_other[potential_NNKs] = [c2]
                    else:
                        binned_other[potential_NNKs].append(c2)

        else:
            low_quality_NNKs += 1
            bad_NNKs.append(potential_NNKs)

    else:
        low_quality_NNKs += 1
        bad_NNKs.append(potential_NNKs)

#TALLYING AND REPORTING # OF BINNED NNKs AND INSERTIONS

total_ins = 0
total_del = 0
total_wt = 0
total_other = 0
for a in all_NNKs:
    if a in binned_by_NNKs:
        ins_values = binned_by_NNKs[a]
        total_ins += len(ins_values)
    if a in binned_deletions:
        del_values = binned_deletions[a]
        total_del += len(del_values)
    if a in binned_wt:
        wt_values = binned_wt[a]
        total_wt += len(wt_values)
    if a in binned_other:
        other_values = binned_other[a]
        total_other = len(other_values)


print('############ SUMMARY ############\n')
print('total n low quality NNKs: ' +str(low_quality_NNKs) + '\n')
print('total n low quality ins: ' +str(low_quality_ins) + '\n')
print('total n reads with no M in ins_window: ' + str(no_Ms) + '\n')
print('total n ins: ' + str(total_ins) + '\n')
print('total n deletions: ' + str(total_del) + '\n')
print('total n wt reads: ' + str(total_wt) + '\n')
print('total n other: ' + str(total_other) + '\n')
print('##################################' + '\n')


#### CALCULATING BIASES AND INSERTION LENGTHS! ###

t_all_NNKs = []
t_binned_by_NNKs = {}     #making new dictionaries wherein NNKs are converted to translated amino acid sequences instead of raw DNA sequences    
t_binned_deletions = {}
t_binned_wt = {}


for NNKs in all_NNKs: #BINNING BY *TRANSLATED* NNKs
    raw_NNKs = Seq(NNKs) 
    t_NNKs = raw_NNKs.translate()
    if str(t_NNKs) not in t_all_NNKs:
        t_all_NNKs.append(str(t_NNKs)) #making a list of all translated NNKs

    if NNKs in binned_by_NNKs.keys(): #BINNING INSERTIONS BY *TRANSLATED* NNKs
        for ins in binned_by_NNKs[NNKs]:
            if str(t_NNKs) not in t_binned_by_NNKs.keys():
                t_binned_by_NNKs[str(t_NNKs)] = [ins]
            else:
                t_binned_by_NNKs[str(t_NNKs)].append(ins)
    
    if NNKs in binned_deletions.keys(): #BINNING DELETIONS BY *TRANSLATED* NNKs
        for d in binned_deletions[NNKs]:
            if str(t_NNKs) not in t_binned_deletions.keys():
                t_binned_deletions[str(t_NNKs)] = [d]
            else:
                t_binned_deletions[str(t_NNKs)].append(d)

    if NNKs in binned_wt.keys(): #BINNING WILDTYPE READS BY *TRANSLATED* NNKs
        for w in binned_wt[NNKs]:
            if str(t_NNKs) not in t_binned_wt.keys():
                t_binned_wt[str(t_NNKs)] = [w]
            else:
                t_binned_wt[str(t_NNKs)].append(w)

count_t_NNKs = 0
TdT_ins_data = []  #Making a list of NNK sequences (i.e. distinct TdT variants) and insertions to organize in a Pandas dataframe
list_NNKs_passing_filters = []

for t_NNKs in t_all_NNKs:
    
    n_ins = 0 #Calculating number of reads for the variant
    n_dels = 0
    n_wt = 0

    if t_NNKs in t_binned_by_NNKs.keys():
        n_ins = len(t_binned_by_NNKs[t_NNKs])       #For a given NNK sequence, calculates the number of reads with insertions
    if t_NNKs in t_binned_deletions.keys():
        n_dels = len(t_binned_deletions[t_NNKs])    #For a given NNK sequence, calculates the number of reads with deletions
    if t_NNKs in t_binned_wt.keys():
        n_wt = len(t_binned_wt[t_NNKs])

    n_reads = n_ins + n_dels + n_wt

    percent_ins = 0 #setting default values that will be replaced if the variant has insertions
    avg_length = 0 

    if t_NNKs in t_binned_by_NNKs.keys():
        lengths_list = []
        bias_list = []  #Making a list of insertions and bias scores to organize in a Pandas dataframe

        for ins in t_binned_by_NNKs[t_NNKs]:    #Calculates the average insertion length for the TdT variant
            ins_length = len(ins)   
            lengths_list.append(ins_length)
        avg_length = round(mean(lengths_list),3)

        percent_ins = 1 #If the TdT variant has insertions but no deletions, the activity fraction is 1
        if t_NNKs in t_binned_deletions.keys(): #Calculates the activity fraction
            percent_ins = round((n_ins / (n_dels + n_ins)),3)    

        if percent_ins < activity_cutoff or n_ins < ins_cutoff:      #Filtering out variants that have low activity or very few reads
            continue

        count_ins = 0
        for ins in t_binned_by_NNKs[t_NNKs]: #Calculating bias score (A/T bias) for each insertion
            count_ins += 1
            length = len(ins)                     
            frac_Gs = ins.count('G')/length                   
            frac_Cs = ins.count('C')/length
            frac_Ts = ins.count('T')/length
            frac_As = ins.count('A')/length
            bias_score = (frac_As + frac_Ts) 
            count = 1   #this will be used to keep track of how many insertions of length L go into the A/T bias calculation 
            single_ins_bias = [str(t_NNKs), length, round((bias_score),3), count]
            bias_list.append(single_ins_bias)

        bias_columns = ['NNKs', 'length', 'bias_score', 'n reads']    #Making a dataframe that will be used to calculate the overall bias of the TdT variant
        bias_dataframe = pd.DataFrame(bias_list, columns = bias_columns)    
        sorted_biases = bias_dataframe.sort_values(by=['NNKs','length'],ascending=[True,True]) 
        d = {'bias_score':'avg_bias_score', 'n reads':'sum_n_reads'}    #Calculating the average bias at each ins length, and adding up n insertions of each length
        avg_biases = sorted_biases.groupby(['length'], as_index=False).agg({'bias_score':'mean', 'n reads':'sum'}).rename(columns=d)
        avg_biases['numerator'] = avg_biases['length']*avg_biases['avg_bias_score']*avg_biases['sum_n_reads']
        avg_biases['denominator'] = avg_biases['length']*avg_biases['sum_n_reads']
        score_numerator = avg_biases['numerator'].sum()
        score_denominator = avg_biases['denominator'].sum()
        weighted_bias = score_numerator/score_denominator #This is the final formula to calculate the weighted bias for the NNK (weighted bias = total A/T bias across all inserted nts)
        
        for ins in t_binned_by_NNKs[t_NNKs]:  #making a larger dataframe with information about every insertion
            length = len(ins)
            ins_counter = 1                        #Will be used to count the total number of insertions
            unique_counter = 1                      #Will be used to count the number of unique insertions                 
            Gs = ins.count('G')                   #Calculating number of each nt in the insertion
            Cs = ins.count('C')
            Ts = ins.count('T')
            As = ins.count('A')
            single_NNK_ins = [str(t_NNKs), str(ins), length, round((Gs/length),3), round((Cs/length),3), round((Ts/length),3), round((As/length),3), avg_length, unique_counter, ins_counter, n_ins, percent_ins, round((weighted_bias),3), n_reads]
            TdT_ins_data.append(single_NNK_ins)

    if percent_ins < activity_cutoff or n_ins < ins_cutoff:      #If the user wants, they can filter out variants that had no insertions 
        continue  

    if t_NNKs not in list_NNKs_passing_filters:
        list_NNKs_passing_filters.append(t_NNKs)
        count_t_NNKs += 1   #Counting the number of NNKs that pass the filters

    if n_ins == 0:
        weighted_bias = 0

    single_NNK_ins = [str(t_NNKs), 'n/a', 0, 0, 0, 0, 0, avg_length, 0, 0, n_ins, percent_ins, round((weighted_bias),3), n_reads] #adding one row to the dataframe specifically to capture variants that had no insertions
    TdT_ins_data.append(single_NNK_ins)

columns = ['NNKs', 'ins', 'length', 'frac_G', 'frac_C', 'frac_T', 'frac_A', 'avg length', 'n unique ins', 'n ins of length L', 'tot n ins', 'percent ins', 'A-T bias score', 'tot n reads']
complete_dataframe = pd.DataFrame(TdT_ins_data, columns = columns)
complete_dataframe.to_csv(input_file + '_COMPLETE_DATAFRAME.csv') #this dataframe gives the bias of every single insertion, plus info about the relevant TdT variant

### --- Making a separate dataframe to report n unique insertions of each length for all NNKs --- ###

df_unique_ins = complete_dataframe.groupby(['NNKs','ins'], as_index=False).agg({'length':'mean','avg length':'mean','n unique ins':'mean','n ins of length L':'sum','tot n ins':'mean','A-T bias score':'mean', 'tot n reads':'mean'})    
df_unique_ins.to_csv(input_file + '_DATAFRAME_UNIQUE_INSERTIONS.csv')

df_count_unique_ins = df_unique_ins.groupby(['NNKs','length'], as_index=False).agg({'avg length':'mean','n unique ins':'sum','n ins of length L':'sum','tot n ins':'mean','A-T bias score':'mean', 'tot n reads':'mean'})
df_count_unique_ins.to_csv(input_file + '_DATAFRAME_UNIQUE_INSERTIONS_GROUPED_BY_LENGTH.csv')


### --- Sorting NNK variants by average insertion length, calculating bias at each length, and calculating n insertions of each length--- ###

grouped_lengths = complete_dataframe.groupby(['NNKs','length'], as_index=False).agg({'frac_G':'mean','frac_C':'mean','frac_T':'mean','frac_A':'mean','avg length':'mean','n ins of length L':'sum','tot n ins':'mean','percent ins':'mean','A-T bias score':'mean','tot n reads':'mean'})
grouped_lengths.to_csv(input_file + '_DATAFRAME_NNKs_lengths_GROUPED.csv')



### --- CALCULATING n NNKs IN THE LIBRARIES --- ###

total_DNA_NNKs = len(all_NNKs)
total_translated_NNKs = len(list_NNKs_passing_filters)

print('Total n NNKs (DNA sequences, good NNK and ins alignments): ' + str(total_DNA_NNKs) +'\n')
print('Total n NNKs (translated, good NNK and ins alignments): ' + str(len(t_all_NNKs)) +'\n')
print('Total n NNKs (translated, after activity and ins filter): ' + str(total_translated_NNKs) + '\n')


### --- CALCULATING DISTRIBUTION OF AVERAGE INSERTION LENGTHS --- ###

collapsed_dataframe = complete_dataframe.groupby('NNKs', as_index=False, sort=False).mean() #making a dataframe that yields 1 row per NNK
dataframe_NNKs = list(collapsed_dataframe['NNKs']) #grabs the NNKs 

ins_lengths = []
for NNKs in list_NNKs_passing_filters:
    grabbed_length = (collapsed_dataframe.loc[collapsed_dataframe['NNKs'] == NNKs, 'avg length']).item()
    ins_lengths.append(grabbed_length)
sns.distplot(ins_lengths,kde=False)
plt.xlabel('average ins length (nt)')
plt.ylabel('n NNKs')
plt.savefig(input_file + '_HISTOGRAM-AVG_LENGTHS_vs_NNKs.pdf')
