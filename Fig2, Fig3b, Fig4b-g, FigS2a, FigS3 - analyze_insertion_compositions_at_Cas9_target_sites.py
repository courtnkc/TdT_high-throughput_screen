'''
This script was used to analyze the data show in Figures 2, 3b, 4b-g, S2a, and S3.
Inputs for this script are demultiplexed fastq files from Illumina sequencing runs, which were quality filtered using the fastq-filter tool (https://pypi.org/project/fastq-filter/)
Demultiplexing was conducted with separate_by_barcodes_modified.py, which can be found at https://github.com/liusynevolab/peCHYRON .
Outputs from this script are average %A/T in the observed insertions, average insertion length, insertion fraction, overall editing efficiency, and a histogram of indel sizes.
'''
import Bio
from Bio.Seq import Seq
import pandas as pd 

#This script reads in fastq files, applies a quality filter to the reads, then profiles edits at the Cas9 target site.

########## USER INPUTS ##########

approx_ins_location_from_5prime = 39       #how far from the beginning of the F read is 1st nt of the insertion? 340 for site3. 280 for CHYRON locus (a.k.a. 680c3).
ins_location_sliding_window = 20    #how many nt upstream and downstream of the expected target site do you want to check for an alignment? 150 for site3. 100 for CHYRON locus.
upstream_seq = 'AGGACGAAACACCGGCCCAGACTGAGCACG'         #CCCTGGCCTGGGTCAATCCTTGGGGCCCAGACTGAGCACG for site3. AGGACGAAACACCGGCCCAGACTGAGCACG for CHYRON locus. What's the sequence immediately upstream of the Cas9 target site? Provide at least 30 nts.
downstream_seq = 'TGAGGGTTAGAGCTAGAAATAGCAAGTTAA'       #TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAGGGC for site3. TGAGGGTTAGAGCTAGAAATAGCAAGTTAA for CHYRON locus. What's the sequence immediately downstream of the Cas9 target site? Provide at least 30 nts.
alignment_length = 15               #how many nts of upstream_seq and downstream_seq do you want to align when searching for the target site?
HD_upstream_seq = 3                 #what's an acceptable Hamming Distance between the expected sequence upstream of the cutsite and the actual sequence? This will be used to check for alignment before extracting the insertion sequences.
HD_downstream_seq = 3               #what's an acceptable Hamming Distance between the expected sequence downstream of the cutsite and the actual sequence? This will be used to identify the end of each insertion sequence.
max_edits = 30                      #what's the maximum number of inserted nts you assume is generated via TdT insertion?
max_del_size = -25

input_files_constant_string = 'TdT4-quality_filtered-' #how are the input fastq files named? This is set up to process multiple samples simultaneously that share a common string in their names but differ by sample numbers appended to the end of the file name.
input_files_sample_numbers = [1, 3, 4, 5, 8, 9, 21, 22]

########## END USER INPUTS ##########


########## DEFINING FUNCTIONS ##########

def calc_H_dist(seq1,seq2): #defining a function to calculate hamming distances between two sequences
    distance = 0
    for a,b in zip(seq1,seq2):
        if a != b:
            distance += 1
    return distance


def grab_insertions(input_seq, input_upstream_seq, input_ins_location, input_downstream_seq, input_HD_downstream_seq, input_HD_upstream_seq, input_max_edits): #determines how long insertion sequences are, then grabs the insertion sequences
    grabbed_ins = [] #this will report 1) whether a good alignment was found and 2) the actual insertion sequence if applicable
    
    upstream_align_position = input_ins_location - alignment_length 
    
    found_alignment = 'false'
    
    actual_upstream = input_seq[upstream_align_position:input_ins_location] #checks alignment upstream of the target site
    if calc_H_dist(actual_upstream, input_upstream_seq) <= input_HD_upstream_seq: #checking if the expected sequence upstream of the edit site is present at the expected position
        ins_location = input_ins_location
        found_alignment = 'true' 

    if found_alignment == 'false': #checking for the expected sequence upstream of the target site, using a scanning window to detect if it is present at a different position than usual
        for i in range(-ins_location_sliding_window, ins_location_sliding_window):
            upstream_align_position = input_ins_location - alignment_length + i
            actual_upstream = input_seq[upstream_align_position:upstream_align_position + alignment_length]
            if calc_H_dist(actual_upstream, input_upstream_seq) <= input_HD_upstream_seq:
                found_alignment = 'true'
                ins_location = input_ins_location + i    
                break

    if found_alignment == 'false': #if an alignment wasn't found above, try aligning the reverse complement of the read
        rev_comp_seq = seq.reverse_complement()
        upstream_align_position = input_ins_location - alignment_length

        actual_upstream = rev_comp_seq[upstream_align_position:input_ins_location] #checks alignment upstream of the target site
        if calc_H_dist(actual_upstream, input_upstream_seq) <= input_HD_upstream_seq: #checking if the expected sequence upstream of the edit site is present at the expected position
            ins_location = input_ins_location
            found_alignment = 'true' 
            input_seq = rev_comp_seq

    if found_alignment == 'false': #checking for the expected sequence upstream of the target site in the REV COMP of the read, using a scanning window to detect if it is present at a different position than usual
        for i in range(-ins_location_sliding_window, ins_location_sliding_window):
            upstream_align_position = input_ins_location - alignment_length + i
            actual_upstream = rev_comp_seq[upstream_align_position:upstream_align_position + alignment_length]
            if calc_H_dist(actual_upstream, input_upstream_seq) <= input_HD_upstream_seq:
                    found_alignment = 'true'
                    ins_location = input_ins_location + i
                    input_seq = rev_comp_seq
                    break
                   
    if found_alignment == 'false':
        grabbed_ins = ['failed alignment', input_seq]

    if found_alignment == 'true':   
        iterate = 'yes'
        end_edit = ins_location #this will hold the location corresponding to the end of the edited target
        nts_inserted = 0

        while iterate == 'yes': #if the locus is edited, start scanning downstream for the expected sequence after the target site. Determine how long the insertion is.
            if nts_inserted == input_max_edits: #a maximum insertion length is set by the user
                iterate = 'no'
                found_alignment = 'false'
                grabbed_ins = ['failed alignment', input_seq]
            elif calc_H_dist(input_seq[end_edit:end_edit+alignment_length], input_downstream_seq) <= input_HD_downstream_seq: #checks if the expected sequence downstream of the target site is present at the expected location for a 20-bp insertion
                iterate = 'no'
            else:
                nts_inserted += 1
                end_edit += 1

        if nts_inserted >= 1 and found_alignment == 'true':
            grabbed_ins_seq = str(input_seq[ins_location: ins_location + nts_inserted]) 
            grabbed_ins = ['good alignment', grabbed_ins_seq]

        if nts_inserted == 0 and found_alignment == 'true':
            grabbed_ins = ['good alignment', ''] #indicates 0nt insertion
            
    return grabbed_ins

########## END OF FUNCTION DEFINITIONS ##########

### GRABBING INSERTION SEQUENCES FROM RAW READS ###

first_alignment_upstream_seq = upstream_seq[-alignment_length:]
first_alignment_downstream_seq = downstream_seq[:alignment_length]

second_alignment_upstream_seq = upstream_seq[-alignment_length - 10: -10]
second_alignment_downstream_seq = downstream_seq[10: 10 + alignment_length]

for sample in input_files_sample_numbers:
    
    list_of_sequences = []

    input_file_name = '{}{}'.format(input_files_constant_string, sample) #this will be used to name the output files
    F_sample_file = open('{}.fastq'.format(input_file_name),'r').read().splitlines() #input file name here (fastq F reads)

    line_counter = 0
    for line in F_sample_file: #extracting the DNA sequences from the fastq input file
        line_counter += 1
        if line_counter %4 == 2:
            if len(line) > approx_ins_location_from_5prime + ins_location_sliding_window:
                list_of_sequences.append(line)
    
    counter = 0
    all_insertion_data = [] #will be used to make a pandas dataframe detailing the alignment and insertion-grabbing results
    for sequence in list_of_sequences:
        counter += 1

        insertion_info = grab_insertions(sequence, first_alignment_upstream_seq, approx_ins_location_from_5prime, first_alignment_downstream_seq, HD_downstream_seq, HD_upstream_seq, max_edits) #determines how long insertion sequences are, then grabs the insertion sequences
        insertion_info.append('1st alignment')

        if insertion_info[0] == 'failed alignment': #if the target site wasn't found in the first attempt at alignment, it may be due to a partial deletion of the cut site. Check for alignment on either side of the cut site.

            insertion_info = grab_insertions(sequence, second_alignment_upstream_seq, approx_ins_location_from_5prime, second_alignment_downstream_seq, HD_downstream_seq, HD_upstream_seq, max_edits+20) #determines how long insertion sequences are, then grabs the insertion sequences
            insertion_info.append('2nd alignment')

        all_insertion_data.append(insertion_info)

    columns = ['description', 'grabbed insertion (or full seq if no insertion was grabbed)', 'which alignment was successful?']
    df_insertion_data = pd.DataFrame(all_insertion_data, columns = columns)
    df_insertion_data.to_csv(input_file_name + '_extracted_insertions.csv')

    
    ### MAKING AN INDEL HISTOGRAM ###
    
    list_of_lengths = []
    for l in range(max_del_size, max_edits + 1):
        list_of_lengths.append(l)

    overall_ins_count = 0 #will calculate total n insertions and deletions, for activity score
    overall_del_count = 0
    overall_wt_count = 0
    total_nts_inserted = 0 #will calculate total nts inserted, so I can get the average insertion length
    list_of_insertions = []

    single_row_indel_count = []
    all_rows_indel_counts = []
    for index, row in df_insertion_data.iterrows():
        indel_length = 0 

        if row['description'] == 'good alignment':

            if row['which alignment was successful?'] == '1st alignment':
                indel_length = len(row['grabbed insertion (or full seq if no insertion was grabbed)'])
                list_of_insertions.append(row['grabbed insertion (or full seq if no insertion was grabbed)'])

            if row['which alignment was successful?'] == '2nd alignment': #adjusting insertion length for 'insertions' that were grabbed after aligning 10 nt upstream and downstream of the cut site
                indel_length = len(row['grabbed insertion (or full seq if no insertion was grabbed)']) - 20       
                if indel_length > 0:
                    trimmed_ins = row['grabbed insertion (or full seq if no insertion was grabbed)'][10:-10]
                    list_of_insertions.append(trimmed_ins)

            if indel_length > 0 and indel_length <= max_edits: #tallying the number of deletions and insertions; this will be used later for insertion fraction calculations

                overall_ins_count += 1
                total_nts_inserted += indel_length

            if indel_length < 0 and indel_length >= -max_edits:
                overall_del_count += 1

            for length in list_of_lengths:
                if length == indel_length:
                    single_row_indel_count.append(1)
                    if length == 0:
                        overall_wt_count += 1
                else:
                    single_row_indel_count.append(0)

            all_rows_indel_counts.append(single_row_indel_count)
            single_row_indel_count = []

    columns = list_of_lengths
    df_indel_data = pd.DataFrame(all_rows_indel_counts, columns = columns)
    df_indel_histogram = df_indel_data.sum(axis = 0)
    df_indel_histogram.to_csv(input_file_name + '_indel_histogram.csv')


    ### CALCULATING THE FRACTION OF EDITED READS THAT CONTAIN AN INSERTION, AND AVERAGE INSERTION LENGTH ###

    if (overall_ins_count + overall_del_count) > 0:
        activity_score = overall_ins_count / (overall_ins_count + overall_del_count)
    else:
        activity_score = 0

    if overall_ins_count > 0:
        avg_ins_length = total_nts_inserted / overall_ins_count
    else:
        avg_ins_length = 0

    total_editing_efficiency = (overall_ins_count + overall_del_count) / (overall_ins_count + overall_del_count + overall_wt_count)


    ### CALCULATING THE COMPOSITION OF EACH INSERTION LENGTH ###
    
    tot_n_Gs = 0
    tot_n_Cs = 0
    tot_n_Ts = 0
    tot_n_As = 0

    bias_list = []
    for length in range(1,max_edits + 1):
        n_ins_of_length_L = 0
        for insertion in list_of_insertions:
            if len(insertion) == length:
                n_ins_of_length_L += 1
                G_count = insertion.count('G')    
                C_count = insertion.count('C')
                T_count = insertion.count('T')
                A_count = insertion.count('A')

                frac_Gs = round(G_count/length , 3)                #Calculating %A/T for each insertion
                frac_Cs = round(C_count/length , 3)
                frac_Ts = round(T_count/length , 3)
                frac_As = round(A_count/length , 3)

                bias_score = (frac_As + frac_Ts) 
                count = 1   #this will be used to keep track of how many insertions of length L go into the weighted bias calculation
                single_ins_bias = [length, frac_Gs, frac_Cs, frac_Ts, frac_As, count, avg_ins_length, activity_score, overall_ins_count, total_editing_efficiency]
                bias_list.append(single_ins_bias)
                
                tot_n_Gs += G_count
                tot_n_Cs += C_count
                tot_n_Ts += T_count
                tot_n_As += A_count

        if len(bias_list) == 0:
            single_ins_bias = [0, 0, 0, 0, 0, 0, avg_ins_length, activity_score, overall_ins_count, total_editing_efficiency]  
            bias_list.append(single_ins_bias)          

    if (tot_n_Gs + tot_n_Cs + tot_n_Ts + tot_n_As) > 0:
        overall_A_T_bias = (tot_n_Ts + tot_n_As) / (tot_n_Gs + tot_n_Cs + tot_n_Ts + tot_n_As)
    else:
        overall_A_T_bias = 0

    for single_ins_data in bias_list:
        single_ins_data.append(overall_A_T_bias)
    
    columns = ['length', 'frac_G', 'frac_C', 'frac_T', 'frac_A', 'n ins of length L', 'avg ins length', 'ins / (ins + del)', 'total n ins', 'editing efficiency', 'overall A/T bias']
    df_bias_breakdown = pd.DataFrame(bias_list, columns = columns)
    
    df_summed_bias_breakdown = df_bias_breakdown.groupby(['length'], as_index=False).agg({'frac_G':'mean','frac_C':'mean','frac_T':'mean','frac_A':'mean','n ins of length L':'sum', 'avg ins length': 'mean', 'ins / (ins + del)': 'mean', 'editing efficiency': 'mean', 'total n ins': 'mean', 'overall A/T bias': 'mean'})
    df_summed_bias_breakdown.to_csv(input_file_name + '_insertion_compositions.csv')
   


        







