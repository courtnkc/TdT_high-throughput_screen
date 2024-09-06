"""
This script was part of the pipeline used to generate data shown in Figures 1c, 1d, 2, 3b, S1, S2b-d, and S3. It also produced the data provided in Tables S1-S4.
The inputs for this script were the raw fastq files from the high-throughput TdT variant assay (in .txt format here).

Given a list of pairs of .fastq files, this script concatenates the forward and reverse reads for each cluster.

The output is a fastq file with the reads merged and associated quality scores.

RUN WITH PYTHON3
"""

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
from Bio.Alphabet import IUPAC
import numpy as np
import multiprocessing as mp
import datetime
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import time
import shutil
import re
local = True

##### USER INPUTS #####

mergeList = [('nR045-L1-P1-ATTACTCG-TATAGCCT-READ1-Sequences.txt','nR045-L1-P1-ATTACTCG-TATAGCCT-READ2-Sequences.txt'),('nR045-L1-PrNotRecog-READ1-Sequences.txt','nR045-L1-PrNotRecog-READ2-Sequences.txt'),('nR045-L2-P1-ATTACTCG-TATAGCCT-READ1-Sequences.txt','nR045-L2-P1-ATTACTCG-TATAGCCT-READ2-Sequences.txt'),('nR045-L2-PrNotRecog-READ1-Sequences.txt','nR045-L2-PrNotRecog-READ2-Sequences.txt')] #list of tuples for all fwd,rvs sequences as .fastq file names that are to be merged
fwdRemove = 0 #number of nucleotides to be removed from the end of the forward reads
rvsRemove = 0 #number of nucleotides to be removed from the beginning of the reverse reads
dirname = r'/pub/courtnkc/courtnkc20032557' #directory where the input data lives

##### END USER INPUTS #####


print(datetime.datetime.now())

#changes current directory to dirname, copies script to that directory with datestamp

if local == True:
    os.chdir(dirname)
    copied_script_name = time.strftime("%Y-%m-%d") + '_' + os.path.basename(__file__)
    shutil.copy(__file__, dirname + os.sep + copied_script_name)

#given two fastq files of NGS reads, returns a new .fastq file containing merged sequences and quality scores
#removes 3rd argument number of characters from the forward reads, removes 4th argument number of characters from the reverse reads

def mergeFqs(fastq1, fastq2, fwdRmv, rvsRmv):
    L = -1
    fqMerge = open('all_TdT_reads_FRmerge_20.5.22.fastq','a') #creates new fastq file to add merged sequences to
    with open(fastq1,'r') as fq1, open(fastq2,'r') as fq2:
        for line1, line2 in zip(fq1,fq2): #simultaneously loop through lines in two fastqs
            L += 1
            if L%4 == 0:
                fqMerge.write(line1)
            elif L%4 == 1: #combine forward read with reverse complement of reverse read, removing specified number of nucleotides from each
                line2Seq = Seq(line2[:-1], IUPAC.unambiguous_dna)
                line2Rcomp = line2Seq.reverse_complement()
                line2Str = str(line2Rcomp)[rvsRmv:]+'\n'
                fqMerge.write(line1[:-(fwdRmv+1)]+str(line2Str))
            elif L%4 == 2:
                fqMerge.write(line1)
            elif L%4 == 3:
                rvsLine2 = line2[::-1]
                fqMerge.write(line1[:-(fwdRmv+1)]+rvsLine2[1+rvsRmv:]+'\n')
    fqMerge.close()
    return

for pair in mergeList:
    print('starting new pair')
    mergeFqs(pair[0],pair[1],fwdRemove,rvsRemove)
print('done')

