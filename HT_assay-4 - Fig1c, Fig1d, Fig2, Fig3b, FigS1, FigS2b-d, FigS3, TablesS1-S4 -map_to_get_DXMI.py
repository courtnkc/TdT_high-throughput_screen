'''
This script was part of the pipeline used to generate data shown in Figures 1c, 1d, 2, 3b, S1, S2b-d, and S3. It also produced the data provided in Tables S1-S4.
Inputs for this script are the demultiplexed reads from "HT_assay-3 - Fig1c, Fig1d, Fig2, Fig3b, FigS1, FigS2b-d, FigS3, TablesS1-S4 -demultiplex.py" and reference sequences specific to each library (library 1: pCKC021_ref_20.5.18.txt, library 2: pCKC022_ref_20.5.18.txt, library 3: pCKC023_ref_20.5.18.txt)
Output is the DNA sequence and corresponding alignment characters (D = deletion, X = mismatch, M = match, I = insertion).
'''

import sys
import subprocess
import glob
import multiprocessing as mp
import datetime
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
import itertools as IT

refseqdir = open("pCKC021_ref_20.5.18.txt", 'r').read().splitlines()[0] #specify reference sequence file in this line. pCKC021_ref_20.5.18.txt for library 1, pCKC022_ref_20.5.18.txt for library 2, pCKC023_ref_20.5.18.txt for library3.

def get_alignments(recordList):
    alignmentOutputsList = [] #list that will be populated with tuples of the form: (sequence, alignment characters), referred to as variables a and b, respectively
    for record in recordList:
        this_seq = str(record.seq)
        proc = subprocess.Popen('./mapp ' + this_seq + '  ' + refseqdir, stdout=subprocess.PIPE, shell=True)
        this_aligned = proc.stdout.read().decode('utf-8')
        alignmentOutputsList.append((this_seq,this_aligned))
    return alignmentOutputsList

def DXMI(filename):

    #Set number of sequences to be processed per processor
    numWorkers = mp.cpu_count()
    seqsPerWorker = 10000 #CHANGE THIS TO ADJUST MULTI-PROCESSING SETUP
    counter = 0

    #parallelize work, breaking up fastq into lists of SeqIO fastq records
    recordIter = SeqIO.parse(open(filename), 'fastq')
    pool = mp.Pool(numWorkers)
    for chunk in iter(lambda: list(IT.islice(recordIter, int(seqsPerWorker*numWorkers))), []): #removes chunks from iterator until only [] remains
        counter += len(chunk)
        print(len(chunk))
        chunk = iter(chunk)
        pieces = list(iter(lambda: (list(IT.islice(chunk, seqsPerWorker))), [])) #removes pieces from iterator to be sent off to worker until iterator is empty
        resultsList = pool.map(get_alignments, pieces) #this variable will be a list of outputs from get_alignments, so it will be a list of lists of (sequence, alignment characters) tuples

        #writing results to output file
        filenameBase = filename.split('.f')[0]
        with open(f'aligned_{filenameBase}.txt','a') as alignOut:
            for result in resultsList:
                for alignment in result:
                    alignOut.write(alignment[0] + '\t' + alignment[1] + '\n')

        #reporting progress
        with open(f'aligned_{filenameBase}_REPORT.txt', 'a') as reportOut:
            reportOut.write(f'results written for {counter} seqs @ {datetime.datetime.now()}\n')

    pool.close()
    pool.join()

if __name__ == '__main__':
    DXMI('all_TdT_reads_PassedFilter_Gs_1_1.fastq')
