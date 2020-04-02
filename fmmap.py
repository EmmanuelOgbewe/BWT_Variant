#fmmap.py
#This is our main code for our basic read alignment. This tool takes
#in two commands index and align. 
# @author Emmanuel Ogbewe

import sys
import os 
import pysam
import helpers
import gzip
import numpy as np
import math 
import readfq as rf

#------Handle input------#
def parseCommandLineArgs(args=[]):
    if 'index' in args:
        # run the index function 
        print("running index function")
        (data,sequenceName) = readFastaFormat()
        index(sequenceName,data)
    elif 'align' in args:
        #run the align function
        print("running align function")
        align()

    else:
        print("Error parsing commanding arguments, invalid command.")

#------Helper Functions-----#
def readFastaFormat():
    result = ''
    with open('data/2019-nCov.fa') as f:
        for line in f.readlines():
            if ">" in line:
                sequenceName = line[1:]
                pass
            else:
                result += line.strip()
    result = result.strip()
    result += "$"
    return (result,sequenceName)

def readGzip():
    with gzip.open("data/reads.fa.gz") as f:
        file_content = ''
        count = 0
        for i in f.readlines():
            if count < 10:
                print(str(i))
                count += 1 
            else:
                break 

    return file_content

def genericWriteToOutput(input,fileName):
    with open(fileName,'w') as f:
        f.write(input)

# align helpers 


#------Start command Functions-------#
"""
Takes two arguments is the path to a FASTA format file and the 
second is the location of an output file where the index will be written.
When invoked in this mode, fmmap will read the input file and will compute
the FM-index of the corresponding string.
retuns index containing the following:
"""
def index(name,s):
    # first build the suffix array 
    # get the BWM and get the first and last column
    name = name 
    length = len(s)
    refString = s
    firstCol = helpers.generateFirstCol(s)
    occurenceTable = helpers.firstOccurTable(firstCol)
    (_,sa) = helpers.suffixArray(s)
    # genericWriteToOutput(str(occurenceTable))
    lastCol = helpers.buildBWT(s)
  
    with open('ref_index','w') as f:
        f.write('name: \n' + name + '\n')
        f.write('reference_string: \n' + refString + '\n')
        f.write('length: \n' + str(length) + '\n')
        f.write('first column: \n' + str(firstCol) + '\n')
        f.write('bwt: \n' + str(lastCol) + '\n')
        f.write('occurence table: \n' + str(occurenceTable) + '\n')
        f.write('suffix array: \n' + str(sa) + '\n')

def align(s):
    # outputFile = 'alignments.sam'
    bwt_index = helpers.BWTIndex()
    
    ninf = float("-inf")
    seed_skip = lambda l: math.floor(l / 5.0)
    gap = 5

    with gzip.open('data/reads.fa.gz', 'rt') as rfile:
        for name, seq, qual in rf.readfq(rfile):
            alignments = []
            read_len = len(seq)
            best_score = ninf
            seed_pos = 0
            skip = seed_skip(read_len)
            for seed_start in range(0, read_len, skip):
                seed_end = min(read_len, seed_start + skip) 

                interval, match_len = bwt_index.get_interval(seq[seed_start:seed_end])

    

#------ End command functions-------#


if __name__ == "__main__":
    #handle command line arguments
   parseCommandLineArgs(sys.argv)
   