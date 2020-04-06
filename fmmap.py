#fmmap.py
#This is our main code for our basic read alignment. This tool takes
#in two commands index and align. 
# @author Emmanuel Ogbewe

import sys
import os 
import pysam
import helpers
import gzip
import math 
import re
import readfq as rf
import timer


#-----Paths sepcified by users ---- #
PATH_FASTA_FILE = ''
PATH_INDEX_OUTPUT = ''

PATH_TO_ALIGNMENTS = ''
PATH_SAM_OUTPUT = ''


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
        #refractor to take in path from command line
        align('data/reads.fa.gz','alignments.sam')

    else:
        print("Error parsing commanding arguments.")
        if 'index' in args or 'align' in args:
            print('Please verify that the path specificed is valid.')

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
def generateCigarTuple(finalCigar):
    tup = ()
    pattern = re.compile(r'([0-9]+)([A-Z]+)')

    for (numbers,letters,) in re.findall(pattern, finalCigar):
        if letters == 'M':
            num = 0
        elif letters == 'I':
            num = 1
        else:
            num = 2

        if not tup == ():
            tup = tup + ((num,int(numbers)),)      
        else:
            tup = ((num,int(numbers)),)
    return tup 

def write_to_sam(output_file,alignment, read_name, query_sequence):
   

    cig = generateCigarTuple(alignment.cigar)
    # print(cig)
    # print(alignment.posInRef)
    # return 
    header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 29882, 'SN':'MN988713.1 Severe acute respiratory syndrome coronavirus 2 isolate 2019-nCoV/USA-IL1/2020, complete genome'}] }

    with pysam.AlignmentFile('assignments.sam', "wb", header=header) as outf:
        a = pysam.AlignedSegment()
        a.query_name = read_name
        a.query_sequence= query_sequence
        a.flag = 0
        a.reference_id = 0
        a.reference_start = alignment.posInRef + 1
        a.mapping_quality = 255
        a.cigar = cig
        # a.next_reference_id = '0'
        # a.next_reference_start= '0'
        a.template_length= len(query_sequence) / 100
        # a.query_qualities = '*'
        # a.tags = (("NM", 1),
        #         ("RG", "L1"))
        outf.write(a)
 


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
    lastCol = helpers.buildBWT(s)
    occurenceTable = helpers.occurenceTable(None,lastCol)
    (_,sa) = helpers.suffixArray(s)
    # genericWriteToOutput(str(occurenceTable))

  
    with open('ref_index','w') as f:
        f.write('name: \n' + name + '\n')
        f.write('reference_string: \n' + refString + '\n')
        f.write('length: \n' + str(length) + '\n')
        f.write('first column: \n' + str(firstCol) + '\n')
        f.write('bwt: \n' + str(lastCol) + '\n')
        f.write('occurence table: \n' + str(occurenceTable) + '\n')
        f.write('suffix array: \n' + str(sa) + '\n')

def align(s,output_file):
    # outputFile = 'alignments.sam'
    bwt_index = helpers.BWTIndex()
    
    ninf = float("-inf")
    seed_skip = lambda l: math.floor(l / 5.0)
    gap = 5
    count = 0
    header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 29882, 'SN':'MN988713.1 Severe acute respiratory syndrome coronavirus 2 isolate 2019-nCoV/USA-IL1/2020, complete genome'}] }

    with pysam.AlignmentFile('assignments.sam', "w", header=header) as outf:
        with gzip.open(s, 'rt') as rfile:
            time = timer.Timer()
            time.start()
            for name, seq, qual in rf.readfq(rfile):
                if count == 40000:
                    time.stop()
                    break
                count += 1
                seq = "CTTCTTAGAGGGAGAAACACTTCCCACAGAAGTGTTAACAGAGGAAGTTGTCTTGAAAACTGGTGATTTACAACCATTAGAACAACCTACTAGTGAAGCT"
                print(count)
                seq = seq.replace("N", "A")
                alignments = []
                read_len = len(seq)
                best_score = ninf
                seed_pos = 0
                skip = seed_skip(read_len)
                for seed_start in range(0, read_len, skip):
                    seed_end = min(read_len, seed_start + skip) 
                    match_len = 0
                    interval = bwt_index.get_interval(seq[seed_start:seed_end])
                if interval == None or interval == 0:
                    continue
                for ref_pos in bwt_index.ref_positions(interval, seed_end, match_len):
                    # print(seq)
                    alignment = bwt_index.fitting_alignment(seq,ref_pos,gap)
                    if alignment.score > best_score:
                        best_score = alignment.score
                        alignments = [alignment]
                    elif alignment.score == best_score:
                        # print(seq)
                        alignments.append(alignment) 
                break
                for al in alignments:
                    cig = generateCigarTuple(al.cigar)
                    a = pysam.AlignedSegment()
                    a.query_name = name
                    a.query_sequence= seq
                    a.flag = 0
                    a.reference_id = 0
                    a.reference_start = alignment.posInRef + 1
                    a.mapping_quality = 255
                    a.cigar = cig
                    # a.next_reference_id = '0'
                    # a.next_reference_start= '0'
                    a.template_length= len(seq) / 100
                    outf.write(a) 
                
                # for a in alignments:
                #     write_to_sam(output_file, a, name,seq)
                
#------ End command functions-------#


if __name__ == "__main__":
    #handle command line arguments
   parseCommandLineArgs(sys.argv)
   