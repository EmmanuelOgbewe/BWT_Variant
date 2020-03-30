#fmmap.py
#This is our main code for our basic read alignment. This tool takes
#in two commands index and align. 
# @author Emmanuel Ogbewe

import sys
import os 
import pysam

#------Handle input and output------#
def parseCommandLineArgs(args=[]):
    if 'index' in args:
        # run the index function 
        print("running index function")
        input = readFastaFormat(args[3])
    elif 'align' in args:
        #run the align function
        print("running align function")
    else:
        print("Error parsing commanding arguments, invalid command.")

#------Helper Functions-----#
def readFastaFormat(file):
    result = ''
    with open(file) as f:
        for line in f.readlines():
            if ">" in line:
                pass
            else:
                result = result + line
    return result


#------Start command Functions-------#
"""
Takes two arguments is the path to a FASTA format file and the 
second is the location of an output file where the index will be written.
When invoked in this mode, fmmap will read the input file and will compute
the FM-index of the corresponding string.
retuns index containing the following:
"""
def index(input,output):
    return 


def align():
    return 
#------ End command functions-------#


if __name__ == "__main__":
    #handle command line arguments
   parseCommandLineArgs(sys.argv)