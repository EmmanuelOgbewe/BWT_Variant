# Bwt_variant 

#### Basic read alignment tool, based on the seed-and-extend using the FM-index for efficient seed finding. ####

Make sure to have the following modules installed if not already:

`pip3 install numpy`

`pip3 install gzip`

`pip3 install pysam`


## To run the index command: 

`python3 fmmap.py fmmap index reference.fa ref_index`

**NOTE**
reference.fa and ref_index must be valid paths.If the following files
are not within the projects directory please make sure you specify the 
full path. 

## To run the align command:

`python3 fmmap.py fmmap align ref_index reads.fa alignments.sam`

**NOTE**
ref_index, reads.fa and alignments.sam must be valid paths.If the following files
are not within the projects directory please make sure you specify the 
full path. 
