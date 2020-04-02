#helpers.py

#==== Suffix array ======#

def getSuffixes(s):
    return sorted([s[i:] for i in range(len(s))])

def printSuffixArray(s):
     suffixArr = suffixArray(s)
     return ', '.join([str(elem) for elem in suffixArr]) 

'''
The suffix array is the list of starting 
positions of these sorted suffixes. 
'''
def suffixArray(s):
    suffixes = getSuffixes(s)
    result = [] 

    for i in suffixes:
        result.append(s.index(i))
    return (suffixes,result)

#==== BWT ======#

def rotations(t):
    tt = t * 2
    return[tt[i:i+len(t)] for i in range(0,len(t))]

def sortMatrix(t):
    return sorted(rotations(t))

def buildBWT(t):
    res = ''.join(map(lambda x:x[-1],sortMatrix(t))) 
    return res 

def generateFirstCol(text):
    suffixes = getSuffixes(text)
    res = [] 
    for suf in suffixes:
        res.append(suf[0])
    return res
    
def firstOccurTable(firstCol):
    """
    creates a lookup table for the first occurence
    """
    firstColumn = firstCol
    table = {} 
    character = ['A','C','G','T']
    for char in character:
        table[char] = firstColumn.index(char)
    return table


def occurenceTable(firstCol):
    #return column of indices for each character
    return {}

def generateLFMapping(t, firstCol, idx):
    transform = t
    match = t[idx]
    characters = firstCol
    numberOfAppearences  = 0 
    firstColAppearences = 0 

    for i in range(0,idx + 1):
        if transform[i] == match:
            numberOfAppearences += 1

    for i in range(len(characters)):
        if characters[i] == match:
            firstColAppearences +=1 
            if(numberOfAppearences == firstColAppearences):
                result = i 
                break
    
    with open("result.txt", "w") as f:
        f.write(str(result))
    return result

#-------BWT Index class -----#

class BWTIndex:

    def __init__(self):
        self.firstCol = ''
        self.bwt = ''
        self.reference = ''
        self.occ = ''

    def processFMIndex(self):
        #populate needed data from ref_index
        firstCol = ''
        bwt = ''
        reference = ''
        occ = ''
        with open('ref_index') as r:
            ''
        return (firstCol,bwt,reference,occ)

    def get_interval(self, s):
        interval = 0 
        match_len = 0
        #perform backward search , return interval and the length of the 
        #query that was matched.

    
        return (interval,match_len)

    def ref_positions(self):
        return []
        





