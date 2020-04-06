#helpers.py
import json
import numpy as np


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
    characters = ['A','C','G','T']
    for char in characters:
        table[char] = firstColumn.index(char)
    return table


def occurenceTable(characterList,lastCol):
    #return column of indices for each character
    characters =  ['A','C','G','T'] if characterList == None else characterList
    result = { }
    row = 0 
    for c1 in range(0,len(lastCol)): 
        for c2 in characters:
            if lastCol[c1] == c2:
                if c2 in result.keys():
                    result[c2].append(result[c2][c1 - 1] + 1)
                else:
                    result[c2] = [1]
            elif c1 > 0:
                if c2 in result.keys():
                    result[c2].append(result[c2][c1 - 1])
            else:
                result[c2] = [0]
    return result 


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

def switchCurrentProcessing(s):
    if 'name:' in s:
        return 'name'
    elif 'reference_string:' in s:
        return 'reference_string'
    elif 'length:' in s:
        return 'length'
    elif 'first column:' in s:
        return 'first column'
    elif 'bwt:' in s:
        return 'bwt'
    elif 'occurence table:' in s:
        return 'occurence table'
    elif 'suffix array:' in s:
        return 'suffix array'
    return None

def processList(s):
    s = s.strip()
    s = s.strip("[").strip("]")
    out = s.replace("'","")
    res = out.split(", ")
    return res

#----Align helpers-----#

#returns cost during comparison
def score(c1,c2):
  if c1 == c2:
    #orig 1,-1, sgap: -1
    return 0
  else:
    return -2

def editDistance(row,col, matrix, directionMatrix, seq,ref):
  sgap = -2

  val_one = score(seq[row - 1], ref[col - 1]) + matrix[row - 1][col - 1] #diagonal
  val_two = matrix[row][col - 1] + sgap #down
  val_three = matrix[row - 1][col] + sgap # left

  an_array = np.array([val_one,val_two,val_three])
  max_index= np.argmax(an_array, axis=0)

  if max_index == 0:
    directionMatrix[row][col] = "D"
  elif max_index == 1:
    directionMatrix[row][col] = "u"
  else:
    directionMatrix[row][col] = "l"

  return an_array[max_index]

def generateAlignmentMatrixes(ref_string,seq):
   
    w,h = len(ref_string) + 1, len(seq) + 1

    #initliaze both matrix 
    valuesMatrix = np.zeros((h,w),dtype=int)
    directionsMatrix = np.zeros((h,w),dtype=str)

    #initliaze first row in matrix
    for i in range(1,len(seq) + 1):
        valuesMatrix[i][0] = i * -1
        directionsMatrix[i][0] = 'u'

    for i in range(1,len(ref_string) + 1):
        #initialize to know direction
        directionsMatrix[0][i] = 'l'
    
    return (valuesMatrix, directionsMatrix)


#backtrack 
def createCigar(directionMatrix,largestIndex):
  cigar = ''
  row, col = len(directionMatrix) - 1, largestIndex

  while(not row == 0):
    #program freezes runs into loop here

    if directionMatrix[row][col] == 'D':
      row = row - 1
      col = col - 1
      cigar = 'M' + cigar
    elif directionMatrix[row][col] == 'u':
      row = row - 1
      cigar = 'I' + cigar
    elif directionMatrix[row][col] == 'l':
      col = col - 1
      cigar = 'D' + cigar

  originalPosInRef = col
#   print(cigar)
#   print("length of cigar" + str(len(cigar)))
  #compress cigar 
  finalCigar = ''
  currentChar = cigar[0]
  count = 0 
    
  for i in range(0,len(cigar)):
    if not cigar[i] == currentChar:
      finalCigar += str(count) + currentChar
      currentChar = cigar[i]
      count = 1
    else:
      count += 1

    if i == len(cigar) - 1:
      if cigar[i] == currentChar:
          finalCigar += str(count) + currentChar

  
#   print(finalCigar)
  return (finalCigar, originalPosInRef)


class Alignment:
    def __init__(self, score,cigar,posInRef):
        self.score = score
        self.cigar = cigar
        self.posInRef = posInRef


class BWTIndex:

    def __init__(self):
        (name,ref,length,firstCol,bwt, occ, suffixArray) = self.processFMIndex()

        self.name = name
        self.length = length
        self.firstCol = firstCol
        self.bwt = bwt
        self.reference = ref
        self.occ = occ
        self.suffixArray = suffixArray
        self.first = firstOccurTable(firstCol)


    def processFMIndex(self):
        #populate needed data from ref_index
        name = ''
        length = ''
        firstCol = ''
        bwt = ''
        reference_string = ''
        occ = ''
        suffixArray = ''
    
        with open('ref_index') as r:
            currProcessing = ''
            for line in r.readlines():
                if not switchCurrentProcessing(line) == None:
                    currProcessing = switchCurrentProcessing(line) 
                    continue 
                
                if not line == '':
                    if currProcessing == 'name':
                        name += line
                    elif currProcessing == 'reference_string':
                        reference_string += line
                    elif currProcessing == 'length':
                        length += line
                    elif currProcessing == 'first column':
                        firstCol += line
                    elif currProcessing == 'occurence table':
                        occ += line
                    elif currProcessing == 'bwt':
                        bwt += line
                    elif currProcessing == 'suffix array':
                        suffixArray += line
        
        name = name.strip()
        reference_string = reference_string.strip()
        length = int(length.strip())
        firstCol = processList(firstCol)
        bwt = bwt.strip()
        json_acceptable_occ = occ.replace("'", "\"")
        occ = json.loads(json_acceptable_occ)
        suffixArray = json.loads(suffixArray)
        print("Finished parsing ref_index")
        return (name,reference_string, length,firstCol, bwt, occ, suffixArray)

    def get_interval(self, s):
        match_len = 0
        #perform backward search , return interval and the length of the 
        #query that was matched.
        occ = self.occ
        first = self.first
        if s[-1] not in first:
            return 0 # character dosen't occur in T 

        left = 1
        right = len(self.reference)
    
        t = len(s) - 2
        while t >= 0 and right > left:
            character = s[t]
            left = first[character] + occ[character][left - 1]
            right = first[character] + occ[character][right - 1]
            t = t - 1
        if left > right:
            return None

        return (left,right)

    def ref_positions(self,interval,seed_end,match_len):
        (l,r) = interval
        result = [] 
        for i in range(l,r):
             result.append(self.suffixArray[i])

        return result

    def fitting_alignment(self,seq,ref_pos, gap):
        # first slice your new reference, the starting position is the ref_pos
        startPosition = ref_pos - gap
        endPosition = ref_pos + (len(seq) + gap)
        if startPosition < 0:
            startPosition = 0
        sl = self.reference[startPosition: endPosition]
        # create matrixes
        (valuesMatrix,dirMatrix) = generateAlignmentMatrixes(sl,seq) 
        print(sl)
        # edit distance
        for row in range(1,len(seq) + 1):
            for col in range(1,len(sl) + 1):
                valuesMatrix[row][col] = editDistance(row,col,valuesMatrix,dirMatrix,seq,sl)

        #get the largest
        largestArr = valuesMatrix[len(valuesMatrix) - 1]
        max_index = np.argmax(largestArr, axis=0)
        largestVal = largestArr[max_index]
        
        #cigar 
        (cigar,posInRef) = createCigar(dirMatrix, max_index)
        #the returned alignment object
        # TODO: consider PosINRef = ref_pos + posInRef
        alignment = Alignment(largestVal,cigar,ref_pos + posInRef)

        return alignment


    
        





