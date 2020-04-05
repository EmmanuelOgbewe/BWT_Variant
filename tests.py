# tests.py
import helpers

# ---- Suffix array ---- #

def testSuffixArray():
    s = 'abaaba$'
    (suffixes,sa) = helpers.suffixArray(s)
    firstCol = helpers.generateFirstCol(s)
    
    occurenceTable = helpers.firstOccurTable(firstCol)
    print(str(occurenceTable))

    assert(firstCol == ['$', 'a', 'a', 'a', 'a', 'b', 'b'])


# ---- BWT ---- #
s = 'abaaba$'

def testOccurenceTable():
    bwt = helpers.buildBWT(s)
    occur = helpers.occurenceTable(['a','b'],bwt)
    print(str(occur))

def testIndexParser():
    bwt_index = helpers.BWTIndex()

def testOneBWTGetInteval():
    bwt_index = helpers.BWTIndex()
    bwt_index.reference = s
    bwt_index.first = {'a': 1, 'b' : 5 }
    bwt_index.occ = {'a' : [1,1,1,2,2,3,4], 'b': [0,1,2,2,2,2,2]}
    bwt_index.suffixArray = helpers.suffixArray(s)[1]
    res = bwt_index.get_interval('aba')
    ref_positions = bwt_index.ref_positions(res,1,1)
    print(str(res))
    # print(ref_positions)
    assert(res == (3,5))


def testTwoBWTGetInteval():
    s = 'mississippi$'
    bwt_index = helpers.BWTIndex()
    bwt_index.reference = s
    bwt_index.first = {'i': 1, 'm' : 5, 'p': 6, 's' : 8 }
    bwt_index.occ = {'i' : [1,1,1,1,1,1,1,2,2,2,3,4], 'm': [0,0,0,0,1,1,1,1,1,1,1,1], 'p': [0,1,1,1,1,1,2,2,2,2,2,2], 's' : [0,0,1,2,2,2,2,2,3,4,4,4]}
    bwt_index.suffixArray = helpers.suffixArray(s)[1]
    res = bwt_index.get_interval('iss')
    ref_positions = bwt_index.ref_positions(res,1,1)
    suffixes = helpers.suffixArray(s)[0]
    print(bwt_index.suffixArray)
    print(suffixes)

    for i in range(res[0],res[1]):
        print(suffixes[i])

    assert(res == (3,5))
    assert(ref_positions == [4,1])

if __name__ == "__main__":
    # testSuffixArray()
    # testOccurenceTable()
    # testIndexParser()
    # testOneBWTGetInteval()
    testTwoBWTGetInteval()
