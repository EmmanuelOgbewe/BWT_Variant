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

if __name__ == "__main__":
    testSuffixArray()