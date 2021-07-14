from nupack import *
from utils import nupackAnalysis

# design complexes

ssStrings = ['(((....)))','((((....))))','(((.....)))','((.((....)).))','((((......))))']

myModel = Model(material='dna', celsius=37,sodium=0.163)

sequences = []
analysisDict = []
for string in ssStrings:
    sequences.append(des(structure = string,model=myModel)[0].replace('U','T'))
    analysisDict.append(nupackAnalysis(sequences[-1], 310, 0.163)[0])

printRecord(sequences)

