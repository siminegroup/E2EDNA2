from nupack import *
from utils import *
from seqfold import dg, fold


def getSeqfoldStructure(sequence,temperature):
    '''
    output the secondary structure for a given sequence at a given condition
    formats - ss string and pair list
    '''
    dg(sequence, temp=temperature)  # get energy of the structure
    # print(round(sum(s.e for s in structs), 2)) # predicted energy of the final structure

    structs = fold(sequence)  # identify structural features
    desc = ["."] * len(sequence)
    pairList = []
    for s in structs:
        pairList.append(s.ij[0])
        pairList[-1]  # list of bound pairs indexed from 1
        if len(s.ij) == 1:
            i, j = s.ij[0]
            desc[i] = "("
            desc[j] = ")"

    ssString = "".join(desc)
    pairList = np.asarray(pairList) + 1

    return ssString, pairList


def getNupackString(sequence,temperature,ionicStrength):
    R = 0.0019872 # ideal gas constant in kcal/mol/K
    gap = R * temperature
    A = Strand(sequence, name='A')
    comp = Complex([A], name='AA')
    set1 = ComplexSet(strands=[A], complexes=SetSpec(max_size=1, include=[comp]))
    model1 = Model(material='dna', celsius = temperature - 273, sodium= ionicStrength)
    results = complex_analysis(set1, model=model1, compute=['pfunc', 'mfe','subopt','pairs'], options={'energy_gap':gap})
    cout = results[comp]
    prob = np.exp(-cout.mfe[0].energy/R/temperature)/float(cout.pfunc) # probability - energy/partition function
    ssString = cout.mfe[0].structure
    # look at suboptimal structures
    subopts = cout.subopt
    subStrings, subProbs, subEns, subStack, subProbs = [[],[],[],[],[]]
    for subopt in subopts:
        #subStrings.append(subopt.structure)
        #subEns.append(subopt.energy)
        #subStack.append(subopt.stack_energy)
        subProbs.append(np.exp(-subopt.energy/R/temperature)/float(cout.pfunc))


    return str(ssString), prob, subopts, subProbs


def getNupackStructure(sequence,temperature,ionicStrength):
    '''
    compute secondary structure using NUPACK potentials
    inputs - sequence, temperature in K, ionic strength (limited range of working ionic strengths)
    '''
    A = Strand(sequence, name='A')
    comp = Complex([A], name='AA')
    set1 = ComplexSet(strands=[A], complexes=SetSpec(max_size=1, include=[comp]))
    model1 = Model(material='dna', celsius = temperature - 273, sodium= ionicStrength)
    results = complex_analysis(set1, model=model1, compute=['pfunc', 'pairs', 'subopt', 'mfe'])
    cout = results[comp]

    '''
    # print the secondary structure strings
    print('\nMFE proxy structures(s) for set')
    for i, s, in enumerate(cout.mfe):
        print('    %2d: %s (%.2f kcal/mol)' % (i, s.structure, s.energy))

    print('\nSuboptimal proxy structures for set')
    for i, s in enumerate(cout.subopt):
        print('    %2d: %s (%.2f kcal/mol)' % (i, s.structure, s.energy))
    '''

    # extract list of most probable base pairs
    pairMat = cout.pairs.to_array()
    pairList = []
    paired = []
    for i in range(len(pairMat)): # this whole thing doesn't seem to work properly to me
        bestMatch = np.argmax(pairMat[i,:])
        if True:#pairMat[i,bestMatch] > 0.5: # if we're relatively confident about this secondary structure feature
            if bestMatch != i: # if the best match is itself, it means it does not bind
                if (not bestMatch in paired) and (not i in paired):# check for duplicates
                    pairList.append([i + 1,bestMatch + 1])
                    paired.append(i)
                    paired.append(bestMatch)

    return cout.mfe[0].structure, np.asarray(pairList) # best structure and pair list


def ssToList(ssString):
    '''
    if for some reason we have a secondary structure string we need to convert to a pair list
    '''
    pairList = []
    paired = []
    for i in range(len(ssString)):
        if ssString[i] == '(': # if it's paired
            counter = 0
            for j in range(1,len(ssString[i:])): # look for the thing paired to it
                if ssString[i+j] == '(':
                    counter += 1
                if ssString[i+j] == ')':
                    if counter > 0:
                        counter -= 1
                    elif counter == 0:
                        # check for duplicates
                        if (not i in paired) and (not i+j in paired):  # check for duplicates
                            paired.append(i)
                            paired.append(i+j)
                            pairList.append([i + 1,i+j + 1]) #make pair list in 1-n basis
                            break

    return pairList



''' test # development
sequences = numbers2letters(np.random.randint(0,4,size=(1000,10))) # we can go even faster I think, by submitting them in batches
strings = []
pairList = []
probs = np.zeros(len(sequences))
nPairs = np.zeros_like(probs)
pFrac = np.zeros_like(probs)
nPins = np.zeros_like(probs)
subProbskT = np.zeros_like(probs) # sum of probabilities within kT
nSubProbskT = np.zeros_like(probs) # number of states within kT
pinSize, pinDist, subopts, subProbs = [[],[],[],[]]
t0 = time.time()
for i in range(len(sequences)):
    out = getNupackString(sequences[i],310,1.0)
    strings.append(out[0])
    probs[i] = out[1]
    subopts.append(out[2])
    subProbs.append(out[3])
    nPairs[i] = out[0].count('(')
    pFrac[i] = 2* nPairs[i] / len(sequences[i])
    pairList.append(ssToList(out[0]))
    subProbskT[i] = np.sum(subProbs[i])
    nSubProbskT[i] = len(subopts[i])

    indA = 0 # hairpin completion index
    for j in range(len(sequences[i])):
        if strings[i][j] == '(':
            indA += 1
        elif strings[i][j] == ')':
            indA -= 1
            if indA == 0: # if we come to the end of a distinct hairpin
                nPins[i] += 1

tf = time.time()
print(tf-t0)
'''

# possible functions on this space
# the secondary structure iteslf - pair list, or ssString
# number of pairs - probs
# pair fraction - pFrac
# structure probability - nPairs
# number of independent hairpins - nPins
# we could also do distance between hairpins
# probability density of cluster of similar subopt structures
