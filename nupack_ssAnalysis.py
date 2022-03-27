# from opendna import opendna
import interfaces
from utils import *
from analysisTools import *

def getSecondaryStructure(aptamerSeq):
    """
    get the secondary structure(s) for a given aptamer's sequence
    using seqfold here - identical features are available in nupack, though results are sometimes different
    :param aptamerSeq: aptamer's sequence
    :return:
    """
    params={}
    params['temperature'] = 298.15       # 25.0C Kevin: used to predict secondary structure and for MD thermostat
    params['ionicStrength'] = 0.138      # Molar: sodium concentration - used to predict secondary structure and add ions to simulation box, must be 1100 M > [Na] > 50 for nupack to run    
    params['[Mg]'] = 0.0005              # Molar: magnesium concentration: 0.2 M > [Mg] > 0 - ONLY applies to NuPack fold - Does NOT add Mg to OpenMM simulation.      
    params['N 2D structures'] = 1

    # printRecord("Getting Secondary Structure(s)")
    # if self.params['secondary structure engine'] == 'seqfold':        
    #     ssString, pairList = getSeqfoldStructure(aptamerSeq, self.params['temperature'])  # seqfold guess
    #     self.ssAnalysis = {'2d string': ssString, 'pair list': pairList, 'displayed 2d string': ssString}
    #     # print('seqfold ssAnalysis: ',self.ssAnalysis)
    #     return pairList
    # elif self.params['secondary structure engine'] == 'NUPACK':

    nup = interfaces.nupack(aptamerSeq, params['temperature'], params['ionicStrength'], params['[Mg]'], stacking='nostacking')  # initialize nupack
    ssAnalysis = nup.run()  # run nupack analysis of possible 2D structures.

    distances = getSecondaryStructureDistance(ssAnalysis['config'])            
    if len(distances) > 1:
        topReps, topProbs, topDists = do2DAgglomerativeClustering(ssAnalysis['config'], ssAnalysis['state prob'], distances)
        topPairLists = []                
        nLists = min(len(topReps), params['N 2D structures'])
        ssAnalysis['displayed 2d string'] = []  # to output the dot-parenthesis notation of selected 2d structure
        for i in range(nLists):
            if (i == 0) or (np.amin(topDists[i, :i]) > 0.05): # only add a config to the list if it's at least 5% different from a previously accepted config
                ssAnalysis['displayed 2d string'].append(configToString(topReps[i].astype(int)))
                topPairLists.append(ssToList(configToString(topReps[i].astype(int))))
        # return topPairLists
        ssAnalysis['topPairLists'] = topPairLists
        return ssAnalysis 

    else:  # if, for some reason, there's only a single plausible structure (unlikely)
        ssAnalysis['displayed 2d string'] = ssAnalysis['2d string']    
        return ssAnalysis


seq = 'ACCTGGGGGAGTATTGCGGAGGAAGGT' # 27nt adenosine-binding aptamer
# seq = 'CCTGGGGGAGTATTGCGGAGGAAGG'
nup_ss = getSecondaryStructure(seq)
pairLists = nup_ss['pair list']

print('how many possible 2D structures: ', len(pairLists))
print('2D structure is: {}'.format(nup_ss['displayed 2d string'][0]))
print(nup_ss['topPairLists'])

