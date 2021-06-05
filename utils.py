'''
utilities
'''
import os
import numpy as np
import MDAnalysis as mda
import glob
import argparse
import matplotlib.pyplot as plt
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB
import MDAnalysis.analysis.nuclinfo as nuclinfo
from MDAnalysis.analysis import distances
import scipy.ndimage as ndimage
import scipy.spatial as spatial
from sklearn.decomposition import PCA
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import calc_bonds


def buildPeptide(sequence):
    '''
    construct a peptide sequence pdb file
    :param sequence:
    :return:
    '''
    #resdict = {"ALA": "A","CYS": "C","ASP": "D","GLU": "E","PHE": "F","GLY": "G","HIS": "H","ILE": "I","LYS": "K","LEU": "L","MET": "M","ASN": "N","PRO": "P","GLN": "Q","ARG": "R","SER": "S","THR": "T","VAL": "V","TRP": "W","TYR": "Y"}
    #PDBdir = "PDBs"
    structure = PeptideBuilder.initialize_res(sequence[0])
    for i in range(1,len(sequence)):
        geo = Geometry.geometry(sequence[i])
        PeptideBuilder.add_residue(structure, geo)

    PeptideBuilder.add_terminal_OXT(structure)

    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save('peptide.pdb')


def combinePDB(file1, file2):
    '''
    combine 2 pdb files into one
    ### actually for this one I think we may defer entirely to the docker
    :param file1:
    :param file2:
    :return:
    '''
    filenames = [file1,file2]
    for file in filenames:
        removeLine(file,'CRYST1')
        removeLine(file,'TITLE')
        removeLine(file,'END')
        if file == 'repStructure.pdb':
            appendLine('repStructure.pdb','TER')

    with open('combined.pdb','w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


def get_input():
    '''
    get the command line in put for the run-num. defaulting to a new run (0)
    :return:
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_num', type=int, default=0)
    cmd_line_input = parser.parse_args()
    run = cmd_line_input.run_num

    return run


def findLine(file, string):
    # return the line number of a given string, indexing from 1
    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.split('\n')
    try:
        lineInd = 1
        for line in text:
            if line == string:
                lineNum = lineInd  # index from 1
            lineInd += 1

        return lineNum
    except:
        raise ValueError("String not found in file!")


def replaceText(file, old_string, new_string):
    # search and replace text in a file, then save the new version
    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.replace(old_string, new_string)
    f = open(file, 'w')
    f.write(text)
    f.close()


def removetwo(target):
    if os.path.exists(target + "_2"):
        os.remove(target)  # delete original
        os.rename(target + "_2", target)  # append updated
    else:
        print("No _2 structure was made!")
        raise ValueError


def copyLine(file, line_number):
    # copy a line of text from a file and return it
    f = open(file, 'r')
    text = f.read()
    f.close()
    return text.split('\n')[line_number - 1]  # copy a line from the file, indexing from 1


def removeLine(file, string):
    '''
    remove every line containing given string from a file
    :param fine:
    :return:
    '''
    f = open(file, 'r')
    text = f.read()
    f.close()
    lines = text.split('\n')
    try:
        for i in range(len(lines)):
            if string in lines[i]:
                lines.pop(i)
    except:
        pass
    text = "\n".join(lines)
    f = open(file, 'w')
    f.write(text)
    f.close()


def addLine(file, string, line_number):
    # add a new line of text to a file
    f = open(file, 'r')
    text = f.read()
    f.close()
    lines = text.split('\n')
    lines.insert(line_number - 1, string)  # replace line ### with the string, indexing from 1
    text = "\n".join(lines)
    f = open(file, 'w')
    f.write(text)
    f.close()


def appendLine(file, string):
    # add a new line of end of a file
    f = open(file, 'r')
    text = f.read()
    f.close()
    lines = text.split('\n')
    lines.insert(len(lines) - 1, string)  # replace line ### with the string, indexing from 1
    text = "\n".join(lines)
    f = open(file, 'w')
    f.write(text)
    f.close()


def writeCheckpoint(text):
    '''
    write some output to the checkpoint file
    :return:
    '''
    f = open('checkpoint.txt', 'a')
    f.write('\n' + text)
    f.close()


def cleanTrajectory(structure,trajectory):
    '''
    remove water, salt from trajectory input
    :param structure: pdb input for initial structure template
    :param trajectory: dcd trajectory file
    :return:
    '''
    u = mda.Universe(structure,trajectory) # load up trajectory
    goodStuff = u.segments[:-2].atoms # cut out salts and solvent
    # write template
    goodStuff.write("clean" + structure)
    # and write trajectory
    with mda.Writer("clean" + trajectory, goodStuff.n_atoms) as W:
        for ts in u.trajectory:
            W.write(goodStuff)


def extractFrame(structure,trajectory,frame):
    '''
    saves a given trajectory frame as a separate pdb file
    :param structure: pdb input for initial structure template
    :param trajectory: dcd trajectory file
    :param frame: frame to be extracted
    :return:
    '''
    u = mda.Universe(structure,trajectory) # load up trajectory
    u.trajectory[frame] # this indexes the trajectory up to the desired frame
    atoms = u.atoms
    atoms.write("repStructure.pdb")


def dnaBasePairDist(u,sequence):
    pairDists = np.zeros((len(u.trajectory),len(sequence),len(sequence))) # matrix of watson-crick distances
    tt = 0
    for ts in u.trajectory:
        for i in range(len(sequence)):
            for j in range(len(sequence)):
                if i > j:
                    pairDists[tt,i,j] = nuclinfo.wc_pair(u,i+1,j+1,seg1='A',seg2='A') # also do WC base-pair distances (can set to only follow secondary structure prediction)
                    pairDists[tt,j,i] = pairDists[tt,i,j]

        tt += 1

    return pairDists


def dnaBaseCOGDist(u,sequence):
    baseDists = np.zeros((len(u.trajectory),len(sequence),len(sequence))) # matrix of watson-crick distances
    tt = 0
    for ts in u.trajectory:
        dnaResidues = u.segments[0]
        posMat = np.zeros((len(sequence), 3))

        for i in range(len(sequence)):
            posMat[i] = dnaResidues.residues[i].atoms.center_of_geometry()

        baseDists[tt, :, :] = distances.distance_array(posMat, posMat, box=u.dimensions)
        tt += 1

    return baseDists

def dnaBaseDihedrals(u,sequence):
    angles = np.zeros((len(u.trajectory),len(sequence),7))
    tt = 0
    for ts in u.trajectory:
        for j in range(2,len(sequence)-1):
            angles[tt,j,:] = nuclinfo.tors(u,seg="A",i=j)

        tt += 1

    return angles


def getPairs(wcDist):
    trajTime = len(wcDist)
    seqLen = wcDist.shape[-1]
    pairedBases = np.zeros((trajTime,seqLen))
    for tt in range(trajTime):
        pairMat = wcDist[tt] + np.eye(seqLen) * 20
        for i in range(wcDist.shape[-1]):
            nearestNeighbour = np.argmin(pairMat[i,:])
            if wcDist[tt,i,nearestNeighbour] < 3.3: # if we're within a hydrogen bond length, count it as a pair
                pairedBases[tt,i] = nearestNeighbour + 1 # indexing from 1

    return pairedBases

def trajectoryPCA(n_components,trajectory):
    '''
    do PCA on some trajectory with pre-selected features
    :param n_components:
    :param trajectory:
    :return: the principal components, their relative contributions, the original trajectory in PC basis
    '''
    pca1 = PCA(n_components=n_components)
    model = pca1.fit(trajectory.reshape(len(trajectory),int(trajectory.shape[-1]*trajectory.shape[-2])))
    components = model.components_
    eigenvalues = pca1.explained_variance_ratio_
    reducedTrajectory = model.transform(trajectory.reshape(len(trajectory),int(trajectory.shape[-1]*trajectory.shape[-2])))

    return components, eigenvalues, reducedTrajectory

def do_kdtree(trajectory,coordinates):
    '''
    returns the index in trajectory closest to 'coordinates'
    :param trajectory:
    :param coordinates:
    :return:
    '''
    mytree = spatial.cKDTree(trajectory)
    dist, index = mytree.query(coordinates)
    return index


def getSecondaryStructureDistance(configs):
    configDistance = np.zeros((len(configs),len(configs)))
    for i in range(len(configs)):
        for j in range(len(configs)):
            configDistance[i,j] = np.sum(configs[i] != configs[j]) / len(configs[0]) # normalized binary distance
    return configDistance


def analyzeSecondaryStructure(pairingTrajectory):
    '''
    given a trajectory of base pairing interactions
    identify the equilibrium structure
    and maybe metastable states
    :param trajectory:
    :return:
    '''
    configs = []
    counter = []
    for tt in range(len(pairingTrajectory)):  # find the equilibrium secondary structure
        if tt == 0:
            configs.append(pairingTrajectory[tt])
            counter.append(0)
        else:
            found = 0
            for i in range(len(configs)):
                if all(pairingTrajectory[tt] == configs[i]):
                    counter[i] += 1  # add one to the population of the state
                    found = 1

            if found == 0:  # if we don't find it, enumerate a new possible state
                configs.append(pairingTrajectory[tt])
                counter.append(0)

    # compute a distance metric between all the seen configs
    configDistance = getSecondaryStructureDistance(configs)

    # identify the most common structures
    # in the future: combine structures with small distances
    # and highlight strongly different structures which nevertheless appear

    bestStructure = configs[np.argmax(counter)]
    return bestStructure

def isolateRepresentativeStructure(trajectory):
    '''
    use PCA to identify collective variables
    identify the most probable structure and save it
    :param trajectory:
    :return:
    '''
    components, eigenvalues, reducedTrajectory = trajectoryPCA(10, trajectory)  # do it with a bunch of components, then re-do it with only the necessary number
    # we want there to be a gap in this spectrum, or at least, to neglect only the small contributions
    n_components = np.sum(np.abs(np.diff(eigenvalues)) > np.mean(np.abs(np.diff(eigenvalues))))  # this is kindof a hack - we can code better logic later
    components, eigenvalues, reducedTrajectory = trajectoryPCA(n_components, trajectory)

    #find minima in the space of collective variables
    try:
        nbins = int(1e5 ** (1 / n_components)) # this is sometimes too large, hence, try-except
        probs, bins = np.histogramdd(reducedTrajectory, bins=nbins)  # multidimensional probability histogram
    except:
        nbins = int(1e3 ** (1 / n_components)) # much smaller - should always work but may impact prediction quality
        probs, bins = np.histogramdd(reducedTrajectory, bins=nbins)  # multidimensional probability histogram
    # smooth it out and take the diverse maxima - compare relative probability
    sigma = 5 # should automate width-setting for this
    smoothProbs = ndimage.gaussian_filter(probs, sigma)
    bestTransformedStructureIndex = np.unravel_index(smoothProbs.argmax(), smoothProbs.shape)

    # Pull representative 3D structures from the original trajectory (rather than back-transforming and generating a new structure from only dihedrals - HARD)
    transformedDimensions = []
    for i in range(len(bins)):
        transformedDimensions.append(bins[i][1:] - np.diff(bins[i][0:2]))

    transformedDimensions = np.asarray(transformedDimensions)
    transformedCoordinates = []
    for i in range(len(transformedDimensions)):
        transformedCoordinates.append(transformedDimensions[i][bestTransformedStructureIndex[i]])  # best structure coordinates in transformed basis

    transformedCoordinates = np.asarray(transformedCoordinates)

    # find structures in the reduced trajectory with these coordinates
    representativeIndex = do_kdtree(reducedTrajectory, transformedCoordinates)

    return representativeIndex


def bindingAnalysis(u,peptide,sequence):
    '''
    analyze the binding of analyte to aptamer by computing relative distances
    :param u:
    :return:
    '''
    # identify base-analyte distances
    if u.segments.n_segments == 2:  # if we have an analyte
        pepNucDists = np.zeros((len(u.trajectory), len(peptide), len(sequence)))  # distances between peptides and nucleotiedes
        tt = 0
        for ts in u.trajectory:
            analyteResidues = u.segments[1]
            baseResidues = u.segments[0]
            posMat1 = np.zeros((len(peptide), 3))
            posMat2 = np.zeros((len(sequence), 3))
            for i in range(len(peptide)):
                posMat1[i] = analyteResidues.residues[i].atoms.center_of_geometry()
            for i in range(len(sequence)):
                posMat2[i] = baseResidues.residues[i].atoms.center_of_geometry()


            pepNucDists[tt, :, :] = distances.distance_array(posMat1, posMat2, box=u.dimensions)
            tt += 1

        # identify contacts
        contactCutoff = 10 # distance within which we say there is a 'contact'
        contacts = []
        for i in range(len(pepNucDists)):
            contacts.append(np.argwhere(pepNucDists[i] < contactCutoff))


    else:
        pepNucDists = [] # return empty
        contacts = []


    return [pepNucDists, contacts]


class atomDistances(AnalysisBase):
    '''
    calculate watson-crick bond lengths
    '''
    def __init__(self, atomgroups, **kwargs):
        '''
        :param atomgroups: a list of atomgroups for which the interatom bond lengths are calculated
        '''
        super(atomDistances,self).__init__(atomgroups[0].universe.trajectory,**kwargs)
        self.atomgroups = atomgroups
        self.ag1 = atomgroups[0]
        self.ag2 = atomgroups[1]

    def _prepare(self):
        self.results = []

    def _single_frame(self):
        distance = calc_bonds(self.ag1.positions,self.ag2.positions,box=self.ag1.dimensions)
        self.results.append(distance)

    def _conclude(self):
        self.results = np.asarray(self.results)


def wcTrajAnalysis(u):
    '''
    use the atomDistances class to calculate the WC base pairing distances between all bases on a sequence
    :param u:
    :return:
    '''