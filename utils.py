"""
utilities
"""
import os
import numpy as np
import MDAnalysis as mda
import mdtraj as md
import glob
import argparse
import matplotlib.pyplot as plt
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB
from MDAnalysis.analysis import distances
import simtk.openmm.app as app
import scipy.ndimage as ndimage
import scipy.spatial as spatial
import sklearn.cluster as cluster
from sklearn.decomposition import PCA
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import calc_bonds
from MDAnalysis.analysis.dihedrals import Dihedral
import time
from nupack import *
from seqfold import dg, fold
from shutil import copyfile, copytree, move


def buildPeptide(peptide):
    """
    construct a peptide sequence pdb file
    :param peptide:
    :return:
    """
    # resdict = {"ALA": "A","CYS": "C","ASP": "D","GLU": "E","PHE": "F","GLY": "G","HIS": "H","ILE": "I","LYS": "K","LEU": "L","MET": "M","ASN": "N","PRO": "P","GLN": "Q","ARG": "R","SER": "S","THR": "T","VAL": "V","TRP": "W","TYR": "Y"}
    structure = PeptideBuilder.initialize_res(peptide[0])
    for i in range(1, len(peptide)):
        geo = Geometry.geometry(peptide[i])
        PeptideBuilder.add_residue(structure, geo)

    # PeptideBuilder.add_terminal_OXT(structure) # OpenMM will not run without this, but LightDock will not run with it. Solution, add terminal oxygen in prepPDB after docking

    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save('peptide.pdb')


class Timer:
    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start


def combinePDB(file1, file2):
    """
    combine 2 pdb files into one
    some special formatting for MDA outputs in particular
    :param file1:
    :param file2:
    :return:
    """
    filenames = [file1, file2]
    for file in filenames:  # remove title, periodic box, endpoints
        removeLine(file, 'CRYST1')
        removeLine(file, 'TITLE')
        removeLine(file, 'END')
        if file == 'repStructure.pdb':
            appendLine('repStructure.pdb', 'TER')

    with open('combined.pdb', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)


def get_input():
    """
    get the command line in put for the run num. defaulting to a new run (0)
    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_num', type=int, default=0)
    parser.add_argument('--sequence', type=str, default='XXX')
    parser.add_argument('--peptide', type=str, default='BBB')
    cmd_line_input = parser.parse_args()
    run = cmd_line_input.run_num
    sequence = cmd_line_input.sequence
    peptide = cmd_line_input.peptide

    return [run, sequence, peptide]


def findLine(file, string):
    # return the line number of a given string in a text file, indexing from 1
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


def readFinalLines(file, lines):
    # return the final N lines of a text file
    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.split('\n')
    finalLines = text[-(lines + 1):]

    return finalLines


def readInitialLines(file, lines):
    # return the first N lines of a text file
    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.split('\n')
    initialLines = text[:lines]

    return initialLines


def replaceText(file, old_string, new_string):
    # search and replace text in a text file, then save the new version
    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.replace(old_string, new_string)
    f = open(file, 'w')
    f.write(text)
    f.close()


def copyLine(file, line_number):
    # copy a line of text from a file and return it
    f = open(file, 'r')
    text = f.read()
    f.close()
    return text.split('\n')[line_number - 1]  # copy a line from the file, indexing from 1


def removeLine(file, string):
    """
    remove every line containing given string from a file
    :param file: the file
    :param string: the string to look for
    :return:
    """
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
    """
    write some output to the checkpoint file
    :return:
    """
    f = open('checkpoint.txt', 'a')
    f.write('\n' + text)
    f.close()


def cleanTrajectory(structure, trajectory):
    """
    remove water, salt from trajectory input
    :param structure: pdb input for initial structure template
    :param trajectory: dcd trajectory file
    :return:
    """
    u = mda.Universe(structure, trajectory)  # load up trajectory
    goodStuff = u.segments[:-2].atoms  # cut out salts and solvent
    goodStuff.write("clean" + structure)  # write topology
    with mda.Writer("clean" + trajectory, goodStuff.n_atoms) as W:  # and write trajectory
        for ts in u.trajectory:  # indexing over the trajectory
            W.write(goodStuff)


# noinspection PyStatementEffect
def extractFrame(structure, trajectory, frame, outFileName):
    """
    saves a given trajectory frame as a separate pdb file
    :param structure: pdb input for initial structure template
    :param trajectory: dcd trajectory file
    :param frame: frame to be extracted
    :param outFileName: the name of the output file
    :return:
    """
    u = mda.Universe(structure, trajectory)  # load up trajectory
    u.trajectory[frame]  # this indexes the trajectory up to the desired frame (weird syntax, I think)
    if u.segments.n_segments > 2:  # if there are more than 2 segments, then there must be solvent and salts (assuming nonzero salt concentration)
        atoms = u.segments[:-2].atoms  # omit solvent and salts
    else:
        atoms = u.atoms
    atoms.write(outFileName)


# noinspection PyTypeChecker
def dnaBaseCOGDist(u):
    """
    given an MDA universe containing ssDNA (segment 1)
    return the trajectory of all the inter-base center-of-geometry distances
    """
    nbases = u.segments[0].residues.n_residues
    baseDists = np.zeros((len(u.trajectory), nbases, nbases))  # matrix of CoG distances
    tt = 0
    for ts in u.trajectory:
        dnaResidues = u.segments[0]
        posMat = np.zeros((nbases, 3))

        for i in range(nbases):
            posMat[i] = dnaResidues.residues[i].atoms.center_of_geometry()

        baseDists[tt, :, :] = distances.distance_array(posMat, posMat, box=u.dimensions)  # fast calculation
        tt += 1

    return baseDists


def getPairs(wcTraj):
    """
    identify bases which are paired according to their WC hydrogen bond lengths
    return the shortest hydrogen bond pair for every base over time
    note - it's possible to have multiple bases paired to one in this algo -
    TO-DO
    post-processing cleanup step where we see what the next closest base is for multi-paired setups
    """
    trajTime = len(wcTraj)
    seqLen = wcTraj.shape[-1]
    pairedBases = np.zeros((trajTime, seqLen)).astype(int)
    problems = []
    for tt in range(trajTime):
        pairMat = wcTraj[tt] + np.eye(seqLen) * 20  # add 20 on the diagonal so it's never counted as the 'nearest neighbour' to itself
        randomIndices = np.arange(wcTraj.shape[-1])
        np.random.shuffle(randomIndices)
        for i in range(wcTraj.shape[-1]): # search for the closest base
            ind = randomIndices[i] # randomize order so that any errors are not correlated
            nearestNeighbour = np.argmin(pairMat[ind, :])
            if nearestNeighbour > ind: # only do the rest if the nearest neighbour is > i - ensures reciprocal algorithm will work
                if pairedBases[tt,nearestNeighbour] == 0: # also always accept the first base and no others - may or may not be optimal but it is clean
                    if wcTraj[tt, ind, nearestNeighbour] < 3.3:  # if we're within a hydrogen bond length, count it as a pair
                        if np.abs(ind-nearestNeighbour) > 2: # also, we cannot pair with a bases directly adjacent
                            pairedBases[tt, ind] = int(nearestNeighbour + 1)  # indexing from 1
                            pairedBases[tt, nearestNeighbour] = int(ind+1) # reciprocally pair this base to base i - indexing from 1
                            # this approach guarantees that sequences are always built from lists of pairs, with no multiple or dangling bonds
                            # possible that we may realize suboptimal structures, but overall this is safe and simple
                            # in principle - we could consider between different 3D configurations, e.g. which pairs are closer, but this I think should work for now
                            # could also minimize error by shrinking the H-bond cutoff range, but we then may reduce pairing too much

        # we then have to clean it, because this does not enforce that bases are mutually paired
        pairingErrors = findPairingErrors(pairedBases[tt], tt)
        if pairingErrors != []:
            problems.append(pairingErrors)  # collect problems timepoint-by-timepoint

    if problems != []:
        print('Encountered error in pairing algorithm!') # should not be possible, currently

    '''
    # now we need to fix the snafu's - or actually nevermind
    for i in len(problems):
        # a base paired to a completely unpaired base is impossible, since the above code is completely symmetric
        # so the only possibility is multiple bases both wanting a third
        # so, we will find bases which are 'wanted'
        # let the closest base keep it
        # second closest base must then go with its next best choice
        config = pairedBases[problems[i][0][0]]
        uniques, counts = np.unique(config, return_counts=True)
        multis = []
        for j in range(len(uniques) - 1): # for nonzero entries, record where we have multiple instances of a given base
            if counts[j+1] > 1:
                multis.append(uniques[j + 1])

        for j in range(len(multis)): # for each multiple-paired base, fix it
            baseInd = multis[j] - 1 # indexing
    '''

    return pairedBases


def pairListToConfig(pairList, sequence):
    '''
    convert a secondary structure for a list of pairs to a 'config'
    which is where every base is listed with its pair-mate (or 0 for unpaired)
    :param pairList:
    :return:
    '''
    config = np.zeros(len(sequence)).astype(int)
    for i in range(len(pairList)):
        config[pairList[i][0] - 1] = pairList[i][1] # indexing
        config[pairList[i][1] - 1] = pairList[i][0]

    return config


# noinspection PyUnresolvedReferences
def trajectoryPCA(n_components, trajectory, transform):
    """
    do PCA on some trajectory array
    :param n_components:
    :param trajectory:
    :return: the principal components, their relative contributions, the original trajectory in PC basis
    """

    if trajectory.ndim == 3:  # convert to 2D, if necessary
        trajectory = trajectory.reshape(len(trajectory), int(trajectory.shape[-1] * trajectory.shape[-2]))

    if n_components > trajectory.shape[-1]:  # if the trajectory is small-dimensional, reduce the number of principal components
        n_components = trajectory.shape[-1]
        print('Warning, attempting PCA on a low-dimensional trajectory!')

    pca1 = PCA(n_components=n_components)
    model = pca1.fit(trajectory)
    components = model.components_
    eigenvalues = pca1.explained_variance_ratio_
    if transform == True:
        reducedTrajectory = model.transform(trajectory)
    else:
        reducedTrajectory = 0

    return components, eigenvalues, reducedTrajectory, pca1


def do_kdtree(trajectory, coordinates):
    """
    returns the index in trajectory closest to 'coordinates'
    :param trajectory:
    :param coordinates:
    :return:
    """
    mytree = spatial.cKDTree(trajectory)
    dist, index = mytree.query(coordinates)
    return index


def getSecondaryStructureDistance(configs):
    """
    normalized binary base pairing distance between sequences of equal length
    :param configs:
    :return:
    """
    configDistance = np.zeros((len(configs), len(configs)))
    nBases = len(configs[0])
    for i in range(len(configs)):
        configDistance[i, :] = np.sum(configs[i] != configs, axis=1) / nBases  # normalized binary distance
    return configDistance


def analyzeSecondaryStructure(pairTraj):
    """
    given a trajectory of base pairing interactions
    identify the equilibrium structure
    and maybe metastable states
    :param trajectory:
    :return:
    """


    # might be useful for clustering
    configs = []
    counter = []
    for tt in range(len(pairTraj)):  # find the equilibrium secondary structure
        if tt == 0:
            configs.append(pairTraj[tt])
            counter.append(0)
        else:
            found = 0
            for i in range(len(configs)):
                if all(pairTraj[tt] == configs[i]):
                    counter[i] += 1  # add one to the population of the state
                    found = 1

            if found == 0:  # if we don't find it, enumerate a new possible state
                configs.append(pairTraj[tt])
                counter.append(1)

    # compute a distance metric between all the seen configs
    bestSingleStructure = configs[np.argmax(counter)]


    if counter[np.argmax(counter)] / np.sum(counter) > 0.5: # if one structure predominates, go with that one
        representativeStructure = bestSingleStructure
    else: # otherwise, clustering, or PCA
        configDistance = getSecondaryStructureDistance(configs)
        agglomerate = cluster.AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='average', compute_full_tree=True, distance_threshold=0.3).fit(configDistance)
        labels = agglomerate.labels_
        nClusters = agglomerate.n_clusters_

        clusters = []
        for i in range(nClusters):
            inds = np.where(labels == i)[0].astype(int)
            clusters.append([configs[j] for j in inds])

        #representativeIndex, reducedTrajectory, eigenvalues = isolateRepresentativeStructure(pairTraj)  # do PCA and multidimensional binning to find free energy minimum
        #representativeStructure = pairTraj[representativeIndex]
        # alternatively we can do  clustering
        # and highlight different structures which nevertheless appear

    return representativeStructure


def doMultiDProbabilityMap(trajectory, range=None, nBins=None):
    """
    multidimensional histogram of system trajectories
    adjust bin size automatically
    """
    # find minima in the space of collective variables
    assert trajectory.ndim == 2

    if nBins is None:
        converged = False
        binList = np.logspace(6, 1, 6)
        ind = 0
        while converged is False:  # if there are too many bins this crashes
            nbins = int(binList[ind] ** (1 / trajectory.shape[-1]))  # this is sometimes too large, hence, try-except
            ind += 1
            probs, bins = np.histogramdd(trajectory, bins=nbins, range=range, density=True)  # multidimensional probability histogram
            converged = True
    else:
        nbins = nBins
        probs, bins = np.histogramdd(trajectory, bins=nBins, range=range, density=True)  # multidimensional probability histogram

    # smooth it out and take the diverse maxima - compare relative probability
    sigma = nbins / 20  # automate sigma to dimension size
    smoothProbs = ndimage.gaussian_filter(probs, sigma)

    return probs, smoothProbs, bins, nbins


def doTrajectoryDimensionalityReduction(trajectory):
    """
    automatically generate dimension-reduced trajectory using PCA
    :param trajectory:
    :return:
    """
    converged = False
    nComponents = 10
    while converged == False:  # add components until some threshold, then take components greater than the average
        components, eigenvalues, reducedTrajectory, pcaModel = trajectoryPCA(nComponents, trajectory, False)  # do it with a bunch of components, then re-do it with only the necessary number
        totVariance = np.sum(eigenvalues)
        if totVariance > 0.85:
            converged = True
        else:
            nComponents += 1

    # we want there to be a gap in this spectrum, or at least, to neglect only the small contributions
    n_components = min(5, np.sum(eigenvalues > np.average(eigenvalues)))  # with the threshold above, this removes some of the variation from having too many components
    components, eigenvalues, reducedTrajectory, pcaModel = trajectoryPCA(n_components, trajectory, True)

    return n_components, reducedTrajectory, pcaModel


def isolateRepresentativeStructure(trajectory):
    """
    use PCA to identify collective variables
    identify the most probable structure and save it
    :param trajectory:
    :return:
    """
    n_components, reducedTrajectory, pcaModel = doTrajectoryDimensionalityReduction(trajectory)
    probs, smoothProbs, bins, nbins = doMultiDProbabilityMap(reducedTrajectory)

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

    return representativeIndex, reducedTrajectory, pcaModel.explained_variance_ratio_


def bindingAnalysis(u, peptide, sequence):
    """
    analyze the binding of analyte to aptamer by computing relative distances
    :param u:
    :return:
    """
    assert u.segments.n_segments == 2
    # identify base-analyte distances
    pepNucDists = peptideAptamerDistances(u, peptide, sequence)
    contacts, nContacts = getPeptideContacts(pepNucDists)
    firstContact = np.nonzero(nContacts[:, 0])[0][0]  # first time when the peptide and aptamer were in close-range contact
    closeContactRatio = np.average(nContacts[firstContact:, 0] > 0)  # amount of time peptide spends in close contact with aptamer
    contactScore = np.average(nContacts[firstContact:, :] / len(peptide))  # per-peptide average contact score, linear average over 8-12 angstrom
    conformationChange = getConformationChange()

    # build directory of outputs
    outDict = {
        'peptide-nucleotide distance': pepNucDists,
        'contact list': contacts,
        '# contacts': nContacts,
        'close contact ratio': closeContactRatio,
        'contact score': contactScore,
        'conformation change': conformationChange
    }

    return outDict


def peptideAptamerDistances(u, peptide, sequence):
    """
    given an MDA universe
    return distances between each peptide and each base
    at every timepoint in the trajectory
    :return: pepNucDists - peptide-nucleicacid distances
    """

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

    return pepNucDists


def getPeptideContacts(pepNucDists):
    """
    get the contacts between peptide and dna aptamer
    return contacts
    return ncontacts (sum of contacts at various ranges
    """
    contactCutoffs = [6, 8, 10, 12]  # Angstroms - distance within which we say there is a 'contact'
    # identify contacts
    contacts = []
    nContacts = np.zeros((len(pepNucDists), len(contactCutoffs)))
    for i in range(len(pepNucDists)):
        contacts.append(np.argwhere(pepNucDists[i] < contactCutoffs[0]))  # loosest cutoff
        for j in range(len(contactCutoffs)):
            nContacts[i, j] = len(np.argwhere(pepNucDists[i] < contactCutoffs[j]))

    return contacts, np.asarray(nContacts)


def getConformationChange():
    """
    compare the pre-complexation (free aptamer) conformation with post-complexation
    NOTE intimately depends on naming conventions for trajectory files!
    """
    try:
        # function to analyze analyte impact on aptamer conformation
        freepdb = 'cleanaptamer.pdb'
        freedcd = 'cleanaptamer.dcd'
        bindpdb = 'cleancomplex.pdb'
        binddcd = 'cleancomplex.dcd'

        freeu = mda.Universe(freepdb, freedcd)
        bindu = mda.Universe(bindpdb, binddcd)

        freeAngles = nucleicDihedrals(freeu)
        bindAngles = nucleicDihedrals(bindu)

        n_components, freeReducedTrajectory, pcaModel = doTrajectoryDimensionalityReduction(freeAngles) # get free aptamer PCA
        bindReducedTrajectory = pcaModel.transform(bindAngles.reshape(len(bindAngles), int(bindAngles.shape[-2] * bindAngles.shape[-1]))) # transform complex trajectory to free aptamer pca basis

        reducedDifferences = np.zeros(n_components)  #
        for i in range(n_components):  # rather than high dimensional probability analysis, we'll instead do differences in averages over dimensions
            # we could instead do this directly with the angles - since it is a bunch of 1D analyses it is cheap. It will then be lossless, but focus less on big movements.
            reducedDifferences[i] = np.average(bindReducedTrajectory[:, i]) / np.average(np.abs(freeReducedTrajectory[:, i])) # the average of all the free trajectories is zero by construction - any nonzero average in the binding trajectory is therefore an aberration

        reducedDifference = np.average(np.abs(reducedDifferences))  # how different is the average binding aptamer conformation from free?
    except:
        print('Missing binding trajectory!')
        reducedDifference = 0

    return reducedDifference


def checkMidTrajectoryBinding(structure, trajectory, peptide, sequence, params, cutoffTime=1):
    """
    check if the analyte has come unbound from the aptamer
    and stayed unbound for a certain amount of time
    :return: True or False
    """
    u = mda.Universe(structure, trajectory)  # load up trajectory
    pepNucDists = peptideAptamerDistances(u, peptide, sequence)
    contacts, nContacts = getPeptideContacts(pepNucDists)
    lastContact = np.nonzero(nContacts[:, 0])[0][-1]  # last time when the peptide and aptamer were in close-range contact

    dt = params['print step']
    detachedTime = (len(nContacts) - lastContact) * dt / 1e3  # in ns, since dt is in ps
    if detachedTime > cutoffTime:  # if we have been unbound longer than the cutoff
        return True
    else:  # otherwise keep going
        return False


class atomDistances(AnalysisBase):
    """
    calculate watson-crick bond lengths
    """

    def __init__(self, atomgroups, **kwargs):
        """
        :param atomgroups: a list of atomgroups for which the interatom bond lengths are calculated
        """
        super(atomDistances, self).__init__(atomgroups[0].universe.trajectory, **kwargs)
        self.atomgroups = atomgroups
        self.ag1 = atomgroups[0]
        self.ag2 = atomgroups[1]

    def _prepare(self):
        self.results = []

    def _single_frame(self):
        distance = calc_bonds(self.ag1.positions, self.ag2.positions, box=self.ag1.dimensions)
        self.results.append(distance)

    def _conclude(self):
        self.results = np.asarray(self.results)


def wcTrajAnalysis(u):
    """
    use the atomDistances class to calculate the WC base pairing distances between all bases on a sequence
    :param u:
    :return:
    """
    n_bases = u.segments[0].residues.n_residues
    atomIndices1 = np.zeros((n_bases, n_bases))
    atomIndices2 = np.zeros_like(atomIndices1)
    # identify relevant atoms for WC distance calculation
    for i in range(1, n_bases + 1):
        for j in range(1, n_bases + 1):
            if u.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DC", "DT", "U", "C", "T", "CYT", "THY", "URA"]:
                a1, a2 = "N3", "N1"
            if u.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DG", "DA", "A", "G", "ADE", "GUA"]:
                a1, a2 = "N1", "N3"
            atoms = u.select_atoms("(resid {0!s} and name {1!s}) or (resid {2!s} and name {3!s}) ".format(i, a1, j, a2))
            atomIndices1[i - 1, j - 1] = atoms[0].id  # bond-mate 1
            if i == j:
                atomIndices2[i - 1, j - 1] = atoms[0].id  # if it's the same base, we want the resulting distance to always be zero
            else:
                atomIndices2[i - 1, j - 1] = atoms[1].id  # bond-mate 2

    # make a flat list of every combination
    bonds = [mda.AtomGroup(atomIndices1.flatten(), u), mda.AtomGroup(atomIndices2.flatten(), u)]
    na = atomDistances(bonds).run()
    traj = na.results.reshape(u.trajectory.n_frames, n_bases, n_bases)

    return traj


def nucleicDihedrals(u):
    """
    use analysis.dihedral to quickly compute dihedral angles for a given DNA sequence
    :param u:
    :return: dihderals
    """
    n_bases = u.segments[0].residues.n_residues
    seg = "A"
    a, b, g, d, e, z, c = [[], [], [], [], [], [], []]
    for i in range(2, n_bases - 1):  # cutoff end bases to ensure we always have 4 atoms for every dihedral unit
        a.append(u.select_atoms(" atom {0!s} {1!s} O3\' ".format(seg, i - 1),
                                " atom {0!s} {1!s} P  ".format(seg, i),
                                " atom {0!s} {1!s} O5\' ".format(seg, i),
                                " atom {0!s} {1!s} C5\' ".format(seg, i)))

        b.append(u.select_atoms(" atom {0!s} {1!s} P    ".format(seg, i),
                                " atom {0!s} {1!s} O5\' ".format(seg, i),
                                " atom {0!s} {1!s} C5\' ".format(seg, i),
                                " atom {0!s} {1!s} C4\' ".format(seg, i)))

        g.append(u.select_atoms(" atom {0!s} {1!s} O5\' ".format(seg, i),
                                " atom {0!s} {1!s} C5\' ".format(seg, i),
                                " atom {0!s} {1!s} C4\' ".format(seg, i),
                                " atom {0!s} {1!s} C3\' ".format(seg, i)))

        d.append(u.select_atoms(" atom {0!s} {1!s} C5\' ".format(seg, i),
                                " atom {0!s} {1!s} C4\' ".format(seg, i),
                                " atom {0!s} {1!s} C3\' ".format(seg, i),
                                " atom {0!s} {1!s} O3\' ".format(seg, i)))

        e.append(u.select_atoms(" atom {0!s} {1!s} C4\' ".format(seg, i),
                                " atom {0!s} {1!s} C3\' ".format(seg, i),
                                " atom {0!s} {1!s} O3\' ".format(seg, i),
                                " atom {0!s} {1!s} P    ".format(seg, i + 1)))

        z.append(u.select_atoms(" atom {0!s} {1!s} C3\' ".format(seg, i),
                                " atom {0!s} {1!s} O3\' ".format(seg, i),
                                " atom {0!s} {1!s} P    ".format(seg, i + 1),
                                " atom {0!s} {1!s} O5\' ".format(seg, i + 1)))

        # c.append(u.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i), # doesn't always get 4 atoms
        #                          " atom {0!s} {1!s} C1\' ".format(seg, i),
        #                          " atom {0!s} {1!s} N9 ".format(seg, i),
        #                          " atom {0!s} {1!s} C4  ".format(seg, i)))

    atomList = [a, b, g, d, e, z]

    dihedrals = np.zeros((u.trajectory.n_frames, len(a), 6))  # initialize combined trajectory
    for i in range(len(atomList)):  # for each type of dihedral, compute the trajectory for all bases
        analysis = Dihedral(atomList[i])
        analysis.run()
        dihedrals[:, :, i] = analysis.angles

    return dihedrals % 360  # convert from -180:180 to 0:360 basis


def numbers2letters(sequences):  # Transforming letters to numbers:
    """
    Converts numerical values to ATGC-format
    :param sequences: numerical DNA sequences to be converted
    :return: DNA sequences in ATGC format
    """
    if type(sequences) != np.ndarray:
        sequences = np.asarray(sequences)

    my_seq = ["" for x in range(len(sequences))]
    row = 0
    for j in range(len(sequences)):
        seq = sequences[j, :]
        assert type(seq) != str, 'Function inputs must be a list of equal length strings'
        for i in range(len(sequences[0])):
            na = seq[i]
            if na == 0:
                my_seq[row] += 'A'
            elif na == 1:
                my_seq[row] += 'T'
            elif na == 2:
                my_seq[row] += 'C'
            elif na == 3:
                my_seq[row] += 'G'
        row += 1
    return my_seq


def killH(structure):
    """
    use simTk modeller to delete all atoms with 'H' in the name
    """
    pdb = app.PDBFile(structure)
    topology = pdb.topology
    positions = pdb.positions
    modeller = app.Modeller(topology, positions)
    modeller.delete(atom for atom in topology.atoms() if "H" in atom.name)
    app.PDBFile.writeFile(modeller.topology, modeller.positions, open(structure.split('.')[0] + '_noH.pdb', 'w'))


def changeSegment(structure, oldSeg, newSeg):
    """
    change the segment ID for all molecule(s) in a pdb file
    """
    replaceText(structure, ' ' + oldSeg + ' ', ' ' + newSeg + ' ')


def addH(structure, pH):
    """
    protonate a given structure
    :param file:
    :return:
    """
    pdb = app.PDBFile(structure)
    topology = pdb.topology
    positions = pdb.positions
    modeller = app.Modeller(topology, positions)
    modeller.addHydrogens(pH=pH)
    app.PDBFile.writeFile(modeller.topology, modeller.positions, open(structure.split('.')[0] + '_H.pdb', 'w'))


def getMoleculeSize(structure):
    """
    use MDAnalysis to determine the cubic xyz dimensions for a given molecule
    return these dimensions
    """
    u = mda.Universe(structure)
    positions = u.atoms.positions
    dimensions = np.zeros(3)
    for i in range(3):
        dimensions[i] = np.ptp(positions[:, i])

    return dimensions


def lightDock(aptamer, analyte, params):
    """
    setup lightdock run, do the run, generate output models, hierarchical clustering and final concatenation of results
    """
    # identify aptamer size
    dimensions = getMoleculeSize(aptamer)
    # compute surface area of rectangular prism with no offset
    surfaceArea = 2 * (dimensions[0] * dimensions[1] + dimensions[1] * dimensions[2] + dimensions[2] * dimensions[1])  # in angstroms
    # compute number of swarms which cover this surface area - swarms are spheres 2 nm in diameter - this is only approximately the right number of swarms, since have no offset, and are using a rectangle
    nSwarms = np.ceil(surfaceArea / (np.pi * 10 ** 2))

    params['swarms'] = int(nSwarms)  # number of glowworm swarms
    params['glowworms'] = 300  # number of glowworms per swarm

    killH(aptamer)  # DNA needs to be deprotonated
    addH(analyte, params['pH'])  # peptide needs to be hydrogenated

    aptamer2 = aptamer.split('.')[0] + "_noH.pdb"
    analyte2 = analyte.split('.')[0] + "_H.pdb"
    changeSegment(analyte2, 'A', 'B')

    # run setup
    os.system(params['ld setup path'] + ' ' + aptamer2 + ' ' + analyte2 + ' -anm -s ' + str(params['swarms']) + ' -g ' + str(params['glowworms']) + ' >> outfiles/lightdockSetup.out')

    # run docking
    os.system(params['ld run path'] + ' setup.json ' + str(params['docking steps']) + ' -s dna >> outfiles/lightdockRun.out')

    # generate docked structures and cluster them
    for i in range(params['swarms']):
        os.chdir('swarm_%d' % i)
        os.system(params['lgd generate path'] + ' ../' + aptamer2 + ' ../' + analyte2 + ' gso_%d' % params['docking steps'] + '.out' + ' %d' % params[
            'glowworms'] + ' > /dev/null 2> /dev/null; >> generate_lightdock.list')  # generate configurations
        os.system(params['lgd cluster path'] + ' gso_%d' % params['docking steps'] + '.out >> cluster_lightdock.list')  # cluster glowworms
        os.chdir('../')

    os.system(params['lgd rank path'] + ' %d' % params['swarms'] + ' %d' % params['docking steps'])  # rank the clustered docking setups

    # generate top structures
    os.system(params['lgd top path'] + ' ' + aptamer2 + ' ' + analyte2 + ' rank_by_scoring.list %d' % params['N docked structures'])
    os.mkdir('top')
    os.system('mv top*.pdb top/')  # collect top structures (clustering currently dubiously working)

    topScores = readInitialLines('rank_by_scoring.list', params['N docked structures'] + 1)[1:]
    for i in range(len(topScores)):
        topScores[i] = float(topScores[i].split(' ')[-1])  # get the last number, which is the score

    return topScores


def appendTrajectory(topology, original, new):
    """
    use mda to combine old and new MD trajectory files
    """
    trajectories = [original, new]
    u = mda.Universe(topology, trajectories)
    with mda.Writer('combinedTraj.dcd', u.atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(u)


def checkTrajPCASlope(topology, trajectory):
    """
    analyze the trajectory to see if it's converged
    """
    converged = 0
    cutoff = 1e-2

    u = mda.Universe(topology, trajectory)

    baseDists = dnaBaseCOGDist(u)  # FAST, base-base center-of-geometry distances
    baseAngles = nucleicDihedrals(u)  # FAST, new, omits 'chi' angle between ribose and base

    mixedTrajectory = np.concatenate((baseDists.reshape(len(baseDists), int(baseDists.shape[-2] * baseDists.shape[-1])), baseAngles.reshape(len(baseAngles), int(baseAngles.shape[-2] * baseAngles.shape[-1]))),
                                     axis=1)  # mix up all our info

    representativeIndex, pcTrajectory, eigenvalues = isolateRepresentativeStructure(mixedTrajectory)  # do dimensionality reduction

    slopes = np.zeros(pcTrajectory.shape[-1])
    for i in range(len(slopes)):
        slopes[i] = np.abs(np.polyfit(np.arange(len(pcTrajectory)), pcTrajectory[:, i], 1)[0])

    combinedSlope = np.sum(slopes * (eigenvalues / np.sum(eigenvalues))) # normalize the components contributions by their eigenvalues

    print('PCA slope average is %.4f' % combinedSlope)

    return combinedSlope


def nupackAnalysis(sequence, temperature, ionicStrength):
    '''
    sequence is DNA FASTA
    temperature in C or K
    ionicStrength in M
    return a bunch of analysis
    '''
    R = 0.0019872  # ideal gas constant in kcal/mol/K
    if temperature > 273: # auto-detect Kelvins
        gap = 2 * R * temperature
        CelsiusTemperature = temperature - 273
    else:
        gap = 2 * R * (temperature + 273) # convert to Kelvin for kT
        CelsiusTemperature = temperature

    A = Strand(sequence, name='A')
    comp = Complex([A], name='AA')
    set1 = ComplexSet(strands=[A], complexes=SetSpec(max_size=1, include=[comp]))
    model1 = Model(material='dna', celsius=CelsiusTemperature, sodium=ionicStrength)
    results = complex_analysis(set1, model=model1, compute=['pfunc', 'mfe', 'subopt', 'pairs'], options={'energy_gap': gap})
    cout = results[comp]
    prob = np.exp(-cout.mfe[0].energy / R / temperature) / float(cout.pfunc)  # probability - energy/partition function
    ssString = cout.mfe[0].structure
    # look at suboptimal structures
    subopts = cout.subopt
    subStrings, subProbs, subEns, subStack, subProbs = [[], [], [], [], []]
    for subopt in subopts:
        subProbs.append(np.exp(-subopt.energy / R / temperature) / float(cout.pfunc))

    ssString = str(ssString)
    pairList = ssToList(ssString)

    # extra analysis
    nPairs = ssString.count('(')  # number of pairs
    pairFrac = 2 * nPairs / len(sequence)  # fraction of paired bases

    nPins = 0  # number of distinct hairpins
    indA = 0
    for j in range(len(sequence)):
        if ssString[j] == '(':
            indA += 1
        elif ssString[j] == ')':
            indA -= 1
            if indA == 0:  # if we come to the end of a distinct hairpin
                nPins += 1

    ssDict = {}  # best structure dictionary
    ssDict['2d string'] = ssString
    ssDict['pair list'] = pairList
    ssDict['num pairs'] = nPairs
    ssDict['pair fact'] = pairFrac
    ssDict['num pins'] = nPins
    ssDict['state prob'] = prob

    return ssDict, [subopts, subProbs]


def getSeqfoldStructure(sequence, temperature):
    """
    output the secondary structure for a given sequence at a given condition
    formats - ss string and pair list
    """
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


def ssToList(ssString):
    """
    if for some reason we have a secondary structure string we need to convert to a pair list
    """
    pairList = []
    paired = []
    for i in range(len(ssString)):
        if ssString[i] == '(':  # if it's paired
            counter = 0
            for j in range(1, len(ssString[i:])):  # look for the thing paired to it
                if ssString[i + j] == '(':
                    counter += 1
                if ssString[i + j] == ')':
                    if counter > 0:
                        counter -= 1
                    elif counter == 0:
                        # check for duplicates
                        if (not i in paired) and (not i + j in paired):  # check for duplicates
                            paired.append(i)
                            paired.append(i + j)
                            pairList.append([i + 1, i + j + 1])  # make pair list in 1-n basis
                            break

    return pairList


def recenterDCD(topology, trajectory):
    """
    topology as pdb
    trajectory as dcd
    creates a new dcd without periodic artifacts
    """
    traj = md.load(trajectory, top = topology)
    traj.image_molecules()
    traj.save(trajectory.split('.')[0] + '_recentered.dcd')


def findPairingErrors(config, tt):
    '''
    given a pair config in the format [x1, x2, x3, x4 ..., xn]
    with x as integers indexed from 1
    identify whether the any elements are paired wth more than one other element
    return a list of errors and the time index
    '''
    pairingErrors = []
    for i in range(len(config)):
        if config[i] != 0:  # if it's paired
            pairMate = int(config[i] - 1)  # indexing - what is 'i' bonded to
            if config[pairMate] != i + 1:  # indexing - is that thing bonded to 'i'?
                pairingErrors.append([tt, i, pairMate])

    return pairingErrors


def configToString(config):
    '''
    convert 2d structure 'config' to dot-bracket string
    :param config:
    :return: string
    '''
    string = '.' * len(config) # start every base as unpaired
    for i in range(len(config)):
        if (config[i] == 0) or (string[i] != '.'): # if the base is unpaired or previously notated, leave it
            pass
        else: # make a pair
            pair = [i, config[i] - 1] # indexing
            stringList = list(string)
            if pair[0] < pair[1]:
                stringList[pair[0]] = '('
                stringList[pair[1]] = ')' # indexing
            elif pair[1] > pair[0]: # possibility for backwards pairs
                stringList[pair[1]] = '('
                stringList[pair[0]] = ')'
            string = ''.join(stringList)
            #print(str(pair) + string) # for debugging


    return string


def secondaryStructureClustering(configs, probs):
    '''
    :configs: list of configurations to cluster
    :return:
    '''

    # ok but I'd like to have a threshold
    # kmeans = cluster.KMeans(n_clusters = 4, random_state=0).fit(dists) # kmeans cluster
    # labels = kmeans.labels_

    # nice, nice - getting somewhere - this one uses a minimum threshold rather than a maximum
    # then algorithm depends on quality of distance metric
    converged = False
    threshold = 0.1  # threshold (on normalized distance metric
    while not converged:
        agglomerate = cluster.AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='average', compute_full_tree=True, distance_threshold=threshold).fit(dists)
        labels = agglomerate.labels_
        nClusters = agglomerate.n_clusters_

        clusters = []
        clusterProbs = []
        clusterDists = []
        probSums = np.zeros(len(np.unique(labels)))
        for i in range(len(np.unique(labels))):
            inds = np.where(labels == i)[0].astype(int)
            clusters.append([configs[j] for j in inds])
            clusterProbs.append([probs[j] for j in inds])
            probSums[i] = np.sum(clusterProbs[-1])  # unnormalized
            clusterDists.append(getSecondaryStructureDistance([configs[j] for j in inds]))

        normedProbSums = probSums / np.sum(probSums)

        if np.average(normedProbSums) < 0.1:  # if the average normed probability of our clusters is less than 10%, we're not converged
            threshold *= 1.1  # if it isn't converged, boost the threshold
        else:
            converged = True

    # find a representative from each cluster
    weightedAvgDistance = []
    clusterReps = []
    for i in range(nClusters):
        weightedAvgDistance.append((clusterDists[i] + np.eye(len(clusterDists[i]))) @ (1 / np.asarray(clusterProbs[i])))  # weight the distance to each config against its probability of being seen
        clusterReps.append(clusters[i][np.argmin(weightedAvgDistance[i])])

    return clusterReps, clusterProbs