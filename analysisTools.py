'''
E2EDNA 2.0 - OpenMM Implementation of E2EDNA !

An automated pipeline for simulating DNA aptamers complexed with target ligands (peptide, DNA, RNA or small molecules).
    
Copyright (C) 2022 Michael Kilgour, Tao Liu and Lena Simine

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

# tools for trajectory analysis
from utils import printRecord
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import scipy.ndimage as ndimage
import scipy.spatial as spatial
import sklearn.cluster as cluster
from sklearn.decomposition import PCA
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import calc_bonds
from MDAnalysis.analysis.dihedrals import Dihedral
from MDAnalysis.analysis.bat import BAT
from seqfold import dg, fold
from openmm import *
from openmm.app import *
import openmm.unit as unit
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
import gzip



# MDA Nucleic + trajectory tools
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

# noinspection PyTypeChecker
def getBaseBaseDistTraj(u):
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

def getPairTraj(wcTraj):
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
        printRecord('Encountered error in pairing algorithm!') # should not be possible, currently

    return pairedBases

def getPepBaseDistTraj(u, peptideSeq, aptamerSeq):
    """
    given an MDA universe
    return distances between each peptide and each base
    at every timepoint in the trajectory
    :return: pepNucDists - peptide-nucleicacid distances
    """

    pepNucDists = np.zeros((len(u.trajectory), len(peptideSeq), len(aptamerSeq)))  # distances between peptides and nucleotiedes
    tt = 0
    for ts in u.trajectory:
        targetResidues = u.segments[1]
        baseResidues = u.segments[0]
        posMat1 = np.zeros((len(peptideSeq), 3))
        posMat2 = np.zeros((len(aptamerSeq), 3))
        for i in range(len(peptideSeq)):
            posMat1[i] = targetResidues.residues[i].atoms.center_of_geometry()
        for i in range(len(aptamerSeq)):
            posMat2[i] = baseResidues.residues[i].atoms.center_of_geometry()

        pepNucDists[tt, :, :] = distances.distance_array(posMat1, posMat2, box=u.dimensions)
        tt += 1

    return pepNucDists

def getPepContactTraj(pepNucDists):
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

def getWCDistTraj(u):
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

def getNucDATraj(u):
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


# Secondary structure analysis utils
def pairListToConfig(pairList, seqLen):
    '''
    convert a secondary structure for a list of pairs to a 'config'
    which is where every base is listed with its pair-mate (or 0 for unpaired)
    :param pairList: list of paired bases
    :param seqLen: length of the DNA sequence
    :return:
    '''
    config = np.zeros(seqLen).astype(int)
    for i in range(len(pairList)):
        config[pairList[i][0] - 1] = pairList[i][1] # indexing
        config[pairList[i][1] - 1] = pairList[i][0]

    return config

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


def getSeqfoldStructure(aptamerSeq, temperature):
    """
    output the secondary structure for a given sequence at a given condition
    formats - ss string and pair list
    """
    dg(aptamerSeq, temp=temperature)  # get energy of the structure
    # printRecord(round(sum(s.e for s in structs), 2)) # predicted energy of the final structure

    structs = fold(aptamerSeq)  # identify structural features
    desc = ["."] * len(aptamerSeq)
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
            #printRecord(str(pair) + string) # for debugging


    return string

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
    configs, counter = enumerate2DConfigurations(pairTraj)
    configDistances = getSecondaryStructureDistance(configs)
    bestSingleStructure = configs[np.argmax(counter)]

    if len(configs) > 1:
        topClusterReps, topClusterProbs, topClusterDists = do2DAgglomerativeClustering(configs, counter, configDistances)
        # compute a distance metric between all the seen configs
        bestClusterRep = topClusterReps[0]
        singleVsClusterDistance = getSecondaryStructureDistance([bestSingleStructure,bestClusterRep])[0,1]
        if singleVsClusterDistance > 0.1:
            printRecord('Best 2D cluster representative is more than 10% different from best single structure')

        return bestClusterRep
    else: # if we only found a single conformation (unlikely except for very short time dynamics)
        return bestSingleStructure

def do2DAgglomerativeClustering(configs, probs, distances, distThreshold=0.1, probThreshold=0.05, minClusters=5):
    '''
    run agglomerative hierarchical clustering on secondary structure configurations
    increase the distance threshold until we reach the desired number and probability density of clusters
    select representative cluster samples, incorporating relative probability of observation
    :return:
    '''
    converged = False
    while not converged:
        agglomerate = cluster.AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='average', compute_full_tree=True, distance_threshold=distThreshold).fit(distances)
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

        nThresholdClusters = np.sum(normedProbSums > probThreshold)
        if (nThresholdClusters >= minClusters) or (nClusters <= minClusters):  # if we have the desired number of clusters above the threshold, we can stop
            converged = True
        else:
            distThreshold *= 1.1  # if it isn't converged, boost the threshold

    # find a representative from each cluster - product of a priori probability and distance from ensemble members
    clusterReps = []
    for i in range(nClusters):
        avgDists = np.average(clusterDists[i], axis=1)
        weightedAvgDists = (1 - avgDists) * clusterProbs[i]
        clusterRepInd = np.argmax(weightedAvgDists)
        clusterReps.append(clusters[i][clusterRepInd])

    # collect and sort the top x representatives
    sortedInds = np.argsort(probSums)
    sortedInds = sortedInds[::-1] # comes ot in reverse order so swap it
    topReps = [clusterReps[j] for j in sortedInds]
    topProbs = normedProbSums[sortedInds]

    clustersToTake = min(nThresholdClusters, nClusters)

    topReps = topReps[:clustersToTake]
    topProbs = topProbs[:clustersToTake]
    topDists = getSecondaryStructureDistance(topReps)

    return topReps, topProbs, topDists

# trajectory analysis functions
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
        printRecord('Warning, attempting PCA on a low-dimensional trajectory!')

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

def enumerate2DConfigurations(pairTraj):
    '''
    explicitly enumerate and count every unique 2D strcutre
    in the pair trajectory
    :param pairTraj:
    :return: config list, counts
    '''
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

    return configs, counter

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
    # converged = False
    # nComponents = 10
    # while converged == False:  # add components until some threshold, then take components greater than the average
    #     components, eigenvalues, reducedTrajectory, pcaModel = trajectoryPCA(nComponents, trajectory, False)  # do it with a bunch of components, then re-do it with only the necessary number
    #     totVariance = np.sum(eigenvalues)
    #     if totVariance > 0.85:
    #         converged = True
    #     else:
    #         nComponents += 1
    # # we want there to be a gap in this spectrum, or at least, to neglect only the small contributions
    # n_components = min(5, np.sum(eigenvalues > np.average(eigenvalues)))  # with the threshold above, this removes some of the variation from having too many components
    # components, eigenvalues, reducedTrajectory, pcaModel = trajectoryPCA(n_components, trajectory, True)
    
    n_components = 5
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
    eigenvalues = pcaModel.explained_variance_ratio_

    for i in range(reducedTrajectory.shape[-1]): # normalize the contribution of each PC by its eigenvalue (importance)
        reducedTrajectory[:,i] = reducedTrajectory[:,i] * eigenvalues[i]

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

    return representativeIndex, reducedTrajectory, eigenvalues

def bindingAnalysis(bindu, freeu, peptideSeq, aptamerSeq):
    """
    analyze the binding of target to aptamer by computing relative distances
    :param u:
    :peptideSeq: peptide sequence
    :aptamerSeq: aptamer sequence
    :return:
    """
    assert bindu.segments.n_segments == 2
    # identify base-target distances
    pepNucDists = getPepBaseDistTraj(bindu, peptideSeq, aptamerSeq)
    contacts, nContacts = getPepContactTraj(pepNucDists)
    if np.nonzero(nContacts[:,0])[0] != []:
        firstContact = np.nonzero(nContacts[:, 0])[0][0]  # first time when the peptide and aptamer were in close-range contact
        closeContactRatio = np.average(nContacts[firstContact:, 0] > 0)  # amount of time peptide spends in close contact with aptamer
        contactScore = np.average(nContacts[firstContact:, :] / len(peptideSeq))  # per-peptide average contact score, linear average over 8-12 angstrom
    else:
        print('Never made contact!')
        closeContactRatio = 0
        contactScore = 0

    conformationChange = getConformationChange(bindu, freeu)

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

def getConformationChange(bindu, freeu):
    """
    compare the pre-complexation (free aptamer) conformation with post-complexation
    NOTE intimately depends on naming conventions for trajectory files!
    """
    # function to analyze target impact on aptamer conformation

    freeAngles = getNucDATraj(freeu)
    bindAngles = getNucDATraj(bindu)

    n_components, freeReducedTrajectory, pcaModel = doTrajectoryDimensionalityReduction(freeAngles) # get free aptamer PCA
    bindReducedTrajectory = pcaModel.transform(bindAngles.reshape(len(bindAngles), int(bindAngles.shape[-2] * bindAngles.shape[-1]))) # transform complex trajectory to free aptamer pca basis

    reducedDifferences = np.zeros(n_components)  #
    for i in range(n_components):  # rather than high dimensional probability analysis, we'll instead do differences in averages over dimensions
        # we could instead do this directly with the angles - since it is a bunch of 1D analyses it is cheap. It will then be lossless, but focus less on big movements.
        reducedDifferences[i] = np.average(bindReducedTrajectory[:, i]) / np.average(np.abs(freeReducedTrajectory[:, i])) # the average of all the free trajectories is zero by construction - any nonzero average in the binding trajectory is therefore an aberration

    reducedDifference = np.average(np.abs(reducedDifferences))  # how different is the average binding aptamer conformation from free?

    return reducedDifference

def checkMidTrajectoryBinding(structure, trajectory, peptideSeq, aptamerSeq, params, cutoffTime=1):
    """
    check if the target has come unbound from the aptamer and stayed unbound for a certain amount of time
    :return: True or False
    """
    u = mda.Universe(structure, trajectory)  # load up trajectory
    pepNucDists = getPepBaseDistTraj(u, peptideSeq, aptamerSeq)
    contacts, nContacts = getPepContactTraj(pepNucDists)
    try:
        lastContact = np.nonzero(nContacts[:, 0])[0][-1]  # last time when the peptide and aptamer were in close-range contact

        dt = params['print_step']
        detachedTime = (len(nContacts) - lastContact) * dt / 1e3  # in ns, since dt is in ps
        if detachedTime > cutoffTime:  # if we have been unbound longer than the cutoff
            return True
        else:  # otherwise keep going
            return False
    except IndexError:
        return False # if we never attached, give up

def checkTrajPCASlope(topology, trajectory, printStep):
    """
    analyze the trajectory to see if it's converged
    """
    converged = 0
    cutoff = 1e-2

    u = mda.Universe(topology, trajectory)

    baseDists = getBaseBaseDistTraj(u)  # FAST, base-base center-of-geometry distances
    baseAngles = getNucDATraj(u)  # FAST, new, omits 'chi' angle between ribose and base

    mixedTrajectory = np.concatenate((baseDists.reshape(len(baseDists), int(baseDists.shape[-2] * baseDists.shape[-1])), baseAngles.reshape(len(baseAngles), int(baseAngles.shape[-2] * baseAngles.shape[-1]))),
                                     axis=1)  # mix up all our info

    representativeIndex, pcTrajectory, eigenvalues = isolateRepresentativeStructure(mixedTrajectory)  # do dimensionality reduction

    slopes = np.zeros(pcTrajectory.shape[-1])
    for i in range(len(slopes)):
        slopes[i] = np.abs(np.polyfit(np.arange(len(pcTrajectory)) * printStep, pcTrajectory[:, i], 1)[0]) # normalize against the time step

    normedSlope = slopes * (eigenvalues / np.sum(eigenvalues)) # normalize the components contributions by their eigenvalues
    combinedSlope = np.linalg.norm(normedSlope)

    printRecord('PCA slope average is %.4f' % combinedSlope)

    return combinedSlope

# Add DeltaGzip analysis:
def round_to_first_dec_multi_10(arr):
    return (np.round(arr, decimals=1)*10).astype(np.int8) # round to the first decimal point


def get_internal_coords(top_file, traj_file):
    u_configs_stateA = mda.Universe(top_file, traj_file, guess_bonds=False)

    protein_backbone_atoms = u_configs_stateA.select_atoms("chainID A and backbone")
    assert len(protein_backbone_atoms.fragments) == 1
    # print(f'Selected atoms are:{np.unique(protein_backbone_atoms.names)}')
    # print(protein_backbone_atoms.n_atoms,'atoms')

    # Pick an atom's Cartesian coordinates to define the translation of the molecule. It has to be a terminal atom: carbonyl oxygen of residue 1
    indx = protein_backbone_atoms.select_atoms('resid 1 and name O').ix[0]
    # print(u_configs_stateA.atoms[indx]) # the returned index is a refernce in the whole universe

    R_backbone_atoms = BAT(ag=protein_backbone_atoms, initial_atom=u_configs_stateA.atoms[indx])
    R_backbone_atoms.run(start=0, stop=5000, verbose=True) # Calculate BAT coordinates
    # print(R_backbone_atoms.results.bat.shape)
    return R_backbone_atoms.results.bat[:,6:] # remove the first 6 DOF, the rest is internal DOF


def compress_each_dof(integer_IC): # (num_frames, num_dof)
    return np.array([len(gzip.compress(np.ndarray.flatten(integer_IC[:,j]), compresslevel=9)) for j in range(integer_IC.shape[1])])


def compress_each_sample(integer_IC, responsive_dof_indx=None): # (num_frames, num_dof)
    if responsive_dof_indx is None:
        return np.array([len(gzip.compress(np.ndarray.flatten(integer_IC[i,:]), compresslevel=9)) for i in range(integer_IC.shape[0])])
    else:
        return np.array([len(gzip.compress(np.ndarray.flatten(integer_IC[i,responsive_dof_indx]), compresslevel=9)) for i in range(integer_IC.shape[0])])


def entropy_DeltaGzip(top_file_bound_state, traj_file_bound_state, top_file_free_state, traj_file_free_state, T):
    # get bound-state internal coords
    internal_coords_bound_state = get_internal_coords(top_file_bound_state, traj_file_bound_state)
    # get free-state internal coords
    internal_coords_free_state = get_internal_coords(top_file_free_state, traj_file_free_state)

    integer_IC_bound_state = round_to_first_dec_multi_10(internal_coords_bound_state)
    integer_IC_free_state = round_to_first_dec_multi_10(internal_coords_free_state)

    # find responsive DOF
    C_each_dof_bound_state = compress_each_dof(integer_IC_bound_state)
    C_each_dof_free_state = compress_each_dof(integer_IC_free_state)
    variable_dof_index = np.where(np.abs(C_each_dof_bound_state - C_each_dof_free_state) > 125)[0] # cutoff = 125 bytes

    # compress the responsive dof
    if len(variable_dof_index) == 0: variable_dof_index = None
    C_dof_each_frame_bound_state = compress_each_sample(integer_IC_bound_state, variable_dof_index) # (num_frames,)
    C_dof_each_frame_free_state = compress_each_sample(integer_IC_free_state, variable_dof_index) # (num_frames,)

    C_responsive_dof_bound_state = np.mean(C_dof_each_frame_bound_state)
    C_responsive_dof_bound_state_err = np.std(C_dof_each_frame_bound_state)
    C_responsive_dof_free_state = np.mean(C_dof_each_frame_free_state) 
    C_responsive_dof_free_state_err = np.std(C_dof_each_frame_free_state)

    kB = 8.314e-3 # kJ/mole/K

    # for a binding process:
    DeltaS_G = kB * np.log(2) * (C_responsive_dof_bound_state - C_responsive_dof_free_state)
    
    return DeltaS_G


def enthalpy_DeltaGzip(top_file, traj_file, lig_sdf_file, r_cutoff=1.0):
    '''
    r_cutoff: nanometer
    '''
    u_traj = mda.Universe(top_file, traj_file) # To retrieve the atomic Cartesian coordinates:
    pdb = PDBFile(top_file)
    ligand_molecule = Molecule.from_file(lig_sdf_file)    

    forcefield = ForceField('amber14/protein.ff14SB.xml') # if water and ions 'amber14/tip3p.xml'
    gaff = GAFFTemplateGenerator(molecules=ligand_molecule, forcefield='gaff-2.11')    
    forcefield.registerTemplateGenerator(gaff.generator) # Register the GAFF template generator to the forcefield object
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=None)
    # system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, ewaldErrorTolerance=5.0e-4, constraints=None, rigidWater=True, hydrogenMass=1.0*amu) # no difference bcuz we only extract the parameters, not run MD
    nonBondedForce = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]

    # find the indices of protein and ligand: 0-indexing
    num_atoms = nonBondedForce.getNumParticles()
    num_protein_atoms = u_traj.select_atoms('(chainID A) or (chainID B)').n_atoms
    num_lig_atoms = u_traj.select_atoms('chainID C').n_atoms
    assert num_atoms == num_protein_atoms + num_lig_atoms
    group_protein = np.arange(num_protein_atoms) # my pdb file lists protein atoms prior to ligand atoms
    group_lig = np.arange(num_lig_atoms) + num_protein_atoms

    #ONE_4PI_EPS0 = 138.935456 # this value is used in OpenMM source code as 1/(4*pi*epsilon0).Unit: (kJ/mole * nm)/e^2
    ## 1eV = 96.487 kJ/mole
    ## epsilon0 = 55.26349406 e^2/(eV * um) from Wikipedia = 0.05526349406 e^2/(eV * nm) = 0.00057276 e^2/(kJ/mole * nm)
    ## 1/(4*np.pi*0.00057276) = 138.9368523394575 (kJ/mole * nm)/e^2
    # consider dielectric constant of water at 30C and 0M NaCl 
    # t=300-273.15; print(87.740 - 0.40008*t+ 9.398*1e-4*t**2 - 1.410*1e-6*t**3)
    ONE_4PI_EPS0 = 138.935456/77.6481
    

    epsilon_pair = np.zeros((len(group_protein), len(group_lig)))
    sigma_pair = np.zeros((len(group_protein), len(group_lig)))
    chargeprod = np.zeros((len(group_protein), len(group_lig)))
    for i, index1 in enumerate(group_protein):
        charge_protein, sigma_protein, epsilon_protein = nonBondedForce.getParticleParameters(index1)
        for j, index2 in enumerate(group_lig):
            charge_lig, sigma_lig, epsilon_lig = nonBondedForce.getParticleParameters(index2)
            
            epsilon_pair[i,j] = np.sqrt(epsilon_protein*epsilon_lig).value_in_unit(kilojoule_per_mole)
            sigma_pair[i,j] = 0.5*(sigma_protein+sigma_lig).value_in_unit(nanometer)
            chargeprod[i,j] = (charge_protein*charge_lig).value_in_unit(elementary_charge**2)
    # Epsilon_pair, sigma_pairm, and chargeprod only need to be calculated once. Next just need to iterate all configurations and compute their pEnergy

    u_traj.trajectory[0]
    sele = u_traj.select_atoms('chainID A or chainID B or chainID C')
    enthalpy_binding = np.zeros(u_traj.trajectory.n_frames)
    for indx, frame in enumerate(tqdm(u_traj.trajectory)):

        for i, index1 in enumerate(group_protein):
            # charge_protein[index1], sigma_protein[index1], epsilon_protein[index1] = nonBondedForce.getParticleParameters(index1)
            # position1 = pdb.getPositions(asNumpy=True)[index1].value_in_unit(nanometer)
            position1 = sele.positions[index1] * 0.1 # nanometer

            for j, index2 in enumerate(group_lig):
                # charge2, sigma2, epsilon2 = nonBondedForce.getParticleParameters(index2)
                # position2 = pdb.getPositions(asNumpy=True)[index2].value_in_unit(nanometer)
                position2 = sele.positions[index2] * 0.1 # nanometer

                # r = Quantity(value=np.linalg.norm(position1-position2), unit=nanometer)
                r = np.linalg.norm(position1-position2) # nanometer
                if r <= r_cutoff: # this is not a great idea because electrostatics decay slowly, need to set much larger cutoff distance than 2.0 nm
                
                # # if index1 in nonBondedForce.getExceptionParameters() and index2 in nonBondedForce.getExceptionParameters():
                # #     p1, p2, chargeProd, sigma_pair, epsilon_pair = nonBondedForce.getExceptionParameters(index1)
                # #     assert p1 == index1
                # #     assert p2 == index2
                # # else:
                # epsilon_pair = np.sqrt(epsilon1*epsilon2).value_in_unit(kilojoule_per_mole)
                # sigma_pair = 0.5*(sigma1+sigma2).value_in_unit(nanometer)
                # chargeprod = (charge1*charge2).value_in_unit(elementary_charge**2)

                    enthalpy_binding[indx] = enthalpy_binding[indx] + 4*epsilon_pair[i,j]*((sigma_pair[i,j]/r)**12 - (sigma_pair[i,j]/r)**6) + ONE_4PI_EPS0*chargeprod[i,j]/r

    return enthalpy_binding #kJ/mole


def dG_DeltaGzip():
    dH_binding = enthalpy_DeltaGzip(top_file_bound_state, traj_file_bound_state, lig_sdf_file, 1.0)
    dS_binding = entropy_DeltaGzip(top_file_bound_state, traj_file_bound_state, top_file_free_state, traj_file_free_state)

    T = 300.0 # Kelvin
    deltaGzip = dH_binding - T*dS_binding

    return deltaGzip
