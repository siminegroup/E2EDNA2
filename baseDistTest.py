import MDAnalysis as mda
import numpy as np
import os
import MDAnalysis.analysis.nuclinfo as nuclinfo
from MDAnalysis.analysis import distances
import time
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import scipy.ndimage as ndimage
import scipy.spatial as spatial
from utils import *


#os.chdir('C:/Users\mikem\Desktop\mmruns\cluster/run5')
#sequence = 'TATGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCAC'

os.chdir('C:/Users\mikem\Desktop\mmruns/run28')
sequence = 'CGCTTTGCG'

structure = 'cleanComplex.pdb'
trajectory = 'cleanTrajectory.dcd'

u = mda.Universe(structure,trajectory)

# SLOW
baseWC = dnaBasePairDist(u, sequence)  # watson-crick base pairing distances (H-bonding)
# FAST
baseDists = dnaBaseCOGDist(u, sequence)  # base-base center-of-geometry distances
# SLOW
baseAngles = dnaBaseDihedrals(u, sequence)  # base-base dihedrals

# 2D structure analysis
pairingTrajectory = getPairs(baseWC)
# for i in range(len(sequence)): # find the equilibrium secondary structure

# 3D structure analysis
representativeIndex = isolateRepresentativeStructure(baseAngles)
# save this structure as a separate file
extractFrame(structure, trajectory, representativeIndex)


# we also need a function which processes the energies spat out by the trajectory thingy

analysisDict = {}  # compile our results

'''
first working version
def pairTest(u,sequence,length,bases):
    pairDists = np.zeros((length,bases,bases)) # matrix of watson-crick distances
    tt = 0
    for ts in u.trajectory[:length]:
        for i in range(bases):
            for j in range(bases):
                if i > j:
                    pairDists[tt,i,j] = nuclinfo.wc_pair(u,i+1,j+1,seg1='A',seg2='A') # also do WC base-pair distances (can set to only follow secondary structure prediction)
                    pairDists[tt,j,i] = pairDists[tt,i,j]
        tt += 1

    return pairDists

def baseDists(u,sequence,length,bases):
    baseDists = np.zeros((length,bases,bases)) # matrix of watson-crick distances
    tt = 0
    for ts in u.trajectory[:length]:
        dnaResidues = u.segments[0]
        posMat = np.zeros((bases, 3))

        for i in range(bases):
            posMat[i] = dnaResidues.residues[i].atoms.center_of_geometry()

        baseDists[tt, :, :] = distances.distance_array(posMat, posMat, box=u.dimensions)
        tt += 1

    return baseDists

def baseDihedrals(u,sequence,length,bases):
    angles = np.zeros((length,bases,7))
    tt = 0
    for ts in u.trajectory[:length]:
        for j in range(2,bases-1):
            angles[tt,j,:] = nuclinfo.tors(u,seg="A",i=j)

        tt += 1

    return angles

def getPairs(wcDist):
    trajTime = len(wcDist)
    seqLen = wcDist.shape[-1]
    pairedBases = np.zeros((trajTime,seqLen))
    for tt in range(trajTime):
        pairMat = wcDist[tt] + np.eye(seqLen) * 20
        for i in range(len(sequence)):
            nearestNeighbour = np.argmin(pairMat[i,:])
            if wcDist[tt,i,nearestNeighbour] < 3.3: # if we're within a hydrogen bond length, count it as a pair
                pairedBases[tt,i] = nearestNeighbour + 1

    return pairedBases

def trajectoryPCA(n_components,trajectory):
    '''
    do PCA on some trajectory with pre-selected features
    :param n_components:
    :param trajectory:
    :return: the principal components, their relative contributions, the original trajectory in PC basis
    '''
    pca1 = PCA(n_components=n_components)
    model = pca1.fit(trajectory.reshape(len(angles),int(angles.shape[-1]*angles.shape[-2])))
    components = model.components_
    eigenvalues = pca1.explained_variance_ratio_
    reducedTrajectory = model.transform(trajectory.reshape(len(angles),int(angles.shape[-1]*angles.shape[-2])))

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


#0) Compute trajectories of the relevant features
sequence = sequence[0:10]
trajTime = 50
wcDist = pairTest(u,sequence,trajTime,len(sequence))
baseDist = baseDists(u,sequence,trajTime,len(sequence))
angles = baseDihedrals(u,sequence,trajTime,len(sequence))

# analysis
#1) identify paired bases
pairedBases = getPairs(wcDist)

#2) identify contacts
#3) identify collective variables # https://mdtraj.org/1.9.4/examples/pca.html
components, eigenvalues, reducedTrajectory = trajectoryPCA(10, angles) # do it with a bunch of components, then re-do it with only the necessary number
# we want there to be a gap in this spectrum, or at least, to neglect only the small contributions
n_components = np.sum(np.abs(np.diff(eigenvalues)) > np.mean(np.abs(np.diff(eigenvalues)))) # this is kindof a hack - we can code better logic later
components, eigenvalues, reducedTrajectory = trajectoryPCA(n_components, angles) # redo it only with the necessary eigenvalues

#4) find minima in the space of collective variables
nbins = int(1e5 ** (1/n_components))
probs, bins = np.histogramdd(reducedTrajectory,bins=nbins) # multidimensional probability histogram
# smooth it out and take the diverse maxima - compare relative probability
sigma = 5
smoothProbs = ndimage.gaussian_filter(probs,sigma)
bestTransformedStructureIndex = np.unravel_index(smoothProbs.argmax(),smoothProbs.shape)

#4b) analyze at the dynamics in the collective basis and see if there is any interesting dynamical information

#5) Pull representative 3D structures from the original trajectory (rather than back-transforming and generating a new structure from only dihedrals - HARD)
transformedDimensions = []
for i in range(len(bins)):
    transformedDimensions.append(bins[i][1:] - np.diff(bins[i][0:2]))

transformedDimensions = np.asarray(transformedDimensions)
transformedCoordinates = []
for i in range(len(transformedDimensions)):
    transformedCoordinates.append(transformedDimensions[i][bestTransformedStructureIndex[i]]) # best structure coordinates in transformed basis

transformedCoordinates = np.asarray(transformedCoordinates)

# find structures in the reduced trajectory with these coordinates
representativeIndex = do_kdtree(reducedTrajectory,transformedCoordinates)

# save this structure as a separate file
extractFrame('cleanComplex.pdb', 'cleanTrajectory.dcd',representativeIndex)
'''