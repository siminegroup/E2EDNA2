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
secondaryStructure = analyzeSecondaryStructure(pairingTrajectory) # find equilibrium secondary structure
predictionError = getSecondaryStructureDistance([secondaryStructure,pairList])# compare to predicted structure

# 3D structure analysis
representativeIndex = isolateRepresentativeStructure(baseAngles)
# save this structure as a separate file
extractFrame(structure, trajectory, representativeIndex)

# we also need a function which processes the energies spat out by the trajectory thingy

analysisDict = {}  # compile our results
