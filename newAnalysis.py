from utils import *
import os
from sklearn.cluster import DBSCAN
from sklearn.cluster import OPTICS
#os.chdir('C:/Users/mikem/Desktop/mmruns/run47') # small , easy
#u = mda.Universe('cleanaptamer.pdb','cleanaptamer.dcd')
#sequence = 'CGCTTTGCG'

os.chdir('C:/Users\mikem\Desktop\mmruns\cluster/run5') # big, hard
sequence = 'TATGCATGTGGGCGACGCAGTGCCCGTGGGATTTACTTGCAC'
u = mda.Universe('cleanComplex.pdb','cleanTrajectory.dcd')



class atomDistances(AnalysisBase):
    '''
    calculate watson-crick bond lengths
    '''
    def __init__(self, atomgroups, **kwargs):
        '''
        :param atomgroups: a list of atomgroups for which the wc lengths are calculated
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

'''
aind = []
bind = []
for i in range(50):
    for j in range(50):
        aind.append(i)
        bind.append(j)

a = mda.AtomGroup(aind,u)
b = mda.AtomGroup(bind,u)

na = atomDistances([a,b]).run()
traj = na.results
'''

def wcTrajAnalysis(u):
    '''
    use the atomDistances class to calculate the WC base pairing distances between all bases on a sequence
    :param u:
    :return:
    '''
    n_bases = u.residues.n_residues
    atomIndices1 = np.zeros((n_bases,n_bases))
    atomIndices2 = np.zeros_like(atomIndices1)
    # identify relevant atoms for WC distance calculation
    for i in range(1,n_bases + 1):
        for j in range(1,n_bases + 1):
            if u.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DC", "DT", "U", "C", "T", "CYT", "THY", "URA"]:
                a1, a2 = "N3", "N1"
            if u.select_atoms(" resid {0!s} ".format(i)).resnames[0] in ["DG", "DA", "A", "G", "ADE", "GUA"]:
                a1, a2 = "N1", "N3"
            atoms = u.select_atoms("(resid {0!s} and name {1!s}) or (resid {2!s} and name {3!s}) ".format(i, a1, j, a2))

            atomIndices1[i-1,j-1]=atoms[0].id # bond-mate 1
            if i == j:
                atomIndices2[i-1,j-1]=atoms[0].id # if it's the same base, we want the resulting distance to always be zero
            else:
                atomIndices2[i-1,j-1]=atoms[1].id # bond-mate 2

    # make a flat list of every combination
    bonds = [mda.AtomGroup(atomIndices1.flatten(),u),mda.AtomGroup(atomIndices2.flatten(),u)]
    na = atomDistances(bonds).run()
    traj = na.results.reshape(u.trajectory.n_frames,n_bases,n_bases)

    return traj


def nucleicDihedrals(u):
    '''
    use analysis.dihedral to quickly compute dihedral angles for a given DNA sequence
    :param u:
    :return: dihderals
    '''
    n_bases = u.residues.n_residues
    seg = "A"
    a,b,g,d,e,z,c = [[],[],[],[],[],[],[]]
    for i in range(2,n_bases-1): # cutoff end bases to ensure we always have 4 atoms for every dihedral unit
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

        #c.append(u.select_atoms(" atom {0!s} {1!s} O4\' ".format(seg, i), # doesn't always get 4 atoms
        #                          " atom {0!s} {1!s} C1\' ".format(seg, i),
        #                          " atom {0!s} {1!s} N9 ".format(seg, i),
        #                          " atom {0!s} {1!s} C4  ".format(seg, i)))

    atomList = [a,b,g,d,e,z]

    dihedrals = np.zeros((u.trajectory.n_frames,len(a),6)) # initialize combined trajectory
    for i in range(len(atomList)): # for each type of dihedral, compute the trajectory for all bases
        analysis = Dihedral(atomList[i])
        analysis.run()
        dihedrals[:,:,i] = analysis.angles

    return dihedrals % 360 # convert from -180:180 to 0:360 basis


baseWC = wcTrajAnalysis(u)

baseAngles = nucleicDihedrals(u)

pairingTrajectory = getPairs(baseWC)

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

optic = OPTICS(min_samples=2,metric='precomputed').fit(configDistance)

avgDist = np.average(configDistance)


configArray = np.asarray(configs)

clusterConfigs = []
clusterDists = []
for i in range(len(np.unique(optic.labels_))-1):
    clusterConfigs.append(configArray[optic.labels_==i])
    clusterDists.append(np.average(getSecondaryStructureDistance(clusterConfigs[i])))

'''

equilibrationOffset = 1000

converged = False
nComponents = 10
while converged == False: # add components until some threshold, then take components greater than the average
    components, eigenvalues0, empty = trajectoryPCA(nComponents, baseAngles[equilibrationOffset:], False)  # do it with a bunch of components, then re-do it with only the necessary number
    totVariance = np.sum(eigenvalues0)
    if totVariance > 0.85:
        converged = True
    else:
        nComponents += 1

# we want there to be a gap in this spectrum, or at least, to neglect only the small contributions
n_components = np.sum(eigenvalues0 > np.average(eigenvalues0)) # with the threshold above, this removes some of the variation from having too many components
components, eigenvalues, reducedTrajectory = trajectoryPCA(n_components, baseAngles[equilibrationOffset:], True)  # do it with a bunch of components, then re-do it with only the necessary number
