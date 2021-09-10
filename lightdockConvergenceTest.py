import os
import numpy as np
from utils import *
from analysisTools import *
import tqdm
import time

'''
script to test convergence of lighdock docking with a small peptide and medium sized aptamer
currently we use swarms = 1 (scaling factor meaning 1x the calculated required number of swarms for a given surface area)
'''

# get relevant paths


params = {}
params['device'] = 'local'  # 'local' or 'cluster'
params['N docked structures'] = 10  # number of docked structures to output from the docker. If running binding, it will go this time (at linear cost)
params['pH'] = 7.4
# paths
if params['device'] == 'local':
    params['workdir'] = '/home/mkilgour/mmruns' #'/mnt/c/Users/mikem/Desktop/mmruns'
    params['mmb dir'] = '/mnt/c/Users/mikem/Desktop/software/Installer.2_14.Linux64'
    params['mmb'] = '/mnt/c/Users/mikem/Desktop/software/Installer.2_14.Linux64/MMB.2_14.Linux64'
    # lightdock python scripts
    params['ld setup path'] = 'ld_scripts/lightdock3_setup.py'
    params['ld run path'] = 'ld_scripts/lightdock3.py'
    params['lgd generate path'] = '../ld_scripts/lgd_generate_conformations.py'
    params['lgd cluster path'] = '../ld_scripts/lgd_cluster_bsas.py'
    params['lg ant path'] = 'ld_scripts/ant_thony.py'
    params['lgd rank path'] = 'ld_scripts/lgd_rank.py'
    params['lgd top path'] = 'ld_scripts/lgd_top.py'

elif params['device'] == 'cluster':
    params['workdir'] = '/home/kilgourm/scratch/mmruns'  # specify your working directory here
    params['mmb dir'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64'  # 'C:/Users/mikem/Desktop/Installer.2_14.Windows/MMB.2_14.exe'
    params['mmb'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/MMB.2_14.Linux64'
    # lightdock python scripts
    params['ld setup path'] = 'python ld_scripts/lightdock3_setup.py'
    params['ld run path'] = 'python ld_scripts/lightdock3.py'
    params['lgd generate path'] = 'python ../ld_scripts/lgd_generate_conformations.py'
    params['lgd cluster path'] = 'python ../ld_scripts/lgd_cluster_bsas.py'
    params['lg ant path'] = 'python ld_scripts/ant_thony.py'
    params['lgd rank path'] = 'python ld_scripts/lgd_rank.py'
    params['lgd top path'] = 'python ld_scripts/lgd_top.py'


# initialize working directory
os.chdir('/mnt/c/Users/mikem/Desktop/mmruns/lightdockConvergenceRuns/run1')

# initialize docker

class ld(): # lightdock
    def __init__(self, aptamer, peptide, params, ind1):
        self.setupPath = params['ld setup path']
        self.runPath = params['ld run path']
        self.genPath = params['lgd generate path']
        self.clusterPath = params['lgd cluster path']
        self.rankPath = params['lgd rank path']
        self.topPath = params['lgd top path']

        self.aptamerPDB = aptamer
        self.peptidePDB = peptide

        self.glowWorms = params['worms']
        self.dockingSteps = params['docking steps']
        self.numTopStructures = params['N docked structures']

        self.pH = params['pH']

        self.ind1 = ind1

    def getSwarmCount(self):
        # identify aptamer size
        dimensions = getMoleculeSize(self.aptamerPDB)
        # compute surface area of rectangular prism with no offset
        surfaceArea = 2 * (dimensions[0] * dimensions[1] + dimensions[1] * dimensions[2] + dimensions[2] * dimensions[1])  # in angstroms
        # compute number of swarms which cover this surface area - swarms are spheres 2 nm in diameter - this is only approximately the right number of swarms, since have no offset, and are using a rectangle
        nSwarms = np.ceil(surfaceArea / (np.pi * 10 ** 2)) * params['swarms']

        self.swarms = int(nSwarms)  # number of glowworm swarms

    def prepPDBs(self):
        killH(self.aptamerPDB)  # DNA needs to be deprotonated
        addH(self.peptidePDB, self.pH)  # peptide needs to be hydrogenated
        self.aptamerPDB2 = self.aptamerPDB.split('.')[0] + "_noH.pdb"
        self.peptidePDB2 = self.peptidePDB.split('.')[0] + "_H.pdb"
        changeSegment(self.peptidePDB2, 'A', 'B')

    def runLightDock(self):
        # run setup
        os.system(self.setupPath + ' ' + self.aptamerPDB2 + ' ' + self.peptidePDB2 + ' -anm -s ' + str(self.swarms) + ' -g ' + str(self.glowWorms) + ' >> outfiles/lightdockSetup.out')

        # run docking
        os.system(self.runPath + ' setup.json ' + str(self.dockingSteps) + ' -s dna >> outfiles/lightdockRun.out')

    def generateAndCluster(self):
        # generate docked structures and cluster them
        for i in range(self.swarms):
            os.chdir('swarm_%d' % i)
            os.system(self.genPath + ' ../' + self.aptamerPDB2 + ' ../' + self.peptidePDB2 + ' gso_%d' % self.dockingSteps + '.out' + ' %d' % self.glowWorms + ' > /dev/null 2> /dev/null; >> generate_lightdock.list')  # generate configurations
            os.system(self.clusterPath + ' gso_%d' % self.dockingSteps + '.out >> cluster_lightdock.list')  # cluster glowworms
            os.chdir('../')

    def rank(self):
        os.system(self.rankPath + ' %d' % self.swarms + ' %d' % self.dockingSteps + ' >> outfiles/ld_rank.out')  # rank the clustered docking setups

    def extractTopStructures(self):
        # generate top structures
        os.system(self.topPath + ' ' + self.aptamerPDB2 + ' ' + self.peptidePDB2 + ' rank_by_scoring.list %d' % self.numTopStructures + ' >> outfiles/ld_top.out')
        os.mkdir('top_%d'%self.ind1)
        os.system('mv top*.pdb top_%d'%self.ind1 + '/')  # collect top structures (clustering currently dubiously working)

    def extractTopScores(self):
        self.topScores = readInitialLines('rank_by_scoring.list', self.numTopStructures + 1)[1:]
        for i in range(len(self.topScores)):
            self.topScores[i] = float(self.topScores[i].split(' ')[-1])  # get the last number, which is the score


    def run(self):
        self.getSwarmCount()
        self.prepPDBs()
        self.runLightDock()
        self.generateAndCluster()
        self.rank()
        self.extractTopStructures()
        self.extractTopScores()
        # cleanup
        dir = 'lightdockOutputs_{}'.format(self.ind1)
        os.mkdir(dir)
        os.system('mv swarm* ' + dir)
        os.system('mv lightdock_* ' + dir)
        os.system('mv setup.json ' + dir)
        os.system('mv *.list ' + dir)
        os.system('mv lightdock.info ' + dir)
        os.system('mv init ' + dir)



# get structures
aptamer = 'cleanComplex.pdb'
peptide = 'peptide.pdb'

# over a loop, do docking


swarms = [0.25, 0.5, 0.75, 1, 1.5, 2]
worms = [10, 25, 50, 100, 200, 300]
steps = [50, 100, 200]

topResults = []

ind = 0
for i in tqdm.tqdm(range(len(swarms))):
    for j in range(len(worms)):
        for k in range(len(steps)):
            t0 = time.time()
            params['docking steps'] = steps[k]  # number of steps for docking simulations
            params['swarms'] = swarms[i]
            params['worms'] = worms[j]
            docker = ld(aptamer,peptide,params,ind)
            docker.run()
            topResults.append(docker.topScores)
            ind += 1
            tf = time.time()
            print('That run took {} seconds'.format(int(tf-t0)))
            np.save('dockingResults', [swarms, worms, steps, topResults])
