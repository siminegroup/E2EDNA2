# import statements
from opendna import *
from utils import *
from simtk.openmm.app import *


'''
script which accepts a DNA sequence and an analyte, and returns the binding affinity

To-Do:
==> auto-equilibration - done, just need a choice of when to sample (easily set!)
==> auto-binding - cut off simulation if the analyte floats away
==> multi-state comparision in 2d and 3d w clustering
==> implicit solvent - ambertools prmtop required
==> outputs & analysis - reduce binding to c-number
==> to cluster
==> check on checkpointing
==> dcd file keeps not existing - problem with WSL, likely won't happen on cluster
==> multi-accuracy system
==> peptide restraints
'''


params = {}
params['device'] = 'local' # 'local' or 'cluster'
params['platform'] = 'CPU' # 'CUDA' or 'CPU'
params['platform precision'] = 'single' # 'single' or 'double' only relevant on 'CUDA' platform

params['explicit run enumeration'] = False # if True, the next run is fresh, in directory 'run%d'%run_num. If false, regular behaviour. Note: ONLY USE THIS FOR FRESH RUNS
if params['device'] == 'cluster':
    params['run num'] = get_input() # option to get run num from command line (default zero)
elif params['device'] == 'local':
    params['run num'] = 0 # manual setting, for 0, do a fresh run, for != 0, pickup on a previous run.

'''
Modes, in order of increasing cost
'2d structure': ssString, pair list and probability
'3d coarse': MMB output, stressed structure, no solvent
'3d smooth': MMB output with short MD relaxation
'coarse dock': best docking scores on coarse MMB structure
'smooth dock': best docking scores on smoothed MMB structure
'free aptamer': evaluate and find representative 3D aptamer structure
'full docking': 'free aptamer' + docking
'full binding': 'full docking' + binding
'''
params['mode'] = 'full binding' # what to do

#Pipeline parameters
params['secondary structure engine'] = 'NUPACK' # 'NUPACK' or 'seqfold' - NUPACK is generally better / more flexible - will become the default/only option
params['equilibration time'] = 0.01 # initial equilibration time in nanoseconds
params['sampling time'] = 0.01 # sampling time in nanoseconds - in auto-sampling, this is the segment-length for each segment
params['auto sampling'] = True # 'True' run sampling until RC's equilibrate, 'False' just run sampling for 'sampling time'
params['time step'] = 2.0 # MD time step in fs
params['print step'] = 1 # MD printout step in ps
params['max autoMD iterations'] = 5 # number of allowable iterations before giving up on auto-sampling - total max simulation length is this * sampling time
params['autoMD convergence cutoff'] = 1e-2 # how small should average of PCA slopes be to count as 'converged'
params['docking steps'] = 100 # number of steps for docking simulations
params['N docked structures'] = 5  # number of docked structures to output from the docker

#OpenMM Parameters
params['force field'] = 'AMBER' # this does nothing
params['water model'] = 'tip3p' # 'tip3p' (runs on amber 14), other explicit models easy to add
params['box offset'] = 1.0 # nanometers
params['barostat interval'] = 25
params['friction'] = 1.0 # 1/picosecond
params['nonbonded method'] = PME
params['nonbonded cutoff'] = 1.0 # nanometers
params['ewald error tolerance'] = 5e-4
params['constraints'] = HBonds
params['rigid water'] = True
params['constraint tolerance'] = 1e-6
params['hydrogen mass'] = 1.0 # in amu

# physical params
params['pressure'] = 1 # atmospheres
params['temperature'] = 310 # Kelvin - used to predict secondary structure and for MD thermostatting
params['ionic strength'] = .163 # mmol - used to predict secondary structure and add ions to simulation box
params['pH'] = 7.4 # simulation will automatically protonate the peptide up to this pH


# paths
if params['device'] == 'local':
    params['workdir'] = '/mnt/c/Users/mikem/Desktop/mmruns'
    params['mmb dir'] = '/mnt/c/Users/mikem/Desktop/software/Installer.2_14.Linux64'
    params['mmb'] = '/mnt/c/Users/mikem/Desktop/software/Installer.2_14.Linux64/MMB.2_14.Linux64'

elif params['device'] == 'cluster':
    params['workdir'] = '/home/kilgourm/scratch/mmruns' # specify your working directory here
    params['mmb dir'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/'#'C:/Users/mikem/Desktop/Installer.2_14.Windows/MMB.2_14.exe'
    params['mmb'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/MMB.2_14.Linux64'

# mmb control files
params['mmb params'] = 'lib/parameters.csv'
params['mmb template'] = 'lib/commands.template.dat'

# lightdock python scripts
params['ld setup path'] = 'lib/lightdock/lightdock3_setup.py'
params['ld run path'] = 'lib/lightdock/lightdock3.py'
params['lgd generate path'] = 'lib/lightdock/lgd_generate_conformations.py'
params['lgd cluster path'] = 'lib/lightdock/lgd_cluster_bsas.py'
params['lg ant path'] = 'lib/lightdock/ant_thony.py'
params['lgd rank path'] = 'lib/lightdock/lgd_rank.py'
params['lgd top path'] = 'lib/lightdock/lgd_top.py'

# structure files
params['analyte pdb'] = 'lib/peptide/peptide.pdb' # optional static analyte - currently not used

'''
==============================================================
'''

if __name__ == '__main__':
    sequence = 'CGCTTTGCG' # random little hairpin #'ACCTGGGGGAGTATTGCGGAGGAAGGT' #ATP binding aptamer
    peptide = 'NNSPRR'#'YQTQTNSPRRAR' or 'False' for DNA analysis only
    opendna = opendna(sequence, peptide, params) # instantiate the class
    opendnaOutput = opendna.run() # retrive binding information (eventually this should become a normalized c-number)
