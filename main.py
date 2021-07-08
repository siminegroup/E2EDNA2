# import statements
from opendna import *
from utils import *
from simtk.openmm.app import *

'''
script which accepts a DNA sequence and an analyte, and returns the binding affinity

Known Issues:
==> On WSL, after using os.rename on .dcd files, os.path and all functions can no longer find it
    -> it's a permissions problem
    -> replace rename with replace -- seems good

To-Do:
==> multi-state comparision in 2d and 3d w clustering
==> implicit solvent - ambertools prmtop required
==> to cluster
==> separate sampling from equilibration
==> rethink checkpointing
==> peptide restraints
==> finish README
    ==> example
==> 2d structure pairing is not exclusive
==> customize MMB temperature
==> add failure output to premature binding cutoff
==> write up .sh script for github
'''

params = {}
params['device'] = 'cluster'  # 'local' or 'cluster'
params['platform'] = 'GPU'  # 'CUDA' or 'CPU'
params['platform precision'] = 'single'  # 'single' or 'double' only relevant on 'CUDA' platform

params['explicit run enumeration'] = False  # if True, the next run is fresh, in directory 'run%d'%run_num. If false, regular behaviour. Note: ONLY USE THIS FOR FRESH RUNS
if params['device'] == 'cluster':
    cmdLineInputs = get_input()
    params['run num'] = cmdLineInputs[0]  # option to get run num from command line (default zero)
    sequence = cmdLineInputs[1] # sequence specified from command line
    peptide = cmdLineInputs[2] # peptide specified from command line
elif params['device'] == 'local':
    params['run num'] = 0  # manual setting, for 0, do a fresh run, for != 0, pickup on a previous run.
    sequence = 'XXX' # placeholder
    peptide = 'BBB' # placeholder

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

params['mode'] = 'full binding'  # what to do

# Pipeline parameters
params['secondary structure engine'] = 'seqfold'  # 'NUPACK' or 'seqfold' - NUPACK is generally better / more flexible - will become the default/only option
params['equilibration time'] = 0.01  # initial equilibration time in nanoseconds
params['sampling time'] = .1  # sampling time in nanoseconds - in auto-sampling, this is the segment-length for each segment
params['auto sampling'] = True  # 'True' run sampling until RC's equilibrate, 'False' just run sampling for 'sampling time'
params['time step'] = 2.0  # MD time step in fs
params['print step'] = 1  # MD printout step in ps
params['max autoMD iterations'] = 3  # number of allowable iterations before giving up on auto-sampling - total max simulation length is this * sampling time
params['autoMD convergence cutoff'] = 1e-2  # how small should average of PCA slopes be to count as 'converged'
params['docking steps'] = 100  # number of steps for docking simulations
params['N docked structures'] = 5  # number of docked structures to output from the docker

# OpenMM Parameters
params['force field'] = 'AMBER'  # this does nothing
params['water model'] = 'tip3p'  # 'tip3p' (runs on amber 14), other explicit models easy to add
params['box offset'] = 1.0  # nanometers
params['barostat interval'] = 25
params['friction'] = 1.0  # 1/picosecond
params['nonbonded method'] = PME
params['nonbonded cutoff'] = 1.0  # nanometers
params['ewald error tolerance'] = 5e-4
params['constraints'] = HBonds
params['rigid water'] = True
params['constraint tolerance'] = 1e-6
params['hydrogen mass'] = 1.0  # in amu

# physical params
params['pressure'] = 1  # atmospheres
params['temperature'] = 310  # Kelvin - used to predict secondary structure and for MD thermostatting
params['ionic strength'] = .163  # mmol - used to predict secondary structure and add ions to simulation box
params['pH'] = 7.4  # simulation will automatically protonate the peptide up to this pH

# paths
if params['device'] == 'local':
    params['workdir'] = '/mnt/c/Users/mikem/Desktop/mmruns'
    params['mmb dir'] = '/mnt/c/Users/mikem/Desktop/software/Installer.2_14.Linux64'
    params['mmb'] = '/mnt/c/Users/mikem/Desktop/software/Installer.2_14.Linux64/MMB.2_14.Linux64'
    # lightdock python scripts
    params['ld setup path'] = 'lightdock/lightdock3_setup.py'
    params['ld run path'] = 'lightdock/lightdock3.py'
    params['lgd generate path'] = 'lightdock/lgd_generate_conformations.py'
    params['lgd cluster path'] = 'lightdock/lgd_cluster_bsas.py'
    params['lg ant path'] = 'lightdock/ant_thony.py'
    params['lgd rank path'] = 'lightdock/lgd_rank.py'
    params['lgd top path'] = 'lightdock/lgd_top.py'

elif params['device'] == 'cluster':
    params['workdir'] = '/home/kilgourm/scratch/mmruns'  # specify your working directory here
    params['mmb dir'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64'  # 'C:/Users/mikem/Desktop/Installer.2_14.Windows/MMB.2_14.exe'
    params['mmb'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/MMB.2_14.Linux64'
    # lightdock python scripts
    params['ld setup path'] = 'python lightdock/lightdock3_setup.py'
    params['ld run path'] = 'python lightdock/lightdock3.py'
    params['lgd generate path'] = 'python lightdock/lgd_generate_conformations.py'
    params['lgd cluster path'] = 'python lightdock/lgd_cluster_bsas.py'
    params['lg ant path'] = 'python lightdock/ant_thony.py'
    params['lgd rank path'] = 'python lightdock/lgd_rank.py'
    params['lgd top path'] = 'python lightdock/lgd_top.py'

# mmb control files
params['mmb params'] = 'lib/mmb/parameters.csv'
params['mmb template'] = 'lib/mmb/commands.template.dat'

# structure files
params['analyte pdb'] = 'lib/peptide/peptide.pdb'  # optional static analyte - currently not used

'''
==============================================================
'''

if __name__ == '__main__':
    if sequence == 'XXX': # if we did not specify a sequence at command line, specify it here
        sequence = 'GCGCTTTTGCGC'  # random little hairpin #'ACCTGGGGGAGTATTGCGGAGGAAGGT' #ATP binding aptamer
    if peptide == 'BBB': # if we did not specify a peptide at command line, specify it here
        peptide = 'NNSPRR'  # 'YQTQTNSPRRAR' or 'False' for DNA analysis only
    opendna = opendna(sequence, peptide, params)  # instantiate the class
    opendnaOutput = opendna.run() # retrive binding information (eventually this should become a normalized c-number)
