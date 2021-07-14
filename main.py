# import statements
from opendna import *
from utils import *
from simtk.openmm.app import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

'''
script which accepts a DNA sequence and an analyte, and returns the binding affinity, and some related things

To-Do:
==> improve equilibration & sampling
==> finish README
    ==> example
==> write automated testing & test other modes
==> cluster auto-zipping script
==> 2d and docking collation and final printout
==> more detailed printouts
==> full trajectory combiner - for printouts
==> do something with trajectory 2D clustering output
==> always run a relaxation after MMB before going to full dynamics (incl. resize water box)

future features

==> identify pep-nuc hydrogen bonds
==> track post-complexation 2D trajectory and compare to pre-complexation, just like the 3D analysis
==> enhanced sampling and/or replica exchange (see openmmtools for implementation)
==> multi-state comparision in 2d and 3d w clustering
==> free energy and/or kd calculation
==> peptide restraints
==> would be nice to recognize local but significant rearrangements rather than just global
==> add nucleoside analytes
    => force field
    => docking
    => analysis & automation
==> implicit solvent - ambertools prmtop required

little things
==> getNucDATraj doesn't wor for terminal 'chi' angles
==> label bare exceptions
'''

params = {}
params['device'] = 'cluster'  # 'local' or 'cluster'
params['platform'] = 'CUDA'  # 'CUDA' or 'CPU'
params['platform precision'] = 'single'  # 'single' or 'double' only relevant on 'CUDA' platform

params['explicit run enumeration'] = False  # if True, the next run is fresh, in directory 'run%d'%run_num. If false, regular behaviour. Note: ONLY USE THIS FOR FRESH RUNS
if params['device'] == 'cluster':
    cmdLineInputs = get_input()
    params['run num'] = cmdLineInputs[0]  # option to get run num from command line (default zero)
    params['sequence'] = cmdLineInputs[1] # sequence specified from command line
    params['peptide'] = cmdLineInputs[2] # peptide specified from command line
elif params['device'] == 'local':
    params['run num'] = 0  # manual setting, for 0, do a fresh run, for != 0, pickup on a previous run.
    params['sequence'] = 'GCGCTTTGCGC' # manually set sequence # ATP aptamer
    params['peptide'] = 'YQTQTNSPRRAR' # manually set peptide

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
params['test mode'] = False # if true, changes params for a short simulation

# Pipeline parameters
params['secondary structure engine'] = 'NUPACK'  # 'NUPACK' or 'seqfold' - NUPACK has many more features and is the only package setup for probability analysis
params['equilibration time'] = 0.1  # initial equilibration time in nanoseconds
params['sampling time'] = 0.1  # sampling time in nanoseconds - in auto-sampling, this is the segment-length for each segment
params['auto sampling'] = False  # 'True' run sampling until RC's equilibrate, 'False' just run sampling for 'sampling time'
params['time step'] = 2.0  # MD time step in fs
params['print step'] = 1  # MD printout step in ps
params['max autoMD iterations'] = 1  # number of allowable iterations before giving up on auto-sampling - total max simulation length is this * sampling time
params['autoMD convergence cutoff'] = 1e-2  # how small should average of PCA slopes be to count as 'converged'
params['docking steps'] = 100  # number of steps for docking simulations
params['N 2D structures'] = 1 # max number of 2D structures to be considered (true number may be smaller depending on clustering)- the cost of this code is roughly linear in this integer
params['N docked structures'] = 1 # number of docked structures to output from the docker. If running binding, it will go this time (at linear cost)
params['fold speed'] = 'normal' # 'quick' 'normal' 'long' # time to spend first fold attempt - faster is cheaper but may not reach correct configuration, particularly for larger aptamers. 'normal' is default

if params['test mode']:
    params['fold speed'] = 'quick'

# physical params
params['pressure'] = 1  # atmospheres
params['temperature'] = 310  # Kelvin - used to predict secondary structure and for MD thermostatting
params['ionic strength'] = .163  # mmol - used to predict secondary structure and add ions to simulation box
params['pH'] = 7.4  # simulation will automatically protonate the peptide up to this pH

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
params['hydrogen mass'] = 1.0  # in amu - we can increase the time if we increase this value


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

# mmb control files
params['mmb params'] = 'lib/mmb/parameters.csv'
params['mmb normal template'] = 'lib/mmb/commands.template.dat'
params['mmb quick template'] = 'lib/mmb/commands.template_quick.dat'
params['mmb long template'] = 'lib/mmb/commands.template_long.dat'


# structure files
params['analyte pdb'] = 'lib/peptide/peptide.pdb'  # optional stasleep 120tic analyte - currently not used

'''
==============================================================
'''


if __name__ == '__main__':
    opendna = opendna(params)  # instantiate the class
    opendnaOutput = opendna.run() # retrive binding information (eventually this should become a normalized c-number)
