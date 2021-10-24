"""
from typing import Dict, Any, Union  # for the dictionary of params
...
params: Dict[str, Union[Union[str, int, bool, float], Any]] = {}  # set up as a dictionary
==> To replace the params = {}
"""
from opendna import *
# from utils import *
# from simtk.openmm.app import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

params = {}

params['device'] = 'local'  # 'local' or 'cluster'
if params['device'] == 'local':  # this is for running MMB. Consult the admin for cluster if running into problems.
    params['local device platform'] = 'macos'  # macos, or linux or windows
params['platform'] = 'CPU'  # 'CUDA' or 'CPU'
params['platform precision'] = 'single'  # 'single' or 'double'. Only relevant on 'CUDA' platform

if params['device'] == 'cluster':
    cmdLineInputs = get_input()  # get input arguments from command lines
    params['run num'] = cmdLineInputs[0]  # option to get run num from command line (default zero)
    params['sequence'] = cmdLineInputs[1]  # DNA aptamer sequence
    params['peptide'] = cmdLineInputs[2]  # target peptide sequence
    params['max walltime'] = cmdLineInputs[3]  # maximum walltime in hours
    # code will adapt maximum sampling steps to finish in less than this time,
    # or give up if this isn't enough time to complete even a minimum run.
elif params['device'] == 'local':
    params['run num'] = 0  # manual setting, for 0, do a fresh run, for != 0, pickup on a previous run
    params['sequence'] = 'AACGCTTTTTGCGTT'  # 'CCCGGGCCCGGG'  # manually set sequence # ATP aptamer
    params['peptide'] = 'YQTQTNSPRRAR'  # 'YRRYRRYRRY'  # 'YQTQTNSPRRAR' # manually set peptide # if no peptide, use ``False``
    params['max walltime'] = 3 * 24  # maximum walltime in hours

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
params['skip MMB'] = False  # if True, it will skip '2d analysis' and 'do MMB'
if params['skip MMB'] is True:
    params['folded initial structure'] = 'foldedSequence_0.pdb'  # if skipping MMB, must provide a folded structure
# if there are >1 folded structures or 2nd structures?

params['mode'] = 'free aptamer'  # 'full binding'  # 'full docking'  #'smooth dock'  #'coarse dock'  #'free aptamer'  # '3d smooth' # 'full binding'  # specify what to do
params['test mode'] = True
params['explicit run enumeration'] = False
params['implicit solvent'] = True  # implicit solvent or explicit solvent
if params['implicit solvent']:
    params['implicit solvent model'] = HCT  # only meaningful if implicit solvent is True
    params['leap template'] = 'leap_template.in'
    params['implicit solvent salt conc'] = 0.163  # molar/unit
    # TODO add more options to params: implicitSolventSaltConc, soluteDielectric, solventDielectric, implicitSolventKappa

# Pipeline parameters
params['secondary structure engine'] = 'NUPACK'  # 'NUPACK' or 'seqfold' - NUPACK has many more features and is the only package set up for probability analysis
params['N 2D structures'] = 1  # 2 # max number of 2D structures to be considered (true number may be smaller depending on clustering)- the cost of this code is roughly linear in this integer # TODO: could expand the candidate size and do something with them
params['fold speed'] = 'normal'  # 'quick', 'normal', 'long' - time to spend on first fold attempt - faster is cheaper but may not reach correct configuration, particularly for larger aptamers. 'normal' is default
params['equilibration time'] = 0.001  # 0.01 # initial equilibration time in nanoseconds
params['smoothing time'] = 0.05  # 1 # ns. MD relax after getting the initial 3D structure from user or MMB before sampling
params['sampling time'] = 0.5  # 1 # sampling time in nanoseconds - in auto-sampling, this is the segment-length for each segment
params['auto sampling'] = False  # 'True': run sampling till RC's equilibrate; 'False': just run sampling for 'sampling time'
params['time step'] = 2.0  # MD time step in fs
params['print step'] = 1  #10 # MD printout step in ps. ns > ps > fs
params['max aptamer sampling iterations'] = 20   # number of allowable iterations before giving on auto-sampling - total max simulation length = this * sampling time
params['max complex sampling iterations'] = 5  # number of iterations for the binding complex
params['autoMD convergence cutoff'] = 1e-2  # how small should average of PCA slopes be to count as 'converged' # TODO: where is the PCA used? to cluster conformations to obtain a representive one? # TODO: another clustering methods
params['docking steps'] = 200  # number of steps for docking simulations
params['N docked structures'] = 1  # 2 # number of docked structures to output from the docker. If running binding, it will go this time (at linear cost) # TODO: "it will go this time"?

if params['test mode']:  # shortcut for debugging
    params['N 2D structures'] = 1  # the clustering algorithm will stop when there are two structures left???
    params['fold speed'] = 'quick'
    params['equilibration time'] = 0.001
    params['smoothing time'] = 0.0001
    params['sampling time'] = 0.0001
    params['auto sampling'] = False
    params['time step'] = 2.0
    params['print step'] = 0.1
    params['max aptamer sampling iterations'] = 2
    params['max complex sampling iterations'] = 2
    params['autoMD convergence cutoff'] = 1e-2
    params['docking steps'] = 10
    params['N docked structures'] = 1
    
# Physical params
params['pressure'] = 1  # atmosphere
params['temperature'] = 310  # Kevin - used to predict secondary structure and for MD thermostat
params['ionic strength'] = 0.163  # Molar - sodium concentration - used to predict secondary structure and add ions to simulation box, must be 1100 M > [Na] > 50 for nupack to run
# TODO: how about adding other ions? Expand the FF as well?
params['[Mg]'] = 0.05  # Molar - magnesium concentration: 0.2 M > [Mg] > 0 - ONLY applies to NuPack fold - Does NOT add Mg to MD simulations
params['pH'] = 7.4  # simulation will automatically protonate the peptide up to this pH. Used in OpenMM for waterBox

# OpenMM params
# Currently, if there is a chk file, automatically resume from the chk file. 
# User needs to specify the chk file here so that it'll be copied to the workdir
params['pick up from chk'] = False
if params['pick up from chk'] is True:
    params['chk file'] = 'foldedSequence_0_amb_processed_state.chk'  # or 'relaxedSequence_0_amb_processed_state.chk'
else:
    params['chk file'] = ""  # empty string
params['force field'] = 'AMBER'  # this does nothing. The force field is specified in __init__ of interfaces.py
params['water model'] = 'tip3p'  # 'tip3p' (runs on Amber 14), other explicit models are also easy to add
params['box offset'] = 1.0  # nanometers
params['barostat interval'] = 25  # NOT USED.
params['friction'] = 1.0  # 1/picoseconds: friction coefficient determines how strongly the system is coupled to the heat bath (OpenMM)
params['nonbonded method'] = PME  # Particle Mesh Ewald: efficient full electrostatics method for use with periodic boundary conditions
params['nonbonded cutoff'] = 1.0  # nanometers
params['ewald error tolerance'] = 5e-4
params['constraints'] = HBonds
params['rigid water'] = True  # By default, OpenMM makes water molecules completely rigid, constraining both their bond lengths and angles. If False, it's good to reduce integration step size to 0.5 fs
params['constraint tolerance'] = 1e-6  # What is this tolerance for? For constraint?
params['hydrogen mass'] = 1.5  # in a.m.u. - we can increase the sampling time if we use heavier hydrogen
params['peptide backbone constraint constant'] = 0  # 10000  # constraint on the peptide's dihedral angles. force constant k.

# Path
if params['device'] == 'local':
    params['workdir'] = '/Users/taoliu/PycharmProjects/myOpenDNA/forCluster/ImSolv/skipMMBLocal/runs'
    params['mmb dir'] = '/Users/taoliu/Desktop/software/Installer.3_0.OSX/lib'
    params['mmb'] = '/Users/taoliu/Desktop/software/Installer.3_0.OSX/bin/MMB'

    # lightdock python scripts: they will be copied to 'ld_scripts' in workdir
    # therefore all addresses are relative to the workdir
    params['ld setup path'] = 'ld_scripts/lightdock3_setup.py'
    params['ld run path'] = 'ld_scripts/lightdock3.py'
    params['lgd generate path'] = '../ld_scripts/lgd_generate_conformations.py'
    params['lgd cluster path'] = '../ld_scripts/lgd_cluster_bsas.py'
    params['lg ant path'] = 'ld_scripts/ant_thony.py'
    params['lgd rank path'] = 'ld_scripts/lgd_rank.py'
    params['lgd top path'] = 'ld_scripts/lgd_top.py'

elif params['device'] == 'cluster':
    params['workdir'] = '/home/taoliu/scratch/implicitSolvent/forOpenDNA/runs'  # specify your working directory here. No / at the end
    params['mmb dir'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64'
    params['mmb'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/MMB.2_14.Linux64'
    # need to tell OS where to find the library files. All MMB files are in the same direcotory.
    # In ~/.bash_profile: export DYLD_LIBRARY_PATH=<mmm dir above>

    # lightdock python scripts: they will be copied to 'ld_scripts' in workdir
    # thereofre all addresses are relative to the workdir
    params['ld setup path'] = 'python ld_scripts/lightdock3_setup.py'
    params['ld run path'] = 'python ld_scripts/lightdock3.py'
    params['lgd generate path'] = 'python ../ld_scripts/lgd_generate_conformations.py'
    params['lgd cluster path'] = 'python ../ld_scripts/lgd_cluster_bsas.py'
    params['lg ant path'] = 'python ld_scripts/ant_thony.py'  # ant? thony?
    params['lgd rank path'] = 'python ld_scripts/lgd_rank.py'
    params['lgd top path'] = 'python ld_scripts/lgd_top.py'

# MMB control files
params['mmb params'] = 'lib/mmb/parameters.csv'
params['mmb normal template'] = 'lib/mmb/commands.template.dat'
params['mmb quick template'] = 'lib/mmb/commands.template_quick.dat'  # fold speed: quick
params['mmb long template'] = 'lib/mmb/commands.template_long.dat'  # fold speed: slow (ie, long)

# structure files: peptide analyte (target)
params['analyte pdb'] = 'lib/peptide/peptide.pdb'  # optional analyte - currently not used

'''
==============================================================
'''
if __name__ == '__main__':
    opendna = opendna(params)  # instantiate the class
    opendnaOutput = opendna.run()  # retrieve binding information (eventually this should become a normalized c-number)
    # TODO: what is a c-number?
