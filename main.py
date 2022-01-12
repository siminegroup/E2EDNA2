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
params['platform'] = 'CPU'  # 'CUDA' or 'CPU'
params['platform precision'] = 'single'  # 'single' or 'double'. Only relevant on 'CUDA' platform

if params['device'] == 'cluster':
    params['local device platform'] = ""  # leave it blank or specify the cluster's OS. This is for running MMB. Consult the admin for cluster if running into problems.
    cmdLineInputs = get_input()  # get input arguments from command lines
                                 # remember to edit the method "get_input()" in utils.py!!!`
    params['run num'] = cmdLineInputs[0]   # option to get run num from command line (default zero)
    params['sequence'] = cmdLineInputs[1]  # DNA aptamer sequence
    params['peptide'] = cmdLineInputs[2]   # target peptide sequence
    params['max walltime'] = cmdLineInputs[3]  # maximum walltime in hours
                                               # After free aptamer MD, "sampling time" might be reduced to ensure bindingDynamics to finish before time out.
                                               # Or give up if this isn't enough time to complete even a minimum run.
    # Physical params
    params['temperature'] = cmdLineInputs[4]     # Kevin - used to predict secondary structure and for MD thermostat
    params['pH'] = cmdLineInputs[5]              # simulation will automatically protonate the peptide up to this pH. Used in OpenMM for waterBox
    params['ionicStrength'] = cmdLineInputs[6]   # Molar - sodium concentration - used to predict secondary structure and add ions to simulation box, must be 1100 M > [Na] > 50 for nupack to run
                                                 # TODO: how about adding other ions? Expand the FF as well?
    params['[Mg]'] = cmdLineInputs[7]            # Molar - magnesium concentration: 0.2 M > [Mg] > 0 - ONLY applies to NuPack fold - Does NOT add Mg to MD simulations
    params['impSolv'] = cmdLineInputs[8]         # this is a string, cannot be used as parameter in prmtop.createSys()!
    params['pressure'] = 1                       # atmosphere

elif params['device'] == 'local':
    params['local device platform'] = 'macos'  # macos, or linux or windows
    params['run num'] = 1  # manual setting, for 0, do a fresh run, for != 0, pickup on a previous run
    params['sequence'] = 'CCTGGGGGAGTATTGCGGAGGAAGG'  # AACGCTTTTTGCGTT # 'CCCGGGCCCGGG'  # manually set sequence # ATP aptamer
    params['peptide'] = False  # 'YQTQTNSPRRAR'  # 'YRRYRRYRRY'  # 'YQTQTNSPRRAR' # manually set peptide # if no peptide, use ``False``
    params['max walltime'] = 3 * 24  # maximum walltime in hours
    # Physical params
    params['pressure'] = 1  # atmosphere
    params['temperature'] = 298 # Kevin - used to predict secondary structure and for MD thermostat
    params['ionicStrength'] = 0.150  # Molar - sodium concentration - used to predict secondary structure and add ions to simulation box, must be 1100 M > [Na] > 50 for nupack to run    
    params['[Mg]'] = 0.005  # Molar - magnesium concentration: 0.2 M > [Mg] > 0 - ONLY applies to NuPack fold - Does NOT add Mg to OpenMM simulation: currently only monovalent ions are supported for explicit solvent.
    params['pH'] = 7.4  # simulation will automatically protonate the peptide up to this pH. Used in OpenMM for waterBox
    params['implicit solvent'] = False  # implicit solvent or explicit solvent
    if params['implicit solvent'] is True:
        params['impSolv'] = 'HCT'  # 'HCT', 'OBC1', 'OBC2', 'GBn' or 'GBn2'
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
params['test mode'] = False
params['explicit run enumeration'] = True  # To resume a previous run from .chk file, use "False" here

# Pipeline parameters
params['secondary structure engine'] = 'NUPACK'  # 'NUPACK' or 'seqfold' - NUPACK has many more features and is the only package set up for probability analysis
params['N 2D structures'] = 1  # 2 # max number of 2D structures to be considered (true number may be smaller depending on clustering)- the cost of this code is roughly linear in this integer # TODO: could expand the candidate size and do something with them
params['fold speed'] = 'normal'  # 'quick', 'normal', 'long' - time to spend on first fold attempt - faster is cheaper but may not reach correct configuration, particularly for larger aptamers. 'normal' is default
params['foldFidelity'] = 0.9  # if folding fidelity < this value, refold; unless the fold speed is 'quick'
params['equilibration time'] = 0.1  # 0.01 # initial equilibration time in nanoseconds
params['smoothing time'] = 1  # ns. MD relax after getting the initial 3D structure from user or MMB before sampling
# smoothing is skipped in this case: after skipping MMB and under "free aptamer" mode: no need to run separate smoothing.
params['skip smoothing'] = False

params['sampling time'] = 100  # sampling time in nanoseconds - in auto-sampling, this is the segment-length for each segment
params['auto sampling'] = False  # 'True': run sampling till RC's equilibrate; 'False': just run sampling for 'sampling time'
params['time step'] = 2.0  # MD time step in fs
params['print step'] = 10 # MD printout step in ps. ns > ps > fs
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

# OpenMM params
# Currently, if there is a chk file, automatically resume from the chk file. 
# User needs to specify the chk file here so that it'll be copied to the workdir
params['pick up from chk'] = False
if params['pick up from chk']:
    params['chk file'] = 'relaxedSequence_0_processed_state.chk'
    # CAUTIOUS: a .chk file created on CUDA platform cannot be run on a CPU platform.
    '''
    A checkpoint contains data that is highly specific to the Context from which it was created.
        It depends on the details of the System, the Platform being used, and the hardware and software
        of the computer it was created on.  If you try to load it on a computer with different hardware,
        or for a System that is different in any way, loading is likely to fail.  Checkpoints created
        with different versions of OpenMM are also often incompatible.  If a checkpoint cannot be loaded,
        that is signaled by throwing an exception.
    '''
    if params['implicit solvent'] is False:
        params['resumed structurePDB'] = 'relaxedSequence_0_processed.pdb'  # everything before _state then + .pdb # Just to provide the structureName
    else:
        params['resumed structurePrmtop'] = 'relaxedSequence_0_processed.top'
        params['resumed structureInpcrd'] = 'relaxedSequence_0_processed.crd'
else:
    params['chk file'] = ""  # empty string <==> not resume a sampling from .chk file

# The following three are only used in explicit solvent
params['force field'] = 'AMBER'  # this does nothing. The force field is specified in __init__ of interfaces.py
params['water model'] = 'tip3p'  # 'tip3p' (runs on Amber 14), other explicit models are also easy to add
params['box offset'] = 1.0  # nanometers

params['barostat interval'] = 25  # NOT USED.
params['friction'] = 1.0  # 1/picoseconds: friction coefficient determines how strongly the system is coupled to the heat bath (OpenMM)
# OpenMM parameters for either explicit or implicit solvent when createSystem()
params['nonbonded method'] = PME  # Particle Mesh Ewald: efficient full electrostatics method for use with periodic boundary conditions to calculate long-range interactions
                                    #  use PME for long-range electrostatics, cutoff for short-range interactions
params['nonbonded cutoff'] = 1.0  # nanometers
params['ewald error tolerance'] = 5e-4  # In implicit solvent: this is the error tolerance to use if nonbondedMethod is Ewald, PME, or LJPME; In explicit solvent, it's "**args": Arbitrary additional keyword arguments
params['constraints'] = HBonds
params['rigid water'] = True  # By default, OpenMM makes water molecules completely rigid, constraining both their bond lengths and angles. If False, it's good to reduce integration step size to 0.5 fs
params['constraint tolerance'] = 1e-6  # What is this tolerance for? For constraint?
params['hydrogen mass'] = 1.5  # in a.m.u. - we can increase the sampling time if we use heavier hydrogen
params['peptide backbone constraint constant'] = 0  # 10000  # constraint on the peptide's dihedral angles. force constant k.

# Specify implicit solvent model
if params['implicit solvent']:
    # Select an implicit solvent model
    if params['impSolv'] == 'HCT':
        params['implicit solvent model'] = HCT
    elif params['impSolv'] == 'OBC1':
        params['implicit solvent model'] = OBC1
    elif params['impSolv'] == 'OBC2':
        params['implicit solvent model'] = OBC2
    elif params['impSolv'] == 'GBn':
        params['implicit solvent model'] = GBn
    elif params['impSolv'] == 'GBn2':
        params['implicit solvent model'] = GBn2
    else:
        raise ValueError('Illegal choice of implicit solvent model. Currently supported: HCT, OBC1, OBC2, GBn or GBn2')
        # sys.exit()

    params['nonbonded method'] = CutoffNonPeriodic  # or NoCutoff
    ''' When building a system in implicit solvent, there's no periodic boundary condition; so cannot use Ewald, PME, or LJPME => either NoCutoff or CutoffNonPeriodic.
    Question: what about params['ewald error tolerance']???? It doesn't matter or will cause conflict? Check source code.
        But can still specify cutoff for electrostatic interactions => params['nonbonded cutoff'] still works
        More on PBC: Periodic boundary conditions are used to avoid surface effects in explicit solvent simulations. Since implicit solvent simulations do not have solvent boundaries (the continuum goes on forever), 
        there is rarely a reason to use periodic boundary conditions in implicit solvent simulations.'''
    params['leap template'] = 'leap_template.in'
    params['implicit solvent salt conc'] = params['ionicStrength']  # molar/unit. Salt conc in solvent: converted to debye length (kappa) using the provided temperature and solventDielectric (see interfaces.py for detailed walk-through)
                                                    # need to be user input
    params['implicit solvent Kappa'] = None  # Debye screening parameter (reciprocal of Debye_length). 1/angstroms. If it's specfied, will ignore "salt conc". 
    params['soluteDielectric'] = 1.0  # default value = 1.0, meaning no screening effect at all.
    params['solventDielectric'] = 78.5  # water: 78.5; This offers a way to model non-aqueous system.

# Path
if params['device'] == 'local':
    params['workdir'] = '/Users/taoliu/PycharmProjects/myOpenDNA/TestImpSolv/mcKeague/runs'
    params['mmb dir'] = '/Users/taoliu/Desktop/software/Installer.3_0.OSX/lib'
    params['mmb']     = '/Users/taoliu/Desktop/software/Installer.3_0.OSX/bin/MMB'  # MMB executable

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
    params['workdir'] = '/home/taoliu/scratch/mcKeague/implicitSolvent/runs'
    params['mmb dir'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64'
    params['mmb'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/MMB.2_14.Linux64'
    # need to tell OS where to find the library files. All MMB files are in the same direcotory.

    # lightdock python scripts: they will be copied to 'ld_scripts' in workdir
    # thereofre all addresses are relative to the workdir
    params['ld setup path'] = 'python ld_scripts/lightdock3_setup.py'
    params['ld run path'] = 'python ld_scripts/lightdock3.py'
    params['lgd generate path'] = 'python ../ld_scripts/lgd_generate_conformations.py'
    params['lgd cluster path'] = 'python ../ld_scripts/lgd_cluster_bsas.py'
    params['lg ant path'] = 'python ld_scripts/ant_thony.py'  # ant? thony?
    params['lgd rank path'] = 'python ld_scripts/lgd_rank.py'
    params['lgd top path'] = 'python ld_scripts/lgd_top.py'

# MMB control files: be sure to download them from this GitHub repo together witht code.
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
