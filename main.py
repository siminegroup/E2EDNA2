#!/usr/bin/env python

'''
E2EDNA 2.0 - OpenMM Implementation of E2EDNA !

An automated pipeline for simulating DNA aptamers complexed with target ligands (peptide, DNA, RNA or small molecules).

Michael Kilgour, Tao Liu, Ilya S. Dementyev, Lena Simine

Copyright 2022 Michael Kilgour, Tao Liu and contributors of E2EDNA 2.0

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
'''
import argparse
import glob
from opendna import *
import os
import shutil
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


def get_args():
    ''' Parse commandline arguments '''

    parser = argparse.ArgumentParser(
        description='E2EDNA: Simulate DNA aptamers complexed with target ligands',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    run_info = parser.add_argument_group('Run Parameters')
    system_info = parser.add_argument_group('System Configuration')
    paths = parser.add_argument_group('Directory Settings')

    run_info.add_argument('-o',
                             '--outdir',
                             metavar='DIR',
                             type=int,
                             required=True,
                             help='Output directory for run created in --workdir')
    run_info.add_argument('-m',
                             '--mode',
                             metavar='MODE',
                             type=str,
                             required=True,
                             help='Run mode')
    run_info.add_argument('-a',
                             '--aptamer',
                             metavar='SEQ',
                             type=str,
                             required=True,
                             help="DNA Aptamer sequence (5'-3')")
    run_info.add_argument('-l',
                             '--ligand',
                             metavar='PDB',
                             type=str,
                             required=True,
                             help='Filename of target ligand')
    run_info.add_argument('-t',
                             '--ligand_type',
                             metavar='TYPE',
                             type=str,
                             default='',
                             help='Target ligand molecule type',
                             choices=['peptide', 'DNA', 'RNA', 'other', ''])
    run_info.add_argument('-s',
                             '--ligand_seq',
                             metavar='SEQ',
                             type=str,
                             default='',
                             help='Target ligand sequence if peptide, DNA, or RNA')
    run_info.add_argument('-f',
                             '--force',
                             action='store_true',
                             help='Overwrite existing --outdir')
    

    system_info.add_argument('-d',
                             '--device',
                             metavar='RUN',
                             type=str,
                             default='local',
                             help='Device configuration',
                             choices=['local', 'cluster'])
    system_info.add_argument('-os',
                             '--system',
                             metavar='OS',
                             type=str,
                             default='macos',
                             help='Operating system',
                             choices=['macos', 'linux', 'WSL'])
    system_info.add_argument('-p',
                             '--platform',
                             metavar='DEV',
                             type=str,
                             default='CPU',
                             help='Processing platform',
                             choices=['CPU', 'CUDA'])

    paths.add_argument('-w',
                          '--workdir',
                          metavar='DIR',
                          type=str,
                          default='./localruns',
                          help='Output directory to store runs')
    paths.add_argument('-md',
                          '--mmb_dir',
                          metavar='DIR',
                          type=str,
                          default='Installer-*/lib',
                          help='MMB library directory')
    paths.add_argument('-mb',
                          '--mmb',
                          metavar='MMB',
                          type=str,
                          default='Installer-*/bin/MMB',
                          help='Path to MMB executable')

    args = parser.parse_args()

    # Check if workdir exists, and create if not
    if not os.path.isdir(args.workdir):
        os.mkdir(args.workdir)

    # Check if output directory for run exists;
    # Either fail or force overwrite if so
    out_dir = os.path.join(args.workdir, 'run' + str(args.outdir))
    if os.path.isdir(out_dir):
        if args.force:
            shutil.rmtree(out_dir)
        else:
            parser.error(f'--outdir already exists at {out_dir}\n'
                        '\t use -f/--force to overwrite')

    # Validate run mode
    valid_modes = ['2d structure', '3d coarse', '3d smooth',
                   'coarse dock', 'smooth dock', 'free aptamer',
                   'full dock', 'full binding']

    if args.mode not in valid_modes:
        parser.error(f'Invalid --mode "{args.mode}". Must be one of: ' +
                     ', '.join(valid_modes) + '.')

    # Read in aptamer sequence if it is in a file
    if os.path.isfile(args.aptamer):
        args.aptamer = open(args.aptamer).read().rstrip()
    
    # If ligand is not 'False', it should be a pdb file
    if not args.ligand == 'False' and not os.path.isfile(args.ligand):
        parser.error(f'Invalid --ligand "{args.ligand}".'
                     'Must be either "False" or a valid PDB file.')

    # If ligand is 'False', these arguments are not used
    if args.ligand == 'False' and args.ligand_type:
        parser.error('--ligand_type is not used if --ligand is "False".')
    if args.ligand == 'False' and args.ligand_type:
        parser.error('--ligand_seq is not used if --ligand is "False".')

    # If ligand is peptide, DNA, or RNA, sequence is required
    if args.ligand_type in ['peptide', 'DNA', 'RNA'] and not args.ligand_seq:
        parser.error(f'--ligand_seq is required if --ligand_type is {args.ligand_type}')

    # Read in ligand sequence if it is in a file
    if os.path.isfile(args.ligand_seq):
        args.ligand_seq = open(args.ligand_seq).read().rstrip()

    # Try to find MMB library directory
    if not glob.glob(args.mmb_dir):
        parser.error(f'MMB library could not be found at {args.mmb_dir}.')
    else:
        found_path = glob.glob(args.mmb_dir)[0]
        if os.getcwd() in found_path: # Already absolute path
            args.mmb_dir = found_path
        else: # Need to make absolute path
            args.mmb_dir = os.path.join(os.getcwd(), found_path)

    # Try to find MMB executable
    if not glob.glob(args.mmb):
        parser.error(f'MMB executable could not be found at {args.mmb}.')
    else:
        found_path = glob.glob(args.mmb)[0]
        if os.getcwd() in found_path: # Already absolute path
            args.mmb = found_path
        else: # Need to make absolute path
            args.mmb = os.path.join(os.getcwd(), found_path)

    return args

args=get_args()
params = {}
# ============================================= Specify your settings within this block for a local test ===================================================
params['device'] = args.device
params['device platform'] = args.system
params['platform'] = args.platform
if params['platform'] == 'CUDA': params['platform precision'] = 'single'  # 'single' or 'double'
if params['device'] == 'local':
    params['workdir'] = args.workdir                # directory manually created to store all future jobs
    params['mmb dir'] = args.mmb_dir      # path to MMB dylib files  # need to tell OS where to find the library files. All MMB files are in the same direcotory.
    params['mmb']     = args.mmb  # path to the MMB executable; MMB is the *name* of the executable here
else:  # params['device'] == 'cluster':
    params['workdir'] = '/home/taoliu/scratch/runs'
    params['mmb dir'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64'
    params['mmb'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/MMB.2_14.Linux64'
params['explicit run enumeration'] = True  # To resume a previous run from .chk file, use ``False`` here
# cmdLineInputs = get_input()  # get input arguments from command lines
params['run num'] = args.outdir
params['mode'] = args.mode
params['aptamerSeq'] = args.aptamer
params['target ligand'] = args.ligand
params['target ligand type'] = args.ligand_type
params['target sequence'] = args.ligand_seq
# params['run num'] = 1  # 0: auto-increase run-num for a fresh run; > 0 AND params['explicit run enumeration'] = True: fresh run; > 0 AND params['explicit run enumeration'] = False: pickup on a previous run;
# params['mode'] = 'full binding'  # specify simulation mode
# params['aptamerSeq'] = 'TAATGTTAATTG'  # manually set DNA aptamer sequence from 5' to 3'
# params['target ligand'] = 'YQTQ.pdb'  # pdb filename of the target ligand. If no target, use 'False', such as in 'free aptamer' mode.
# params['target ligand type'] = 'peptide'  # 'peptide' or 'DNA' or 'RNA' or 'other'; This is ignored if no target.
# params['target sequence'] = 'YQTQTNSPRRAR'  # empty string, unless target ligand has sequence.
params['example target pdb'] = 'lib/peptide/peptide.pdb'   # an example of target ligand: a peptide, used when no given target but want to do docking
params['example peptide sequence'] = 'YQTQTNSPRRAR'        # YQTQ.pdb
''' params['mode'] can be:
    '2d structure': ssString, pair list and probability
    '3d coarse': MMB output, stressed structure, no solvent
    '3d smooth': MMB output with short MD relaxation
    'coarse dock': best docking scores on coarse MMB structure
    'smooth dock': best docking scores on smoothed MMB structure
    'free aptamer': evaluate and find representative 3D aptamer structure
    'full dock': 'free aptamer' + docking
    'full binding': 'full dock' + binding
'''
# ================================== Modifying anything below this block is NOT necessary, especically for a test run =======================================


# **************************************************************  Default setting ****************************************************************************
params['test mode'] = True                       # a quick test: short MD sampling and docking (if any)
params['secondary structure engine'] = 'NUPACK'  # 'NUPACK' or 'seqfold' - NUPACK has many more features and is the only package set up for probability analysis
params['foldFidelity'] = 0.9                     # MMB: if folding fidelity < this value, refold; unless the fold speed is 'quick'

params['implicit solvent'] = False       # ``False``: explicit solvent
params['auto sampling'] = False          # ``False``: just run sampling for params['sampling time']; ``True``: run sampling till RC's equilibrate; 
params['skip MMB'] = False               # ``True``: it will skip '2d analysis' and 'do MMB', then proceed to MD smoothing or MD sampling or dock (if "coarse dock")
params['skip smoothing'] = True          # ``True``: from MMB folding straight to MD sampling. Eg, if skipping MMB in 'free aptamer' mode
# The following two "pick up" control: if True, will ignore "skip MMB" and "skip smoothing"
params['pick up from freeAptamerChk'] = False  # ``True``: resume a MD of free aptamer, must provide a .chk file below. Skip everything before.
params['pick up from complexChk'] = False      # ``True``: resume a MD of aptamer-ligand, must provide a .chk file below. Skip everything before. Not complete
# Can only pick up from freeAptamerChk or complexChk. It makes no sense to resume both free aptamer and aptamer-ligand dynamics

params['max walltime'] = 24          # hours
params['pressure'] = 1               # 1 standard atmosphere
params['temperature'] = 298          # Kevin: used to predict secondary structure and for MD thermostat
params['ionicStrength'] = 0.150      # Molar: sodium concentration - used to predict secondary structure and add ions to simulation box, must be 1100 M > [Na] > 50 for nupack to run    
params['[Mg]'] = 0.005               # Molar: magnesium concentration: 0.2 M > [Mg] > 0 - ONLY applies to NuPack fold - Does NOT add Mg to OpenMM simulation.
params['pH'] = 7.4                   # simulation will automatically protonate the target, such as a peptide, up to this pH. Used in OpenMM for waterBox

if params['implicit solvent'] is True: 
    params['impSolv'] = 'HCT'  # 'HCT', 'OBC1', 'OBC2', 'GBn' or 'GBn2'
    params['DNA force field'] = 'DNA.OL15'  # 'DNA.OL15' or 'DNA.bsc1': used by free aptamer MD sampling in implicit solvent.
                                            # For complex MD sampling: if target is peptide or RNA, will add "leaprc.protein.ff14SB" or "source leaprc.RNA.OL3" to leap input file.
params['build a peptide as target ligand'] = False
params['peptide backbone constraint constant'] = 0  # if target ligand is a peptide, we can choose to put constraint on the peptide's dihedral angles. force constant k. 
# Ask Ilya: What's its unit? Only works if the target ligand is a peptide?
# ************************************************************************************************************************************************************

if params['skip MMB'] is True: params['folded initial structure'] = 'foldedAptamer_0.pdb'  
    # if skipping MMB, must provide a folded structure
    # what if >1 folded structures or 2nd structures?

if params['test mode'] is True:          # shortcut for debugging or running the code for the first time (validation)
    params['N 2D structures'] = 1        # the clustering algorithm will stop when there are two structures left???
    params['fold speed'] = 'quick'      # 'quick'
    params['equilibration time'] = 0.0001  # ns. 50 steps
    params['smoothing time'] = 0.001       # ns. 500 steps
    params['sampling time'] = 0.002        # ns. 1000 steps
    params['time step'] = 2.0              # fs
    params['print step'] = 0.05            # ps. print out every 25 steps
    params['max aptamer sampling iterations'] = 2
    params['max complex sampling iterations'] = 2
    params['autoMD convergence cutoff'] = 1e-2
    params['docking steps'] = 10
    params['N docked structures'] = 1
else:
    params['N 2D structures'] = 1    # 2 # max number of 2D structures to be considered (true number may be smaller depending on clustering) - the cost of this code is roughly linear in this integer 
                                     # TODO: could expand the candidate size and do something with them
    params['fold speed'] = 'normal'  # 'quick', 'normal', 'long' - time to spend on first fold attempt - faster is cheaper but may not reach correct configuration, particularly for larger aptamers. 'normal' is default
    params['equilibration time'] = 0.1  # 0.01 # initial equilibration time in nanoseconds
    params['smoothing time'] = 1        # ns. MD relax after getting the initial 3D structure from user or MMB before sampling
    params['sampling time'] = 100       # sampling time in nanoseconds - in auto-sampling, this is the segment-length for each segment
    params['time step'] = 2.0           # MD time step in fs
    params['print step'] = 10           # MD printout step in ps. ns > ps > fs
    params['max aptamer sampling iterations'] = 20   # number of allowable iterations before giving on auto-sampling - total max simulation length = this * sampling time
    params['max complex sampling iterations'] = 5    # number of iterations for the binding complex
    params['autoMD convergence cutoff'] = 1e-2       # how small should average of PCA slopes be to count as 'converged' 
                                                     # TODO: where is the PCA used? to cluster conformations to obtain a representive one? 
                                                     # TODO: another clustering methods
    params['docking steps'] = 200      # number of steps for docking simulations
    params['N docked structures'] = 1  # 2 # number of docked structures to output from the docker. If running binding, it will go this time (at linear cost) # TODO: "it will go this time"?

# OpenMM params
# If wants to resume a simulation, user must specify the chk file, which will be copied to the workdir
if (params['pick up from freeAptamerChk'] is True) or (params['pick up from complexChk'] is True):
    params['chk file'] = 'relaxedAptamer_0_processed_state.chk'
    # CAUTIOUS: a .chk file created on CUDA platform cannot be run on a CPU platform.
    ''' A checkpoint contains data that is highly specific to the Context from which it was created.
        It depends on the details of the System, the Platform being used, and the hardware and software
        of the computer it was created on.  If you try to load it on a computer with different hardware,
        or for a System that is different in any way, loading is likely to fail.  Checkpoints created
        with different versions of OpenMM are also often incompatible.  If a checkpoint cannot be loaded,
        that is signaled by throwing an exception.
    '''
    if params['implicit solvent'] is False:  # in explicit solvent
        params['resumed structurePDB'] = 'relaxedAptamer_0_processed.pdb'  # everything before _state then + .pdb # Just to provide the structureName
    else:
        params['resumed structurePrmtop'] = 'relaxedAptamer_0_processed.top'  # not complete: resume from top and crd file in implicit solvent
        params['resumed structureInpcrd'] = 'relaxedAptamer_0_processed.crd'
else:
    params['chk file'] = ""  # empty string <==> not resume a sampling from .chk file
# Can only pick up from freeAptamerChk or complexChk. It makes no sense to resume both free aptamer and aptamer-ligand dynamics
# if params['pick up from complexChk'] is True:
#     params['chk file'] = 'complex_1_2_processed_state.chk'
#     if params['implicit solvent'] is False:  # in explicit solvent
#         params['resumed structurePDB'] = 'complex_1_2_processed.pdb'  # everything before _state then + .pdb # Just to provide the structureName
#     else:
#         params['resumed structurePrmtop'] = 'complex_1_2_processed.top'  # not complete: resume from top and crd file in implicit solvent
#         params['resumed structureInpcrd'] = 'complex_1_2_processed.crd'
# else:
#     params['chk file'] = ""

# For guidance on adjustments, check out: http://docs.openmm.org/latest/userguide/application.html
params['force field'] = 'amber14-all'    # other choices: amber14/protein.ff14SB, amber14/DNA.OL15, amber14/RNA.OL3, amber14/lipid17, etc.
params['water model'] = 'amber14/tip3p'  # other choices: amber14/tip3pfb, amber14/tip4pew, amber14/spce, etc.
''' For complete list of available force fields that are bundled with OpenMM: 
    http://docs.openmm.org/latest/userguide/application/02_running_sims.html?highlight=forcefield#force-fields    
    The force field is specified in __init__ of interfaces.py
'''
params['box offset'] = 1.0  # nanometers
# The three params above are only used in explicit solvent.
params['barostat interval'] = 25  # NOT USED for a constant volume simulation. Useful when using a barostat in a constant pressure simulation.
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

# Specify implicit solvent model
if params['implicit solvent'] is True: # Select an implicit solvent model
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
    params['implicit solvent Kappa'] = None  # Debye screening parameter (reciprocal of Debye_length). 1/angstroms. If it's specfied, will ignore "salt conc". 
    params['soluteDielectric'] = 1.0         # default value = 1.0, meaning no screening effect at all.
    params['solventDielectric'] = 78.5       # water: 78.5; This offers a way to model non-aqueous system.

# Path: no need to change the following paths.
# MMB: Be sure to download ``lib/mmb`` folder from the GitHub repo together witht code.
params['mmb params'] = 'lib/mmb/parameters.csv'                       # parameter for the MMB executable
params['mmb normal template'] = 'lib/mmb/commands.template.dat'       # MMB folding protocol template No.1
params['mmb quick template'] = 'lib/mmb/commands.template_quick.dat'  # MMB folding protocol template No.2
params['mmb long template'] = 'lib/mmb/commands.template_long.dat'    # MMB folding protocol template No.3
# LightDock: these python scripts are available in the virtual environment's "bin" directory. 
# Once the environment is activated, these scripts are available at any path.
params['ld setup']     = 'lightdock3_setup.py'
params['ld run']       = 'lightdock3.py'
params['lgd generate'] = 'lgd_generate_conformations.py'
params['lgd cluster']  = 'lgd_cluster_bsas.py'
params['lgd rank']     = 'lgd_rank.py'
params['lgd top']      = 'lgd_top.py'
# params['lg ant']       = 'ant_thony.py'

# ================================ main function entrance ==============================
if __name__ == '__main__':
    opendna = opendna(params)  # instantiate the class
    opendnaOutput = opendna.run()  # retrieve binding information (eventually this should become a normalized c-number)    
    printRecord('\nPipeline successfully stopped.')
