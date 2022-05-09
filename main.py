#!/usr/bin/env python

'''
E2EDNA 2.0 - OpenMM Implementation of E2EDNA !

An automated pipeline for simulating DNA aptamers complexed with target ligands (peptide, DNA, RNA or small molecules).

Copyright 2022 Michael Kilgour, Tao Liu and Lena Simine

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
import os
import platform
import shutil
from opendna import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)



parser = argparse.ArgumentParser(
    description='E2EDNA: Simulate DNA aptamers complexed with target ligands',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

compute_info = parser.add_argument_group('Compute Platform Configuration')
paths = parser.add_argument_group('Directory Settings')
run_info = parser.add_argument_group('Run Parameters')


# =========================== YAML configuration file ===========================
parser.add_argument('-yaml',
                         '--yaml_config',
                         metavar='\b', # now to show metavar in usage message
                         type=str,
                         default='simu_config.yaml',
                         required=True,
                         help='A YAML configuration file that can specify all the arguments.')
parser.add_argument('-f',
                         '--force',
                         action='store_true',
                         required=False,
                         help='Overwrite existing --run_num')

# =========================== Compute platform information ===========================
compute_info.add_argument('-d',
                         '--device',
                         metavar='\b',
                         type=str,
                         default='local',
                         help='Device configuration',
                         choices=['local', 'cluster'])
compute_info.add_argument('-os',
                         '--operating_system',
                         metavar='\b',
                         type=str,
                         default='macos',
                         help='Operating system',
                         choices=['macos', 'linux', 'WSL'])
compute_info.add_argument('-p',
                         '--platform',
                         metavar='\b',
                         type=str,
                         default='CPU',
                         help='Processing platform',
                         choices=['CPU', 'CUDA'])
compute_info.add_argument('--CUDA_precision', 
                         metavar='\b',
                         type=str,
                         default='single',
                         help='Precision of CUDA, if used',
                         choices=['single', 'double'])

# =========================== Directory information ===========================
paths.add_argument('-w',
                    '--workdir',
                    metavar='DIR',
                    type=str,
                    default='./localruns',
                    help='Working directory to store individual output runs')
paths.add_argument('-mbdir',
                    '--mmb_dir',
                    metavar='\b',
                    type=str,
                    default='Installer*/lib',
                    help='MMB library directory')
paths.add_argument('-mb',
                    '--mmb',
                    metavar='\b',
                    type=str,
                    default='Installer*/bin/MMB*',
                    help='Path to MMB executable')


# =========================== Parameters to run pipeline ===========================
# "required" is set to be False as default, so that it is possible to only set it up in yaml config file

run_info.add_argument('--test_mode',
                         metavar='\b',
                         type=str,
                         default='Yes',
                         help='Rapidly run through pipeline using default test parameters',
                         choices=['Yes','No'])
run_info.add_argument('--explicit_run_enumeration',
                         metavar='\b',
                         type=str,
                         default='Yes',
                         help="'Yes': start a new run_num; 'No': output will be written to the an existing run_num folder.",
                         choices=['Yes','No'])
run_info.add_argument('-r',
                         '--run_num',
                         metavar='\b',
                         type=int,
                         default=1,
                         help='Run number. Output will be written to {--workdir}/run{--run_num}/')


run_info.add_argument('-m',
                         '--mode',
                         metavar='\b',
                         type=str,
                         help='Run mode',
                         choices=['2d structure', '3d coarse', '3d smooth',
                                  'coarse dock', 'smooth dock', 'free aptamer',
                                  'full dock', 'full binding'])
''' --mode:
    '2d structure': ssString, pair list and probability
    '3d coarse': MMB output, stressed structure, no solvent
    '3d smooth': MMB output with short MD relaxation
    'coarse dock': best docking scores on coarse MMB structure
    'smooth dock': best docking scores on smoothed MMB structure
    'free aptamer': evaluate and find representative 3D aptamer structure
    'full dock': 'free aptamer' + docking
    'full binding': 'full dock' + binding
'''

run_info.add_argument('-a',
                         '--aptamer_seq',
                         metavar='\b',
                         type=str,                         
                         help="DNA Aptamer sequence (5'->3')")
run_info.add_argument('-l',
                         '--ligand',
                         metavar='\b',
                         type=str,
                         default=None,
                         help="Name of PDB file for ligand structure; None if not to have ligand")
run_info.add_argument('-lt',
                         '--ligand_type',
                         metavar='\b',
                         type=str,
                         default=None,
                         help='Type of ligand molecule',
                         choices=['peptide', 'DNA', 'RNA', 'other'])
run_info.add_argument('-ls',
                         '--ligand_seq',
                         metavar='\b',
                         type=str,
                         default=None,
                         help='Ligand sequence if peptide, DNA, or RNA')
# Examples:
run_info.add_argument('--example_target_pdb',
                         metavar='\b',
                         type=str,
                         default='examples/peptide_ligand.pdb',
                         help='An example peptide ligand included in E2EDNA package: used when wish to test docking')
run_info.add_argument('--example_peptide_seq',
                         metavar='\b',
                         type=str,
                         default='YQTQTNSPRRAR',
                         help='The sequence of the example peptide ligand')



# # NUPACK and MMB modules:
run_info.add_argument('--skip_MMB',
                         metavar='\b',
                         type=str,
                         default='No',
                         help="'Yes': skip both 2D structure analysis and MMB folding, and start with a known --init_structure",
                         choices=['Yes','No'])
run_info.add_argument('-init',
                         '--init_structure',
                         metavar='\b',
                         type=str,
                         default=None,
                         help='Name of PDB file if starting pipeline on a DNA aptamer with known structure')

run_info.add_argument('--secondary_structure_engine',
                         metavar='\b',
                         type=str,
                         default='NUPACK',
                         help='Pipeline module that is used to predict secondary structures')
run_info.add_argument('--N_2D_structures',
                         metavar='\b',
                         type=int,
                         default=1,
                         help='Number of predicted secondary structures')
run_info.add_argument('--Mg_conc',
                         metavar='\b',
                         type=float,
                         default=0.0,
                         help='Magnesium molar concentration used in NUPACK: [0, 0.2], default = 0.0')

run_info.add_argument('--fold_fidelity',
                         metavar='\b',
                         type=float,
                         default=0.9,
                         help="Refold in MMB if score < fidelity (default = 0.9) unless the fold speed is 'quick'")
run_info.add_argument('--fold_speed',
                         metavar='\b',
                         type=str,
                         default='normal',
                         help='MMB folding speed',
                         choices=['quick','normal','slow'])
run_info.add_argument('--mmb_normal_template',
                         metavar='\b',
                         type=str,
                         default='lib/mmb/commands.template.dat',
                         help='Path to MMB folding protocol of normal speed')
run_info.add_argument('--mmb_quick_template',
                         metavar='\b',
                         type=str,
                         default='lib/mmb/commands.template_quick.dat',
                         help='Path to MMB folding protocol of quick speed')
run_info.add_argument('--mmb_slow_template',
                         metavar='\b',
                         type=str,
                         default='lib/mmb/commands.template_long.dat',
                         help='Path to MMB folding protocol of slow speed')
run_info.add_argument('--mmb_params',
                         metavar='\b',
                         type=str,
                         default='lib/mmb/parameters.csv',
                         help='Path to parameter file bundled with MMB package')




# OpenMM module: 
# # resume a MD sampling
run_info.add_argument('-pk',
                         '--pickup',
                         metavar='\b',
                         type=str,
                         default='No',
                         help='Whether the run is to resume MD sampling of an unfinished run or an old run.',
                         choices=['Yes', 'No'])
run_info.add_argument('--pickup_from_freeAptamerChk',
                         metavar='\b',
                         type=str,
                         default='No',
                         help='Resume MD sampling of free aptamer: skip everything before',
                         choices=['Yes','No'])
run_info.add_argument('--pickup_from_complexChk',
                         metavar='\b',
                         type=str,
                         default='No',
                         help='Resume MD sampling of aptamer-ligand: skip everything before',
                         choices=['Yes','No']) # Not complete
run_info.add_argument('--chk_file',
                         metavar='\b',
                         type=str,
                         default=None,
                         help='Name of checkpoint file for resuming MD sampling (file is placed in the current run{--run_num}')
run_info.add_argument('--pickup_pdb',
                         metavar='\b',
                         type=str,
                         default=None,
                         help='PDB file (topology + coordinates) for resuming MD sampling in explicit solvent (file is placed in the current run{--run_num}')
# not complete: resume from top and crd file in implicit solvent
# In opendna.py:
# outputDict['free aptamer results {}'.format(self.i)] = self.freeAptamerDynamics(self.params['pickup_pdb'], implicitSolvent=self.params['implicit_solvent'])
run_info.add_argument('--pickup_top',
                         metavar='\b',
                         type=str,
                         default=None,
                         help='Topology file if to resume MD sampling in Amber implicit sovlent') # '.top' 
# run_info.add_argument('--pickup_crd',
#                          metavar='\b',
#                          type=str,
#                          default=None,
#                          help='Coordinate file if to resume MD sampling in Amber implicit sovlent') # '.crd'

# # 
run_info.add_argument('--pressure',
                         metavar='\b',
                         type=float,
                         default=1.0,
                         help='Pressure in the unit of atm: default = 1.0')
run_info.add_argument('--temperature',
                         metavar='\b',
                         type=float,
                         default=298.0,
                         help='Temperature in Kelvin: default = 298.0')
run_info.add_argument('--ionicStrength',
                         metavar='\b',
                         type=float,
                         default=0.1,
                         help='Sodium molar concentration (could be used by NUPACK and OpenMM): default = 0.1') # could be used in implicit solvent as well
run_info.add_argument('--pH',
                         metavar='\b',
                         type=float,
                         default=7.4,
                         help='Could be used by OpenMM: Default = 7.4')
# run_info.add_argument('--barostat_frequency',
#                          metavar='\b',
#                          type=int,
#                          default=25, # unit: number of steps
#                          help='Useful when using a barostat in a constant pressure simulation')


# # MD sampling
run_info.add_argument('--auto_sampling',
                         metavar='\b',
                         type=str,
                         default='No',
                         help="'Yes': run MD sampling till convergence",
                         choices=['Yes','No'])
run_info.add_argument('--autoMD_convergence_cutoff',
                         metavar='\b',
                         type=float,
                         default=1.0e-2,
                         help='Convergence cutoff in {--auto_sampling}: default = 1.0e-2')
run_info.add_argument('--max_aptamer_sampling_iter',
                         metavar='\b',
                         type=int,
                         default=20,
                         help='Max number of iterations for free aptamer MD sampling if --auto_sampling')
run_info.add_argument('--max_complex_sampling_iter',
                         metavar='\b',
                         type=int,
                         default=5,
                         help='Max number of iterations for aptamer-ligand MD sampling if --auto_sampling')
run_info.add_argument('--max_walltime',
                         metavar='\b',
                         type=float,
                         default=24.0,
                         help='Walltime in hours to check runtime: default = 24.0')

run_info.add_argument('--skip_smoothing',
                         metavar='\b',
                         type=str,
                         default='Yes',
                         help="Default is 'Yes': no short MD relaxation before MD sampling",
                         choices=['Yes','No'])
run_info.add_argument('--equilibration_time',
                         metavar='\b',
                         type=float,
                         default=0.1,
                         help='Equilibration time in nanoseconds after energy minimization and before MD sampling')
run_info.add_argument('--smoothing_time',
                         metavar='\b',
                         type=float,
                         default=None,
                         help='Time in nanoseconds for short MD relaxation before MD sampling, if any')
run_info.add_argument('--sampling_time',
                         metavar='\b',
                         type=float,
                         default=1.0,
                         help='MD sampling time in nanoseconds')
run_info.add_argument('--time_step',
                         metavar='\b',
                         type=float,
                         default=2.0,
                         help='time step in femtoseconds in MD sampling: default = 2.0')
run_info.add_argument('--print_step',
                         metavar='\b',
                         type=float,
                         default=10.0,
                         help='Printout step in picoseconds in MD sampling: default = 10.0')


run_info.add_argument('--force_field',
                         metavar='\b',
                         type=str,
                         default='amber14-all',
                         help='Force field used in OpenMM') 
                         # amber14/protein.ff14SB, amber14/DNA.OL15,amber14/RNA.OL3, amber14/lipid17, etc.
                         # http://docs.openmm.org/latest/userguide/application/02_running_sims.html#force-fields
run_info.add_argument('--hydrogen_mass',
                         metavar='\b',
                         type=float,
                         default=1.5,
                         help='Default = 1.5 and unit is amu')


run_info.add_argument('--water_model',
                         metavar='\b',
                         type=str,
                         default='amber14/tip3p',
                         help='Explicit water solvent model used in OpenMM') # amber14/tip3pfb,amber14/tip4pew,amber14/spce,etc.
run_info.add_argument('--box_offset',
                         metavar='\b',
                         type=float,
                         default=1.0,
                         help='Buffering offset in nanometers on explicit solvent box (if used): default = 1.0')
run_info.add_argument('--constraints',
                         metavar='\b',
                         type=str,
                         default='None',
                         help="Specify which bond angles and/or lengths should be implemented with constraints: default = 'None'",
                         choices=['None','HBonds','AllBonds','HAngles'])
run_info.add_argument('--constraint_tolerance',
                         metavar='\b',
                         type=float,
                         default=1.0e-6,
                         help='Distance tolerance for constraint in OpenMM integrator: default = 1.0e-6')
run_info.add_argument('--rigid_water',
                         metavar='\b',
                         type=str,
                         default='Yes',
                         help='Whether to make water molecules completely rigid (bond lengths and angles)',
                         choices=['Yes','No'])
run_info.add_argument('--nonbonded_method',
                         metavar='\b',
                         type=str,
                         default='NoCutoff',
                         help='Type of nonbonded interactions',
                         choices=['Ewald','PME','LJPME','CutoffPeriodic','CutoffNonPeriodic','NoCutoff'])
                        # In implicit solvent: only 'CutoffPeriodic', 'CutoffNonPeriodic' or 'NoCutoff'
run_info.add_argument('--nonbonded_cutoff',
                         metavar='\b',
                         type=float,
                         default=1.0,
                         help='The cutoff distance in nanometers to use for nonbonded interactions: default = 1.0')
run_info.add_argument('--ewald_error_tolerance',
                         metavar='\b',
                         type=float,
                         default=5.0e-4,
                         help="Error tolerance if --nonbonded_method is Ewald, PME, or LJPME: default = 5.0e-4")
run_info.add_argument('--friction',
                         metavar='\b',
                         type=float,
                         default=1.0,
                         help='Friction coefficient in unit of 1/ps, used in Langevin integrator: default = 1.0')


run_info.add_argument('--implicit_solvent',
                         metavar='\b',
                         type=str,
                         default='No',
                         help='Whether to use Amber implicit solvent',
                         choices=['Yes','No'])
run_info.add_argument('--impSolv',
                         metavar='\b',
                         type=str,
                         default=None,
                         help="Specify Amber GB implicit solvent model if needed",
                         choices=['HCT','OBC1','OBC2','GBn','GBn2'])
run_info.add_argument('--soluteDielectric',
                         metavar='\b',
                         type=float,
                         default=1.0,
                         help='The solute dielectric constant to use in the implicit solvent model: default = 1.0')
run_info.add_argument('--solventDielectric',
                         metavar='\b',
                         type=float,
                         default=78.5,
                         help='The solvent dielectric constant to use in the implicit solvent model: default = 78.5')
run_info.add_argument('--leap_template',
                         metavar='\b',
                         type=str,
                         default='leap_template.in',
                         help='A script for running LEap program in Ambertools21, provided in E2EDNA package.')
run_info.add_argument('--DNA_force_field',
                         metavar='\b',
                         type=str,
                         default='DNA.OL15',
                         help='Force field for DNA used by LEap program',
                         choices=['DNA.OL15','DNA.bsc1'])


# Docking module:
run_info.add_argument('--docking_steps',
                         metavar='\b',
                         type=int,
                         default=200,
                         help='Number of steps for docking simulations')
run_info.add_argument('--N_docked_structures',
                         metavar='\b',
                         type=int,
                         default=1,
                         help='Number of docked structures output from the docker')


# Binary arguments: 'Yes' or 'No': more readable
# Then convert these binary arguments into boolean after parse_args(). See below.

# numerical arguments: integer
# numerical arguments: float
# string arguments:

# Finished reading arguments from command line input
params = parser.parse_args()

# Convert 10 'Yes'/'No' binary arguments into boolean:
params.test_mode                  = params.test_mode=='Yes'
params.explicit_run_enumeration   = params.explicit_run_enumeration=='Yes'
params.skip_MMB                   = params.skip_MMB=='Yes'
params.pickup                     = params.pickup=='Yes'
params.pickup_from_freeAptamerChk = params.pickup_from_complexChk=='Yes'
params.pickup_from_complexChk     = params.pickup_from_complexChk=='Yes'
params.auto_sampling              = params.auto_sampling=='Yes'
params.skip_smoothing             = params.skip_smoothing=='Yes'
params.rigid_water                = params.rigid_water=='Yes'
params.implicit_solvent           = params.implicit_solvent=='Yes'

params = vars(params) # convert from namespace to dictionary

# load up configuraton yaml file
if params.yaml_config is not None:    
    yaml_path = Path(params.yaml_config)
    assert yaml_path.exists()
    assert yaml_path.suffix in {".yaml", ".yml"}
    with yaml_path.open("r") as f:
        yaml_config = yaml.safe_load(f) # return a dict


# Merge parameter settings from command line and config yaml
params.update(yaml_config) # use yaml file to overwrite command line inputs (especially the default values in params). 
# yaml file does not have to set up all arguments


# Check if workdir exists, and create if not
if not os.path.isdir(params['workdir']):
    os.mkdir(params['workdir'])

# Check if output directory for run exists;
# Either fail or force overwrite if so
out_dir = os.path.join(params['workdir'], 'run' + str(params['run_num']))
if os.path.isdir(out_dir):
    if params['force']:
        shutil.rmtree(out_dir) # Delete an entire directory tree
    elif params['pickup'] is False:
        parser.error(f"--run_num already exists at {out_dir}\n"
                    '\t use -f/--force to overwrite')
    # The third scenario: resume the run in a certain output directory (ie, pickup is True), do nothing.

# Read in aptamer sequence if it is in a file
if os.path.isfile(params['aptamer_seq']):
    params['aptamer_seq'] = open(params['aptamer_seq']).read().rstrip()

# If ligand is not None, it should be a pdb file
if (params['ligand'] is not None) and (not os.path.isfile(params['ligand'])):
    parser.error(f"Invalid --ligand '{params['ligand']}'. Must be either None or a valid PDB filename.")

# If ligand is None, these two arguments are not used: ligand_type and ligand_seq
if (params['ligand'] is None) and (params['ligand_type'] is not None):
    parser.error('--ligand_type should not used if --ligand is None.')
if (params['ligand'] is None) and (params['ligand_seq'] is not None):
    parser.error('--ligand_seq should not used if --ligand is None.')

# If ligand_type is peptide, DNA, or RNA, its sequence is required
if (params['ligand_type'] in ['peptide', 'DNA', 'RNA']) and (params['ligand_seq'] is None):
    parser.error(f"--ligand_seq is required if --ligand_type is {params['ligand_type']}")
# Read in ligand sequence if it is in a file
if os.path.isfile(params['ligand_seq']):
    params['ligand_seq'] = open(params['ligand_seq']).read().rstrip()

if params['skip_MMB'] is False:
    # Try to find MMB library directory
    if not glob.glob(params['mmb_dir']):
        parser.error(f"MMB library could not be found at {params['mmb_dir']}.")
    else:
        found_path = glob.glob(params['mmb_dir'])[0]
        if os.getcwd() in found_path: # Already absolute path
            params['mmb_dir'] = found_path
        else: # Need to make absolute path
            params['mmb_dir'] = os.path.join(os.getcwd(), found_path)
    # Try to find MMB executable
    if not glob.glob(params['mmb']):
        parser.error(f"MMB executable could not be found at {params['mmb']}.")
    else:
        found_path = glob.glob(params['mmb'])[0]
        if os.getcwd() in found_path: # Already absolute path
            params['mmb'] = found_path
        else: # Need to make absolute path
            params['mmb'] = os.path.join(os.getcwd(), found_path)
# If skip_MMB: make sure init_structure is provided
else: 
    if params['init_structure'] is None:
        parser.error(f"--init_structure is required if --skip_MMB is {params['skip_MMB']}")

if params['pickup'] is True:
    # only one MD step is resumed
    if (params['pickup_from_freeAptamerChk'] is True) and (params['pickup_from_complexChk'] is True):
        parser.error(f"--pickup_from_freeAptamerChk and --pickup_from_complexChk cannot be 'Yes' together if --pickup is 'Yes'.")
    if (params['pickup_from_freeAptamerChk'] is False) and (params['pickup_from_complexChk'] is False):
        parser.error(f"One of --pickup_from_freeAptamerChk and --pickup_from_complexChk should be 'Yes' if --pickup is 'Yes'.")
    # chk_file cannot be None
    if params['chk_file'] is None:
        parser.error(f"--chk_file is required if --pickup is 'Yes'.")
    # pdb or top file cannot be None
    if params['implicit_solvent'] is False:
        if params['pickup_pdb'] is None:
            parser.error(f"--pickup_pdb is required if --pickup is 'Yes' in explicit solvent.")
    else: # not complete: resume from top and crd file in implicit solvent
        if params['pickup_top'] is None:
            parser.error(f"--pickup_top is required if --pickup is 'Yes' in implicit solvent.")
# If not to resume, all pickup flags should be 'No' and no file should be provided
else:
    if (params['pickup_from_freeAptamerChk'] is True) or (params['pickup_from_complexChk'] is True):
        parser.error(f"--pickup_from_freeAptamerChk and --pickup_from_complexChk should both be 'No' if --pickup is 'No'.")
    
    if (params['chk_file'] is not None) or (params['pickup_pdb'] is not None) or (params['pickup_top'] is not None):
        parser.error(f"--chk_file, --pickup_pdb and --pickup_top should be None --pickup is 'No'.")

# params = params

# When OpenMM is required:
if params['mode'] in ['3d smooth','smooth dock','free aptamer','full dock','full binding']:

    # Convert some string variables to corresponding objects in OpenMM
    if params['constraints'] == 'HBonds':
        params['constraints'] = HBonds
    elif params['constraints'] == 'HAngles':
        params['constraints'] = HAngles
    elif params['constraints'] == 'AllBonds':
        params['constraints'] = AllBonds
    elif params['constraints'] == 'None':
        params['constraints'] = None

    if params['implicit_solvent'] is False:
        # In explicit solvent: more choices for nonbonded method
        if params['nonbonded_method'] == 'NoCutoff':
            params['nonbonded_method'] = NoCutoff
        elif params['nonbonded_method'] == 'CutoffNonPeriodic':
            params['nonbonded_method'] = CutoffNonPeriodic
        elif params['nonbonded_method'] == 'CutoffPeriodic':
            params['nonbonded_method'] = CutoffPeriodic
        elif params['nonbonded_method'] == 'LJPME':
            params['nonbonded_method'] = LJPME
        elif params['nonbonded_method'] == 'Ewald':
            params['nonbonded_method'] = Ewald
        elif params['nonbonded_method'] == 'PME':
            params['nonbonded_method'] = PME
        # Particle Mesh Ewald: efficient full electrostatics method for 
        # use with periodic boundary conditions to calculate long-range interactions
        else:
            raise ValueError("Illegal choice of nonbonded_method in explicit solvent model. "
                "Choose from ['Ewald','PME','LJPME','CutoffPeriodic','CutoffNonPeriodic','NoCutoff'].")
    # implicit solvent model
    else:
        params['implicit_solvent_salt_conc'] = params['ionicStrength']  
        '''Salt conc in solvent: converted to debye length (kappa) using the provided temperature and solventDielectric 
        (see interfaces.py for detailed walk-through)'''

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

        # In implicit solvent: only CutoffPeriodic or CutoffNonPeriodic or NoCutoff
        if params['nonbonded_method'] == 'NoCutoff':
            params['nonbonded_method'] = NoCutoff
        elif params['nonbonded_method'] == 'CutoffNonPeriodic':
            params['nonbonded_method'] = CutoffNonPeriodic
        elif params['nonbonded_method'] == 'CutoffPeriodic':
            params['nonbonded_method'] = CutoffPeriodic
        else:
            raise ValueError("Illegal choice of nonbonded_method in implicit solvent model. "
                "Choose 'CutoffPeriodic', 'CutoffNonPeriodic' or 'NoCutoff'.")
        ''' When building a system in implicit solvent, no periodic boundary condition => cannot use Ewald, PME, or LJPME 
                => only NoCutoff or CutoffNonPeriodic or CutoffPeriodic.        
            PBC is used to avoid surface effects in explicit solvent and implicit solvent does not have solvent boundaries (the continuum goes on forever).
        ''' 



# Default values for shortcut during debugging or testing code
if params['test_mode'] is True:
    params['skip_smoothing'] = True
    params['N_2D_structures'] = 1
    params['fold_speed'] = 'quick'
    params['equilibration_time'] = 0.0001  # ns. 50 steps
    params['sampling_time'] = 0.002        # ns. 1000 steps
    params['time_step'] = 2.0              # fs
    params['print_step'] = 0.05            # ps. print out every 25 steps
    params['max_aptamer_sampling_iter'] = 2
    params['max_complex_sampling_iter'] = 2
    params['autoMD_convergence_cutoff'] = 1e-2
    params['docking_steps'] = 10
    params['N_docked_structures'] = 1


# Default settings of lightdock
params['ld_setup']     = 'lightdock3_setup.py'
params['ld_run']       = 'lightdock3.py'
params['lgd_generate'] = 'lgd_generate_conformations.py'
params['lgd_cluster']  = 'lgd_cluster_bsas.py'
params['lgd_rank']     = 'lgd_rank.py'
params['lgd_top']      = 'lgd_top.py'


# ================================ main function entrance ==============================
if __name__ == '__main__':
    opendna = opendna(params)  # instantiate the class
    opendnaOutput = opendna.run()  # retrieve binding information (eventually this should become a normalized c-number)    
    printRecord('\nPipeline successfully stopped.')

