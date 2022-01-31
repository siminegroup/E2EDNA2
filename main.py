from opendna import *
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

params = {}
# ============================================= Specify your settings within this block for a local test ===================================================
params['device'] = 'local'           # 'local' or 'cluster'
params['device platform'] = 'macos'  # 'macos' or 'linux' or 'WSL' (Windows Subsystem for Linux). Not supporting pure Windows OS (due to NUPACK)
params['platform'] = 'CPU'           # 'CPU' or 'CUDA'
if params['platform'] == 'CUDA': params['platform precision'] = 'single'  # 'single' or 'double'
if params['device'] == 'local':
    params['workdir'] = '/path-to-E2EDNA2-main/localruns'                  # directory manually created to store all future jobs
    params['mmb dir'] = '/path-to-E2EDNA2-main/Installer.3_0.OSX/lib'      # path to MMB dylib files  # need to tell OS where to find the library files. All MMB files are in the same direcotory.
    params['mmb']     = '/path-to-E2EDNA2-main/Installer.3_0.OSX/bin/MMB'  # path to the MMB executable; Note here MMB is the name of the executable file
else:  # params['device'] == 'cluster':
    params['workdir'] = '/home/taoliu/scratch/runs'
    params['mmb dir'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64'
    params['mmb'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/MMB.2_14.Linux64'
params['explicit run enumeration'] = True  # To resume a previous run from .chk file, use ``False`` here
cmdLineInputs = get_input()  # get input arguments from command lines
params['run num'] = cmdLineInputs[0]
params['mode'] = cmdLineInputs[1]
params['aptamerSeq'] = cmdLineInputs[2]
params['target ligand'] = cmdLineInputs[3]
params['target ligand type'] = cmdLineInputs[4]
params['target sequence'] = cmdLineInputs[5]
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
    params['DNA force field'] = 'DNA.OL15'  # 'DNA.OL15' or 'DNA.bsc1'. For free aptamer MD sampling if implicit solvent.
                                            # If target is peptide or RNA, will also add "leaprc.protein.ff14SB" or "source leaprc.RNA.OL3" to the "leap.in"
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
if params['pick up from freeAptamerChk'] is True:
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
if params['pick up from complexChk'] is True:
    params['chk file'] = 'complex_1_2_processed_state.chk'
    if params['implicit solvent'] is False:  # in explicit solvent
        params['resumed structurePDB'] = 'complex_1_2_processed.pdb'  # everything before _state then + .pdb # Just to provide the structureName
    else:
        params['resumed structurePrmtop'] = 'complex_1_2_processed.top'  # not complete: resume from top and crd file in implicit solvent
        params['resumed structureInpcrd'] = 'complex_1_2_processed.crd'
else:
    params['chk file'] = ""

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
# LightDock: these python scripts will be copied to 'ld_scripts' in workdir, therefore all addresses are relative to the workdir
if params['device'] == 'local':        
    params['ld setup path'] = 'python ld_scripts/lightdock3_setup.py'
    params['ld run path'] = 'python ld_scripts/lightdock3.py'
    params['lgd generate path'] = 'python ../ld_scripts/lgd_generate_conformations.py'
    params['lgd cluster path'] = 'python ../ld_scripts/lgd_cluster_bsas.py'
    params['lg ant path'] = 'python ld_scripts/ant_thony.py'
    params['lgd rank path'] = 'python ld_scripts/lgd_rank.py'
    params['lgd top path'] = 'python ld_scripts/lgd_top.py'
elif params['device'] == 'cluster':    
    params['ld setup path'] = 'python ld_scripts/lightdock3_setup.py'
    params['ld run path'] = 'python ld_scripts/lightdock3.py'
    params['lgd generate path'] = 'python ../ld_scripts/lgd_generate_conformations.py'
    params['lgd cluster path'] = 'python ../ld_scripts/lgd_cluster_bsas.py'
    params['lg ant path'] = 'python ld_scripts/ant_thony.py'  # ant? thony?
    params['lgd rank path'] = 'python ld_scripts/lgd_rank.py'
    params['lgd top path'] = 'python ld_scripts/lgd_top.py'

# ================================ main function entrance ==============================
if __name__ == '__main__':
    opendna = opendna(params)  # instantiate the class
    opendnaOutput = opendna.run()  # retrieve binding information (eventually this should become a normalized c-number)    
    printRecord('\nPipeline successfully stopped.')
