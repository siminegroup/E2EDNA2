# import statements
from opendna import *
from utils import *
from simtk.openmm.app import *


'''
script which accepts a DNA sequence and an analyte, and returns the binding affinity

To-Do:
==> testing
==> docking
==> add back PCA to 2nd structure, instead of clustering
==> auto-equilibration
    ==> need to implement checkpointing
==> multi-state comparision in 2d and 3d
==> implicit solvent - amber10 codes don't agree
==> sub seqfold for NUPACK
==> add docker
'''


params = {}
params['device'] = 'local' # 'local' or 'cluster'
params['platform'] = 'CPU' # 'CUDA' or 'CPU'
params['platform precision'] = 'single'

if params['device'] == 'cluster':
    params['run num'] = get_input()
elif params['device'] == 'local':
    params['run num'] = 0 # manual setting, for 0, do a fresh run, for != 0, pickup on a previous run.

# Simulation parameters
params['force field'] = 'AMBER' # this does nothing
params['water model'] = 'tip3p' # 'tip3p' (runs on amber 14), 'implicit' (runs on amber 10 - not working)
params['equilibration time'] = 0.001 # equilibration time in nanoseconds
params['sampling time'] = 0.01 # sampling time in nanoseconds
params['auto sampling'] = False # 'True' run sampling until RC's equilibrate + 'sampling time', 'False' just run sampling for 'sampling time'
params['time step'] = 3.0 # in fs
params['print step'] = 1 # printout step in ps - want there to be more than 2 and less than 100 total frames in any trajectory

params['box offset'] = 1.0 # nanometers
params['barostat interval'] = 25
params['friction'] = 1.0 # 1/picosecond
params['nonbonded method'] = PME
params['nonbonded cutoff'] = 1.0 # nanometers
params['ewald error tolerance'] = 5e-4
params['constraints'] = HBonds
params['rigid water'] = True
params['constraint tolerance'] = 1e-6
params['hydrogen mass'] = 4.0 # in amu

# physical params
params['pressure'] = 1 # atmospheres
params['temperature'] = 310 # Kelvin
params['ionic strength'] = .163 # mmol
params['pH'] = 7.4 # simulation will automatically protonate the peptide up to this pH

# paths
if params['device'] == 'local':
    params['workdir'] = '/mnt/c/Users/mikem/Desktop/mmruns' #'C:/Users\mikem\Desktop/mmruns'
    params['mmb dir'] = '/mnt/c/Users/mikem/Desktop/software/Installer.2_14.Linux64'#'C:/Users/mikem/Desktop/Installer.2_14.Windows/MMB.2_14.exe'
    params['mmb'] = '/mnt/c/Users/mikem/Desktop/software/Installer.2_14.Linux64/MMB.2_14.Linux64'#'C:/Users/mikem/Desktop/Installer.2_14.Windows/MMB.2_14.exe'
    params['mmb params'] = 'lib/parameters.csv'
    params['mmb template'] = 'lib/commands.template.dat'
elif params['device'] == 'cluster':
    params['workdir'] = '/home/kilgourm/scratch/mmruns' # specify your working directory here
    params['mmb dir'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/'#'C:/Users/mikem/Desktop/Installer.2_14.Windows/MMB.2_14.exe'
    params['mmb'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/MMB.2_14.Linux64'
    params['mmb params'] = 'lib/parameters.csv'
    params['mmb template'] = 'lib/commands.template.dat'

# structure files
params['analyte pdb'] = 'lib/peptide/peptide.pdb' # optional - currently not used


'''
==============================================================
'''

if __name__ == '__main__':
    sequence = 'CGTTTCG' #'ACCTGGGGGAGTATTGCGGAGGAAGGT' #ATP binding aptamer
    peptide = False #'YQT'#'YQTQTNSPRRAR'
    opendna = opendna(sequence,peptide, params)
    opendnaOutput = opendna.run() # retrieve binding score and center-of-mass time-series

    #os.chdir('C:/Users\mikem\Desktop/tinkerruns\clusterTests/fullRuns/run36')
    #comProfile, bindingScore = evaluateBinding('complex_sampling.arc') # normally done inside the run loop
    '''
    timeEq, potEq, kinEq = getEnergy()
    timeSa, potSa, kinSa = getEnergy()

    outputs = {}
    outputs['time equil'] = timeEq
    outputs['pot equil'] = potEq
    outputs['kin equil'] = kinEq
    outputs['time sampling'] = timeSa
    outputs['pot sampling'] = potSa
    outputs['kin sampling'] = kinSa
    outputs['binding score'] = bindingScore
    outputs['com profile'] = comProfile
    outputs['params'] = params
    np.save('bindingOutputs', outputs) # unpack with load then outputs.item()
    '''

