# import statements
from opendna import *
from utils import *
from simtk.openmm.app import *


'''
script which accepts a DNA sequence and an analyte, and returns the binding affinity

To-Do:
==> auto-equilibration - done, just need a choice of when to sample (easily set!)
==> multi-state comparision in 2d and 3d w clustering
==> implicit solvent - ambertools prmtop required
==> outputs & analysis - reduce binding to c-number
==> to cluster
'''


params = {}
params['device'] = 'local' # 'local' or 'cluster'
params['platform'] = 'CPU' # 'CUDA' or 'CPU'
params['platform precision'] = 'single' # only relevant for 'CUDA'

if params['device'] == 'cluster':
    params['run num'] = get_input()
elif params['device'] == 'local':
    params['run num'] = 0 # manual setting, for 0, do a fresh run, for != 0, pickup on a previous run.

# Simulation parameters
params['secondary structure engine'] = 'NUPACK' # 'NUPACK' or 'seqfold' - NUPACK is generally better / more flexible
params['force field'] = 'AMBER' # this does nothing
params['water model'] = 'tip3p' # 'tip3p' (runs on amber 14), other explicit models easy to add
params['equilibration time'] = 0.01 # initial equilibration time in nanoseconds
params['sampling time'] = 0.01 # sampling time in nanoseconds - in auto-sampling, this is the segment-length for each segment
params['auto sampling'] = True # NON FUNCTIONAL 'True' run sampling until RC's equilibrate + 'sampling time', 'False' just run sampling for 'sampling time'
params['time step'] = 2.0 # in fs
params['print step'] = 1 # printout step in ps

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
params['N docked structures'] = 5  # number of docked structures to output from the docker

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
    params['setup path'] = '/home/mkilgour/miniconda3/bin/lightdock3_setup.py'
    params['lightdock path'] = '/home/mkilgour/miniconda3/bin/lightdock3.py'
    params['lgd generate path'] = '/home/mkilgour/miniconda3/bin/lgd_generate_conformations.py'
    params['lgd cluster path'] = '/home/mkilgour/miniconda3/bin/lgd_cluster_bsas.py'
    params['anthony py'] = '/home/mkilgour/miniconda3/bin/ant_thony.py'
    params['lgd rank'] = '/home/mkilgour/miniconda3/bin/lgd_rank.py'
    params['lgd top'] = '/home/mkilgour/miniconda3/bin/lgd_top.py'

elif params['device'] == 'cluster':
    params['workdir'] = '/home/kilgourm/scratch/mmruns' # specify your working directory here
    params['mmb dir'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/'#'C:/Users/mikem/Desktop/Installer.2_14.Windows/MMB.2_14.exe'
    params['mmb'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/MMB.2_14.Linux64'
    params['mmb params'] = 'lib/parameters.csv'
    params['mmb template'] = 'lib/commands.template.dat'
    params['setup path'] = '~/lightDock_scripts/lightdock3_setup.py'
    params['lightdock path'] = '~/lightDock_scripts/lightdock3.py'
    params['lgd generate path'] = '~/lightDock_scripts/lgd_generate_conformations.py'
    params['lgd cluster path'] = '~/lightDock_scripts/lgd_cluster_bsas.py'
    params['anthony py'] = '~/lightDock_scripts/ant_thony.py'
    params['lgd rank'] = '~/lightDock_scripts/lgd_rank.py'
    params['lgd top'] = '~/lightDock_scripts/lgd_top.py'

# structure files
params['analyte pdb'] = 'lib/peptide/peptide.pdb' # optional static analyte - currently not used


'''
==============================================================
'''

if __name__ == '__main__':
    sequence = 'GCCCTTTCGGA' #'ACCTGGGGGAGTATTGCGGAGGAAGGT' #ATP binding aptamer
    peptide = 'NNSPRR'#'YQTQTNSPRRAR'
    opendna = opendna(sequence,peptide, params)
    opendnaOutput = opendna.run() # retrieve binding score and center-of-mass time-series
    np.save('opendnaOutput',opendnaOutput) # save outputs
