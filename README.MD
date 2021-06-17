# E2EDNA - OpenMM Implementation

Michael Kilgour, Tao Liu, Lena Simine
mjakilgour gmail com

## INSTALLATION
* Get the repository
```
it clone git@github.com:InfluenceFunctional/OpenDNA.git WHAT_YOU_WANT_IT_TO_BE_CALLED
```
* Setup python environment with the following packages
```
biopython==1.78
matplotlib==3.3.2
MDAnalysis==1.0.0
numpy==1.19.2
pandas==1.1.3
ParmEd==3.2.0
Pillow==7.2.0
pyparsing==2.4.7
scipy==1.4.1
OpenMM=7.5.0
PeptideBuilder==1.1.0
scikit-learn

```
* Install appropriate MacroMoleculeBuilder for your system from [SimTK here](https://simtk.org/projects/rnatoolbox) - old and new versions have been tested and seem to work
* Create working directory where runs will take place
* Set paths in main.py. Note, there is a toggle for running on a 'local' machine vs 'cluster', with distinct paths.
```
* params['mmb'] --> set to the MacroMoleculeBuilder executable on your system
* params['workdir'] --> the working directory created in the previous step
* params['analyte'] --> path to the analyte structure in .pdb format
* various .key files --> parameters may be specific to your analyte if custom force-field parameters are required. Adjust keyfiles accordingly, leaving the keywords on top as-is.
```
	

# __ none of this below is current __
## RUNNING A JOB
Before running a job, ensure installation, analyte parameterization and simulation parameters have been set.
1. Set 'run num' in main.py to zero for a fresh run, which will create a new run directory in your workdir. If you want to pick-up on a previous run, put its number here. 
   - Using argparse, one can directly set this via command line from a submission script rather than editing main.py itself.
2. Set simulation parameters in main.py, and select the mode to run in as params['simulation type']. 'free aptamer' runs the aptamer by itself, 'binding' runs a representative structure of the aptamer complexed with the analyte, and 'analysis' redoes reaction coordinate analysis given pre-existing trajectories.
4. Set the DNA sequence as a FASTA string to the 'sequence' variable at the bottom of main.py.
5. Run the code, e.g., for run_num = 0
```
python main.py --run_num=0
```
Inputs and outputs will go to /workdir/run0. In this directory, /outfiles contains outputs for debugging the various processes within the script.

In general, one would run first a 'free aptamer' simulation to retrieve a representaitive 3D structure and confirm the secondary structure stability. Then, one would change the operating mode to 'binding' to assess aptamer-analyte binding.
6. Simulation parameters and outputs including energy traces are saved to 'e2ednaOutputs.npy', which can be reloaded as
   - outputs = np.load('e2ednaOutputs.npy',allow_pickle=True).item()
   

## EXAMPLE
Here is an example of executing this code on a minimal system.
### System
For speed and simplicity, we will look at a tiny hairpin, GCTTTGC, and asses it's binding with UTP(-4) in PBS.

### Setup
We will use our default settings, and track the distances between bases which are expected to pair (1-7 and 2-6) as our reaction coordinates. 
The specific atoms to follow in these bases must be determined manually by users.
```
# Simulation parameters
params['reaction coordinates'] = [[26,213],[53,185]]
params['force field'] = 'AMOEBA'
params['target'] = 'UTP-4' # 
params['minimization gradrms'] = 0.05
params['equilibration time'] = 0.01 # equilibration time in nanoseconds
params['sampling time'] = 0.01 # sampling time in nanoseconds
params['time step'] = 2.0 # in fs
params['num charges'] = 4 # number of positive counterions (Na+) to add to simulation box to counter the analyte.
params['NaCl concentration'] = 163 # in mmol
params['box offset'] = 5 # in angstroms 
params['print step'] = 0.1 # printout step in ps 
params['heavy hydrogen'] = False 
params['outisde secondary structure'] = False 

# analyte placement - only relevant if 'simulation type' is 'binding'
params['analyte position'] = 'random' # 'random' - will be at least 2 nm from the aptamer in one of 8 cubic directions, 'manual' places the aptamer center at the given coordinates relative to origin
if params['analyte position'] == 'manual':
    params['analyte coordinates'] = [10, 10, 10] # places analyte [X, Y, Z] angstroms from the origin

# AMOEBA parameters
params['vdW radius'] = 12 # van der Waals radius
params['polarization version'] = 'OPT4' # OPT3 or OPT4 - 4 is slower but more accurate
params['polar eps'] = 1e-5 # induced dipole convergence parameter default 1e-6
params['polar predict'] = True # whether to run with polar-predict or not
```

### Structure evaluation

### Binding evaluation