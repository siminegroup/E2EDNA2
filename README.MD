# E2EDNA - OpenMM Implementation - E2EDNA 2.0!

### A fully-automated code for DNA aptamer 2D and 3D structure analysis, and binding analysis with small molecules and short peptides.

**Michael Kilgour, Tao Liu, Ilya S. Dementyev, Lena Simine**

_mjakilgour gmail com_

For original version of E2EDNA: J. Chem. Inf. Model. 2021, 61, 9, 4139–4144 (https://doi.org/10.1021/acs.jcim.1c00696)

## Installation
* Get the repository
```
git clone git@github.com:InfluenceFunctional/OpenDNA.git WHAT_YOU_WANT_IT_TO_BE_CALLED
```
* Set up a Python environment with the following packages
  * Note to Simine Group Members running on Compute Canada clusters: skip this section and go directly to paths. Your virtual environment is automatically setup via the instructions in [``EXAMPLE_SUB_cluster.sh``](https://github.com/InfluenceFunctional/OpenDNA/blob/TL-cluster-code/EXAMPLE_SUB_cluster.sh) in the [``TL-cluster-code``](https://github.com/InfluenceFunctional/OpenDNA/tree/TL-cluster-code) branch. Modify this bash script accordingly for other types of computing clusters.
  * For local runs, follow the [``EXAMPLE_SUB.sh``](https://github.com/InfluenceFunctional/OpenDNA/blob/main/EXAMPLE_SUB.sh) in the [``main``](https://github.com/InfluenceFunctional/OpenDNA) branch to install required packages and their dependences.
  * In particular, as to the ``MacroMoleculeBuilder (MMB)`` module:
    * Download appropriate version for your system from [SimTK](https://simtk.org/projects/rnatoolbox): old and new versions have been tested and seem to work.
  * As shown in the [``requirements.txt``](https://github.com/InfluenceFunctional/OpenDNA/blob/main/requirements.txt), the required dependences are as follows. Use ``pip install -r requirements.txt`` to install them all in once:
```
biopython==1.78
certifi==2020.12.5
cycler==0.10.0
decorator==4.4.2
GridDataFormats==0.5.0
gsd==1.9.3
joblib==1.0.1
kiwisolver==1.3.1
matplotlib==3.3.2
mmtf-python==1.1.2
mock==4.0.3
msgpack==0.6.2
networkx==2.5
pandas==1.1.3
ParmEd==3.2.0
Pillow==7.2.0
pyparsing==2.4.7
python-dateutil==2.8.1
pytz==2021.1
scipy==1.4.1
six==1.15.0
tqdm==4.59.0
scikit_learn==0.24.2
Cython==0.29.23
mdtraj==1.9.5
```
* Create a working directory where runs will take place on your system.
  * Each new run will create a serially-numbered directory in this folder containing all the relevant outputs.
* If running on a cluster, customize your submission shell script (example: [``EXAMPLE_SUB_cluster.sh``](https://github.com/InfluenceFunctional/OpenDNA/blob/TL-cluster-code/EXAMPLE_SUB_cluster.sh)). 
* Specify the file paths in [``main.py``](https://github.com/InfluenceFunctional/OpenDNA/blob/main/main.py):
  * Note that there is a toggle for running on a ``local`` machine vs. a computing ``cluster``.
  * ``params['workdir']``: set to the working directory you created in the previous step.
  * ``params['mmb dir']``: set to the MMB directory for its library files. Usually a name like ``.../Installer.#_##.Linux(or OSX)/lib``, depending on the software version.
  * ``params['mmb']``: set to the MMB executable on your system. Usually a name like ``.../Installer.#_##.Linux(or OSX)/bin/MMB``, depending on the software version.
* Not necessary to change the following paths, especially doing a test run:
  * ``params['mmb params']``: path to the MMB parameter file in [``lib/mmb/``](https://github.com/InfluenceFunctional/OpenDNA/tree/main/lib/mmb) of this repository. Modifying this file is highly not recommended. **Be sure to downlaod the entire ``lib/`` together with code.**
  * ``params['mmb normal template']``, ``params['mmb quick template']`` and ``params['mmb long template']``: paths to 3 template protocols of running MMB (``normal``, ``quick``, and ``long`` modes) in [``lib/mmb/``](https://github.com/InfluenceFunctional/OpenDNA/tree/main/lib/mmb) of this repository. They do not need to be changed unless to create your own MMB folding protocol.
  * ``params['ld setup path']``, ``params['ld run path']``, etc: paths to various [``LightDock``](https://lightdock.org/) Python scripts, providing they are in a working directory. No need to change these paths. The paths are used to automatically copy the scripts to ``ld_scripts/`` within a working directory from [``lib/lightdock/``](https://github.com/InfluenceFunctional/OpenDNA/tree/main/lib/lightdock).
    * Alternatively, these scripts automatically come with an installation of ``pip install lightdock``. Depending on your installation, you may need to add ``python`` in the paths of these scripts, as shown in the paths for a ``cluster`` run.


## Modes
OpenDNA takes in a DNA aptamer sequence in FASTA format, and optionally a short peptide or other small molecule, and returns details of the aptamer structure and binding behaviour.
This code implements several distinct analysis modes so users may customize the level of computational cost and accuracy.

* `2d structure` &rarr; returns NUPACK or seqfold analysis of aptamer secondary structure. Very fast, O(<1s). If using NUPACK, includes probability of observing a certain fold and of suboptimal folds within kT of the minimum.
* `3d coarse` &rarr; returns MMB fold of the best secondary structure. Fast O(5-30 mins). Results in a strained 3D structure which obeys base pairing rules and certain stacking interactions.
* `3d smooth` &rarr; identical to '3d coarse', with a short MD relaxation in solvent. ~Less than double the cost of '3d coarse' depending on relaxation time.
* `coarse dock` &rarr; uses the 3D structure from '3d coarse' as the initial condition for a LightDock simulation, and returns best docking configurations and scores. Depending on docking parameters, adds O(5-30mins) to '3d coarse'.
* `smooth dock` &rarr; identical to 'coarse dock', instead using the relaxed structure from '3d smooth'. Similar cost. 
* `free aptamer` &rarr; fold the aptamer in MMB and run extended MD sampling to identify a representative, equilibrated 2D and 3D structure. Slow O(hours).
* `full dock` &rarr; Return best docking configurations and scores from a LightDock run using the fully-equilibrated aptamer structure 'free aptamer'. Similar cost (LightDock is relatively cheap)
* `full binding` &rarr; Same steps as 'full docking', with follow-up extended MD simulation of the best binding configuration. Slowest O(hours).


## Running a job

### Quickstart
* Set 'params' in [``main.py``](https://github.com/InfluenceFunctional/OpenDNA/blob/main/main.py), as indicated in "Installation".
* Run the bash script ``automate_tests.sh`` to test all 8 modes automatically.
* Alternatively, a single run can be carried out by run_num, mode, aptamer sequence, and ligand's structural file. For example,
```
python main.py --run_num=1 --mode='free aptamer' --aptamerSeq='TAATGTTAATTG' --ligand='YQYQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR'

#--ligand='False'        # if no ligand. --ligandType and --ligandSeq will be ignored.
#--ligandType='peptide'  # or 'DNA' or 'RNA' or 'other'. Assuming 'other' ligand can be described by Amber14 force field.
#--ligandSeq=''          # if no sequence. For instance, when ligandType is 'other'

```

# __ work in progress__

### Physical Parameters

Default force field is AMBER 14. Other AMBER fields and explicit water models are trivial to implement. Implicit water requires moving to building systems from AMBER prmtop files. CHARMM may also be easily implemented, but hasn't been tested. AMOEBA 2013 parameters do not include nucleic acids, and AMOEBABIO18 parameters are not implemented in OpenMM.  
```
* params['force field'] = 'AMBER'
* params['water model'] = 'tip3p'
```

Default parameters here - for guidance on adjustments start [here](http://docs.openmm.org/latest/userguide/application.html).

```
params['box offset'] = 1.0 # nanometers
params['barostat interval'] = 25
params['friction'] = 1.0 # 1/picosecond
params['nonbonded method'] = PME
params['nonbonded cutoff'] = 1.0 # nanometers
params['ewald error tolerance'] = 5e-4
params['constraints'] = HBonds
params['rigid water'] = True
params['constraint tolerance'] = 1e-6
params['pressure'] = 1 
```

Increasing hydrogen mass e.g., to 4 AMU enables longer time-steps up to ~3-4 fs. See documentation for details.

```
params['hydrogen mass'] = 1.0 # in amu
```

Temperature, pH and ionic strength are taken into account for 2D folding in NUPACK, ion concentration in MD simulation, and protonation of molecules for MD (safest near 7-7.4).

```
params['temperature'] = 310 # Kelvin - used to predict secondary structure and for MD thermostatting
params['ionic strength'] = .163 # mmol - used to predict secondary structure and add ions to simulation box
params['pH'] = 7.4 # simulation will automatically protonate the peptide up to this pH
```

The peptide backbone constraint constant is the constant used to constrain backbone dihedrals. 
A minimum of 10000, as it is currently set, is recommended for good constraints (deviations < 5° were always seen with this value).
For more info, please read README_CONSTRAINTS.md.

```
params['peptide backbone constraint constant'] = 10000
```
#### Implicit Solvent
```
params['implicit solvent'] = True
if params['implicit solvent']:
    params['implicit solvent model'] = OBC1  # only meaningful if implicit solvent is True
    params['leap template'] = 'leap_template.in'
    # TODO add more options to params: implicitSolventSaltConc, soluteDielectric, solventDielectric, implicitSolventKappa
```

#### Starting with a folded DNA aptamer structure (instead of just a FASTA sequence)
```
params['skip MMB'] = True  # it will skip '2d analysis' and 'do MMB'
if params['skip MMB'] is True:
    params['folded initial structure'] = 'foldedSequence_0.pdb'  # if wishing to skip MMB, must provide a folded structure
```

