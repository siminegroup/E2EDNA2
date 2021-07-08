# E2EDNA - OpenMM Implementation - OpenDNA!

### A fully-automated code for DNA aptamer 2D and 3D structure analysis, and binding analysis with small molecules and short peptides.

**Michael Kilgour, Tao Liu, Lena Simine**

_mjakilgour gmail com_

## INSTALLATION
* Get the repository
```
git clone git@github.com:InfluenceFunctional/OpenDNA.git WHAT_YOU_WANT_IT_TO_BE_CALLED
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
scikit-learn==0.24.2
Cython==0.29.23
LightDock==0.9.0
nupack==4.0.0.23
mdtraj==1.9.5
```
* Download appropriate MacroMoleculeBuilder for your system from [SimTK here](https://simtk.org/projects/rnatoolbox) - old and new versions have been tested and seem to work. We'll need the 'Installer' directory and MMB executable. 
* Create working directory where runs will take place on your system. 
  Each new run will create a serially-numbered directory in this folder containing all the relevant outputs.
* Set paths in main.py and submission shell script, if relevant (example enclosed). Note, there is a toggle for running on a 'local' machine vs 'cluster', with distinct paths for e.g., developing on a local vs cluster environment.
```
* params['mmb'] --> set to the MacroMoleculeBuilder executable on your system
* params['mmb dir'] --> set to the MMB directory, usually a name like /Installer.#_##.Linux
- params['mmb template'] and params['mmb params'] --> included in this repository, should not have to be changed
* params['ld path'] --> paths to various lightdock python scripts (included in this package, alternatively automatically come with a lightdock pip installation)
* params['workdir'] --> the working directory you created in the previous step
* various .key files --> parameters may be specific to your analyte if custom force-field parameters are required. Adjust keyfiles accordingly, leaving the keywords on top as-is.
```

## Modes
OpenDNA takes in a DNA aptamer sequence in FASTA format, and optionally a short peptide or other small molecule, and returns details of the aptamer structure and binding behaviour.
This code implements several distinct analysis modes so users may customize the level of computational cost and accuracy.

* `2d structure` &rarr; returns NUPACK or seqfold analysis of aptamer secondary structure. Very fast, O(<1s). If using NUPACK, includes probability of observing a certain fold and of suboptimal folds within kT of the minimum.
* `3d coarse` &rarr; returns MMB fold of the best secondary structure. Fast O(5-30 mins). Results in a strained 3D structure which obeys base pairing rules and certain stacking interactions.
* `3d smooth` &rarr; identical to '3d coarse', with a short MD relaxation in solvent. ~Less than double the cost of '3d coarse' depending on relaxation time.
* `coarse dock` &rarr; uses the 3D structure from '3d coarse' as the initial condition for a LightDock simulation, and returns best docking configurations and scores. Depending on docking parameters, adds O(5-30mins) to '3d coarse'.
* `smooth dock` &rarr; identical to 'coarse dock', instead using the relaxed structure from '3d smooth'. Similar cost. 
* `free aptamer` &rarr; fold the aptamer in MMB and run extended MD sampling to identify a representative, equilibrated 2D and 3D structure. Slow O(hours).
* `full docking` &rarr; Return best docking configurations and scores from a LightDock run using the fully-equilibrated aptamer structure 'free aptamer'. Similar cost (LightDock is relatively cheap)
* `full binding` &rarr; Same steps as 'full docking', with follow-up extended MD simulation of the best binding configuration. Slowest O(hours).

# __ work in progress__


## Running a job

### Quickstart
* Set 'params' in main.py
* Set sequence = 'ATGC' where ATGC is your desired FASTA sequence
* Set peptide = 'YYYY' where YYYY is your desired peptide (restricted structures TBD)
* Run main, with desired run num. Zero for a fresh run, nonzero for either picking up on a prior run or explicitly enumerated new run (TEST)

`python main.py --run-num=0` In cluster mode, or set manually in local mode

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
