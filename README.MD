# Documentation
## E2EDNA 2.0 - OpenMM Implementation of E2EDNA !

### An automated pipeline for simulation of DNA aptamers complexed with small molecules and short peptides.

**Michael Kilgour, Tao Liu, Ilya S. Dementyev, Lena Simine**

_mjakilgour gmail com_

For original version of E2EDNA: J. Chem. Inf. Model. 2021, 61, 9, 4139–4144 (https://doi.org/10.1021/acs.jcim.1c00696)

## Installation
This installation path has been tested on macOS and it relies on conda and pip package managers. 
1. Download in any convenient way the E2EDNA package from this repository.
2. Register and download NUPACK from http://www.nupack.org/downloads, note the path to ~/nupack-###/source/package directory
3. In the e2edna directory please modify 2 lines the local_installation.sh script: anaconda3 or miniconda3 path, and the path to nupack (see step 2)
4. From e2edna folder run either local_installation.sh or macos_installation.sh (if using the latter please adjust the last line in the .yml file)
5. Register and download MMB from https://simtk.org/projects/rnatoolbox . Place the Installer### folder into the e2edna folder. NB: Do not specify DYLD_LIBRARY_PATH against the recommendations of the MMB installation guide. This is to avoid interference with the OpenMM module.
6. Update 3 paths in main.py: 
       params['workdir'] = '~/localruns'   # working directory                                                                                                     
       params['mmb dir'] = '~/e2edna/Installer.###/lib'      # path to MMB dylib files
       params['mmb']     = '~/e2edna/Installer.###/bin/MMB-executable' #path to MMB executable    
       
## Running a job

### Quickstart
* Set 'params' in [``main.py``](https://github.com/InfluenceFunctional/OpenDNA/blob/main/main.py), as indicated in "Installation".
* Run the bash script [``automate_tests.sh``](https://github.com/InfluenceFunctional/OpenDNA/blob/main/automate_tests.sh) to test all 8 modes automatically.
* Alternatively, a single run can be carried out by run_num, mode, aptamer sequence, and ligand's structural file. For example,
```
python main.py --run_num=1 --mode='free aptamer' --aptamerSeq='TAATGTTAATTG' --ligand='False' --ligandType='' --ligandSeq=''
python main.py --run_num=2 --mode='full dock' --aptamerSeq='TAATGTTAATTG' --ligand='YQTQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR'
    
#--ligand='False'        # if no ligand. --ligandType and --ligandSeq will be ignored.
#--ligandType='peptide'  # or 'DNA' or 'RNA' or 'other'. Assuming 'other' ligand can be described by Amber14 force field.
#--ligandSeq=''          # if no sequence. For instance, when ligandType is 'other'

```


## Functionality: Eight different modes of operation
E2EDNA 2.0 takes in a DNA aptamer sequence in FASTA format, and optionally a short peptide or other small molecule, and returns details of the aptamer structure and binding behaviour.
This code implements several distinct analysis modes so users may customize the level of computational cost and accuracy.

* `2d structure` &rarr; returns NUPACK or seqfold analysis of aptamer secondary structure. Very fast, O(<1s). If using NUPACK, includes probability of observing a certain fold and of suboptimal folds within kT of the minimum.
* `3d coarse` &rarr; returns MMB fold of the best secondary structure. Fast O(5-30 mins). Results in a strained 3D structure which obeys base pairing rules and certain stacking interactions.
* `3d smooth` &rarr; identical to '3d coarse', with a short MD relaxation in solvent. ~Less than double the cost of '3d coarse' depending on relaxation time.
* `coarse dock` &rarr; uses the 3D structure from '3d coarse' as the initial condition for a LightDock simulation, and returns best docking configurations and scores. Depending on docking parameters, adds O(5-30mins) to '3d coarse'.
* `smooth dock` &rarr; identical to 'coarse dock', instead using the relaxed structure from '3d smooth'. Similar cost. 
* `free aptamer` &rarr; fold the aptamer in MMB and run extended MD sampling to identify a representative, equilibrated 2D and 3D structure. Slow O(hours).
* `full dock` &rarr; Return best docking configurations and scores from a LightDock run using the fully-equilibrated aptamer structure 'free aptamer'. Similar cost (LightDock is relatively cheap)
* `full binding` &rarr; Same steps as 'full docking', with follow-up extended MD simulation of the best binding configuration. Slowest O(hours).

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

