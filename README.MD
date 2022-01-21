# Documentation
## E2EDNA 2.0 - OpenMM Implementation of E2EDNA !

### An automated pipeline for simulation of DNA aptamers complexed with small molecules and short peptides.

**Michael Kilgour, Tao Liu, Ilya S. Dementyev, Lena Simine**

_mjakilgour gmail com_

For original version of E2EDNA: J. Chem. Inf. Model. 2021, 61, 9, 4139–4144 (https://doi.org/10.1021/acs.jcim.1c00696)

## Installation
This installation path has been tested on macOS and it relies on conda and pip package managers. 
1. Download the E2EDNA 2.0 package from this repository.
2. Register and download NUPACK from http://www.nupack.org/downloads, you will need the path to ~/nupack-###/source/package directory
3. In the E2EDNA2 directory please modify the [``macos_installation.sh``](https://github.com/InfluenceFunctional/E2EDNA2/blob/main/macos_installation.sh) script: update the path to nupack (see step 2)
4. From E2EDNA2 folder run [``macos_installation.sh``](https://github.com/InfluenceFunctional/E2EDNA2/blob/main/macos_installation.sh). 
   * Caveat: in case `conda activate e2edna` command gives an error or if after the script finishes e2edna enviroment has not been activated, please either replace the activation command with with `source activate /path-to-env/e2edna` 
   * OR alternatively copy paste commands from the script without modifications to command line and run one by one, this will go around the unconfigured shell issue.
6. Register and download MMB from https://simtk.org/projects/rnatoolbox . Place the Installer### folder into the e2edna folder. NB: Do not specify `DYLD_LIBRARY_PATH` against the recommendations of the MMB installation guide. This is to avoid interference with the OpenMM module.
7. Update 3 paths in main.py:
```
      params['workdir'] = '/path-to-e2edna/localruns'                         # working directory   
       
      params['mmb dir'] = '/path-to-e2edna/e2edna/Installer.###/lib'          # path to MMB dylib files
       
      params['mmb']     = '/path-to-e2edna/Installer.###/bin/MMB-executable'  # path to MMB executable    
```
## Running a job

### Quickstart
* Set 'params' in [``main.py``](https://github.com/InfluenceFunctional/E2EDNA2/blob/main/main.py), as indicated in "Installation".
* Run the bash script [``automate_tests.sh``](https://github.com/InfluenceFunctional/E2EDNA2/blob/main/automate_tests.sh) to test all 8 modes automatically.
* Alternatively, a single run can be carried out by run_num, mode, aptamer sequence, and ligand's structural file. For example,
```
python main.py --run_num=1 --mode='free aptamer' --aptamerSeq='TAATGTTAATTG' --ligand='False' --ligandType='' --ligandSeq=''
python main.py --run_num=2 --mode='full dock' --aptamerSeq='TAATGTTAATTG' --ligand='YQTQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR'
    
# --ligand='False'        # if no ligand. --ligandType and --ligandSeq will be ignored.
# --ligandType='peptide'  # or 'DNA' or 'RNA' or 'other'. Assuming 'other' ligand can be described by Amber14 force field.
# --ligandSeq=''          # if no sequence. For instance, when ligandType is 'other'

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

# Test run: inputs and outcomes
Running this script [``automate_tests.sh``](https://github.com/InfluenceFunctional/E2EDNA2/blob/main/automate_tests.sh) will automatically run simple very light simulations of all 8 modes.
Here we explain what outputs to look for and what success looks like.

* Mode 1:`2d structure` 
Input: fasta sequence, e.g, CGCGCGCGCGCGC

Outputs:

Success evaluation: observe the dot-parenthesis representation for 2d structure, e.g., ..(...)..

* Mode 2:`3d coarse`

Input: ‘3d unrefined’,  fasta sequence, e.g, CGCGCGCGCGCGC

Outputs: 

Success evaluation: Visualize foldedAptamer_0.pdb in VMD or PyMOL


* Mode 3:`3d coarse`

Input: fasta sequence, e.g, CGCGCGCGCGCGC 

Outputs: 

Success evaluation:

* Mode 4:`coarse dock`

Input: fasta sequence, e.g, CGCGCGCGCGCGC; PDB of target ligand, e.g., ‘target.pdb’

Outputs: 

Success evaluation:

* Mode 5:`smooth dock`
Input: fasta sequence, e.g, CGCGCGCGCGCGC; PDB of target ligand, e.g., ‘target.pdb’

Outputs: 

Success evaluation:


* Mode 6: `free aptamer`
Given a DNA sequence, its secondary structure will be predicted and represented by a contact map of dot-parenthesis notation. Under the guidance of the predicted secondary structure, the sequence will then be folded into an initial three dimensional conformation. Last step is to run a molecular-dynamics simulation to sample its conformational space and find out the representative conformation from the MD trajectory. (if we ask contact predictor for >1 ssStructure)

Input: fasta sequence, e.g, CGCGCGCGCGCGC

Modifications to the code: 
set params[‘mode’] = ‘free aptamer’
params['sequence'] =’CGCGCGCGCGCGC’

Outputs: 
Secondary structure prediction: such as ((....))....((.(...).)).. in “record.txt”
MMB folded structure: “foldedAptamer_0.pdb”

MD simulation: 
Binary trajectory: “clean_foldedAptamer_0_processed_complete_trajectory.dcd”
Topology: “clean_foldedAptamer_0_processed.pdb”
Representative conformation: “repStructure_0.pdb”

Success evaluation: 
The DCD trajectory file is generated, and file “log.txt” shows that the MD sampling of free aptamer is 100% complete.
Visualize MD trajectory of free aptamer using the topology and the binary trajectory file.
Visualize representative conformation of the DNA aptamer.

* Mode 7: `full dock`
Given a DNA sequence, its secondary structure will be predicted and represented by a contact map of dot-parenthesis notation. Under the guidance of the predicted secondary structure, the sequence will then be folded into an initial three dimensional conformation. Next is to run a molecular-dynamics simulation to sample its conformational space and find out the representative conformation from the MD trajectory. Finally, the representative structure will be docked by a target ligand of interest (its structure must be provided as a PDB file).

Input: fasta sequence, e.g, CGCGCGCGCGCGC; PDB of target ligand, e.g., ‘target.pdb’
Modifications to the code: 
set params[‘mode’] = ‘full docking’ 
params['sequence'] =’CGCGCGCGCGCGC’
params[‘target’] = ‘target.pdb’  # need to update the code for this.

Outputs: 
Secondary structure prediction: such as ((....))....((.(...).)).. in “record.txt”
MMB folded structure: “foldedAptamer_0.pdb”

MD simulation: 
Binary trajectory: “clean_foldedAptamer_0_processed_complete_trajectory.dcd”
Topology: “clean_foldedAptamer_0_processed.pdb”
Representative conformation: “repStructure_0.pdb”
Docking: Aptamer-ligand complex structure: “top_1.pdb”. Docking score is in “record.txt”.

Success evaluation: 
The DCD trajectory file is generated, and file “log.txt” shows that the MD sampling of free aptamer is 100% complete.
Visualize MD trajectory of free aptamer using the topology and the binary trajectory file.
Visualize representative conformation of the DNA aptamer.
Visualize aptamer-ligand complex structure.

* Mode 8: `full binding`
Given a DNA sequence, its secondary structure will be predicted and represented by a contact map of dot-parenthesis notation. Under the guidance of the predicted secondary structure, the sequence will then be folded into an initial three dimensional conformation. Next is to run a molecular-dynamics simulation to sample its conformational space and find out the representative conformation from the MD trajectory. The representative structure will be docked by a target ligand of interest (its structure must be provided as a PDB file). Finally, the aptamer-ligand complex molecule will be sampled by MD simulation to investigate its dynamics.

Input: fasta sequence, e.g, CGCGCGCGCGCGC; PDB of target ligand, e.g., ‘target.pdb’
Modifications: 
set params[‘mode’] = ‘full binding’
params['sequence'] =’CGCGCGCGCGCGC’
params[‘target’] = ‘target.pdb’  # need to update the code for this.

Outputs: 
Secondary structure prediction: such as ((....))....((.(...).)).. in “record.txt”
MMB folded structure: “foldedAptamer_0.pdb”
MD simulation of free aptamer: 
Binary trajectory: “clean_foldedAptamer_0_processed_complete_trajectory.dcd”
Topology: “clean_foldedAptamer_0_processed.pdb”
Representative conformation: “repStructure_0.pdb”
Docking: Aptamer-ligand complex structure: “top_1.pdb”. Docking score is in “record.txt”.
MD simulation of aptamer-ligand complex:
Binary trajectory: “clean_complex_0_0_processed_trajectory.dcd”
Topology: “clean_complex_0_0_processed.pdb”

Success evaluation: 
File “log.txt” shows that the MD sampling of free aptamer is 100% complete and the DCD trajectory file is generated.
Visualize MD trajectory of free aptamer using the topology and the binary trajectory file.
Visualize representative conformation of the DNA aptamer.
Visualize aptamer-ligand complex structure
The DCD trajectory file is generated, and file “log_complex.txt” shows that the MD sampling of aptamer-ligand is 100% complete.
Visualize MD trajectory of aptamer-ligand using its binary and topolog file. It is worth noting that the aptamer might seem far apart from the target ligand, which could be a result of the periodic boundary condition. Should we correct it or leave user to do it?

MD simulation might stop at the onset with “Particle coordinate is nan” error. It could be due to the energy minimization being too aggressive so tha the coordinate gets out of boundary, then integrator cannot work on those non-sense coordinate values. In this case, re-run the pipeline.

MMB folding could take a while if multiple refolding takes place for any tricky sequence.



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

