# Documentation
## E2EDNA 2.0 - OpenMM Implementation of E2EDNA !

### An automated pipeline for simulating DNA aptamers complexed with target ligands (peptide, DNA, RNA or small molecules).

**Michael Kilgour, Tao Liu, Ilya S. Dementyev, Lena Simine**

* For original version of E2EDNA: J. Chem. Inf. Model. 2021, 61, 9, 4139–4144 (https://doi.org/10.1021/acs.jcim.1c00696) and https://github.com/InfluenceFunctional/E2EDNA

* Interested in contributing to developing E2EDNA? Check out how to contribute [here](https://github.com/siminegroup/E2EDNA2/blob/JOSS/CONTRIBUTING.md).
## 1. Installation
This installation path has been tested on macOS and it relies on conda and pip package managers. 
1. Download the E2EDNA 2.0 package from this repository.
2. Register and download NUPACK package from http://www.nupack.org/downloads, find out the path to the downloaded package, such as ``~/Downloads/nupack-4.0.0.27/package``.
3. In the E2EDNA2 directory please modify the [macos_installation.sh](https://github.com/siminegroup/E2EDNA2/blob/JOSS/macos_installation.sh) script: update the path to NUPACK installation package (see step 2).
4. From the E2EDNA2 directory run ``source macos_installation.sh`` in command line to install a virtual environment "e2edna" and required dependences. Note that (e2edna) environment should be activated when the script finishes. 
   * Caveat: If the scripts fails to activate the environment automatically, this is likely because `conda activate e2edna` command gives an error such as "Your shell has not been properly configured to use 'conda activate'".
   * To fix this issue, please, first remove the installed environment by ``conda remove --name e2edna --all`` in command line, then either replace the conda activation command in the script with `source activate /path-to-envs/e2edna`, e.g., `source activate ~/miniconda3/envs/e2edna` and run the script again 
   * *OR* run the commands one by one interactively in your terminal. 
5. Register and download MMB from https://simtk.org/projects/rnatoolbox. Place the "Installer.*version.OS*" folder, such as "Installer.3_0.OSX" into the E2EDNA2 directory. N.B.: Do not specify `DYLD_LIBRARY_PATH` for macOS against the recommendations of the MMB installation guide. This is to avoid interference with the OpenMM module.
6. Update 3 paths in [main.py](https://github.com/siminegroup/E2EDNA2/blob/JOSS/main.py):
```
params['workdir'] = '/path-to-E2EDNA2-JOSS/localruns'                                # directory manually created to store all future jobs
params['mmb dir'] = '/path-to-E2EDNA2-JOSS/Installer.version.OS/lib'                 # path to MMB dylib files
params['mmb']     = '/path-to-E2EDNA2-JOSS/Installer.version.OS/bin/MMB-executable'  # path to MMB executable
```
## 2. Running a job
### Quickstart
* Set 3 paths in `params`, as indicated in step6 of [Installation](https://github.com/siminegroup/E2EDNA2/tree/JOSS#1-installation) section. Leave the other keys of `params` unchanged.
* Run the bash script [automate_tests.sh](https://github.com/siminegroup/E2EDNA2/blob/JOSS/automate_tests.sh) to test run all 8 modes automatically. Check out [Test runs](https://github.com/siminegroup/E2EDNA2/tree/JOSS#4-test-runs-inputs-and-outcomes) section for details. 
* All 8 modes take about 15 minutes to run on a macbook pro laptop and the outputs are saved in the work directory. We chose very simple systems for the the test-run and where ever docking is performed the docking configuration may or may not be found due to the limitations of the chosen system. Please note that a failure to find a docked configuration is not a failure of E2EDNA. 
* Alternatively, a single run can be carried out by run_num, mode, aptamer sequence, and ligand's structural file. For example,
```
python main.py --run_num=1 --mode='free aptamer' --aptamerSeq='TAATGTTAATTG' --ligand='False' --ligandType='' --ligandSeq=''
python main.py --run_num=2 --mode='full dock' --aptamerSeq='TAATGTTAATTG' --ligand='YQTQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR'
    
# --ligand='False'        # if no ligand. --ligandType and --ligandSeq will be ignored.
# --ligandType='peptide'  # or 'DNA' or 'RNA' or 'other'. Assuming 'other' ligand can be described by Amber14 force field.
# --ligandSeq=''          # if no sequence. For instance, when ligandType is 'other'
```
### Customizing parameters
* parameters passed as command line inputs:
  * `--run_num`: is used to name the folder for each job.
  * `--aptamerSeq`: DNA aptamer sequence (5'->3').
  * `--ligand`: PDB filename of target ligand; If no ligand, `--ligand='False'` and the following two inputs can be empty strings.
  * `--ligandType`: 'peptide' or 'DNA' or 'RNA' or 'other', assuming 'other' ligand can be described by force field used in MD simulation (default: Amber14).
  * `--ligandSeq`: target ligand's sequence (5'->3') if target ligand is a peptide or DNA or RNA. If `--ligandType='other'`, this input is empty string.
* A handful of parameters can be adjusted in [main.py](https://github.com/siminegroup/E2EDNA2/blob/JOSS/main.py) to allow for flexible control of the simulation conditions, including device platform (CPU or CUDA), solvent ionic strength, temperature, pH, implicit solvent, simulation force field, etc.

## 3. Functionality: Eight different modes of operation
E2EDNA 2.0 takes in a DNA aptamer sequence (5' to 3'), and optionally a short peptide or other small molecule, and returns details of the aptamer structure and binding behaviour.
This code implements several distinct operation modes so users may customize the level of computational cost and accuracy.

* `2d structure` &rarr; returns NUPACK analysis of aptamer secondary structure. Very fast, O(<1s). If using NUPACK, includes probability of observing a certain fold and of suboptimal folds within kT of the minimum.
* `3d coarse` &rarr; returns MMB fold of the best secondary structure. Fast, O(5-30 mins). Results in a strained 3D structure which obeys base pairing rules and certain stacking interactions.
* `3d smooth` &rarr; identical to '3d coarse', with a short MD relaxation in solvent. ~Less than double the cost of '3d coarse' depending on relaxation time.
* `coarse dock` &rarr; uses the 3D structure from '3d coarse' as the initial condition for a LightDock simulation, and returns best docking configurations and scores. Depending on docking parameters, adds O(5-30mins) to '3d coarse'.
* `smooth dock` &rarr; identical to 'coarse dock', instead using the relaxed structure from '3d smooth'. Similar cost to 'coarse dock'.
* `free aptamer` &rarr; fold the aptamer in MMB and run extended MD sampling to identify a representative, equilibrated 2D and 3D structure. Slow, O(hours).
* `full dock` &rarr; Return best docking configurations and scores from a LightDock run using the fully-equilibrated aptamer structure 'free aptamer'. Similar cost (LightDock is relatively cheap)
* `full binding` &rarr; Same steps as 'full dock', with follow-up extended MD simulation of the best binding configuration. Slowest, O(hours).
## 4. Test runs: inputs and outcomes
Running this script [automate_tests.sh](https://github.com/siminegroup/E2EDNA2/blob/JOSS/automate_tests.sh) will automatically run light tests of all 8 modes.
Here we explain what outputs to look for and what success looks like for each mode.

1. `--mode='2d structure'`

* **Key inputs**: DNA aptamer sequence

* **Outputs**: predicted secondary structure in record.txt

* **Success evaluation**: observe the dot-bracket notion for secondary structure, such as .(((....))).

2. `--mode='3d coarse'`

* **Key inputs**: DNA aptamer sequence

* **Outputs**: predicted secondary structure in record.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb

* **Success evaluation**: visualize MMB-folded aptamer structure in software like VMD or PyMOL 


3. `--mode='3d smooth'`

* **Key inputs**: DNA aptamer sequence

* **Outputs**: predicted secondary structure in record.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; Short MD relaxation trajectory of free aptamer: foldedAptamer_0_processed_trajectory.dcd and clean_foldedAptamer_0_processed_trajectory.dcd (without solvent and ions); Relaxed aptamer structure: relaxedAptamer_0.pdb

* **Success evaluation**: simulation logfile, MDlog_freeAptamerSmoothing.txt, indicates 100% completion; visualize the relaxation trajectory and relaxed structure in software like VMD or PyMOL 

4. `--mode='coarse dock'`

* **Key inputs**: DNA aptamer sequence; PDB filename of target ligand

* **Outputs**: predicted secondary structure in record.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; MMB-folded aptamer docked by target ligand: complex_0_0.pdb (if docking happened)

* **Success evaluation**: visualize the docked structure, if docking happened, in software like VMD or PyMOL 

5. `--mode='smooth dock'`

* **Key inputs**: DNA aptamer sequence; PDB filename of target ligand

* **Outputs**: predicted secondary structure in record.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; Short MD relaxation trajectory of free aptamer: foldedAptamer_0_processed_trajectory.dcd and clean_foldedAptamer_0_processed_trajectory.dcd (without solvent and ions); Relaxed aptamer structure: relaxedAptamer_0.pdb; Relaxed aptamer docked by target ligand: complex_0_0.pdb (if docking happened)

* **Success evaluation**: visualize the docked structure, if docking happened, in software like VMD or PyMOL 


6. `--mode='free aptamer'`

* **Key inputs**: DNA aptamer sequence

* **Outputs**: predicted secondary structure in record.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; Long MD sampling trajectory of free aptamer: foldedAptamer_0_processed_complete_trajectory.dcd and clean_foldedAptamer_0_processed_complete_trajectory.dcd (without solvent and ions); Representative structure of free aptamer: repStructure_0.pdb

* **Success evaluation**: simulation logfile, MDlog_freeAptamerSampling.txt, indicates 100% completion; visualize the sampling trajectory and representative structure of free aptamer in software like VMD or PyMOL 

7. `--mode='full dock'`

* **Key inputs**: DNA aptamer sequence; PDB filename of target ligand
* **Outputs**: predicted secondary structure in record.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; Long MD sampling trajectory of free aptamer: foldedAptamer_0_processed_complete_trajectory.dcd and clean_foldedAptamer_0_processed_complete_trajectory.dcd (without solvent and ions); Representative structure of free aptamer: repStructure_0.pdb; Representative aptamer docked by target ligand: complex_0_0.pdb (if docking happened)
* **Success evaluation**: visualize the docked structure, if docking happened, in software like VMD or PyMOL 

8. `--mode='full binding'`

* **Key inputs**: DNA aptamer sequence; PDB filename of target ligand

* **Outputs**: predicted secondary structure in record.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; Long MD sampling trajectory of free aptamer: foldedAptamer_0_processed_complete_trajectory.dcd and clean_foldedAptamer_0_processed_complete_trajectory.dcd (without solvent and ions); Representative structure of free aptamer: repStructure_0.pdb; Representative aptamer docked by target ligand: complex_0_0.pdb (if docking happened); Long MD sampling trajectory of aptamer-ligand complex: complex_0_0_processed_complete_trajectory.dcd and clean_complex_0_0_processed_complete_trajectory.dcd (without solvent and ions)

* **Success evaluation**: simulation logfile, MDlog_complexSampling.txt, indicates 100% completion; visualize the sampling trajectory of aptamer-ligand complex in software like VMD or PyMOL 

