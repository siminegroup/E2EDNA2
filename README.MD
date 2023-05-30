<!-- # Documentation -->
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04182/status.svg)](https://doi.org/10.21105/joss.04182)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6546661.svg)](https://doi.org/10.5281/zenodo.6546661)
[![GitHub release](https://img.shields.io/github/v/release/siminegroup/E2EDNA2)](https://github.com/siminegroup/E2EDNA2/releases)

## E2EDNA 2.0 - OpenMM Implementation of E2EDNA !

### An automated pipeline for simulating DNA aptamers complexed with target ligands (peptide, DNA, RNA or small molecules).

* To view the original version of E2EDNA, refer to its [GitHub repo](https://github.com/InfluenceFunctional/E2EDNA) and [DOI](https://doi.org/10.1021/acs.jcim.1c00696). <!-- J. Chem. Inf. Model. 2021, 61, 9, 4139–4144 -->

* Interested in contributing to developing E2EDNA? Check out how to contribute [here](https://github.com/siminegroup/E2EDNA2/blob/JOSS/CONTRIBUTING.md).

## Reference
If you use this code in any future publications, please cite our work using `Kilgour et al., (2022). E2EDNA 2.0: Python Pipeline for Simulating DNA Aptamers with Ligands. Journal of Open Source Software, 7(73), 4182`

E2EDNA pipeline makes use of several other open-sourced software packages, therefore please be mindful of citing them as well:
- [NUPACK](http://www.nupack.org/home/references)
- [MMB](https://simtk.org/plugins/publications/index.php?type=group&id=359&pluginname=publications)
- [OpenMM](http://docs.openmm.org/7.7.0/userguide/introduction.html#referencing-openmm)
- [LightDock](https://lightdock.org/reference/) 
<!-- - or [LightDock_GitHub](https://github.com/lightdock/lightdock#2-reference) -->

## Table of contents
- [Installation](#1-installation)
- [Usage](#2-usage)
- [Running a job](#3-running-a-job)
  - [Using the scripts for automated tests](#using-the-scripts-for-automated-tests)
  - [Running on your own data with customized input arguments](#running-on-your-own-data-with-customized-input-arguments)
  - [Check results in an output directory](#check-results-in-an-output-directory)
- [Functionality of eight different operation modes](#4-functionality-of-eight-different-operation-modes)
- [Automated test runs](#5-automated-test-runs)


## 1. Installation
1. Download the E2EDNA 2.0 package from this repository.
2. Locate `macos_installation.sh` in the downloaded E2EDNA2 codebase directory. Then at the codebase directory, run `$ source macos_installation.sh` in command line to create a conda virtual environment named `e2edna` and install required dependences. The `e2edna` environment should be activated when the installation script finishes, which means a string '(e2edna)' should show up at the beginning of the command line prompt.
   * If the script fails to activate the environment automatically, this is likely because `$ conda activate e2edna` command in the script gives an error such as `Your shell has not been properly configured to use 'conda activate'`. 
   * If so, manually run `$ source activate <path_to_e2edna_conda_environment>` to activate the environment. To help find the path, run `$ conda info -e` to list all conda environments and their paths on your computer.
3. As the message indicates at the end of installation process, if you wish to execute E2EDNA pipeline with a DNA aptamer sequence rather than its 3D structure, please register and download MMB from https://simtk.org/projects/rnatoolbox. Then copy or move the downloaded MMB folder to the codebase directory and remember to fill `MMB-related paths` section in the configuration file `simu_config.yaml`
<!-- N.B.: Our experience recommends macOS users not specify `DYLD_LIBRARY_PATH` against the MMB installation guide, to avoid interference with the OpenMM module. No longer need to worry, it is taken care of in the interfaces.py-->


## 2. Usage
The usage and help statements can be accessed with the `-h/--help` flags:
```
(e2edna)$ ./main.py --help
usage: main.py [-h] -yaml [-ow] [-d] [-os] [-p] [--CUDA_precision] [-w DIR] [-mbdir] [-mb] [--quick_check_mode] [-r] [-m] [-a]
               [-l] [-lt] [-ls] [--example_target_pdb] [--example_peptide_seq] [--skip_MMB] [-init] [--secondary_structure_engine]
               [--N_2D_structures] [--Mg_conc] [--fold_fidelity] [--fold_speed] [--mmb_normal_template] [--mmb_quick_template]
               [--mmb_slow_template] [--mmb_params] [-pk] [--pickup_from_freeAptamerChk] [--pickup_from_complexChk] [--chk_file]
               [--pickup_pdb] [--pressure] [--temperature] [--ionicStrength] [--pH] [--auto_sampling] [--autoMD_convergence_cutoff]
               [--max_aptamer_sampling_iter] [--max_walltime] [--skip_smoothing] [--equilibration_time] [--smoothing_time]
               [--aptamer_sampling_time] [--complex_sampling_time] [--time_step] [--print_step] [--force_field] [--hydrogen_mass]
               [--water_model] [--box_offset] [--constraints] [--constraint_tolerance] [--rigid_water] [--nonbonded_method]
               [--nonbonded_cutoff] [--ewald_error_tolerance] [--friction] [--implicit_solvent] [--implicit_solvent_model]
               [--soluteDielectric] [--solventDielectric] [--implicit_solvent_Kappa] [--leap_template] [--DNA_force_field]
               [--docking_steps] [--N_docked_structures]

E2EDNA: Simulate DNA aptamers complexed with target ligands

optional arguments:
  -h, --help            show this help message and exit
  -yaml, --yaml_config 
                        A YAML configuration file that can specify all the arguments (default: simu_config.yaml)
  -ow, --overwrite      Overwrite existing --run_num (default: False)

Compute Platform Configuration:
  -d, --device      Device configuration (default: local)
  -os, --operating_system 
                        Operating system (default: macos)
  -p, --platform    Processing platform (default: CPU)
  --CUDA_precision    Precision of CUDA, if used (default: single)

Directory Settings:
  -w DIR, --workdir DIR
                        Working directory to store individual output runs (default: ./localruns)
  -mbdir, --mmb_dir 
                        MMB library directory (default: None)
  -mb, --mmb        Path to MMB executable (default: None)

Run Parameters:
  --quick_check_mode  Rapidly run a certain mode for quick check using default test parameters (default: Yes)
  -r, --run_num     Run number. Output will be written to {--workdir}/run{--run_num} (default: 1)
  -m, --mode        Run mode (default: None)
  -a, --aptamer_seq 
                        DNA Aptamer sequence (5'->3') (default: None)
  -l, --ligand      Name of PDB file for ligand structure; None if not to have ligand (default: None)
  -lt, --ligand_type 
                        Type of ligand molecule (default: None)
  -ls, --ligand_seq 
                        Ligand sequence if peptide, DNA, or RNA (default: None)
  --example_target_pdb 
                        An example peptide ligand included in E2EDNA package: used when wish to test docking (default:
                        examples/example_peptide_ligand.pdb)
  --example_peptide_seq 
                        The sequence of the example peptide ligand (default: YQTQTNSPRRAR)
  --skip_MMB          If `Yes`: skip both 2D structure analysis and MMB folding, and start with a known --init_structure (default: No)
  -init, --init_structure 
                        Name of PDB file if starting pipeline on a DNA aptamer with known structure (default: None)
  --secondary_structure_engine 
                        Pipeline module that is used to predict secondary structures (default: NUPACK)
  --N_2D_structures   Number of predicted secondary structures (default: 1)
  --Mg_conc           Magnesium molar concentration used in NUPACK: [0, 0.2] (default: 0.0)
  --fold_fidelity     Refold in MMB if score < `fold_fidelity` unless the `fold_speed` is `quick` (default: 0.9)
  --fold_speed        MMB folding speed (default: normal)
  --mmb_normal_template 
                        Path to MMB folding protocol of normal speed (default: lib/mmb/commands.template.dat)
  --mmb_quick_template 
                        Path to MMB folding protocol of quick speed (default: lib/mmb/commands.template_quick.dat)
  --mmb_slow_template 
                        Path to MMB folding protocol of slow speed (default: lib/mmb/commands.template_long.dat)
  --mmb_params        Path to parameter file bundled with MMB package (default: lib/mmb/parameters.csv)
  -pk, --pickup     Whether the run is to resume MD sampling of an unfinished run or an old run (default: No)
  --pickup_from_freeAptamerChk 
                        Resume MD sampling of free aptamer: skip everything before it (default: No)
  --pickup_from_complexChk 
                        Resume MD sampling of aptamer-ligand: skip everything before it (default: No)
  --chk_file          Name of checkpoint file for resuming MD sampling, format: <path>/<filename>.chk (default: None)
  --pickup_pdb        PDB file (topology+coordinates) for resuming MD sampling in explicit solvent, format: <path>/<filename>.pdb (default: None)
  --pressure          Pressure in the unit of atm (default: 1.0)
  --temperature       Temperature in Kelvin (default: 298.0)
  --ionicStrength     Sodium molar concentration (could be used by NUPACK and OpenMM) (default: 0.1)
  --pH                Could be used by OpenMM (default: 7.4)
  --auto_sampling     If `Yes`: run MD sampling till convergence, currently only feasible in free aptamer sampling (default: No)
  --autoMD_convergence_cutoff 
                        Convergence cutoff if doing auto_sampling (default: 0.01)
  --max_aptamer_sampling_iter 
                        Max number of iterations for free aptamer MD sampling if doing auto_sampling (default: 20)
  --max_walltime      Walltime in hours to check runtime (default: 24.0)
  --skip_smoothing    If `Yes`: no short MD relaxation before MD sampling (default: Yes)
  --equilibration_time 
                        Equilibration time in nanoseconds after energy minimization and before MD sampling (default: 0.1)
  --smoothing_time    Time in nanoseconds for short MD relaxation before MD sampling, if any (default: None)
  --aptamer_sampling_time 
                        MD sampling time in nanoseconds for free aptamer dynamics (default: None)
  --complex_sampling_time 
                        MD sampling time in nanoseconds for aptamer-ligand complex dynamics (default: None)
  --time_step         time step in femtoseconds in MD sampling (default: 2.0)
  --print_step        Printout step in picoseconds in MD sampling (default: 10.0)
  --force_field       Force field used in OpenMM (default: amber14-all)
  --hydrogen_mass     Unit is amu (default: 1.5)
  --water_model       Explicit water solvent model used in OpenMM (default: amber14/tip3p)
  --box_offset        Buffering offset in nanometers on solvent box, if using explicit solvent (default: 1.0)
  --constraints       Specify which bond angles and/or lengths should be implemented with constraints (default: None)
  --constraint_tolerance 
                        Distance tolerance for constraint in OpenMM integrator (default: 1e-06)
  --rigid_water       Whether to make water molecules completely rigid at bond lengths and angles (default: Yes)
  --nonbonded_method  Type of nonbonded interactions (default: NoCutoff)
  --nonbonded_cutoff  The cutoff distance in nanometers to use for nonbonded interactions (default: 1.0)
  --ewald_error_tolerance 
                        Error tolerance if `nonbonded_method` is `Ewald`, `PME`, or `LJPME` (default: 0.0005)
  --friction          Friction coefficient in unit of 1/ps, used in Langevin integrator (default: 1.0)
  --implicit_solvent  Whether to use an Amber GB implicit solvent model (default: No)
  --implicit_solvent_model 
                        Specify an Amber GB implicit solvent model if needed (default: None)
  --soluteDielectric  The solute dielectric constant to use in the implicit solvent model (default: 1.0)
  --solventDielectric 
                        The solvent dielectric constant to use in the implicit solvent model (default: 78.5)
  --implicit_solvent_Kappa 
                        Debye screening parameter; If specified by user, OpenMM will ignore {--ionicStrength} in implicit solvent. (default: None)
  --leap_template     A script for running LEap program in Ambertools21, provided in E2EDNA package. (default: leap_template.in)
  --DNA_force_field   Force field for DNA used by LEap program (default: DNA.OL15)
  --docking_steps     Number of steps for docking simulations (default: 10)
  --N_docked_structures 
                        Number of docked structures output from the docker (default: 1)
```

## 3. Running a job
### Using the scripts for automated tests
* `examples/automated_tests/` folder provides `.sh` and `.yaml` for automated testing. Check out `how_to_run_automated_tests.txt` in the folder for instructions. Each set of automated tests takes about 15 minutes to complete on a macbook pro laptop. We chose a simple DNA aptamer system for the "test runs" purpose.
<!-- and whereever docking is performed the docking configuration may or may not be found due to the limitations of the chosen system. Please note that a failure to find a docked configuration is not a failure of E2EDNA.  -->

### Running on your own data with customized input arguments
All the input arguments listed in [Usage](#2-usage) can be customized either in command line or a `.yaml` configuration file, or both. The configuration `.yaml` is designed to be superior to command line inputs. Therefore, if an argument is specified in both `.yaml` and command line, the one in command line will be ignored.

A single run can be carried out by specifying all necessary input arguments in a configuration file:
```
(e2edna)$ ./main.py --yaml_config=simu_config.yaml
```
Alternatively, part or all of those necessary input arguments can also be passed into the pipeline via command line while having been commented out in the `.yaml` file. See `simu_config_automated_tests.yaml` and `automated_tests.sh` in `examples/automated_tests/` as an example.

Below list the part of parameters in [Usage](#2-usage) which call for particular attention, due to their values have limited choices, for example.
* Compute Platform Configuration
  * `-d/--device`: running device; either `local` or `cluster`
  * `-os/--operating_system`: operating sysmtem; `macos` or `linux` or `WSL`
  * `-p/--platform`: processing platform; either `CPU` or `CUDA`
  * `--CUDA_precision`: precision of CUDA if used; either `single` or `double`; Default is `single`

* Directory Settings
  * `-w/--workdir`: directory to write results for each run; Default is `./localruns`
  * `-md/--mmb_dir`: path to MMB library directory. Both absolute and relative (to the codebase directory) paths are accepted.
  <!-- Wildcards can be used to describe location, for instance the default is "Installer-/lib" which will match both "Installer-3.0-Ubuntu18.04/lib" and "Installer.3_0.OSX/lib".  -->
  * `-mb/--mmb`: path to MMB executable. Both absolute and relative (to the codebase directory) paths are accepted.
  <!-- Wild cards can be used to describe location. -->

* Run Parameters    
  * `-m/--mode`: mode of operation; Must be one of the modes described in [Functionality](#4-functionality-of-eight-different-operation-modes), ie, `'2d structure'`, `'3d coarse'`, `'3d smooth'`, `'coarse dock'`, `'smooth dock'`, `'free aptamer'`, `'full dock'`, `'full binding'`
  * `-a/--aptamerSeq`: DNA aptamer sequence (5'->3'). A string made of case sensitive letters from {A, G, C, T} only.
  <!-- or the name of a readable text file containing the sequence. -->
  * `-l/--ligand`: PDB filename of target ligand; If no ligand, `--ligand` and the following two inputs (`--ligand_type` and `--ligand_seq`) should be left off.
  * `-lt/--ligand_type`: `peptide`, `DNA`, `RNA`, or `other`, assuming `other` ligand can be described by force field used in MD simulation (default is Amber14).
  * `-ls/--ligand_seq`: a string of target ligand's sequence if `--ligand_type` is a `peptide` or `DNA` or `RNA`. If `--ligand_type=other`, do not use `--ligand_seq` flag.
  <!-- or the name of a readable text file containing the sequence.  -->
  * `--fold_speed`: three choices; `quick`, `normal` or `slow`
  * `--force_field` and `--water_model`: plenty of choices, check out options of [forcefields in OpenMM/7.7.0](http://docs.openmm.org/7.7.0/userguide/application/02_running_sims.html#force-fields).
  * `--constraints`: four choices; `HBonds`, `AllBonds`, `HAngles` or `None`.
  * `--nonbonded_method`: at most six choices; `Ewald`, `PME`, `LJPME`, `CutoffPeriodic`, `CutoffNonPeriodic` or `NoCutoff`; Only the last three if using an Amber GB implicit solvent model.
  * `--implicit_solvent_model`: five choices; `HCT`, `OBC1`, `OBC2`, `GBn` or `GBn2`
  * `--DNA_force_field`: two choices; `DNA.OL15` or `DNA.bsc1`

### Check results in an output directory
Output of a single run will be written to `{--workdir}/run{--run_num}` directory. 

Only the main outputs remain in the `run{--run_num}/` folder, including a log file named `run_output_log.txt`. All intermediate or temporary files are grouped to subfolders, such as `run1/md_aptamer_sampling_runfiles_0`, based on which module they were generated from. In each `run_output_log.txt`, generation statements of those main output files are marked with `>>>` and one can selectively print out:
```
$ cat run_output_log.txt | grep '>>>'
>>> Predicted 2D structure #0                   : .(((....))).
>>> MMB folded the aptamer and generated folded structure: foldedAptamer_0.pdb
>>> Generated aptamer-ligand complex structure: complex_0_0.pdb
$
```

## 4. Functionality of eight different operation modes
The pipeline could implement several distinct operation modes so users may customize the level of computational cost and accuracy.

* `'2d structure'` &rarr; returns NUPACK analysis of aptamer secondary structure. Very fast, O(<1s). If using NUPACK, includes probability of observing a certain fold and of suboptimal folds within kT of the minimum.
* `'3d coarse'` &rarr; returns MMB fold of the best secondary structure. Fast, O(5-30 mins). Results in a strained 3D structure which obeys base pairing rules and certain stacking interactions.
* `'3d smooth'` &rarr; identical to `'3d coarse'`, with a short MD relaxation in solvent. About less than double the cost of `'3d coarse'` depending on relaxation time.
* `'coarse dock'` &rarr; uses the 3D structure from `'3d coarse'` as the initial condition for a LightDock simulation, and returns best docking configurations and scores. Depending on docking parameters, adds O(5-30mins) to `'3d coarse'`.
* `'smooth dock'` &rarr; identical to 'coarse dock', instead using the relaxed structure from `'3d smooth'`. Similar cost to `'coarse dock'`.
* `'free aptamer'` &rarr; fold the aptamer in MMB and run extended MD sampling to identify a representative, equilibrated 2D and 3D structure. Slow, O(hours).
* `'full dock'` &rarr; Return best docking configurations and scores from a LightDock run using the fully-equilibrated aptamer structure `'free aptamer'`. Similar cost (LightDock is relatively cheap)
* `'full binding'` &rarr; Same steps as `'full dock'`, with follow-up extended MD simulation of the best binding configuration. Slowest, O(hours).

## 5. Automated test runs
Running the scripts of automated tests mentioned in [Running a job](#3-running-a-job) will automatically run light tests of 8 modes. Here we explain the inputs, what outputs to look for, and what a successful run should look like for each mode. 

1. `--mode='2d structure'`

* **Key inputs**: DNA aptamer sequence

* **Outputs**: predicted secondary structure in run_output_log.txt

* **Success evaluation**: observe the dot-bracket notion for secondary structure, such as .(((....))).

2. `--mode='3d coarse'`

* **Key inputs**: DNA aptamer sequence

* **Outputs**: predicted secondary structure in run_output_log.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb

* **Success evaluation**: visualize MMB-folded aptamer structure in software like VMD or PyMOL 


3. `--mode='3d smooth'`

* **Key inputs**: DNA aptamer sequence

* **Outputs**: predicted secondary structure in run_output_log.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; Short MD relaxation trajectory of free aptamer: foldedAptamer_0_processed_trajectory.dcd and clean_foldedAptamer_0_processed_trajectory.dcd (without solvent and ions); Relaxed aptamer structure: relaxedAptamer_0.pdb

* **Success evaluation**: simulation logfile, MDlog_freeAptamerSmoothing.txt, indicates 100% completion; visualize the relaxation trajectory and relaxed structure in software like VMD or PyMOL 

4. `--mode='coarse dock'`

* **Key inputs**: DNA aptamer sequence; PDB filename of target ligand

* **Outputs**: predicted secondary structure in run_output_log.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; MMB-folded aptamer docked by target ligand: complex_0_0.pdb (if docking happened)

* **Success evaluation**: visualize the docked structure, if docking happened, in software like VMD or PyMOL 

5. `--mode='smooth dock'`

* **Key inputs**: DNA aptamer sequence; PDB filename of target ligand

* **Outputs**: predicted secondary structure in run_output_log.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; Short MD relaxation trajectory of free aptamer: foldedAptamer_0_processed_trajectory.dcd and clean_foldedAptamer_0_processed_trajectory.dcd (without solvent and ions); Relaxed aptamer structure: relaxedAptamer_0.pdb; Relaxed aptamer docked by target ligand: complex_0_0.pdb (if docking happened)

* **Success evaluation**: visualize the docked structure, if docking happened, in software like VMD or PyMOL 


6. `--mode='free aptamer'`

* **Key inputs**: DNA aptamer sequence

* **Outputs**: predicted secondary structure in run_output_log.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; Long MD sampling trajectory of free aptamer: foldedAptamer_0_processed_complete_trajectory.dcd and clean_foldedAptamer_0_processed_complete_trajectory.dcd (without solvent and ions); Representative structure of free aptamer: repStructure_0.pdb

* **Success evaluation**: simulation logfile, MDlog_freeAptamerSampling.txt, indicates 100% completion; visualize the sampling trajectory and representative structure of free aptamer in software like VMD or PyMOL 

7. `--mode='full dock'`

* **Key inputs**: DNA aptamer sequence; PDB filename of target ligand
* **Outputs**: predicted secondary structure in run_output_log.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; Long MD sampling trajectory of free aptamer: foldedAptamer_0_processed_complete_trajectory.dcd and clean_foldedAptamer_0_processed_complete_trajectory.dcd (without solvent and ions); Representative structure of free aptamer: repStructure_0.pdb; Representative aptamer docked by target ligand: complex_0_0.pdb (if docking happened)
* **Success evaluation**: visualize the docked structure, if docking happened, in software like VMD or PyMOL 

8. `--mode='full binding'`

* **Key inputs**: DNA aptamer sequence; PDB filename of target ligand

* **Outputs**: predicted secondary structure in run_output_log.txt; MMB-folded aptamer structure: foldedAptamer_0.pdb; Long MD sampling trajectory of free aptamer: foldedAptamer_0_processed_complete_trajectory.dcd and clean_foldedAptamer_0_processed_complete_trajectory.dcd (without solvent and ions); Representative structure of free aptamer: repStructure_0.pdb; Representative aptamer docked by target ligand: complex_0_0.pdb (if docking happened); Long MD sampling trajectory of aptamer-ligand complex: complex_0_0_processed_complete_trajectory.dcd and clean_complex_0_0_processed_complete_trajectory.dcd (without solvent and ions)

* **Success evaluation**: simulation logfile, MDlog_complexSampling.txt, indicates 100% completion; visualize the sampling trajectory of aptamer-ligand complex in software like VMD or PyMOL 
