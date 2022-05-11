#!/bin/bash

# Assuming running on a local macos machine and planning to use NUPACK and MMB modules.

printf "Note: make sure to have copied automated_tests.sh and simu_config_automated_tests.yaml from examples/automated_tests/ to the codebase directory which contains main.py etc.\n\n"

printf "All the following tests start the pipeline with DNA sequence, therefore require NUPACK and MMB modules to be available.\n"
printf "The \"macos_installation.sh\" script provides instructions on setting up NUPACK and MMB.\n"


printf "\nStart automating tests one by one...\n"
printf "===================================="

printf "\nTESTING MODE #1: '2d structure'\n"
python ./main.py --yaml_config=simu_config_automated_tests.yaml --run_num=1 --mode='2d structure' --aptamer_seq='TAATGTTAATTG' #> run1.out
printf "\nEND OF TEST #1. Results are saved to folder './localruns/run1'. Check out the log file 'run_output_log.txt' there.\n"
printf "===================================="

printf "\nTESTING MODE #2: '3d coarse'\n"
python ./main.py --yaml_config=simu_config_automated_tests.yaml --run_num=2 --mode='3d coarse' --aptamer_seq='TAATGTTAATTG'  #> run2.out
printf "\nEND OF TEST #2. Results are saved to folder './localruns/run2'. Check out the log file 'run_output_log.txt' there.\n"
printf "===================================="

printf "\nTESTING MODE #3: '3d smooth'\n"
python ./main.py --yaml_config=simu_config_automated_tests.yaml --run_num=3 --mode='3d smooth' --aptamer_seq='TAATGTTAATTG'  --smoothing_time=0.001 #> run3.out
printf "\nEND OF TEST #3. Results are saved to folder './localruns/run3'. Check out the log file 'run_output_log.txt' there.\n"
printf "===================================="

printf "\nTESTING MODE #4: 'coarse dock'\n"
python ./main.py --yaml_config=simu_config_automated_tests.yaml --run_num=4 --mode='coarse dock' --aptamer_seq='TAATGTTAATTG' --ligand='examples/example_peptide_ligand.pdb' --ligand_type='peptide' --ligand_seq='YQTQTNSPRRAR' #> run4.out
printf "\nEND OF TEST #4. Results are saved to folder './localruns/run4'. Check out the log file 'run_output_log.txt' there.\n"
printf "===================================="

printf "\nTESTING MODE #5: 'smooth dock'\n"
python ./main.py --yaml_config=simu_config_automated_tests.yaml --run_num=5 --mode='smooth dock' --aptamer_seq='TAATGTTAATTG' --smoothing_time=0.001 --ligand='examples/example_peptide_ligand.pdb' --ligand_type='peptide' --ligand_seq='YQTQTNSPRRAR' #> run5.out
printf "\nEND OF TEST #5. Results are saved to folder './localruns/run5'. Check out the log file 'run_output_log.txt' there.\n"
printf "===================================="

printf "\nTESTING MODE #6: 'free aptamer'\n"
python ./main.py --yaml_config=simu_config_automated_tests.yaml --run_num=6 --mode='free aptamer' --aptamer_seq='TAATGTTAATTG' --aptamer_sampling_time=0.001 #> run6.out
printf "\nEND OF TEST #6. Results are saved to folder './localruns/run6'. Check out the log file 'run_output_log.txt' there.\n"
printf "===================================="

printf "\nTESTING MODE #7: 'full dock'\n"
python ./main.py --yaml_config=simu_config_automated_tests.yaml --run_num=7 --mode='full dock' --aptamer_seq='TAATGTTAATTG' --aptamer_sampling_time=0.001 --ligand='examples/example_peptide_ligand.pdb' --ligand_type='peptide' --ligand_seq='YQTQTNSPRRAR' #> run7.out
printf "\nEND OF TEST #7. Results are saved to folder './localruns/run7'. Check out the log file 'run_output_log.txt' there.\n"
printf "===================================="

printf "\nTESTING MODE #8: 'full binding'\n"
python ./main.py --yaml_config=simu_config_automated_tests.yaml --run_num=8 --mode='full binding' --aptamer_seq='TAATGTTAATTG' --aptamer_sampling_time=0.001 --ligand='examples/example_peptide_ligand.pdb' --ligand_type='peptide' --ligand_seq='YQTQTNSPRRAR' --complex_sampling_time=0.001 #> run8.out
printf "\nEND OF TEST #8. Results are saved to folder './localruns/run8'. Check out the log file 'run_output_log.txt' there.\n"
printf "====================================\n"
