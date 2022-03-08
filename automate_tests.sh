#!/bin/bash

# Assumption: lcoal, macos; user has defined work dir, MMB path in main.py and provided target ligand's pdb etc.
# in short, it is ready to run and test different modes

printf "Start automating tests one by one...\n"
printf "===================================="

printf "\nTESTING MODE #1: '2d structure'\n"
python main.py --outdir 1 --mode '2d structure' --aptamer 'TAATGTTAATTG' --ligand 'False' --ligand_type '' --ligand_seq '' #> run1.out
printf "\nEND OF TEST #1. Results are saved to folder \"run1\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "===================================="

printf "\nTESTING MODE #2: '3d coarse'\n"
python main.py --outdir 2 --mode '3d coarse' --aptamer 'TAATGTTAATTG' --ligand 'False' --ligand_type '' --ligand_seq '' #> run2.out
printf "\nEND OF TEST #2. Results are saved to folder \"run2\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_0.pdb\n"
printf "===================================="

printf "\nTESTING MODE #3: '3d smooth'\n"
python main.py --outdir 3 --mode '3d smooth' --aptamer 'TAATGTTAATTG' --ligand 'False' --ligand_type '' --ligand_seq '' #> run3.out
printf "\nEND OF TEST #3. Results are saved to folder \"run3\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_0.pdb\n"
printf "\tShort MD relaxation trajectory of free aptamer: foldedAptamer_0_processed_trajectory.dcd and clean_foldedAptamer_0_processed_trajectory.dcd (without solvent and ions)\n"
printf "\tRelaxed aptamer structure: relaxedAptamer_0.pdb\n"
printf "===================================="

printf "\nTESTING MODE #4: 'coarse dock'\n"
python main.py --outdir 4 --mode 'coarse dock' --aptamer 'TAATGTTAATTG' --ligand 'YQTQ.pdb' --ligand_type 'peptide' --ligand_seq 'YQTQTNSPRRAR' #> run4.out
printf "\nEND OF TEST #4. Results are saved to folder \"run4\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_0.pdb\n"
printf "\tMMB-folded aptamer docked by target ligand: complex_0_0.pdb (if docking happened)\n"
printf "===================================="

printf "\nTESTING MODE #5: 'smooth dock'\n"
python main.py --outdir 5 --mode 'smooth dock' --aptamer 'TAATGTTAATTG' --ligand 'YQTQ.pdb' --ligand_type 'peptide' --ligand_seq 'YQTQTNSPRRAR' #> run5.out
printf "\nEND OF TEST #5. Results are saved to folder \"run5\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_0.pdb\n"
printf "\tShort MD relaxation trajectory of free aptamer: foldedAptamer_0_processed_trajectory.dcd and clean_foldedAptamer_0_processed_trajectory.dcd (without solvent and ions)\n"
printf "\tRelaxed aptamer structure: relaxedAptamer_0.pdb\n"
printf "\tRelaxed aptamer docked by target ligand: complex_0_0.pdb (if docking happened)\n"
printf "===================================="

printf "\nTESTING MODE #6: 'free aptamer'\n"
python main.py --outdir 6 --mode 'free aptamer' --aptamer 'TAATGTTAATTG' --ligand 'False' --ligand_type '' --ligand_seq '' #> run6.out
printf "\nEND OF TEST #6. Results are saved to folder \"run6\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_0.pdb\n"
printf "\tLong MD sampling trajectory of free aptamer: foldedAptamer_0_processed_complete_trajectory.dcd and clean_foldedAptamer_0_processed_complete_trajectory.dcd (without solvent and ions)\n"
printf "\tRepresentative structure of free aptamer: repStructure_0.pdb\n"
printf "===================================="

printf "\nTESTING MODE #7: 'full dock'\n"
python main.py --outdir 7 --mode 'full dock' --aptamer 'TAATGTTAATTG' --ligand 'YQTQ.pdb' --ligand_type 'peptide' --ligand_seq 'YQTQTNSPRRAR' #> run7.out
printf "\nEND OF TEST #7. Results are saved to folder \"run7\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_0.pdb\n"
printf "\tLong MD sampling trajectory of free aptamer: foldedAptamer_0_processed_complete_trajectory.dcd and clean_foldedAptamer_0_processed_complete_trajectory.dcd (without solvent and ions)\n"
printf "\tRepresentative structure of free aptamer: repStructure_0.pdb\n"
printf "\tRepresentative aptamer docked by target ligand: complex_0_0.pdb (if docking happened)\n"
printf "===================================="

printf "\nTESTING MODE #8: 'full binding'\n"
python main.py --outdir 8 --mode 'full binding' --aptamer 'TAATGTTAATTG' --ligand 'YQTQ.pdb' --ligand_type 'peptide' --ligand_seq 'YQTQTNSPRRAR' #> run8.out
printf "\nEND OF TEST #8. Results are saved to folder \"run8\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_0.pdb\n"
printf "\tLong MD sampling trajectory of free aptamer: foldedAptamer_0_processed_complete_trajectory.dcd and clean_foldedAptamer_0_processed_complete_trajectory.dcd (without solvent and ions)\n"
printf "\tRepresentative structure of free aptamer: repStructure_0.pdb\n"
printf "\tRepresentative aptamer docked by target ligand: complex_0_0.pdb (if docking happened)\n"
printf "\tLong MD sampling trajectory of aptamer-ligand complex: complex_0_0_processed_complete_trajectory.dcd and clean_complex_0_0_processed_complete_trajectory.dcd (without solvent and ions)\n"
printf "====================================\n"
