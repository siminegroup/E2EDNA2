# What should we set as command inputs to automate the tests?
# assumption: lcoal, macos; user has defined work dir, MMB path in main.py and provided target ligand's pdb etc.
# in short, it is ready to run and test different modes

printf "Start automating tests one by one...\n"
printf "===================================="

printf "\nTESTING MODE #1: '2d structure'\n"
python main.py --run_num=1 --mode='2d structure' --aptamerSeq='TAATGTTAATTG' --ligand='False' --ligandType='' --ligandSeq='' > run1.out
printf "\nEND OF TEST #1. Results are saved to folder \"run1\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "===================================="

printf "\nTESTING MODE #2: '3d coarse'\n"
python main.py --run_num=2 --mode='3d coarse' --aptamerSeq='TAATGTTAATTG' --ligand='False' --ligandType='' --ligandSeq='' > run2.out
printf "\nEND OF TEST #2. Results are saved to folder \"run2\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_#.pdb\n"
printf "===================================="

printf "\nTESTING MODE #3: '3d smooth'\n"
python main.py --run_num=3 --mode='3d smooth' --aptamerSeq='TAATGTTAATTG' --ligand='False' --ligandType='' --ligandSeq='' > run3.out
printf "\nEND OF TEST #3. Results are saved to folder \"run3\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_#.pdb\n"
printf "\tShort MD relaxation trajectory of free aptamer: foldedAptamer_#_processed_trajectory.dcd or clean_foldedAptamer_#_processed_trajectory.dcd (without solvent and ions)\n"
printf "\tRelaxed aptamer structure: relaxedAptamer_#.pdb\n"
printf "===================================="

printf "\nTESTING MODE #4: 'coarse dock'\n"
python main.py --run_num=4 --mode='coarse dock' --aptamerSeq='TAATGTTAATTG' --ligand='YQTQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR' > run4.out
printf "\nEND OF TEST #4. Results are saved to folder \"run4\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_#.pdb\n"
printf "\tMMB-folded aptamer docked by target ligand: complex_#_#.pdb (if any docking)\n"
printf "===================================="

printf "\nTESTING MODE #5: 'smooth dock'\n"
python main.py --run_num=5 --mode='smooth dock' --aptamerSeq='TAATGTTAATTG' --ligand='YQTQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR' > run5.out
printf "\nEND OF TEST #5. Results are saved to folder \"run5\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_#.pdb\n"
printf "\tShort MD relaxation trajectory of free aptamer: foldedAptamer_#_processed_trajectory.dcd or clean_foldedAptamer_#_processed_trajectory.dcd (without solvent and ions)\n"
printf "\tRelaxed aptamer structure: relaxedAptamer_#.pdb\n"
printf "\tRelaxed aptamer docked by target ligand: complex_#_#.pdb (if any docking)\n"
printf "===================================="

printf "\nTESTING MODE #6: 'free aptamer'\n"
python main.py --run_num=6 --mode='free aptamer' --aptamerSeq='TAATGTTAATTG' --ligand='False' --ligandType='' --ligandSeq='' > run6.out
printf "\nEND OF TEST #6. Results are saved to folder \"run6\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_#.pdb\n"
printf "\tLong MD sampling trajectory of free aptamer: foldedAptamer_#_processed_complete_trajectory.dcd or clean_foldedAptamer_#_processed_complete_trajectory.dcd (without solvent and ions)\n"
printf "\tRepresentative structure of free aptamer: repStructure_#.pdb\n"
printf "===================================="

printf "\nTESTING MODE #7: 'full dock'\n"
python main.py --run_num=7 --mode='full dock' --aptamerSeq='TAATGTTAATTG' --ligand='YQTQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR' > run7.out
printf "\nEND OF TEST #7. Results are saved to folder \"run7\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_#.pdb\n"
printf "\tLong MD sampling trajectory of free aptamer: foldedAptamer_#_processed_complete_trajectory.dcd or clean_foldedAptamer_#_processed_complete_trajectory.dcd (without solvent and ions)\n"
printf "\tRepresentative structure of free aptamer: repStructure_#.pdb\n"
printf "\tRepresentative aptamer docked by target ligand: complex_#_#.pdb (if any docking)\n"
printf "===================================="

printf "\nTESTING MODE #8: 'full binding'\n"
python main.py --run_num=8 --mode='full binding' --aptamerSeq='TAATGTTAATTG' --ligand='YQTQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR' > run8.out
printf "\nEND OF TEST #8. Results are saved to folder \"run8\", where:\n"
printf "\t2d structure: in record.txt\n"
printf "\tMMB-folded aptamer structure: foldedAptamer_#.pdb\n"
printf "\tLong MD sampling trajectory of free aptamer: foldedAptamer_#_processed_complete_trajectory.dcd or clean_foldedAptamer_#_processed_complete_trajectory.dcd (without solvent and ions)\n"
printf "\tRepresentative structure of free aptamer: repStructure_#.pdb\n"
printf "\tRepresentative aptamer docked by target ligand: complex_#_#.pdb (if any docking)\n"
printf "\tLong MD sampling trajectory of aptamer-ligand complex: complex_#_#_processed_complete_trajectory.dcd or clean_complex_#_#_processed_complete_trajectory.dcd (without solvent and ions)\n"
printf "====================================\n"
