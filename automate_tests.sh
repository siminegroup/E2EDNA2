# What should we set as command inputs to automate the tests?
# assumption: lcoal, macos; user has defined work dir, MMB path in main.py and provided target ligand's pdb etc.
# in short, it is ready to run and test different modes

printf "Starting the automated testing...\n"
printf "===================================="

printf "\nTESTING MODE #1: '2d structure'\n"
python main.py --run_num=1 --mode='2d structure' --aptamerSeq='TAATGTTAATTG' --ligand='False' --ligandType='' --ligandSeq='' > run1.out
printf "\nEND OF TEST #1: should see info of 2d structure in \"record.txt\" of working dir \"run1\"\n"
printf "===================================="

# printf "\nTESTING MODE #2: '3d coarse'\n"
# python main.py --run_num=2 --mode='3d coarse' --aptamerSeq='TAATGTTAATTG' --ligand='False' --ligandType='' --ligandSeq='' > run2.out
# printf "\nEND OF TEST #2: should see info of 2d structure + MMB folded structure in \"run2\"\n"
# printf "===================================="

printf "\nTESTING MODE #3: '3d smooth'\n"
python main.py --run_num=3 --mode='3d smooth' --aptamerSeq='TAATGTTAATTG' --ligand='False' --ligandType='' --ligandSeq='' > run3.out
printf "\nEND OF TEST #3: should see info of 2d structure + MMB folded structure + short MD relaxation in \"run3\"\n"
printf "===================================="

printf "\nTESTING MODE #4: 'coarse dock'\n"
python main.py --run_num=4 --mode='coarse dock' --aptamerSeq='TAATGTTAATTG' --ligand='YQTQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR' > run4.out
printf "\nEND OF TEST #4: should see info of 2d structure + MMB folded structure + docking by provided target ligand in \"run4\"\n"
printf "===================================="

printf "\nTESTING MODE #5: 'smooth dock'\n"
python main.py --run_num=5 --mode='smooth dock' --aptamerSeq='TAATGTTAATTG' --ligand='YQTQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR' > run5.out
printf "\nEND OF TEST #5: should see info of 2d structure + MMB folded structure + short MD relaxation + docked by provided target ligand in \"run5\"\n"
printf "===================================="

# printf "\nTESTING MODE #6: 'free aptamer'\n"
# python main.py --run_num=6 --mode='free aptamer' --aptamerSeq='TAATGTTAATTG' --ligand='False' --ligandType='' --ligandSeq='' > run6.out
# printf "\nEND OF TEST #6: should see info of 2d structure + MMB folded structure + long MD sampling in \"run6\"\n"
# printf "===================================="

# printf "\nTESTING MODE #7: 'full dock'\n"
# python main.py --run_num=7 --mode='full dock' --aptamerSeq='TAATGTTAATTG' --ligand='YQTQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR' > run7.out
# printf "\nEND OF TEST #7: should see info of 2d structure + MMB folded structure + long MD sampling + docked by provided target ligand in \"run7\"\n"
# printf "===================================="

# printf "\nTESTING MODE #8: 'full binding'\n"
# python main.py --run_num=8 --mode='full binding' --aptamerSeq='TAATGTTAATTG' --ligand='YQTQ.pdb' --ligandType='peptide' --ligandSeq='YQTQTNSPRRAR' > run8.out
# printf "\nEND OF TEST #8: should see info of 2d structure + MMB folded structure + long MD sampling + docked by provided target ligand + long MD sampling in \"run8\"\n"
# printf "====================================\n"
