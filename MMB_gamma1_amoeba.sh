#!/bin/bash
#SBATCH --account=def-simine
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=0-02:00
#SBATCH --mail-user=tao.liu7@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --error="error-%j_gamma1_amoeba.out"
#SBATCH --output="slurm-%j_gamma1_amoeba.out"

# put your MMB path here. Don't use ~
export LD_LIBRARY_PATH="/home/taoliu/projects/def-simine/programs/MMB/Installer.2_14.Linux64"
# # More careful:
#if [ -z "$LD_LIBRARY_PATH" ]; then
#	export LD_LIBRARY_PATH="/home/taoliu/projects/def-simine/programs/MMB/Installer.2_14.Linux64"
#else
#	export LD_LIBRARY_PATH="/home/taoliu/projects/def-simine/programs/MMB/Installer.2_14.Linux64:$LD_LIBRARY_PATH"
#fi

module load python/3.8.10
module load StdEnv/2020  gcc/9.3.0  cuda/11.4  openmpi/4.0.3
module load openmm/7.7.0  # To use AMOEBA2018 FF

# To load AmberTools21 (from https://docs.computecanada.ca/wiki/AMBER) and use tleap to generate .prmtop file
#module load StdEnv/2020 gcc/9.3.0 cuda/11.4 openmpi/4.0.3 
#module load ambertools/21
#source $EBROOTAMBERTOOLS/amber.sh
# Or load AmberTools20 (GPU version):
#module load StdEnv/2020  gcc/9.3.0  cuda/11.4  openmpi/4.0.3
module load amber/20.12-20.15
source $EBROOTAMBER/amber.sh
#The two options above require Python/3.8; if have to use Python/3.7, use the one below and downgrade chosen version of OpenMM to 7.5.0 (due to cuda/11.0)
#module load StdEnv/2020  gcc/9.3.0  cuda/11.0  openmpi/4.0.3
#module load amber/20.9-20.15

virtualenv --no-download $SLURM_TMPDIR/env  # need to load a Python module before running virtualenv
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt
pip install --no-index MDAnalysis==2.0.0
pip install --no-index lightdock==0.9.0
pip install --no-index pdbfixer==1.8 #"import openmm...." -> compatible with OpenMM/7.7.0; In pdbfixer/1.7: import simtk.openmm.... -> compatible with OpenMM/7.5.0
pip install --no-index nupack==4.0.0.27

#pip install --no-index numpy==1.20.2
#pip install /cvmfs/soft.computecanada.ca/custom/python/wheelhouse/gentoo/generic/numpy-1.20.2+computecanada-cp37-cp37m-linux_x86_64.whl
#ambertools21 and ambertools20 seems to contain numpy already: numpy-1.21.2
#pip install --no-index --upgrade numpy
# mdtraj is no longer used in utils and replaced by MDAnalysis (2022-04-05).
# Without this upgrade, numpy tends to give error for some packages such as mdtraj on beluga (2021-10-08)

pip install ~/projects/def-simine/programs/opendnaWheels/seqfold-0.7.7-py3-none-any.whl
pip install ~/projects/def-simine/programs/opendnaWheels/PeptideBuilder-1.1.0-py3-none-any.whl 
# https://pypi.org/project/PeptideBuilder/1.1.0/
# https://www.wheelodex.org/projects/peptidebuilder/
#pip install ~/projects/def-simine/programs/opendnaWheels/nupack-4.0.0.23-cp37-cp37m-linux_x86_64.whl
#python -c "import simtk.openmm" #no longer needed for OpenMM/7.7.0

date -u
# friction = 1.0
python ./main.py --run_num=1 --mode='free aptamer' --aptamerSeq='ACCTGGGGGAGTATTGCGGAGGAAGGT' --ligand='False' --ligandType='' --ligandSeq='' --friction=1.0
# python ./main.py --run_num=1 --sequence='AACGCTTTTTGCGTA' --peptide='A' --walltime=10 --temperature=310 --pH=7.4 --ionicStrength=0.163 --Mg=0.05 --impSolv='HCT'

date -u
