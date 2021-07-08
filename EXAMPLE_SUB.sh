#!/bin/bash
#SBATCH --account=ACCOUNT
#SBATCH --gres=gpu:v100:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0-4:00
#SBATCH --mail-user=YOUR EMAIL
#SBATCH --mail-type=END

export LD_LIBRARY_PATH= # put your MMB path here

module load python/3.7
module load StdEnv/2020
module load gcc/9.3.0
module load openmpi/4.0.3
module load cuda/11.0
module load python openmm

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index MDAnalysis
pip install --no-index lightdock
pip install --no-index -r requirements.txt
pip install numpy==1.20.2
pip install ~/wheels/seqfold-0.7.7-py3-none-any.whl
pip install ~/wheels/PeptideBuilder-1.1.0-py3-none-any.whl 
python -c "import simtk.openmm"

python ./main.py --run_num=0 --sequence='AACGCTTTTTGCGTT' --peptide='YQTQTNSPRRAR'
