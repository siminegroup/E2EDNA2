#!/bin/bash
# Installation for a local run

# Create a conda environment ``myDNAEnv`` with python 3.7.11
conda create -n myDNAEnv python=3.7.11
# Activate the conda envs (using miniconda3 as an example)
source activate ~/miniconda3/envs/myDNAEnv

# Make sure to install all the required packages and their dependence
pip install --upgrade pip
pip install -r requirements.txt

# MDAnalysis: https://www.mdanalysis.org/pages/installation_quick_start/
conda config --add channels conda-forge
conda install mdanalysis

# OpenMM: http://docs.openmm.org/7.5.0/userguide/application.html#installing-openmm
conda install -c conda-forge openmm

# AmberTools21: https://ambermd.org/GetAmber.php#ambertools 
conda install -c conda-forge ambertools=21 compilers

# LightDock: # https://lightdock.org/
pip install lightdock  

# Seqfold: https://pypi.org/project/seqfold/
pip install seqfold

# NUPACK: download wheel from here: http://www.nupack.org/downloads
pip install <downloaded wheel for your OS>

# PeptideBuilder: # https://pypi.org/project/PeptideBuilder/1.1.0/ # https://www.wheelodex.org/projects/peptidebuilder/
pip install PeptideBuilder

# MMB: see the README on the main branch.
