#!/bin/bash
# For local run 

# Make sure to install all the required packages and their dependence
pip install --upgrade pip
pip install -r requirements.txt

# MDAnalysis: https://www.mdanalysis.org/pages/installation_quick_start/
conda config --add channels conda-forge
conda install mdanalysis
# LightDock: # https://lightdock.org/
pip install lightdock  

# Seqfold: https://pypi.org/project/seqfold/
pip install seqfold
# NUPACK: download wheel from here: http://www.nupack.org/downloads
pip install <downloaded wheel for your OS>
# PeptideBuilder: # https://pypi.org/project/PeptideBuilder/1.1.0/ # https://www.wheelodex.org/projects/peptidebuilder/
pip install PeptideBuilder

date -u
python ./main.py
date -u
