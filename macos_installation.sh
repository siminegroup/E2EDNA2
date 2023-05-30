#!/bin/bash

echo -e "#########################################################################################################\n"
echo -e "  Welcome to the ultimate DNA simulator E2EDNA!"
echo -e "  This script is to create a conda virtual environment and install necessary Python packages for E2EDNA\n"
echo -e "#########################################################################################################\n"

echo -e "\nCREATING E2EDNA ENV"
conda env create -f e2edna-env.yml

echo -e "\nACTIVATING E2EDNA ENV"
conda activate e2edna
# source activate ~/miniconda3/envs/e2edna

echo -e "\n\nBefore ending this installation:"
echo -e "\nNote1: If the e2edna environment was not successfully activated,
	which means the string '(e2edna)' is not shown at the beginning of your terminal prompt,
	run this activation command in the command line:"
echo -e "               $ source activate <path_to_e2edna_conda_environment>\n"
echo -e "	To find the path above, list all conda environments and their paths on your computer by running:"
echo -e "               $ conda info -e\n"

echo -e "\nNote2: If you wish to start E2EDNA pipeline from DNA aptamer sequence rather than its 3D structure,"
echo -e "	please download MMB from https://simtk.org/projects/rnatoolbox, respectively.\n"
echo -e "	Then copy or move the downloaded MMB folder to the codebase directory (containing scripts like main.py), 
	and fill 'MMB-related paths' in the yaml configuration file.\n"

echo -e "\nGOOD LUCK!\n"
