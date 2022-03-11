#!/bin/bash

echo Welcome to the ultimate DNA simulator E2EDNA!

echo CREATING E2EDNA ENV
conda env create -f e2edna-env.yml

echo ACTIVATING E2EDNA ENV
conda activate e2edna
# source activate ~/miniconda3/envs/e2edna

echo INSTALLING NUPACK
python -m pip install -U nupack -f ~/Downloads/nupack-4.0.0.27/package

echo THE FINAL STEP IS TO DOWNLOAD MMB from https://simtk.org/projects/rnatoolbox and set it up
echo GOOD LUCK!
