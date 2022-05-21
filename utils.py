'''
E2EDNA 2.0 - OpenMM Implementation of E2EDNA !

An automated pipeline for simulating DNA aptamers complexed with target ligands (peptide, DNA, RNA or small molecules).
    
Copyright (C) 2022 Michael Kilgour, Tao Liu and Lena Simine

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

"""
Utitilies
"""
from openmm.app import *
import openmm.unit as unit

import argparse
import os
import csv
import numpy as np
import time

import Bio.PDB  # biopython
# import mdtraj as md  # mdtraj calls "simtk.openmm", hence "Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead."
import MDAnalysis as mda

from pdbfixer import PDBFixer # pdbfixer package has been installed by conda
from collections import Counter

from PeptideBuilder import Geometry  # pip install PeptideBuilder
import PeptideBuilder


# def recenterDCD(topology, trajectory):
#     """
#     topology as pdb
#     trajectory as dcd
#     creates a new dcd without periodic artifacts
#     """
#     traj = md.load(trajectory, top = topology)
#     traj.image_molecules()
#     traj.save(trajectory.split('.')[0] + '_recentered.dcd')
# # MDAnalysis provides something similar: https://userguide.mdanalysis.org/1.1.1/examples/transformations/center_protein_in_box.html
# # In the future, implement it using MDAnalysis, because mdtraj-1.9.7 still calls "simtk.openmm" which is a deprecated package.
# # In fact, MDAnalysis also interacts with openmm package, unfortunately, also by calling "simtk.openmm". So the depreceation warning persists but not shows up in the beginning of the command line output.


class Timer:
    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start


def printRecord(statement, filepath=""):
    """
    print a string to command line output and a text file
    :param statement:
    :return:
    """
    print(statement)
    if os.path.exists(filepath + 'run_output_log.txt'):
        with open(filepath + 'run_output_log.txt', 'a') as file:
            file.write('\n' + statement)
    else:
        with open(filepath + 'run_output_log.txt', 'w') as file:
            file.write('\n' + statement)


def prepPDB(file, boxOffset, pH, ionicStrength, MMBCORRECTION=False, waterBox=True):
    """
    Soak pdb file in water box
    :param file:
    :param boxOffset:
    :param pH:
    :param ionicStrength:
    :param MMBCORRECTION: if the input pdb file is an MMB output, we need to apply a correction, since MMB is a little weird formatting-wise
                           https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=359&t=13397&p=0&start=0&view=&sid=bc6c1b9005122914ec7d572999ba945b
    :param waterBox:
    :return:
    """

    if MMBCORRECTION:
        replaceText(file, '*', "'")  # due to a bug in this version of MMB - structures are encoded improperly - this fixes it

    fixer = PDBFixer(filename=file)
    padding, boxSize, boxVectors = None, None, None
    geompadding = float(boxOffset) * unit.nanometer  # TODO what does unit.nanometer do

    boxMode = 'cubic'  # TODO toggle for box type - look at openmm-setup source code for other box types

    if boxMode == 'cubic':
        padding = geompadding  # for cubic box
    elif boxMode == 'rectangular prism':
        # or we can make a rectangular prism which (maybe) cuts off sides of the cube
        u = mda.Universe(file)
        coords = u.atoms.positions
        xrange = np.ptp(coords[:, 0])  # get maximum dimension
        yrange = np.ptp(coords[:, 1])
        zrange = np.ptp(coords[:, 2])
        maxsize = max([xrange, yrange, zrange])
        # minimum dimension is half the longest TODO why?
        # also convert to nanometer
        xrange = max([xrange, maxsize / 2]) / 10
        yrange = max([yrange, maxsize / 2]) / 10
        zrange = max([zrange, maxsize / 2]) / 10

        # TODO may also need an EWALD offset
        xrange = xrange + 2 * boxOffset
        yrange = yrange + 2 * boxOffset
        zrange = zrange + 2 * boxOffset

        boxSize = [xrange, yrange, zrange] * unit.nanometer  # for rectangular prism

    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()  # may need to optimize bonding here

    fixer.addMissingHydrogens(pH=pH)  # add missing hydrogens
    # TODO: Don't understand these above yet 
    
    if waterBox == True:
        ionicStrength = float(ionicStrength) * unit.molar
        positiveIon = 'Na+'  # params['positiveion']+'+'
        negativeIon = 'Cl-'  # params['negativeion']+'-'
        fixer.addSolvent(boxSize, padding, boxVectors, positiveIon, negativeIon, ionicStrength)
    
    PDBFile.writeFile(fixer.topology, fixer.positions, open(file.split('.pdb')[0] + '_processed.pdb', 'w'))


def replaceText(file, old_string, new_string):
    # search and replace string in a text file, then save the new version
    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.replace(old_string, new_string)
    f = open(file, 'w')
    f.write(text)
    f.close()


def findLine(file, string):
    """
    return the line number of a given string in a text file, indexing from 1
    :param file:
    :param string:
    :return lineNum:
    """
    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.split('\n')
    try:
        lineInd = 1
        for line in text:
            if line == string:
                lineNum = lineInd  # index from 1
            lineInd += 1
        return lineNum
    except:
        raise ValueError("String not found in the file!")


def addLine(file, string, line_number):
    """
    add a new line of text to a file
    :param file:
    :param string:
    :param line_number:
    :return:
    """
    f = open(file, 'r')
    text = f.read()
    f.close()
    lines = text.split('\n')
    lines.insert(line_number - 1, string)
    text = "\n".join(lines)
    f = open(file, 'w')
    f.write(text)
    f.close()


def cleanTrajectory(structure, trajectory):
    """
    Remove water, salt from trajectory
    :param structure:
    :param trajectory:
    :return:
    """
    structureName  = os.path.basename(structure)
    trajectoryName = os.path.basename(trajectory)
    structureDir  = os.path.dirname(structure)
    trajectoryDir = os.path.dirname(trajectory)
    
    # cleaned files are in the same directory as the uncleaned ones.
    clean_structure  = os.path.join(structureDir, 'clean_'+structureName)
    clean_trajectory = os.path.join(trajectoryDir, 'clean_'+trajectoryName)

    u = mda.Universe(structure, trajectory)
    # TODO: if u.segments.n_segments > 2:  # if > 2 segments, then there must be solvent and salts (assuming nonzero salt concentration)
    goodStuff = u.segments[:-2].atoms  # cut out salts and solvents

    goodStuff.write(clean_structure)  # write topology
    with mda.Writer(clean_trajectory, goodStuff.n_atoms) as W:
        for ts in u.trajectory:  # indexing over the trajectory
            W.write(goodStuff)

def cleanPDB(structure):
    """
    Remove water, salt from PDB
    :param structure:    
    :return:
    """
    u = mda.Universe(structure)    
    goodStuff = u.segments[:-2].atoms  # cut out salts and solvents
    goodStuff.write("clean_" + structure)  # write topology

def extractFrame(structure, trajectory, frame, outFileName):
    """
    Save a given trajectory frame as a separate pdb file
    :param structure:
    :param trajectory:
    :param frame:
    :param outFileName:
    :return:
    """
    u = mda.Universe(structure, trajectory)  # load up the whole trajectory
    u.trajectory[frame]  # Michael: this indexes the trajectory up to the desired frame (weird syntax, I think)  # Me: this takes out the desired frame
    if u.segments.n_segments > 2:  # if there are more than 2 segments, then there must be solvent and salts (assuming nonzero salt concentration)
        atoms = u.segments[:-2].atoms  # omit solvent and salts
    else:
        atoms = u.atoms
    atoms.write(outFileName)
    if frame == -1:
        printRecord('MDAnalysis: save the last frame into: {}'.format(outFileName))


def mode(lst):
    """
    Returns the statistical mode(s) of a given list, only used for checking backbone_dihedrals
    Source: https://stackabuse.com/calculating-mean-median-and-mode-in-python/
    :param lst:
    :return mode_list, a list of the mode(s) found in the input list:
    """
    c = Counter(lst)
    mode_list = [k for k, v in c.items() if v == c.most_common(1)[0][1]]
    return mode_list


def killH(structure):
    """
    Use simTk modeller to delete all atoms with "H" in the name
    :param structure:
    :return:
    """
    pdb = PDBFile(structure)
    topology = pdb.topology
    positions = pdb.positions
    modeller = Modeller(topology, positions)
    modeller.delete(atom for atom in topology.atoms() if "H" in atom.name)
    PDBFile.writeFile(modeller.topology, modeller.positions, open(structure.split('.')[0] + '_noH.pdb', 'w'))


def addH(structure, pH):
    """
    Protonate a given structure
    :param structure:
    :param pH:
    :return:
    """
    pdb = PDBFile(structure)
    topology = pdb.topology
    positions = pdb.positions
    modeller = Modeller(topology, positions)
    modeller.addHydrogens(pH=pH)
    PDBFile.writeFile(modeller.topology, modeller.positions, open(structure.split('.')[0] + '_H.pdb', 'w'))


def changeSegment(structure, oldSeg, newSeg):
    """
    Change the segment ID for all molecule(s) in a pdb file
    :param structure:
    :param oldSeg:
    :param newSeg:
    :return:
    """
    replaceText(structure, ' ' + oldSeg + ' ', ' ' + newSeg + ' ')


def readInitialLines(file, lines):
    # Return the first N lines of a text file
    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.split('\n')
    initialLines = text[:lines]

    return initialLines


def readFinalLines(file, lines):
    # Return the final N lines of a text file
    f = open(file, 'r')
    text = f.read()
    f.close()
    text = text.split('\n')
    finalLines = text[-(lines + 1):]  # TODO: this will return "lines+1" rows from the end of the file

    return finalLines


def appendTrajectory(topology, original, new):
    """
    Use mda to combine old and new MD trajectories into one nice video
    """
    trajectories = [original, new]
    u = mda.Universe(topology, trajectories)
    with mda.Writer('combinedTraj.dcd', u.atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(u)


def removeLine(file, string):
    """
    remove every line containing given string from a file
    :param file: the file
    :param string: the string to look for
    :return:
    """
    f = open(file, 'r')
    text = f.read()
    f.close()
    lines = text.split('\n')
    try:
        for i in range(len(lines)):
            if string in lines[i]:
                lines.pop(i)
    except:
        pass
    text = "\n".join(lines)
    f = open(file, 'w')
    f.write(text)
    f.close()


def appendLine(file, string):
    # add a new line of end of a file
    f = open(file, 'r')
    text = f.read()
    f.close()
    lines = text.split('\n')
    lines.insert(len(lines) - 1, string)  # replace line ### with the string, indexing from 1
    text = "\n".join(lines)
    f = open(file, 'w')
    f.write(text)
    f.close()
