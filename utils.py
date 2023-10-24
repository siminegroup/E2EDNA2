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
import os, re
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
from shutil import copyfile, copytree

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


def prepPDB_with_Mg(file, boxOffset, MgCl2_molar, NaCl_molar, MMBCORRECTION=False, waterBox=True):
    """
    Soak pdb file in water box, and add Mg2+ ions
    :param file: file to be processed
    :param boxOffset: in unit of nanometer
    :param NaCl_molar: in unit of molar
    :param MgCl2_molar: in unit of molar
    :param MMBCORRECTION: if the input pdb file is an MMB output, we need to apply a correction, since MMB is a little weird formatting-wise
                           https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=359&t=13397&p=0&start=0&view=&sid=bc6c1b9005122914ec7d572999ba945b
    :param waterBox:
    :return:
    """

    if MMBCORRECTION:
        replaceText(file, '*', "'")  # due to a bug in this version of MMB - structures are encoded improperly - this fixes it

    fixer = PDBFixer(filename=file)
    padding, boxSize, boxVectors = None, None, None
    boxMode = 'cubic'  # TODO toggle for box type - look at openmm-setup source code for other box types

    if boxMode == 'cubic':
        padding = float(boxOffset) * unit.nanometer  # for cubic box
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

    # fixer.findMissingResidues()
    # fixer.findMissingAtoms()
    # fixer.addMissingAtoms()  # may need to optimize bonding here
    # fixer.addMissingHydrogens(pH=pH)  # add missing hydrogens, but not really working for theophylline or DNA though even at pH 1.0 lol
    
    # save the CONECT records, if any, as they will be removed once solvent is added
    my_file = open(file, 'r')
    my_file_lines = my_file.readlines()
    my_file.close()
    conect_records = []
    for line in my_file_lines:
        if line.find('CONECT') == 0: conect_records.append(line) # if finding 'CONECT' as the beginning of the line

    has_conect_records = True if len(conect_records) > 0 else False
    if has_conect_records:
        conect_records_string = ''.join(conect_records)
        if conect_records_string[-1] == '\n': # if any, remove the last '\n' because need to append an 'End' 
            conect_records_string = conect_records_string[:-1]
        # with open('CONECT_records.txt', 'w') as conect_records_file:
        #     conect_records_file.write(conect_records_string)
        # Q: when there is a ligand and a "TER" between DNA and ligand residues, does it offset CONECT records by 1??
        # Good news: MDAnalysis removes all the TER records, so we will add CONECT after using MDAnalysis to center the system


    if waterBox == True:
        MgCl2_concentration = float(MgCl2_molar) * unit.molar
        NaCl_concentration = float(NaCl_molar) * unit.molar
        positiveIon = 'Na+'  # params['positiveion']+'+'
        negativeIon = 'Cl-'  # params['negativeion']+'-'

        # fixer.addSolvent(boxSize, padding, boxVectors, positiveIon, negativeIon, ionicStrength)
            # ''' PDBFixer actually calls openmm.app.modeller.Modeller class: https://github.com/openmm/pdbfixer/blob/master/pdbfixer/pdbfixer.py#L1061
            # Modeller.addSolvent() needs a forcefield to determine van der Waals radii and atomic charges. The default forcefield used in PDBfixer is amber14 and amber14/tip3p
            # If using another water model, could choose to use Modeller.addSolvent() directly;
            # Or, if the water model is similar to tip3p, then the downstream local energy minimization can correct for the small differences between the models. 
            # This is how they explained in the Modeller API why they picked tip3p as the default water forcefield
            # '''
        # Add NaCl to the same [MgCl2], so that later we can replace all Na+ with Mg2+
        # here, add ions without neutralizing the solute!
        modeller = Modeller(fixer.topology, fixer.positions)
        fixer_forcefield = fixer._createForceField(fixer.topology, water=True) # default: using amber14-all and amber14/tip3p models
        # if system contains ligand, maybe add GAFFTemplateGenerator to the forcefield as usual, then ligand can be recognized as well.
        modeller.addSolvent(fixer_forcefield, padding=padding, positiveIon=positiveIon, negativeIon=negativeIon, ionicStrength=MgCl2_concentration, neutralize=False)
        chains = list(modeller.topology.chains())
        if len(chains) == 1:
            chains[0].id = 'A'
        else:
            chains[-1].id = chr(ord(chains[-2].id)+1)
        # fixer.topology = modeller.topology
        # fixer.positions = modeller.positions

        modeller.deleteWater() # remove all water
        modeller.delete(atom for atom in modeller.topology.atoms() if atom.name == 'Cl') # remove all Cl-
        # Now the only ions left are Na+, same concentraion as desired Mg2+:
        system_Na_dry_pdb = file.split('.pdb')[0] + '_Na_dry.pdb'
        PDBFile.writeFile(modeller.topology, modeller.positions, open(system_Na_dry_pdb, 'w'))
        # Replace all Na+ with Mg2+ in pdb file: here assuming "Na" and "NA" only occur to sodium records. 
        # A possible exception may be the ligand residue name: it may contain "NA", then need to replace the whole word
        # But, still need to be mindful of this step
        replaceText(system_Na_dry_pdb, 'Na', 'Mg', swap_whole_word=False)
        replaceText(system_Na_dry_pdb, 'NA', 'MG', swap_whole_word=True)
        # change the pdb filename 
        system_Mg_dry_pdb = file.split('.pdb')[0] + '_Mg_dry.pdb'
        os.rename(system_Na_dry_pdb, system_Mg_dry_pdb)

        if has_conect_records:
            # Add back the CONECT records
            f = open(system_Mg_dry_pdb, 'r')
            string = f.read()
            f.close()
            # remove 'END'
            end_loc = string.find('END') # find where the 'END' appears
            # replace everything after END with the CONECT records
            new_string = string[:end_loc] + conect_records_string + '\nEND\n'
            with open(system_Mg_dry_pdb, 'w') as new_file:
                new_file.write(new_string)

        # add NaCl
        fixer = PDBFixer(system_Mg_dry_pdb) # recall the system now carries several Mg2+
        fixer.addSolvent(positiveIon=positiveIon, negativeIon=negativeIon, ionicStrength=NaCl_concentration)
        # here we do not specify box size at all, in order to use the existing topology's box vectors, so that the box size does not change
        # Some Na+ and/or Cl- will be added first to neutralize the system including Mg2+, then (Na+, Cl-) pairs will be added to give the requested [NaCl]

    solvated_file = file.split('.pdb')[0] + '_processed_not_centered.pdb'
    PDBFile.writeFile(fixer.topology, fixer.positions, open(solvated_file, 'w'))

    # center the solvated molecule in the simulation box
    u = mda.Universe(solvated_file)
    dim = u.dimensions[:3] # retrieve the box dimension (the last three are the box angles)
    box_center = dim / 2
    all_atoms = u.select_atoms("all")
    # print(all_atoms.center_of_mass(pbc=True)) # translate the molecules so that their centers are within the same periodic box
    # print(all_atoms.center_of_mass(pbc=False))
    COM = all_atoms.center_of_mass(pbc=False)
    u.atoms.translate(box_center - COM)
    centered_file = file.split('.pdb')[0] + '_processed.pdb'
    all_atoms.write(centered_file)

    # Good news: MDAnalysis removes all the TER records, so we will add CONECT after using MDAnalysis to center the system
    if has_conect_records:
        # Add back the CONECT records
        f = open(centered_file, 'r')
        string = f.read()
        f.close()
        # remove 'END'
        end_loc = string.find('END') # find where the 'END' appears
        # replace everything after END with the CONECT records
        new_string = string[:end_loc] + conect_records_string + '\nEND\n'
        with open(centered_file, 'w') as new_file:
            new_file.write(new_string)

    
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
    geompadding = float(boxOffset) * unit.nanometer

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

    # fixer.findMissingResidues()
    # fixer.findMissingAtoms()
    # fixer.addMissingAtoms()  # may need to optimize bonding here

    # fixer.addMissingHydrogens(pH=pH)  # add missing hydrogens, but not really working for theophylline or DNA though even at pH 1.0 lol
    
    if waterBox == True:
        ionicStrength = float(ionicStrength) * unit.molar
        positiveIon = 'Na+'  # params['positiveion']+'+'
        negativeIon = 'Cl-'  # params['negativeion']+'-'
        fixer.addSolvent(boxSize, padding, boxVectors, positiveIon, negativeIon, ionicStrength)
        ''' PDBFixer actually calls openmm.app.modeller.Modeller class: https://github.com/openmm/pdbfixer/blob/master/pdbfixer/pdbfixer.py#L1061
        Modeller.addSolvent() needs a forcefield to determine van der Waals radii and atomic charges. The default forcefield used in PDBfixer is amber14 and amber14/tip3p
        If using another water model, could choose to use Modeller.addSolvent() directly;
        Or, if the water model is similar to tip3p, then the downstream local energy minimization can correct for the small differences between the models. 
        This is how they explained in the Modeller API why they picked tip3p as the default water forcefield
        '''
    
    solvated_file = file.split('.pdb')[0] + '_processed_not_centered.pdb'
    PDBFile.writeFile(fixer.topology, fixer.positions, open(solvated_file, 'w'))

    # center the solvated molecule in the simulation box
    u = mda.Universe(solvated_file)
    dim = u.dimensions[:3] # retrieve the box dimension (the last three are the box angles)
    box_center = dim / 2
    all_atoms = u.select_atoms("all")
    # print(all_atoms.center_of_mass(pbc=True)) # translate the molecules so that their centers are within the same periodic box
    # print(all_atoms.center_of_mass(pbc=False))
    COM = all_atoms.center_of_mass(pbc=False)
    u.atoms.translate(box_center - COM)
    centered_file = file.split('.pdb')[0] + '_processed.pdb'
    all_atoms.write(centered_file)
    

def replaceText(file, old_string, new_string, swap_whole_word=False):
    # search and replace string in a text file, then save the new version
    f = open(file, 'r')
    text = f.read() # as string
    f.close()
    if swap_whole_word:
        text = re.sub(r"\b%s\b" % old_string, new_string, text)
    else:
        text = text.replace(old_string, new_string)
        # text = re.sub(old_string, new_string, text)  # another way
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
    # if u.segments.n_segments > 2:  # if > 2 segments, then there must be solvent and salts (assuming nonzero salt concentration)
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
