"""
Unused imports:
import os
import re
import numpy as np
import MDAnalysis as mda
"""
import sys
import glob
from shutil import copyfile, copytree

import interfaces
from utils import *
from analysisTools import *

# noinspection PyPep8Naming  #meaning??
'''
To-do:
1. Implicit solvent    
    -- add more options to the implicit solvent branch: let user specify more detail of this mode    
    -- Other than: MDSmoothing, autoMD, omm class: need to specify implicit solvent in other methods? e.g., analyzeTrajectory
Doing:

Done:
1. skip MMB locally.
'''


class opendna:
    def __init__(self, params):
        self.workDir = ""  # rather use an empty string, not an empty list
        self.params = params
        self.sequence = self.params['sequence']
        self.peptide = self.params['peptide']
        self.pdbDict = {}
        self.dcdDict = {}  # output file format is .dcd???

        self.actionDict = {}
        self.getActionDict()  # specify the actions based on the selected mode
        if self.actionDict['make workdir']:
            self.setup()  # if we don't need a workdir & MMB files (eg, give a 3D structure), don't make one.

        self.i = int(-1)
        self.j = int(-1)

    def getActionDict(self):
        """
        generate a binary sequence of 'to-do's', given user's input
        :return:
        """
        '''
        Modes, in order of increasing cost
        '2d structure': ssString, pair list and probability
        '3d coarse': MMB output, stressed structure, no solvent
        '3d smooth': MMB output with short MD relaxation
        'coarse dock': best docking scores on coarse MMB structure
        'smooth dock': best docking scores on smoothed MMB structure
        'free aptamer': evaluate and find representative 3D aptamer structure
        'full docking': 'free aptamer' + docking
        'full binding': 'full docking' + binding
        '''

        if self.params['mode'] == '2d structure':
            self.actionDict['make workdir'] = False
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = False
            self.actionDict['do smoothing'] = False
            self.actionDict['get equil repStructure'] = False
            self.actionDict['do docking'] = False
            self.actionDict['do binding'] = False

        elif self.params['mode'] == '3d coarse':
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = False
            self.actionDict['get equil repStructure'] = False
            self.actionDict['do docking'] = False
            self.actionDict['do binding'] = False

        elif self.params['mode'] == '3d smooth':
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = True
            self.actionDict['get equil repStructure'] = False
            self.actionDict['do docking'] = False
            self.actionDict['do binding'] = False

        elif self.params['mode'] == 'coarse dock':  # dock the structure out of MMB. Is it useful by any chance?
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = False
            self.actionDict['get equil repStructure'] = False
            self.actionDict['do docking'] = True
            self.actionDict['do binding'] = False

        elif self.params['mode'] == 'smooth dock':  # no MD sampling for equilibrated structure. Useful?
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = True
            self.actionDict['get equil repStructure'] = False
            self.actionDict['do docking'] = True
            self.actionDict['do binding'] = False

        elif self.params['mode'] == 'free aptamer':
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = True
            self.actionDict['get equil repStructure'] = True
            self.actionDict['do docking'] = False
            self.actionDict['do binding'] = False

        elif self.params['mode'] == 'full docking':
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = True
            self.actionDict['get equil repStructure'] = True
            self.actionDict['do docking'] = True
            self.actionDict['do binding'] = False

        elif self.params['mode'] == 'full binding':
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = True
            self.actionDict['get equil repStructure'] = True
            self.actionDict['do docking'] = True
            self.actionDict['do binding'] = True        

    def setup(self):
        """
        set up a working directory
        copy in relevant xyz and keyfiles
        move to the working directory
        :return:
        """

        if (self.params['explicit run enumeration'] is True) or (self.params['run num'] == 0):  # make a new directory
            if self.params['run num'] == 0:  # auto-increment the max run_num folder
                self.makeNewWorkingDirectory()
            else:  # create a working directory defined by the user
                self.workDir = self.params['workdir'] + '/run%d' % self.params['run num']
                os.mkdir(self.workDir)
        else:
            self.workDir = self.params['workdir'] + '/' + 'run%d' % self.params['run num']

        # copy relevant files to the workDir
        os.mkdir(self.workDir + '/outfiles')
        if self.params['implicit solvent'] is True:
            copyfile(self.params['leap template'], self.workDir + '/leap_template.in')
        # copy structure files
        copyfile(self.params['analyte pdb'], self.workDir + '/analyte.pdb')  # ie, target of the aptamer


        # If we have a folded structure to start with, then skip MMB:
        if self.params['skip MMB'] is True:  # otherwise, 'do 2d analysis' and 'do MMB' are always True
            self.actionDict['do 2d analysis'] = False  
            self.actionDict['do MMB'] = False
            printRecord('Skipping 2D analysis and MMB folding. Starting with an folded structure.', self.workDir+'/')
            
            if not os.path.isfile('./' + self.params['folded initial structure']):
            # if not os.path.exists(self.params['folded initial structure']):
                printRecord('Designated folded structure does not exist! Terminating the pipeline.', self.workDir+'/')
                self.terminateRun()  # maybe do it via "raise an error"?
            else:
                printRecord('Copying the folded structure: {} to workdir.'.format(self.params['folded initial structure']), self.workDir+'/')
                copyfile(self.params['folded initial structure'], self.workDir + '/foldedSequence_0.pdb')
        else:
            # copy MMB params and command script
            copyfile(self.params['mmb params'], self.workDir + '/parameters.csv')
            copyfile(self.params['mmb normal template'], self.workDir + '/commands.template.dat')
            copyfile(self.params['mmb quick template'], self.workDir + '/commands.template_quick.dat')
            copyfile(self.params['mmb long template'], self.workDir + '/commands.template_long.dat')        
               
        if self.params['pick up from chk'] is True:
            # everything before MD sampling is skipped
            self.actionDict['do 2d analysis'] = False  
            self.actionDict['do MMB'] = False            
            self.actionDict['do smoothing'] = False
        
            printRecord('Copying the checkpoint file: {} to workdir.'.format(self.params['chk file']), self.workDir+'/')
            if os.path.exists(self.params['chk file']):
                copyfile(self.params['chk file'], self.workDir + '/' + self.params['chk file'])
            else:
                printRecord('Designated checkpoint file does not exist! Terminating the pipeline.', self.workDir+'/')
                self.terminateRun()

        # copy lightdock scripts
        if self.actionDict['do docking'] is True:
            copytree('lib/lightdock', self.workDir + '/ld_scripts')  # if destination dir does not already exist, create one then copy the whole source dir

        # copy csv file (dihedral restraints) if restraints are turned on
        if self.params['peptide backbone constraint constant'] != 0:
            copyfile('backbone_dihedrals.csv', self.workDir + '/backbone_dihedrals.csv')

        # move to working dir
        os.chdir(self.workDir)
        printRecord('Simulating {} with {}'.format(self.sequence, self.peptide))

        if (self.params['skip MMB'] is False) and (self.params['device'] == 'local'):
            if self.params['local device platform'] == 'linux':
                os.environ["LD_LIBRARY_PATH"] = self.params['mmb dir']  
                # Michael: export the path to the MMB library - only necessary in WSL environments. No need on Windows machine
            else:
                pass 
              # Tao: do NOT use os.environ["DYLD_LIBRARY_PATH"] for macos machine! It wouldn't help MMB locate the library AND it would confuse OpenMM python package.
              # Tao: see the MMB class on how to run MMB on macos machine.

    def makeNewWorkingDirectory(self):
        """
        make a new working directory: not overlapping previously existing directories
        :return:
        """
        workdirs = glob.glob(self.params['workdir'] + '/' + 'run*')  # check for prior working directories
        if len(workdirs) > 0:
            prev_runs = []
            for i in range(len(workdirs)):
                prev_runs.append(int(workdirs[i].split('run')[-1]))

            prev_max = max(prev_runs)
            self.workDir = self.params['workdir'] + '/run%d' % (prev_max + 1)
            os.mkdir(self.workDir)
            printRecord('Starting Fresh Run %d' % (prev_max + 1))
        else:
            self.workDir = self.params['workdir'] + '/' + 'run1'
            printRecord('Starting Fresh Run 1')
            os.mkdir(self.workDir)

    def run(self):
        """
        run the end-to-end simulation pipeline for the chosen mode
        consult checkpoints to not repeat prior steps
        :return:
        """

        # outputDict = {'params': self.params}
        outputDict = {}
        outputDict['params'] = self.params
        np.save('opendnaOutput', outputDict)  # Save an array to a binary file in NumPy ``.npy`` format.

        if self.actionDict['do 2d analysis']:   # get secondary structure
            self.pairLists = self.getSecondaryStructure(self.sequence)
            outputDict['2d analysis'] = self.ssAnalysis
            np.save('opendnaOutput', outputDict)  # save 2d structure results

            printRecord('Running over %d' % len(self.pairLists) + ' possible 2D structures')
            num_2dSS = len(self.pairLists)

        elif self.params['pick up from chk'] is True:  
            printRecord('Skip all the steps before MD sampling to resume previous sampling from .chk file')
            num_2dSS = 1  # quick and dirty

        else:  # just skipped nupack or seqfold
            printRecord('Starting with an existing folded strcuture.')
            num_2dSS = 1  # quick and dirty

        for self.i in range(num_2dSS):  # loop over all possible secondary structures
            
            if self.actionDict['do 2d analysis'] is True:  # self.ssAnalysis only exists if we "do 2d analysis"
                printRecord('2D structure #{} is {}'.format(self.i, self.ssAnalysis['2d string'][self.i]))
                self.pairList = np.asarray(self.pairLists[self.i])  # be careful!!!: .pairList vs. .pairLists

            if self.actionDict['do MMB']:  # fold 2D into 3D
                self.foldSequence(self.sequence, self.pairList)
            elif self.params['pick up from chk'] is False:  # start with a folded initial structure: skipped MMB but will MD smooth
                self.pdbDict['folded sequence {}'.format(self.i)] = self.params['folded initial structure']
                self.pdbDict['representative aptamer {}'.format(self.i)] = self.params['folded initial structure']  # in "coarse dock" mode.

            if self.actionDict['do smoothing']:
                if self.params['skip MMB'] is False:
                    self.MDSmoothing(self.pdbDict['mmb folded sequence {}'.format(self.i)], relaxationTime=self.params['smoothing time'], implicitSolvent=self.params['implicit solvent'])  # relax for xx nanoseconds
                else:
                    self.MDSmoothing(self.pdbDict['folded sequence {}'.format(self.i)], relaxationTime=self.params['smoothing time'], implicitSolvent=self.params['implicit solvent'])  # relax for xx nanoseconds

            if self.actionDict['get equil repStructure']:  # definitely did smoothing if want an equil structure
                outputDict['free aptamer results {}'.format(self.i)] = self.runFreeAptamer(self.pdbDict['relaxed sequence {}'.format(self.i)], implicitSolvent=self.params['implicit solvent'])
                np.save('opendnaOutput', outputDict)  # save outputs


            if self.actionDict['do docking'] and (self.peptide is not False):  # find docking configuration for the complexed structure
                # coarse dock: no smoothing;
                # smooth dock: smooth + dock
                # full dock: smooth + dock + equil structure
                # full binding: smooth + dock + equil structure + sampling dynamics
                outputDict['dock scores {}'.format(self.i)] = self.dock(self.pdbDict['representative aptamer {}'.format(self.i)], 'peptide.pdb')
                # pdbDict['representative aptamer {}' is defined at MMb folding, MD smoothing and runFreeAptamer
                # TODO: does lightdock also support Amber implicit solvent model?
                np.save('opendnaOutput', outputDict)  # save outputs

            # N docked structures are specified by user: how many docked structures do we want to investigate
            printRecord('Running over %d' % self.params['N docked structures'] + ' docked structures')

            for self.j in range(self.params['N docked structures']):  # loop over docking configurations for a given secondary structure
                printRecord('Docked structure #{}'.format(self.j))
                if self.actionDict['do binding']:  # run MD on the complexed structure
                    outputDict['binding results {} {}'.format(self.i, self.j)] = self.bindingDynamics(self.pdbDict['binding complex {} {}'.format(self.i, int(self.j))], implicitSolvent=self.params['implicit solvent'])
                    # TODO why need int(self.j)? Why sometimes %d % string, but sometimes {}.format?

                    np.save('opendnaOutput', outputDict)

        return outputDict

    # ======================================================================================
    # ======================================================================================
    # ====================== supporting functions are followed (part1) =====================
    # ======================================================================================
    # ======================================================================================
    def getSecondaryStructure(self, sequence):
        """
        get the secondary structure(s) for a given sequence
        using seqfold here - identical features are available in nupack, though results are sometimes different TODO: really??????
        :param sequence:
        :return:
        """
        printRecord("Getting Secondary Structure(s)")
        if self.params['secondary structure engine'] == 'seqfold':
            ssString, pairList = getSeqfoldStructure(sequence, self.params['temperature'])  # seqfold guess
            self.ssAnalysis = [ssString, pairList]
            # TODO what are in ssString and pairList?
            return pairList

        elif self.params['secondary structure engine'] == 'NUPACK':
            nup = interfaces.nupack(sequence, self.params['temperature'], self.params['ionic strength'], self.params['[Mg]'])  # initialize nupack
            self.ssAnalysis = nup.run()  # run nupack analysis of possible 2D structures.

            distances = getSecondaryStructureDistance(self.ssAnalysis['config'])
            if len(distances) > 1:
                topReps, topProbs, topDists = do2DAgglomerativeClustering(self.ssAnalysis['config'], self.ssAnalysis['state prob'], distances)
                topPairLists = []
                nLists = min(len(topReps), self.params['N 2D structures'])
                for i in range(nLists):
                    if (i == 0) or (np.amin(topDists[i, :i]) > 0.05): # only add a config to the list if it's at least 5% different from a previously accepted config
                        topPairLists.append(ssToList(configToString(topReps[i].astype(int))))
                return topPairLists
            else:  # if, for some reason, there's only a single plausible structure (unlikely)
                return self.ssAnalysis['pair list']

    def foldSequence(self, sequence, pairList):
        """
        generate a coarsely folded 3D structure for the aptamer
        :param sequence: the ssDNA aptamer sequence
        :param pairList: list of binding base pairs
        :return:
        """
        # write pair list as fictitious forces to the MMB command file
        printRecord("Folding Sequence. Fold speed={}".format(self.params['fold speed']))        

        mmb = interfaces.mmb(sequence, pairList, self.params, self.i)
        MDAfoldFidelity, MMBfoldFidelity = mmb.run()
        foldFidelity = MMBfoldFidelity

        printRecord('Initial fold fidelity = %.3f' % foldFidelity)
        attempts = 1
        if self.params['fold speed'] != 'quick':  # don't rerun if we're quick folding
            intervalLength = 5  # initial interval length
            while (foldFidelity <= 0.9) and (attempts < 5):  # if it didn't fold properly, try again with a longer annealing time - up to XX times
                self.params['fold speed'] = 'long'
                printRecord("Refolding Sequence")
                mmb = interfaces.mmb(sequence, pairList, self.params, self.i, intervalLength)  # extra intervalLength argument                
                foldFidelity = mmb.run()
                printRecord('Subsequent fold fidelity = %.3f' % foldFidelity)
                attempts += 1
                intervalLength += 5  # on repeated attempts, try harder
                # TODO: what does intervalLength do in nupack???

        self.pdbDict['mmb folded sequence {}'.format(self.i)] = mmb.foldedSequence  # mmb.foldedSequence = foldedSequence_{}.pdb: defined in the mmb.run()
        os.system('mv commands.run* ' + mmb.fileDump)  # mmb.fileDump is a directory for intermediate files during running MMB
        printRecord("Folded Sequence")

    def MDSmoothing(self, structure, relaxationTime=0.01, implicitSolvent=False):  # default relaxation time is 0.01 ns
        """
        Run a short MD sampling to relax the coarse MMB structure
        Relaxation time in nanoseconds, print time in picoseconds
        :param structure:
        :param relaxationTime: does not seem to be used???
        :param implicitSolvent:
        :return:
        """
        structureName = structure.split('.')[0]
        printRecord('\nRunning relaxation (smoothing)')
        if implicitSolvent is False:
            # set up periodic box and condition: pH and ionic strength => protons, ions and their concentrations
            prepPDB(structure, self.params['box offset'], self.params['pH'], self.params['ionic strength'], MMBCORRECTION=True, waterBox=True)

            print('Done preparing files with waterbox. Start openmm.')

        else:  # prepare prmtop and crd file using LEap in ambertools
            printRecord('Implicit solvent: running LEap to generate .prmtop and .crd for folded aptamer...')
            #os.system('pdb4amber {}.pdb > {}_amb_processed.pdb 2> {}_pdb4amber_out.log'.format(structureName, structureName, structureName))
            copyfile(structure, structureName + '_amb_processed.pdb')
            copyfile('./leap_template.in', 'leap.in')
            replaceText('leap.in', 'myDNASEQ', structureName)
            os.system('tleap -f leap.in > leap.out')
            os.system('tail -1 leap.out')  # show last line
            # printRecord(readFinalLines('leap.out', 1))  # show last line. problematic
            structureName += '_amb'  # after pdb4amber and saveAmberParm, the file name became structureName_amb_processed.pdb/top/crd

        processedStructure = structureName + '_processed.pdb'
        processedStructureTrajectory = structureName + '_processed_trajectory.dcd'
        omm = interfaces.omm(structure=processedStructure, params=self.params, simTime=self.params['smoothing time'], implicitSolvent=implicitSolvent)
        self.ns_per_day = omm.doMD()  # run MD in OpenMM framework

        printRecord('Pre-relaxation simulation speed %.1f' % self.ns_per_day + 'ns/day')  # print out sampling speed

        if implicitSolvent is False:
            cleanTrajectory(processedStructure, processedStructureTrajectory)  # remove water and salt from trajectory
        else:  # no water or salt to remove
            copyfile(processedStructure, 'clean_' + processedStructure)  # TODO no cleaning for now
            copyfile(processedStructureTrajectory, 'clean_' + processedStructureTrajectory)  # TODO no cleaning for now

        self.dcdDict['relaxed sequence {}'.format(self.i)] = 'clean_' + processedStructureTrajectory        
        self.pdbDict['relaxed sequence {}'.format(self.i)] = 'relaxedSequence_{}.pdb'.format(self.i)  # specify the file name for final frame, ie, relaxed sequence
        
        # extractFrame(processedStructure, processedStructureTrajectory, -1, self.pdbDict['relaxed sequence {}'.format(self.i)])  # extract final frame        
        # # Current warning from MDA: UserWarning: Unit cell dimensions not found. CRYST1 record set to unitary values. I think: MDA cannot find unit cell dimension info. Can we add info when using implicit solvent?            

        # Another possible way to extract the last frame is to use OpenMM:            
        omm.extractLastFrame(self.pdbDict['relaxed sequence {}'.format(self.i)])
        # # testing
        # extractFrame(processedStructure, processedStructureTrajectory, -1, "MDA_relaxed.pdb")
        # omm.extractLastFrame("OpenMM_relaxed.pdb")

        self.pdbDict['representative aptamer {}'.format(self.i)] = 'relaxedSequence_{}.pdb'.format(self.i)  # in "smooth dock" mode

    def runFreeAptamer(self, aptamer, implicitSolvent=False):
        """
        Run MD sampling for free aptamer and relevant analysis
        :param implicitSolvent:
        :param aptamer:
        :return:
        """
        structureName = aptamer.split('.')[0]  # e.g., aptamer: "relaxedSequence_0.pdb"
        printRecord('\nRunning free aptamer dynamics')
        if implicitSolvent is False:
            # set up periodic box and condition: pH and ionic strength => protons, ions and their concentrations
            prepPDB(aptamer, self.params['box offset'], self.params['pH'], self.params['ionic strength'], MMBCORRECTION=True, waterBox=True)
        else:  # prepare prmtop and crd file using LEap in ambertools
            printRecord('Implicit solvent: running LEap to generate .prmtop and .crd for relaxed aptamer...')
            # os.system('pdb4amber {}.pdb > {}_amb_processed.pdb 2> {}_pdb4amber_out.log'.format(structureName, structureName, structureName))
            copyfile(aptamer, structureName + '_amb_processed.pdb')
            copyfile('./leap_template.in', 'leap.in')
            replaceText('leap.in', 'myDNASEQ', structureName)
            os.system('tleap -f leap.in > leap.out')
            os.system('tail -1 leap.out')  # show last line
            # printRecord(readFinalLines('leap.out', 1))  # show last line. problematic            
            structureName += '_amb'  # after pdb4amber and saveAmberParm, the file name became structureName_amb_processed.pdb/top/crd

        processedAptamer = structureName + '_processed.pdb'
        processedAptamerTrajectory = structureName + '_processed_complete_trajectory.dcd'  # this is output file of autoMD

        self.autoMD(structure=processedAptamer, binding=False, implicitSolvent=implicitSolvent)  # run MD sampling till converged to equilibrium sampling of RC's

        printRecord('Free aptamer simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print out sampling speed
        self.checkRuntime()
        print("\nChecked time. Try cleaning")

        if implicitSolvent is False:
            cleanTrajectory(processedAptamer, processedAptamerTrajectory)  # clean up trajectory for later use. by doing what?
            print("Cleaned traj. Start analyzing.")    
        else:  # no water or salt to remove
            copyfile(processedAptamer, 'clean_' + processedAptamer)  # TODO no cleaning for now
            copyfile(processedAptamerTrajectory, 'clean_' + processedAptamerTrajectory)  # TODO no cleaning for now
            print("No cleaning in implicit solvent, just copied traj. Start analyzing.")
            
        self.dcdDict['sampled aptamer {}'.format(self.i)] = 'clean_' + processedAptamerTrajectory
        self.pdbDict['sampled aptamer {}'.format(self.i)] = 'clean_' + processedAptamer
        aptamerDict = self.analyzeTrajectory(self.pdbDict['sampled aptamer {}'.format(self.i)], self.dcdDict['sampled aptamer {}'.format(self.i)])
        # Within analyzeTrajectory, the last step is also to save an representative frame. We can also replace it using OpenMM??
            # self.omm.extractLastFrame('repStructure_%d' % self.i + '.pdb', representativeIndex) # need more scripting to complete it
        # aptamerDict = {}
        # TODO: analyzeTraj --> getNucDAtraj --> Dihedral: raise ValueError("All AtomGroups must contain 4 atoms")        
        self.pdbDict['representative aptamer {}'.format(self.i)] = 'repStructure_{}.pdb'.format(self.i)

        printRecord('Free aptamer sampling complete')

        return aptamerDict

    def dock(self, aptamer, peptide):
        """
        Use LightDock to run docking and isolate good complex structures
        :param aptamer:
        :param peptide:
        :return:
        """
        printRecord('Docking')
        if bool(self.params['peptide backbone constraint constant']):  # it constant != 0, bool=True
            printRecord('Peptide will be constrained on their dihidral angles')

        buildPeptide(self.peptide, customAngles=bool(self.params['peptide backbone constraint constant']))
        ld = interfaces.ld(aptamer, peptide, self.params, self.i)  # ld is a new class, therefore need to pass in this class's params: self.params
        ld.run()
        topScores = ld.topScores

        for i in range(self.params['N docked structures']):
            copyfile('top_%d' % self.i + '/top_%d' % int(i+1) + '.pdb', 'complex_{}_{}.pdb'.format(self.i, int(i)))
            # TODO: why int
            self.pdbDict['binding complex {} {}'.format(self.i, int(i))] = 'complex_{}_{}.pdb'.format(self.i, int(i))

        printRecord('Docking complete, best scores are {}'.format(topScores))
        return topScores

    def bindingDynamics(self, complex, implicitSolvent=False):
        """
        Run MD sampling for a complex structure and relevant analysis
        :param implicitSolvent:
        :param complex:
        :return:
        """
        # structureName = complex.split('.')[0]
        printRecord('Running Binding Simulation')
        # set up periodic box and condition: pH and ionic strength => protons, ions and their concentrations
        prepPDB(complex, self.params['box offset'], self.params['pH'], self.params['ionic strength'], MMBCORRECTION=True, waterBox=True)
        processedComplex = complex.split('.')[0] + '_processed.pdb'
        processedComplexTrajectory = processedComplex.split('.')[0] + '_complete_trajectory.dcd'  # this is output file of autoMD

        self.autoMD(processedComplex, binding=True)

        printRecord('Complex simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print out sampling speed
        self.checkRuntime()

        cleanTrajectory(processedComplex, processedComplexTrajectory)

        self.dcdDict['sampled complex {} {}'.format(self.i, self.j)] = 'clean_' + processedComplexTrajectory
        self.pdbDict['sampled complex {} {}'.format(self.i, self.j)] = 'clean_' + processedComplex

        bindingDict = self.analyzeBinding(self.pdbDict['sampled complex {} {}'.format(self.i, self.j)],
                                          self.dcdDict['sampled complex {} {}'.format(self.i, self.j)],
                                          self.pdbDict['sampled aptamer {}'.format(self.i)],
                                          self.dcdDict['sampled aptamer {}'.format(self.i)])
        # print findings
        printRecord('Binding Results: Contact Persistence = {:.2f}, Contact Score = {:.2f}, Conformation Change = {:.2f}'.format(bindingDict['close contact ratio'], bindingDict['contact score'], bindingDict['conformation change']))
        # TODO: why {:.2f} instead of {.2f}?
        printRecord('Complex sampling complete')

        return bindingDict

    def autoMD(self, structure, binding=False, implicitSolvent=False):
        """
        Run MD till either for a set amount of time or till the dynamics 'converge'
        Optionally, sample after we reach convergence ("equilibrium") -- not implemented yet (?)
        :param structure:
        :param binding:
        :return:
        """
        if binding:  # whether it's complex or free aptamer
            maxIter = self.params['max complex sampling iterations']
        else:
            maxIter = self.params['max aptamer sampling iterations']
        cutoff = self.params['autoMD convergence cutoff']

        structureName = structure.split('.')[0]  # e.g., structure: relaxedSequence_0_amb_processed.pdb
        if self.params['auto sampling'] is False:  # just run MD for the given sampling time
            self.analyteUnbound = False
            omm = interfaces.omm(structure=structure, params=self.params, implicitSolvent=implicitSolvent)
            self.ns_per_day = omm.doMD()  # run MD in OpenMM framework
            print('Generated:', structureName + '_trajectory.dcd')
            os.replace(structureName + '_trajectory.dcd', structureName + "_complete_trajectory.dcd")
            print('Replaced ^ with:', structureName + '_complete_trajectory.dcd')

        elif self.params['auto sampling'] is True:  # run MD till convergence (equilibrium)
            converged = False
            iter = 0
            self.analyteUnbound = False

            while (converged is False) and (iter < maxIter):
                iter += 1
                omm = interfaces.omm(structure=structure, params=self.params, implicitSolvent=implicitSolvent)
                self.ns_per_day = omm.doMD()

                if iter > 1:  # if we have multiple trajectory segments, combine them
                    appendTrajectory(structure, structureName + '_trajectory-1.dcd', structureName + '_trajectory.dcd')
                    os.replace('combinedTraj.dcd', structureName + '_trajectory-1.dcd')
                else:
                    os.replace(structureName + '_trajectory.dcd', structureName + '_trajectory-1.dcd')  # in case we need to combine two trajectories

                combinedSlope = checkTrajPCASlope(structure, structureName + '_trajectory-1.dcd', self.params['print step'])
                # TODO what is the slope and what it for?
                if binding:
                    self.analyteUnbound = checkMidTrajectoryBinding(structure, structureName + '_trajectory-1.dcd', self.peptide, self.sequence, self.params, cutoffTime=1)
                    if self.analyteUnbound:
                        printRecord('Analyte came unbound!')

                if (combinedSlope < cutoff) or (self.analyteUnbound is True):
                    converged = True

            os.replace(structureName + '_trajectory-1.dcd', structureName + '_complete_trajectory.dcd')

    def analyzeTrajectory(self, structure, trajectory):
        """
        Analyze trajectory for aptamer fold
        :param trajectory:
        :return:
        """
        u = mda.Universe(structure, trajectory)

        # extract distance info through the trajectory
        wcTraj = getWCDistTraj(u)  # watson-crick base pairing distances (H-bonding)
        baseDistTraj = getBaseBaseDistTraj(u)  # FAST, base-base center-of-geometry distances
        nucleicAnglesTraj = getNucDATraj(u)  # FAST, new, omits 'chi' angle between ribose and base

        # 2D structure analysis
        pairTraj = getPairTraj(wcTraj)
        secondaryStructure = analyzeSecondaryStructure(pairTraj)  # find equilibrium secondary structure
        if self.params['skip MMB'] is False:
            printRecord('Predicted 2D structure :' + self.ssAnalysis['2d string'][self.i])  # if we did not fold the structure, we probably do not know the secondary structure.
        printRecord('Actual 2D structure    :' + configToString(secondaryStructure))

        # 3D structure analysis
        mixedTrajectory = np.concatenate((baseDistTraj.reshape(len(baseDistTraj), int(baseDistTraj.shape[-2] * baseDistTraj.shape[-1])), nucleicAnglesTraj.reshape(len(nucleicAnglesTraj), int(nucleicAnglesTraj.shape[-2] * nucleicAnglesTraj.shape[-1]))), axis=1)  # mix up all our info
        representativeIndex, pcTrajectory, eigenvalues = isolateRepresentativeStructure(mixedTrajectory)

        # save this structure a separate file
        extractFrame(structure, trajectory, representativeIndex, 'repStructure_%d' % self.i + '.pdb')

        # we also need  function which processes the energies spat out by the trajectory logs
        analysisDict = {}  # compile our results
        analysisDict['2D structure'] = secondaryStructure
        analysisDict['RC trajectories'] = pcTrajectory
        analysisDict['representative structure index'] = representativeIndex

        return analysisDict

    def analyzeBinding(self, bindStructure, bindTrajectory, freeStrcuture, freeTrajectory):
        """
        Analyze trajectory for aptamer fold + analyte binding information
        :param bindStructure:
        :param bindTrajectory:
        :param freeStrcuture:
        :param freeTrajectory:
        :return:
        """
        bindu = mda.Universe(bindStructure, bindTrajectory)
        freeu = mda.Universe(freeStrcuture, freeTrajectory)
        bindingDict = bindingAnalysis(bindu, freeu, self.peptide, self.sequence)  # look for contacts between analyte and aptamer
        if self.analyteUnbound:
            bindingDict['analyte came unbound'] = True
        else:
            bindingDict['analyte came Unbound'] = False

        return bindingDict

    # =======================================================
    # ====== supporting functions are followed (part2) ======
    # =======================================================

    def checkRuntime(self):
        """
        Estimate the maximum runtime for this code, given a ns/day calculation from existing work done
        :return:
        """
        projectedMaxRuntime = self.computeRuntime()
        buffer = 0.8  # ideally would like a buffer of 20% of the walltime
        doReduction = 0
        # if the code is going to take too long, shorten the sampling time until it fits
        while projectedMaxRuntime > (buffer * self.params['max walltime']):  # max walltime is in hours
            doReduction = 1
            self.params['sampling time'] *= 0.95
            projectedMaxRuntime = self.computeRuntime()

        if doReduction == 1:
            if self.params['sampling time'] < 0.1:
                printRecord('Sampling time reduced to 100 ps - this run is too expensive, and will now terminate')
                self.terminateRun()
            else:
                printRecord('Sampling time chunk reduced to {} ns'.format(self.params['sampling time']))

    def computeRuntime(self):
        dt = self.params['sampling time']  # ns in MD sampling after getting (smoothed) 3D structure
        ssStructures = self.params['N 2D structures']
        dockingStructures = self.params['N docked structures']
        aptamerSteps = self.params['max aptamer sampling iterations']
        complexSteps = self.params['max complex sampling iterations']
        maxNS = dt * (ssStructures * aptamerSteps + ssStructures * dockingStructures * complexSteps)  # maximum possible sampling time in nanoseconds
        projectedMaxRuntime = maxNS / self.ns_per_day * 24  # how many hours needed to complete our maximum sampling time

        return projectedMaxRuntime

    def terminateRun(self):
        """ for some reason, the run needs to end """
        sys.exit()
