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

import sys
import glob
import shutil
from shutil import copyfile, copytree

import interfaces
from utils import *
from analysisTools import *

# noinspection PyPep8Naming  # ignore inspetion based on python PEP8 Naming convention


class opendna:
    def __init__(self, params):
        self.workDir = ""  # rather use an empty string, not an empty list
        self.params = params
        self.aptamerSeq = self.params['aptamer_seq']        
        self.targetPDB = self.params['ligand']        # either pdb filename (eg, "/path/peptide.pdb") or None
        self.targetType = self.params['ligand_type']  # 'peptide' or 'DNA' or 'RNA' or 'other' or None. Assume now the target is a peptide
        self.targetSeq = self.params['ligand_seq']    # if target is not a peptide/DNA/RNA, targetSeq is None
        #analyzeBinding method is commented out: no analysis in this simulator.
        
        self.pdbDict = {}
        self.dcdDict = {}
        self.actionDict = {}
        self.getActionDict()  # specify the actions based on the selected mode
        
        # modify actions if necessary
        if self.params['skip_MMB'] is True:  
            self.actionDict['do MMB'] = False
            if self.params['mode'] != '2d structure':
                self.actionDict['do 2d analysis'] = False  
        
        if self.params['pickup'] is True:
            printRecord('Copying the checkpoint file: {} to run folder: {}.'.format(self.params['chk_file'],self.workDir), self.workDir+'/')
            printRecord('Copying the topology file: {} to run folder: {}.'.format(self.params['pickup_pdb'],self.workDir), self.workDir+'/')

            # These conditions actually have been checked in main.py, here is for redundancy
            if os.path.exists(self.params['chk_file']) and os.path.exists(self.params['pickup_pdb']):
                # user provides a .chk file stored at params['chk_file']
                copyfile(self.params['chk_file'], os.path.join(self.workDir, os.path.basename(self.params['chk_file'])))
                self.params['chk_file'] = os.path.basename(self.params['chk_file']) # get its filename
                # user provides a .pdb file stored at params['pickup_pdb']
                copyfile(self.params['pickup_pdb'], os.path.join(self.workDir, os.path.basename(self.params['pickup_pdb'])))
                self.params['pickup_pdb'] = os.path.basename(self.params['pickup_pdb']) # get its filename
            else:
                printRecord('Designated checkpoint file does not exist! Terminating the pipeline.', self.workDir+'/')
                self.terminateRun()

            if self.params['pickup_from_freeAptamerChk'] is True:
                # everything before "freeAptamerDynamics" is skipped
                self.actionDict['do 2d analysis'] = False
                self.actionDict['do MMB'] = False            
                self.actionDict['do smoothing'] = False
                self.actionDict['get equil repStructure'] = True  # making sure to run "freeAptamerDynamics"
            elif self.params['pickup_from_complexChk'] is True:
                # everything before "complexDynamics" is skipped
                self.actionDict['do 2d analysis'] = False
                self.actionDict['do MMB'] = False            
                self.actionDict['do smoothing'] = False            
                self.actionDict['get equil repStructure'] = False
                self.actionDict['do docking'] = False
                self.actionDict['do binding'] = True  # making sure to run "complexDynamics"

        if self.actionDict['make workdir']:
            self.setup()  # if we don't need a workdir & MMB files (eg, give a 3D structure), don't make one.

        self.i = int(-1)
        self.j = int(-1)

    def getActionDict(self):
        """
        generate a binary sequence of 'to-do's', given user's input into self.actionDict
        :return: none
        """
        if self.params['mode'] == '2d structure':
            self.actionDict['make workdir'] = True
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
            self.params['skip_smoothing'] = False

        elif self.params['mode'] == 'coarse dock':  # directly dock the structure generated in MMB. Is it useful by any chance?
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = False
            self.actionDict['get equil repStructure'] = False
            self.actionDict['do docking'] = True
            self.actionDict['do binding'] = False

        elif self.params['mode'] == 'smooth dock':  # dock with no MD sampling. Useful?
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = True
            self.actionDict['get equil repStructure'] = False
            self.actionDict['do docking'] = True
            self.actionDict['do binding'] = False
            self.params['skip_smoothing'] = False

        elif self.params['mode'] == 'free aptamer':
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = True
            self.actionDict['get equil repStructure'] = True
            self.actionDict['do docking'] = False
            self.actionDict['do binding'] = False

        elif self.params['mode'] == 'full dock':
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

    def autoMakeNewWorkDir(self):
        """
        make a new working directory: increment --run_num by 1 => not overlapping previously existing directories
        :return:
        """
        workdirs = glob.glob(os.path.join(self.params['workdir'], 'run*'))  # check for prior working directories
        if len(workdirs) > 0: # there are exisiting run{run_num}}]
            prev_runs = []
            for i in range(len(workdirs)):
                prev_runs.append(int(workdirs[i].split('run')[-1]))

            prev_max = max(prev_runs)
            self.workDir = os.path.join(self.params['workdir'], 'run' + str(prev_max + 1)) # increment by 1 => {workdir}/run{run_num+1}
            os.mkdir(self.workDir)
            printRecord('Starting Fresh Run %d' % (prev_max + 1), self.workDir+'/')
        else:
            self.workDir = os.path.join(self.params['workdir'], 'run1')
            os.mkdir(self.workDir)
            printRecord('Starting Fresh Run 1', self.workDir+'/')

    def setup(self):
        """
        set up a working directory, copy in relevant xyz and keyfiles, move to the working directory
        :return:
        """
        # Make a new directory
        if self.params['run_num'] == 0:  # auto-increment the max run_num folder
            self.autoMakeNewWorkDir()
        else:  # create a working directory defined by the user
            self.workDir = os.path.join(self.params['workdir'], 'run' + str(self.params['run_num'])) #{workdir}/run{run_num}
            os.mkdir(self.workDir)
            printRecord('Starting Fresh Run {}'.format(self.params['run_num']), self.workDir+'/')

        # copy relevant files to the workDir
        os.mkdir(os.path.join(self.workDir, 'outfiles_of_pipeline_modules'))
                
        # If we skip_MMB because we have a folded structure to start with
        if self.params['skip_MMB'] is True:            
            printRecord('Skipping 2D analysis and MMB folding. Starting with a folded structure.', self.workDir+'/')            
            if not os.path.isfile(self.params['init_structure']):
            # if not os.path.exists(self.params['init_structure']):
                printRecord('Designated folded structure does not exist! Terminating the pipeline.', self.workDir+'/')
                self.terminateRun()  # maybe: do it via "raise an error"
            else:
                printRecord('Copying the folded structure: {} to run folder: {}.'.format(self.params['init_structure'],self.workDir), self.workDir+'/')
                init_structure_in_workDir = os.path.basename(self.params['init_structure']) # remove the direcotry info
                copyfile(self.params['init_structure'], os.path.join(self.workDir, init_structure_in_workDir)) # copy to workDir
                self.params['init_structure'] = init_structure_in_workDir
        elif (self.params['mode'] != '2d structure') and (self.params['pickup'] is False):
            # copy MMB params and command script
            copyfile(self.params['mmb_params'], os.path.join(self.workDir, 'parameters.csv'))
            copyfile(self.params['mmb_normal_template'], os.path.join(self.workDir, 'commands.template.dat'))
            copyfile(self.params['mmb_quick_template'], os.path.join(self.workDir, 'commands.template_quick.dat'))
            copyfile(self.params['mmb_slow_template'], os.path.join(self.workDir, 'commands.template_long.dat'))
        
        # if do docking, need targetPDB
        if (self.actionDict['do docking'] is True):  
            if self.targetPDB is None:  # if user does not provide a target ligand pdb, use our provided example target ligand: a peptide.
                printRecord('Use the peptide ligand in examples/...', self.workDir+'/')
                self.targetPDB = 'example_peptide_ligand.pdb' # give it a name
                copyfile(self.params['example_target_pdb'], os.path.join(self.workDir, self.targetPDB))  # ie, target of the aptamer
                self.targetSeq = self.params['example_peptide_seq']  # self.targetSeq is not empty because we are dealing with a peptide.
            else:  # user provides a targetPDB stored at self.targetPDB (ie, params['ligand'])
                copyfile(self.targetPDB, os.path.join(self.workDir, os.path.basename(self.targetPDB)))
                # get its filename:
                self.targetPDB = os.path.basename(self.targetPDB)
                self.params['ligand'] = self.targetPDB
            # From now on, "self.targetPDB is None" is equivalent to "'do docking' is False".

        if self.params['implicit_solvent'] is True:
            copyfile(self.params['leap_template'], self.workDir + '/leap_template.in')

        # move to working dir
        os.chdir(self.workDir)
        
        # specify the DNA_force_field for free aptamer MD sampling if implicit solvent. Assume only dealing with a DNA aptamer.
        if self.params['implicit_solvent'] is True:
            replaceText('leap_template.in', 'DNA_FF', self.params['DNA_force_field'])
        
        # Print prompt: free aptamer or aptamer + target ligand
        printRecord('Simulation mode: {}'.format(self.params['mode']))
        if (self.actionDict['do docking'] is False):  # no given target ligand, nor do docking
            printRecord('Simulating free aptamer: {}'.format(self.aptamerSeq))
        else:
            printRecord('Simulating {} with {}'.format(self.aptamerSeq, self.targetPDB))

        if (self.params['skip_MMB'] is False) and (self.params['device'] == 'local'):
            if (self.params['operating_system'] == 'linux') or (self.params['operating_system'] == 'WSL'):
                os.environ["LD_LIBRARY_PATH"] = self.params['mmb_dir']  # Add MMB library's path to environment variable
                # When running MMB on 'local':
                    # Windows: no action needed. But this pipeline does not support Windows OS, due to NUPACK.
                    # WSL or linux: explicitly add MMB library path.
                    # macos: no action needed here.                                
                # warning: do NOT use os.environ["DYLD_LIBRARY_PATH"] for macos! It does not help MMB locate the library AND it confuses OpenMM.
                # tip: refer to 'mmb' class in interfaces.py on how to run MMB on macos machine.


    def run(self):
        """
        run the end-to-end simulation pipeline for the chosen mode
        consult checkpoints to not repeat prior steps
        :return:
        """
        # outputDict = {'params': self.params}
        outputDict = {}
        outputDict['params'] = self.params
        np.save('final_output_dict', outputDict)  # Save an array to a binary file in NumPy ``.npy`` format.
        
        if self.actionDict['do 2d analysis']:  # get secondary structure
            self.pairLists = self.getSecondaryStructure(self.aptamerSeq)
            outputDict['2d analysis'] = self.ssAnalysis
            np.save('final_output_dict', outputDict)  # save 2d structure results
            printRecord('Running over %d' % len(self.pairLists) + ' possible 2D structures.')
            num_2dSS = len(self.pairLists)
        
        else:  # there are three cases where 2d analysis is not carried out.
            if self.params['pickup_from_complexChk'] is True:
                printRecord('Resume a sampling of aptamer-ligand complex from .chk_file, therefore skip all the steps before it.')
            elif self.params['pickup_from_freeAptamerChk'] is True:  
                printRecord('Resume a sampling of free aptamer from .chk_file, therefore skip all the steps before it.')                
            else:
                printRecord('Starting with an existing folded strcuture, therefore skip 2d analysis and MMB folding.')
            num_2dSS = 1  # quick and dirty

        for self.i in range(num_2dSS):  # loop over all possible secondary structures
            if self.actionDict['do 2d analysis'] is True:  # self.ssAnalysis only exists if we "do 2d analysis"
                printRecord('>>> Predicted 2D structure #{}                   : {}'.format(self.i, self.ssAnalysis['displayed 2d string'][self.i]))
                self.pairList = np.asarray(self.pairLists[self.i])  # be careful!!!: .pairList vs. .pairLists                                                

            if self.actionDict['do MMB']:  # fold 2D into 3D
                self.foldSequence(self.aptamerSeq, self.pairList)
            elif self.params['skip_MMB'] is True:   
                # situation here: no MMB + not to resume a MD of either free aptamer or aptamer-ligand complex: ie, (self.params['pickup_from_freeAptamerChk'] is False) and (self.params['pickup_from_complexChk'] is False)
                # meaning of the condition: just skip_MMB and preceed to next step: MD smoothing or MD sampling or dock (if "coarse dock")
                self.pdbDict['folded aptamer {}'.format(self.i)] = self.params['init_structure'] # filen has been copied to workDir and only store the file's basename
                self.pdbDict['representative aptamer {}'.format(self.i)] = self.params['init_structure']  # for "coarse dock" mode.            

            if self.params['skip_smoothing'] is True:
                self.actionDict['do smoothing'] = False
                if self.params['mode'] != '2d structure':
                    printRecord('\nNo relaxation (smoothing) of the folded aptamer.')
                self.pdbDict['relaxed aptamer {}'.format(self.i)] = 'foldedAptamer_0.pdb'
            elif self.actionDict['do smoothing']:  # just to double check
                aptamer_smooth_dir = 'md_aptamer_relaxation_runfiles_%d' % self.i
                os.mkdir(aptamer_smooth_dir)

                if self.params['skip_MMB'] is False:
                    self.MDSmoothing(runfilesDir=aptamer_smooth_dir, structure=self.pdbDict['mmb folded aptamer {}'.format(self.i)], relaxationTime=self.params['smoothing_time'], implicitSolvent=self.params['implicit_solvent'])  # relax for xx nanoseconds
                else:  # if given a folded structure, skip_MMB and directly use the provided structure
                    self.MDSmoothing(runfilesDir=aptamer_smooth_dir, structure=self.pdbDict['folded aptamer {}'.format(self.i)], relaxationTime=self.params['smoothing_time'], implicitSolvent=self.params['implicit_solvent'])  # relax for xx nanoseconds

            if self.actionDict['get equil repStructure']:
                aptamer_sampling_dir = 'md_aptamer_sampling_runfiles_%d' % self.i
                os.mkdir(aptamer_sampling_dir)
                if self.params['pickup_from_freeAptamerChk'] is False:
                    outputDict['free aptamer results {}'.format(self.i)] = self.freeAptamerDynamics(runfilesDir=aptamer_sampling_dir, aptamerPDB=self.pdbDict['relaxed aptamer {}'.format(self.i)], simTime=self.params['aptamer_sampling_time'], implicitSolvent=self.params['implicit_solvent'])
                else:  # user must make sure they have provided a .chk_file
                    if self.params['implicit_solvent'] is True:
                        raise ValueError("Resuming sampling is currently supported for explicit solvent run only, because implicit solvent run is fast.")
                    else:
                        outputDict['free aptamer results {}'.format(self.i)] = self.freeAptamerDynamics(runfilesDir=aptamer_sampling_dir, aptamerPDB=self.params['pickup_pdb'], simTime=self.params['aptamer_sampling_time'], implicitSolvent=self.params['implicit_solvent'])
                np.save('final_output_dict', outputDict)  # save outputs

            if (self.actionDict['do docking'] is True) and (self.targetPDB is not None):  # find docking configuration for the complexed structure
                # coarse dock: no smoothing;
                # smooth dock: smooth + dock;
                # full dock: smooth + dock + equil structure;
                # full binding: smooth + dock + equil structure + sampling dynamics;                
                outputDict['dock scores {}'.format(self.i)] = self.dock(self.pdbDict['representative aptamer {}'.format(self.i)], self.targetPDB)
                # pdbDict['representative aptamer {}' is defined at MMB folding, MD smoothing or freeAptamerDynamics                
                np.save('final_output_dict', outputDict)  # save outputs
                # N_docked_structures are specified by user: how many docked structures do we want to investigate
                # num_docked_structure = self.params['N_docked_structures']
                num_docked_structure = len(self.topDockingScores)
                
            if self.params['pickup_from_complexChk'] is True:  # it means we did not dock but we have a complex structure.
                num_docked_structure = 1
                printRecord('Resuming the sampling of 1 docked aptamer-ligand structure')

            if self.actionDict['do binding'] is True:  # run MD on the aptamer-ligand structure
                printRecord('\nSampling %d' % num_docked_structure + ' docked structures...')                
                for self.j in range(num_docked_structure):  # loop over docking configurations for a given secondary structure
                    printRecord('\nDocked structure #{}:'.format(self.j))

                    complex_sampling_dir = 'md_complex_sampling_runfiles_{}_{}'.format(self.i, self.j)
                    os.mkdir(complex_sampling_dir)
                    if self.params['pickup_from_complexChk'] is False:
                        outputDict['binding results {} {}'.format(self.i, self.j)] = self.complexDynamics(runfilesDir=complex_sampling_dir, complexPDB=self.pdbDict['binding complex {} {}'.format(self.i, int(self.j))], simTime=self.params['complex_sampling_time'], implicitSolvent=self.params['implicit_solvent'])
                    else:  # user must make sure they have provided a .chk_file
                        if self.params['implicit_solvent'] is True:
                            raise ValueError("Resuming sampling is currently supported for explicit solvent run only, because implicit solvent run is fast.")
                        else:
                            outputDict['binding results {} {}'.format(self.i, self.j)] = self.complexDynamics(runfilesDir=complex_sampling_dir, complexPDB=self.params['pickup_pdb'], simTime=self.params['complex_sampling_time'], implicitSolvent=self.params['implicit_solvent'])
                    np.save('final_output_dict', outputDict)
            
        printRecord('\nAll results are saved to output folder: ' + self.workDir + ', including:')
        printRecord('       a log file: run_output_log.txt')
        printRecord('       a NumPy file of main outputs: final_output_dict.npy')
        printRecord('\nExiting pipeline...')

        return outputDict

    # ======================================================================================
    # ======================================================================================
    # ====================== supporting functions are followed (part1) =====================
    # ======================================================================================
    # ======================================================================================
    def getSecondaryStructure(self, aptamerSeq):
        """
        get the secondary structure(s) for a given aptamer's sequence
        using seqfold here - identical features are available in nupack, though results are sometimes different
        :param aptamerSeq: aptamer's sequence
        :return:
        """
        printRecord("\nPredicting Secondary Structure(s)")
        if self.params['secondary_structure_engine'] == 'seqfold':
            ssString, pairList = getSeqfoldStructure(aptamerSeq, self.params['temperature'])  # seqfold guess
            self.ssAnalysis = {'2d string': ssString, 'pair list': pairList, 'displayed 2d string': ssString}
            # print('seqfold ssAnalysis: ',self.ssAnalysis)
            return pairList
        elif self.params['secondary_structure_engine'] == 'NUPACK':
            nup = interfaces.nupack(aptamerSeq, self.params['temperature'], self.params['ionicStrength'], self.params['Mg_conc'])  # initialize nupack
            self.ssAnalysis = nup.run()  # run nupack analysis of possible 2D structures.

            distances = getSecondaryStructureDistance(self.ssAnalysis['config'])            
            if len(distances) > 1:
                topReps, topProbs, topDists = do2DAgglomerativeClustering(self.ssAnalysis['config'], self.ssAnalysis['state prob'], distances)
                topPairLists = []                
                nLists = min(len(topReps), self.params['N_2D_structures'])
                self.ssAnalysis['displayed 2d string'] = []  # to output the dot-parenthesis notation of selected 2d structure
                for i in range(nLists):
                    if (i == 0) or (np.amin(topDists[i, :i]) > 0.05): # only add a config to the list if it's at least 5% different from a previously accepted config
                        self.ssAnalysis['displayed 2d string'].append(configToString(topReps[i].astype(int)))
                        topPairLists.append(ssToList(configToString(topReps[i].astype(int))))
                return topPairLists

            else:  # if, for some reason, there's only a single plausible structure (unlikely)
                self.ssAnalysis['displayed 2d string'] = self.ssAnalysis['2d string']
                return self.ssAnalysis['pair list']

    def foldSequence(self, aptamerSeq, pairList):
        """
        generate a coarsely folded 3D structure for the aptamer
        :param aptamerSeq: the ssDNA aptamer's sequence
        :param pairList: list of binding base pairs
        :return:
        """
        # write pair list as fictitious forces to the MMB command file
        printRecord("\nFolding Aptamer from Sequence. Fold speed = {}.".format(self.params['fold_speed']))        

        mmb = interfaces.mmb(aptamerSeq, pairList, self.params, self.i)
        MDAfoldFidelity, MMBfoldFidelity = mmb.run()
        foldFidelity = MMBfoldFidelity
        printRecord('Initial fold fidelity = %.3f' % foldFidelity)
        printRecord('Initial fold fidelity = %.3f (from MDAnalysis)' % MDAfoldFidelity)
        attempts = 1
        # refold if folding fidelity < threshold
        if self.params['fold_speed'] != 'quick':  # do not rerun if we're quick folding
            intervalLength = 5  # initial interval length
            while (foldFidelity <= self.params['fold_fidelity']) and (attempts < 5):  # if it didn't fold properly, try again with a longer annealing time - up to XX times
                self.params['fold_speed'] = 'slow'
                
                printRecord("Refolding Aptamer from Sequence")
                mmb = interfaces.mmb(aptamerSeq, pairList, self.params, self.i, intervalLength)  # extra intervalLength argument                
                MDAfoldFidelity, MMBfoldFidelity = mmb.run()
                foldFidelity = MMBfoldFidelity
                printRecord('Subsequent fold fidelity = %.3f' % foldFidelity)
                printRecord('Subsequent MDA fold fidelity = %.3f' % MDAfoldFidelity)
                attempts += 1
                intervalLength += 5  # on repeated attempts, try harder

        self.pdbDict['mmb folded aptamer {}'.format(self.i)] = mmb.foldedAptamerSeq  # mmb.foldedAptamerSeq = foldedAptamer_{}.pdb: defined in the mmb.run()
        
        # mmb.fileDump is a directory for intermediate files during running MMB
        os.system('mv commands.run* ' + mmb.fileDump)
        os.system('mv parameters.csv ' + mmb.fileDump)
        os.system('mv commands.template.dat ' + mmb.fileDump)
        os.system('mv commands.template_quick.dat ' + mmb.fileDump)
        os.system('mv commands.template_long.dat ' + mmb.fileDump)
    
        printRecord(">>> MMB folded the aptamer and generated folded structure: {}".format(self.pdbDict['mmb folded aptamer {}'.format(self.i)]))
        
        # For "coarse dock" mode:
        self.pdbDict['representative aptamer {}'.format(self.i)] = self.pdbDict['mmb folded aptamer {}'.format(self.i)]  # ie, 'relaxedAptamer_{}.pdb'.format(self.i)

    def MDSmoothing(self, runfilesDir, structure, relaxationTime=0.01, implicitSolvent=False):  # default relaxation time is 0.01 ns
        """
        Run a short MD sampling to relax the coarse MMB structure
        Relaxation time in nanoseconds, print time in picoseconds
        :param runfilesDir: a directory to store runfiles generated during MD, eg, "md_aptamer_relaxation_runfiles_0"
        :param structure:
        :param relaxationTime: params['smoothing_time']
        :param implicitSolvent:
        :return:
        """
        structureName = structure.split('.')[0] # structureName should have no directory info and just the file basename
        printRecord('\nRunning Short Relaxation (Smoothing) of Free Aptamer')
        if implicitSolvent is False:
            # set up periodic box and condition: pH and ionic strength => protons, ions and their concentrations
            prepPDB(structure, self.params['box_offset'], self.params['pH'], self.params['ionicStrength'], MMBCORRECTION=True, waterBox=True)
            printRecord('Done preparing files with waterbox. Start OpenMM.')
        else:  # prepare prmtop and crd file using LEap in ambertools
            printRecord('Implicit solvent: running LEap to generate .prmtop and .crd for folded aptamer...')
            #os.system('pdb4amber {}.pdb > {}_amb_processed.pdb 2> {}_pdb4amber_out.log'.format(structureName, structureName, structureName))
            copyfile(structure, structureName + '_amb_processed.pdb')
            copyfile('./leap_template.in', 'leap_MDSmoothing.in')
            replaceText('leap_MDSmoothing.in', 'mySTRUCTURE', structureName)
            os.system('tleap -f leap_MDSmoothing.in > leap_MDSmoothing.out')
            os.rename('leap.log', 'leap_MDSmoothing.log')
            os.system('mv leap_MDSmoothing.in outfiles_of_pipeline_modules/')
            os.system('mv leap_MDSmoothing.out outfiles_of_pipeline_modules/')
            os.system('mv leap_MDSmoothing.log outfiles_of_pipeline_modules/')
            # os.system('tail -1 leap_MDSmoothing.out')  # show last line
            # printRecord(readFinalLines('leap_MDSmoothing.out', 1))  # show last line. problematic
            structureName += '_amb'  # after pdb4amber and saveAmberParm, will generate structureName_amb_processed.pdb/top/crd

        processedStructure = structureName + '_processed.pdb'
        processedStructureTrajectory = structureName + '_processed_trajectory.dcd'
        omm = interfaces.omm(runfilesDir=runfilesDir, structurePDB=processedStructure, params=self.params, simTime=relaxationTime, binding=False, implicitSolvent=implicitSolvent)
        self.ns_per_day = omm.doMD()  # run MD in OpenMM framework

        printRecord('Short relaxation (MD smoothing): simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print out sampling speed
        printRecord("Try cleaning...")

        stru_in_runfiles = os.path.join(runfilesDir, processedStructure)
        traj_in_runfiles = os.path.join(runfilesDir, processedStructureTrajectory)
        # move relaxed pdb used by OpenMM to runfilesDir (original trajectory was already written to there)        
        os.replace(processedStructure, stru_in_runfiles)

        if implicitSolvent is False:
            cleanTrajectory(stru_in_runfiles, traj_in_runfiles)  # remove water and salt from trajectory. Prepend 'clean_' to files' basenames
            printRecord("Cleaned short MD relaxation trajectory of free aptamer.") # cleaned files are in the same directory as the uncleaned ones.
            printRecord('>>> Generated: clean_' + processedStructureTrajectory + ' in the folder ./' + runfilesDir)
            self.dcdDict['relaxed aptamer {}'.format(self.i)] = 'clean_' + processedStructureTrajectory
        else:  # no water or salt to remove
            # copyfile(processedStructure, 'clean_' + processedStructure)
            # copyfile(processedStructureTrajectory, 'clean_' + processedStructureTrajectory)
            printRecord("No need for cleaning in implicit solvent.")
            printRecord('>>> Generated short MD trajectory of free aptamer: ' + processedStructureTrajectory + ' in the folder ./' + runfilesDir)
            self.dcdDict['relaxed aptamer {}'.format(self.i)] = processedStructureTrajectory

        # extractFrame(processedStructure, processedStructureTrajectory, -1, self.pdbDict['relaxed aptamer {}'.format(self.i)])  # extract final frame and omit solvent and salts
        # # Current warning from MDA: UserWarning: Unit cell dimensions not found. CRYST1 record set to unitary values. MDA cannot find unit cell dimension info. Can we add info when using implicit solvent?
        # Another possible way to extract the last frame is to use OpenMM:        
        printRecord("Extracting the last frame of the relaxation trajectory as the relaxed structure...")
        relaxedAptamerPDB = 'relaxedAptamer_{}.pdb'.format(self.i)  # specify the file name for the extracted final frame
        omm.extractLastFrame(lastFrameFileName=relaxedAptamerPDB)
        cleanPDB(relaxedAptamerPDB)  # remove ion and generate the following clean_relaxedAptamer.pdb
        clean_relaxedAptamerPDB = 'clean_' + relaxedAptamerPDB
        os.system('mv ' + clean_relaxedAptamerPDB + ' ' + relaxedAptamerPDB)  # overwrite the relaxedAptamer file with the one without ions and solvent
        # printRecord('>>> OpenMM: saved the last frame into: ' + relaxedAptamerPDB)
        printRecord('>>> MD smoothing generated relaxed structure: ' + relaxedAptamerPDB)
        self.pdbDict['relaxed aptamer {}'.format(self.i)] = relaxedAptamerPDB
        
        # in "smooth dock" mode:
        self.pdbDict['representative aptamer {}'.format(self.i)] = self.pdbDict['relaxed aptamer {}'.format(self.i)]  # ie, 'relaxedAptamer_{}.pdb'.format(self.i)

    def freeAptamerDynamics(self, runfilesDir, aptamerPDB, simTime, implicitSolvent=False):
        """
        Run MD sampling for free aptamer and relevant analysis
        :param runfilesDir: a directory to store runfiles generated during MD, eg, "md_aptamer_sampling_runfiles_0"
        :param implicitSolvent:
        :param aptamerPDB: PDB file name of the aptamer to be simulated. It contains no directory info
        :return:
        """
        structureName = aptamerPDB.split('.')[0]  # aptamerPDB: eg "foldedAptamer_0.pdb" or "foldedAptamer_0_processed.pdb" (if resuming sampling)

        if self.params['pickup_from_freeAptamerChk'] is True:
            printRecord('\nResuming a previous free aptamer dynamics')
            self.autoMD(runfilesDir=runfilesDir, structurePDB=aptamerPDB, simTime=simTime, binding=False, implicitSolvent=implicitSolvent)  # In this case, relaxedAptamer_0_processed.pdb (due to resuming)
            processedAptamer = aptamerPDB
            processedAptamerTrajectory = structureName + '_complete_trajectory.dcd'  # this is output file of autoMD
        else:
            printRecord('\nRunning Fresh Free Aptamer Dynamics')
            if implicitSolvent is False:
                # set up periodic box and condition: pH and ionic strength => protons, ions and their concentrations
                prepPDB(aptamerPDB, self.params['box_offset'], self.params['pH'], self.params['ionicStrength'], MMBCORRECTION=True, waterBox=True)
                printRecord('Done preparing files with waterbox. Start OpenMM.')
            else:  # prepare prmtop and crd file using LEap in ambertools
                printRecord('Implicit solvent: running LEap to generate .prmtop and .crd for aptamer...')
                # os.system('pdb4amber {}.pdb > {}_amb_processed.pdb 2> {}_pdb4amber_out.log'.format(structureName, structureName, structureName))
                copyfile(aptamerPDB, structureName + '_amb_processed.pdb')
                copyfile('./leap_template.in', 'leap_freeAptamerMD.in')
                replaceText('leap_freeAptamerMD.in', 'mySTRUCTURE', structureName)
                os.system('tleap -f leap_freeAptamerMD.in > leap_freeAptamerMD.out')
                os.rename('leap.log', 'leap_freeAptamerMD.log')
                os.system('mv leap_freeAptamerMD.in outfiles_of_pipeline_modules/')
                os.system('mv leap_freeAptamerMD.out outfiles_of_pipeline_modules/')
                os.system('mv leap_freeAptamerMD.log outfiles_of_pipeline_modules/')
                # os.system('tail -1 leap_freeAptamerMD.out')  # show last line
                # printRecord(readFinalLines('leap_freeAptamerMD.out', 1))  # show last line. problematic            
                structureName += '_amb'  # after pdb4amber and saveAmberParm, will generate structureName_amb_processed.pdb/top/crd

            processedAptamer = structureName + '_processed.pdb'
            self.autoMD(runfilesDir=runfilesDir, structurePDB=processedAptamer, simTime=simTime, binding=False, implicitSolvent=implicitSolvent)  # run MD sampling either for a set amount of time or till converged to equilibrium sampling of RC's
            processedAptamerTrajectory = structureName + '_processed_complete_trajectory.dcd'  # this is output file of autoMD

        printRecord('Free aptamer MD sampling: simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print out sampling speed
        # self.checkRuntime()
        # printRecord("Checked time. Try cleaning...")
        printRecord("Try cleaning...")

        stru_in_runfiles = os.path.join(runfilesDir, processedAptamer)
        traj_in_runfiles = os.path.join(runfilesDir, processedAptamerTrajectory)        
        # move processed pdb used by OpenMM to runfilesDir (original trajectory was already written to there)
        os.replace(processedAptamer, stru_in_runfiles)

        if implicitSolvent is False:
            cleanTrajectory(stru_in_runfiles, traj_in_runfiles)  # clean up trajectory for later use. remove water and salt from trajectory. Prepend 'clean_' to files' basenames.
            printRecord("Cleaned MD sampling trajectory of free aptamer.")
            printRecord('>>> Generated: clean_' + processedAptamerTrajectory + ' in the folder ./' + runfilesDir)
            self.dcdDict['sampled aptamer {}'.format(self.i)] = 'clean_' + processedAptamerTrajectory
            self.pdbDict['sampled aptamer {}'.format(self.i)] = 'clean_' + processedAptamer
        else:  # no water or salt to remove
            # copyfile(processedAptamer, 'clean_' + processedAptamer)  # TODO no cleaning for now
            # copyfile(processedAptamerTrajectory, 'clean_' + processedAptamerTrajectory)  # TODO no cleaning for now
            printRecord("No need for cleaning in implicit solvent.")
            printRecord('>>> Generated MD sampling trajectory of free aptamer: ' + processedAptamerTrajectory + ' in the folder ./' + runfilesDir)
            self.dcdDict['sampled aptamer {}'.format(self.i)] = processedAptamerTrajectory
            self.pdbDict['sampled aptamer {}'.format(self.i)] = processedAptamer
        
        printRecord("Analyzing the MD sampling trajectory to look for free aptamer's representative configuration...")
        aptamerDict = self.analyzeTrajectory(os.path.join(runfilesDir, self.pdbDict['sampled aptamer {}'.format(self.i)]), os.path.join(runfilesDir, self.dcdDict['sampled aptamer {}'.format(self.i)]))
        # old issue: analyzeTraj --> getNucDAtraj --> Dihedral: raise ValueError("All AtomGroups must contain 4 atoms")
        self.pdbDict['representative aptamer {}'.format(self.i)] = 'repStructure_{}.pdb'.format(self.i)        
        printRecord('>>> Generated representative structure of free aptamer from the trajectory: {}.'.format(self.pdbDict['representative aptamer {}'.format(self.i)]))
        printRecord('Free aptamer MD sampling complete.')

        return aptamerDict

    def dock(self, aptamerPDB, targetPDB):
        """
        Use LightDock to run docking and isolate good complex structures
        :param aptamer: pdb file of aptamer after MD simulation
        :param targetPDB: pdb file of the target ligand. Either user provides pdb of target, or use the example pdb we provide in examples/
        :return:
        """
        printRecord('\nDocking Ligand with Aptamer')

        ld = interfaces.ld(aptamerPDB, targetPDB, self.params, self.i)  # ld is a new class, therefore need to pass in this class's params: self.params
        ld.run()
        topScores = ld.topScores  # empty list or contains float value of docking score

        for i in range(len(topScores)):
            complexFileName = 'complex_{}_{}.pdb'.format(self.i, int(i))

            # copy file such as lightdock_topscores_0/top_1.pdb as complex_{}_{}.pdb
            copyfile('lightdock_topscores_%d'%self.i + '/top_%d'%int(i+1) + '.pdb', complexFileName)
            self.pdbDict['binding complex {} {}'.format(self.i, int(i))] = complexFileName
            printRecord('>>> Generated aptamer-ligand complex structure: ' + complexFileName)

        self.topDockingScores = topScores
        if len(topScores) == 0:
            printRecord('Docking complete, no docking scores are obtained.')
            return 'Not found'
        elif len(topScores) < self.params['N_docked_structures']:
            printRecord('Number of identified docked structures is less than what you are looking for!')
            printRecord('Docking complete, best score(s): {}.'.format(topScores))
            return topScores
        else:
            printRecord('Docking complete, best score(s): {}.'.format(topScores))
            return topScores

    def complexDynamics(self, runfilesDir, complexPDB, simTime, implicitSolvent=False):
        """
        Run MD sampling for a complex structure and relevant analysis
        :param runfilesDir: a directory to store runfiles generated during MD, eg, 'md_complex_sampling_runfiles_0_0'
        :param implicitSolvent:
        :param complexPDB: complex structure's pdb like 'complex_1_2.pdb'. It contains no directory info
        :return:
        """    
        structureName = complexPDB.split('.')[0]  # complexPDB: eg "complex_1_2.pdb" or "complex_1_2_processed.pdb" (if resuming sampling)

        if self.params['pickup_from_complexChk'] is True:
            printRecord('\nResuming a previous dynamics of aptamer-ligand complex')
            self.autoMD(runfilesDir=runfilesDir, structurePDB=complexPDB, simTime=simTime, binding=True, implicitSolvent=implicitSolvent)  # In this case, complex_{}_{}_processed.pdb (due to resuming)
            processedComplex = complexPDB
            processedComplexTrajectory = structureName + '_complete_trajectory.dcd'  # this is output file of autoMD
        else:
            printRecord('\nRunning Fresh Dynamics of Aptamer-Ligand Complex')
            if implicitSolvent is False:
                # set up periodic box and condition: pH and ionic strength => protons, ions and their concentrations
                prepPDB(complexPDB, self.params['box_offset'], self.params['pH'], self.params['ionicStrength'], MMBCORRECTION=True, waterBox=True)
                printRecord('Done preparing files with waterbox. Start OpenMM.')
            else:  # prepare prmtop and crd file using LEap in ambertools
                printRecord('Implicit solvent: running LEap to generate .prmtop and .crd for aptamer-ligand complex...')
                # os.system('pdb4amber {}.pdb > {}_amb_processed.pdb 2> {}_pdb4amber_out.log'.format(structureName, structureName, structureName))
                copyfile(complexPDB, structureName + '_amb_processed.pdb')
                copyfile('./leap_template.in', 'leap_complexMD.in') 
                replaceText('leap_complexMD.in', 'mySTRUCTURE', structureName)
                # If target is peptide or RNA, will also add "leaprc.protein.ff14SB" or "source leaprc.RNA.OL3" to the leap_complexMD.in
                if self.targetType == 'peptide':  # 'peptide' or 'DNA' or 'RNA' or 'other'
                    replaceText('leap_complexMD.in', '# TARGET_IS_PEPTIDE ', '')
                elif self.targetType == 'RNA':
                    replaceText('leap_complexMD.in', '# TARGET_IS_RNA ', '')
                os.system('tleap -f leap_complexMD.in > leap_complexMD.out')
                os.rename('leap.log', 'leap_complexMD.log')
                os.system('mv leap_complexMD.in outfiles_of_pipeline_modules/')
                os.system('mv leap_complexMD.out outfiles_of_pipeline_modules/')
                os.system('mv leap_complexMD.log outfiles_of_pipeline_modules/')
                # os.system('tail -1 leap_complexMD.out')  # show last line
                # printRecord(readFinalLines('leap_complexMD.out', 1))  # show last line. problematic            
                structureName += '_amb'  # after pdb4amber and saveAmberParm, will generate structureName_amb_processed.pdb/top/crd

            processedComplex = structureName + '_processed.pdb'          
            self.autoMD(runfilesDir=runfilesDir, structurePDB=processedComplex, simTime=simTime, binding=True, implicitSolvent=implicitSolvent)  # run MD sampling either for a set amount of time or till converged to equilibrium sampling of RC's
            processedComplexTrajectory = structureName + '_processed_complete_trajectory.dcd'  # this is output file of autoMD          

        printRecord('Aptamer-ligand complex MD sampling: simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print out sampling speed
        # self.checkRuntime()
        # printRecord("Checked time. Try cleaning...")
        printRecord("Try cleaning...")

        stru_in_runfiles = os.path.join(runfilesDir, processedComplex)
        traj_in_runfiles = os.path.join(runfilesDir, processedComplexTrajectory)        
        # move processed pdb used by OpenMM to runfilesDir (original trajectory was already written to there)
        os.replace(processedComplex, stru_in_runfiles)

        if implicitSolvent is False:
            cleanTrajectory(stru_in_runfiles, traj_in_runfiles)  # clean up trajectory for later use. remove water and salt from trajectory. Prepend 'clean_' to files' basenames
            printRecord("Cleaned MD sampling trajectory of aptamer-ligand complex.")
            printRecord('>>> Generated: clean_' + processedComplexTrajectory + ' in the folder ./' + runfilesDir)
            self.dcdDict['sampled complex {} {}'.format(self.i, self.j)] = 'clean_' + processedComplexTrajectory
            self.pdbDict['sampled complex {} {}'.format(self.i, self.j)] = 'clean_' + processedComplex
        else:  # no water or salt to remove
            # copyfile(processedComplex, 'clean_' + processedComplex)  # TODO no cleaning for now
            # copyfile(processedComplexTrajectory, 'clean_' + processedComplexTrajectory)  # TODO no cleaning for now
            printRecord("No need for cleaning in implicit solvent.")
            printRecord('>>> Generated MD sampling trajectory of aptamer-ligand complex: ' + processedComplexTrajectory + ' in the folder ./' + runfilesDir)
            self.dcdDict['sampled complex {} {}'.format(self.i, self.j)] = processedComplexTrajectory
            self.pdbDict['sampled complex {} {}'.format(self.i, self.j)] = processedComplex

        # # Not carry out any analysis for current simulator. 
        # bindingDict = self.analyzeBinding(self.pdbDict['sampled complex {} {}'.format(self.i, self.j)],
        #                                   self.dcdDict['sampled complex {} {}'.format(self.i, self.j)],
        #                                   self.pdbDict['sampled aptamer {}'.format(self.i)],
        #                                   self.dcdDict['sampled aptamer {}'.format(self.i)])
        # # print findings
        # printRecord('Binding Results: Contact Persistence = {:.2f}, Contact Score = {:.2f}, Conformation Change = {:.2f}'.format(bindingDict['close contact ratio'], bindingDict['contact score'], bindingDict['conformation change']))        
        printRecord('Aptamer-ligand complex MD sampling complete.')
        bindingDict = {}
        return bindingDict

    def autoMD(self, runfilesDir, structurePDB, simTime, binding=False, implicitSolvent=False):
        """
        Run MD either for a set amount of time or till the dynamics 'converge'
        :param runfilesDir: a directory to store runfiles generated during MD, eg, "md_aptamer_sampling_runfiles_0" or "md_complex_sampling_runfiles_0"
        :param structurePDB: pdb file name of the structure; Could be free aptamer or aptamer-ligand complex. It contains no directory info
        :param binding: ``True`` - complex; ``False``: free aptamer
        :return:
        """        

        cutoff = self.params['autoMD_convergence_cutoff']

        structureName = structurePDB.split('.')[0]  # e.g., free aptamer: relaxedAptamer_0_amb_processed.pdb (implicit solvent) or relaxedAptamer_0_processed.pdb (explicit solvent)
        if (self.params['auto_sampling'] is False) or (binding is True):  # just run MD for the given sampling time
            # aptamer-ligand complex: no auto_sampling
            self.targetUnbound = False
            omm = interfaces.omm(runfilesDir=runfilesDir, structurePDB=structurePDB, params=self.params, simTime=simTime, binding=binding, implicitSolvent=implicitSolvent)
            self.ns_per_day = omm.doMD()  # run MD in OpenMM framework            
            
            # print('Generated:', structureName + '_trajectory.dcd')
            traj          = os.path.join(runfilesDir, structureName+'_trajectory.dcd')
            traj_complete = os.path.join(runfilesDir, structureName+'_complete_trajectory.dcd')
            os.replace(traj, traj_complete)
            # print('Replaced ^ with:', structureName + '_complete_trajectory.dcd')
            printRecord('>>> Generated: ' + structureName + '_complete_trajectory.dcd (auto_sampling was OFF) in the folder ./' + runfilesDir)

        elif (self.params['auto_sampling'] is True) and (binding is False):  # run MD till convergence: applicable to free aptamer sampling only for now
            maxIter = self.params['max_aptamer_sampling_iter']
            converged = False
            iter = 0
            # self.targetUnbound = False

            traj_1 = os.path.join(runfilesDir, structureName+'_trajectory-1.dcd')
            traj   = os.path.join(runfilesDir, structureName+'_trajectory.dcd')
            traj_complete = os.path.join(runfilesDir, structureName+'_complete_trajectory.dcd')

            while (converged is False) and (iter < maxIter):
                iter += 1
                omm = interfaces.omm(runfilesDir=runfilesDir, structurePDB=structurePDB, params=self.params, simTime=simTime, binding=False, implicitSolvent=implicitSolvent)
                self.ns_per_day = omm.doMD()

                if iter > 1:  # if we have multiple trajectory segments, combine them
                    appendTrajectory(structurePDB, traj_1, traj)
                    # os.replace('combinedTraj.dcd', structureName + '_trajectory-1.dcd')
                    os.replace('combinedTraj.dcd', traj_1)
                else:
                    # os.replace(structureName + '_trajectory.dcd', structureName + '_trajectory-1.dcd')  # in case we need to combine two trajectories
                    os.replace(traj, traj_1)  # in case we need to combine two trajectories

                # combinedSlope = checkTrajPCASlope(structurePDB, structureName + '_trajectory-1.dcd', self.params['print_step'])
                combinedSlope = checkTrajPCASlope(structurePDB, traj_1, self.params['print_step'])
                
                # # No binding analysis
                # if binding is True:
                #     if self.targetType == 'peptide':
                #         self.targetUnbound = checkMidTrajectoryBinding(structurePDB, structureName + '_trajectory-1.dcd', self.targetSeq, self.aptamerSeq, self.params, cutoffTime=1)
                #         if self.targetUnbound:
                #             printRecord('Peptide ligand came unbound!')
                # if (combinedSlope < cutoff) or (self.targetUnbound is True):
                #     converged = True
                if (combinedSlope < cutoff): converged = True

            # os.replace(structureName + '_trajectory-1.dcd', structureName + '_complete_trajectory.dcd')
            os.replace(traj_1, traj_complete)
            printRecord('>>> Generated: ' + structureName + '_complete_trajectory.dcd (auto_sampling was ON) in the folder ./' + runfilesDir)

    def analyzeTrajectory(self, structure, trajectory):
        """
        Analyze trajectory of free aptamer (MD sampling)
        :param structure: topology file of the free aptamer, eg, its PDB file.
        :param trajectory: free aptamer MD trajectory
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

        # 3D structure analysis
        mixedTrajectory = np.concatenate((baseDistTraj.reshape(len(baseDistTraj), int(baseDistTraj.shape[-2] * baseDistTraj.shape[-1])), nucleicAnglesTraj.reshape(len(nucleicAnglesTraj), int(nucleicAnglesTraj.shape[-2] * nucleicAnglesTraj.shape[-1]))), axis=1)  # mix up all our info
        representativeIndex, pcTrajectory, eigenvalues = isolateRepresentativeStructure(mixedTrajectory)

        # save this structure a separate file
        extractFrame(structure, trajectory, representativeIndex, 'repStructure_%d' % self.i + '.pdb')
        if self.actionDict['do 2d analysis'] is True:
            printRecord('Predicted 2D structure                      :' + self.ssAnalysis['displayed 2d string'][self.i])  
            # if we did not run 2d analysis, we probably do not know the secondary structure.        
        printRecord('Equilibrium 2D structure from the trajectory:' + configToString(secondaryStructure))

        # we also need function which processes the energies spat out by the trajectory logs
        analysisDict = {}  # compile our results
        analysisDict['2D structure'] = secondaryStructure
        analysisDict['RC trajectories'] = pcTrajectory
        analysisDict['representative structure index'] = representativeIndex

        return analysisDict

    # def analyzeBinding(self, bindStructure, bindTrajectory, freeStructure, freeTrajectory):
    #     """
    #     Analyze trajectory of both free aptamer and aptamer-target complex to find out about binding information
    #     :param bindStructure: complex structure's PDB
    #     :param bindTrajectory: complex structure's MD trajectory
    #     :param freeStructure: free aptamer structure's PDB
    #     :param freeTrajectory: free aptamer structure's MD trajectory
    #     :return:
    #     """
    #     bindu = mda.Universe(bindStructure, bindTrajectory)
    #     freeu = mda.Universe(freeStructure, freeTrajectory)

    #     if self.targetType == 'peptide':
    #         bindingDict = bindingAnalysis(bindu, freeu, self.targetSeq, self.aptamerSeq)  # look for contacts between target ligand and aptamer.
    #         if self.targetUnbound:
    #             bindingDict['target came unbound'] = True
    #         else:
    #             bindingDict['target came Unbound'] = False

    #         return bindingDict

    # =======================================================
    # ====== supporting functions are followed (part2) ======
    # =======================================================
    
    # def checkRuntime(self):
    #     """
    #     Estimate the maximum runtime for this code, given a ns/day calculation from existing work done
    #     :return:
    #     """
    #     projectedMaxRuntime = self.computeRuntime()
    #     buffer = 0.9  # ideally would like a buffer of 20% of the walltime
    #     doReduction = 0
    #     # if the code is going to take too long, shorten the sampling time until it fits
    #     while projectedMaxRuntime > (buffer * self.params['max_walltime']):  # max_walltime is in hours
    #         doReduction = 1
    #         self.params['sampling_time'] *= 0.95
    #         projectedMaxRuntime = self.computeRuntime()

    #     if doReduction == 1:
    #         if self.params['sampling_time'] < 0.1:
    #             printRecord('Sampling time reduced to 100 ps - this run is too expensive, and will now terminate')
    #             self.terminateRun()
    #         else:
    #             printRecord('Binding sampling time chunk reduced to {} ns'.format(self.params['sampling_time']))

    # def computeRuntime(self):
    #     dt = self.params['sampling_time']  # ns in MD sampling after getting (smoothed) 3D structure
    #     ssStructures = self.params['N_2D_structures']
    #     dockingStructures = self.params['N_docked_structures']
    #     aptamerSteps = self.params['max_aptamer_sampling_iter']
    #     maxNS = dt * (ssStructures * aptamerSteps + ssStructures * dockingStructures)  # maximum possible sampling_time in nanoseconds
    #     # complexSteps = self.params['max_complex_sampling_iter']
    #     # maxNS = dt * (ssStructures * aptamerSteps + ssStructures * dockingStructures * complexSteps)  # maximum possible sampling_time in nanoseconds
    #     projectedMaxRuntime = maxNS / self.ns_per_day * 24  # how many hours needed to complete our maximum sampling_time

    #     return projectedMaxRuntime
    
    def terminateRun(self):
        """ for some reason, the run needs to end """
        sys.exit()
