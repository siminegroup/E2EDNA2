import sys
import glob
from shutil import copyfile, copytree

import interfaces
from utils import *
from analysisTools import *

# noinspection PyPep8Naming  #meaning??


class opendna:
    def __init__(self, params):
        self.workDir = ""  # rather use an empty string, not an empty list
        self.params = params
        self.aptamerSeq = self.params['aptamerSeq']        
        self.targetPDB = self.params['target ligand']        # either pdb filename (eg, "peptide.pdb") or boolean ``False``
        self.targetType = self.params['target ligand type']  # 'peptide' or 'DNA' or 'RNA' or 'other'. Assume now the target is a peptide
        self.targetSeq = self.params['target sequence']
        # empty string, if target is not a peptide/DNA/RNA: how to deal with the "analyzeBinding"? No analysis in this simulator.
        
        self.pdbDict = {}
        self.dcdDict = {}
        self.actionDict = {}
        self.getActionDict()  # specify the actions based on the selected mode
        # modify actions based on three toggles.
        if self.params['skip MMB'] is True:  
            self.actionDict['do 2d analysis'] = False  
            self.actionDict['do MMB'] = False        
        if self.params['pick up from freeAptamerChk'] is True:
            # everything before "freeAptamerDynamics" is skipped
            self.actionDict['do 2d analysis'] = False
            self.actionDict['do MMB'] = False            
            self.actionDict['do smoothing'] = False
            self.actionDict['get equil repStructure'] = True  # making sure to run "freeAptamerDynamics"
            # # Not used for now:
            # printRecord('Copying the checkpoint file: {} to workdir.'.format(self.params['chk file']), self.workDir+'/')
            # if os.path.exists(self.params['chk file']):
            #     copyfile(self.params['chk file'], self.workDir + '/' + self.params['chk file'])
            # else:
            #     printRecord('Designated checkpoint file does not exist! Terminating the pipeline.', self.workDir+'/')
            #     self.terminateRun()
        elif self.params['pick up from complexChk'] is True:
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
            self.params['skip smoothing'] = False

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
            self.params['skip smoothing'] = False

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

    def setup(self):
        """
        set up a working directory, copy in relevant xyz and keyfiles, move to the working directory
        :return:
        """
        if (self.params['explicit run enumeration'] is True) or (self.params['run num'] == 0):  # make a new directory
            if self.params['run num'] == 0:  # auto-increment the max run_num folder
                self.makeNewWorkingDirectory()
            else:  # create a working directory defined by the user
                self.workDir = self.params['workdir'] + '/run%d' % self.params['run num']
                os.mkdir(self.workDir)
                printRecord('Starting Fresh Run {}'.format(self.params['run num']), self.workDir+'/')
        else:
            self.workDir = self.params['workdir'] + '/' + 'run%d' % self.params['run num']

        # copy relevant files to the workDir
        os.mkdir(self.workDir + '/outfiles')
                
        # If we skip MMB because we have a folded structure to start with
        if self.params['skip MMB'] is True:            
            printRecord('Skipping 2D analysis and MMB folding. Starting with an folded structure.', self.workDir+'/')            
            if not os.path.isfile('./' + self.params['folded initial structure']):
            # if not os.path.exists(self.params['folded initial structure']):
                printRecord('Designated folded structure does not exist! Terminating the pipeline.', self.workDir+'/')
                self.terminateRun()  # maybe do it via "raise an error"?
            else:
                printRecord('Copying the folded structure: {} to workdir.'.format(self.params['folded initial structure']), self.workDir+'/')
                copyfile(self.params['folded initial structure'], self.workDir + '/foldedAptamer_0.pdb')
        else:
            # copy MMB params and command script
            copyfile(self.params['mmb params'], self.workDir + '/parameters.csv')
            copyfile(self.params['mmb normal template'], self.workDir + '/commands.template.dat')
            copyfile(self.params['mmb quick template'], self.workDir + '/commands.template_quick.dat')
            copyfile(self.params['mmb long template'], self.workDir + '/commands.template_long.dat')        
        
        # if do docking, need targetPDB
        if (self.actionDict['do docking'] is True):  
            if self.targetPDB == 'False':  # if user does not provide a target ligand pdb, use our provided example target ligand: a peptide.
                self.targetPDB = 'example_target_peptide.pdb'             # self.targetPDB is no longer False.
                self.targetSeq = self.params['example peptide sequence']  # self.targetSeq is not empty because we are dealing with a peptide.
                copyfile(self.params['example target pdb'], self.workDir + '/' + self.targetPDB)  # ie, target of the aptamer
            else:  # user provides a targetPDB
                copyfile(self.targetPDB, self.workDir + '/' + self.targetPDB)
            # From now on, "self.targetPDB == 'False'" is equivalent to "'do docking' is False".
            # copy lightdock scripts        
            copytree('lib/lightdock', self.workDir + '/ld_scripts')  # if destination dir does not already exist, create one then copy the whole source dir

        # copy csv file (dihedral restraints) if target ligand is a peptide and constraint is enabled
        if self.params['peptide backbone constraint constant'] != 0:
            copyfile('backbone_dihedrals.csv', self.workDir + '/backbone_dihedrals.csv')

        if self.params['implicit solvent'] is True:
            copyfile(self.params['leap template'], self.workDir + '/leap_template.in')

        # move to working dir
        os.chdir(self.workDir)
        
        # specify the DNA force field for free aptamer MD sampling if implicit solvent. Assume only dealing with a DNA aptamer.
        if self.params['implicit solvent'] is True:
            replaceText('leap_template.in', 'DNA_FF', self.params['DNA force field'])
        
        # Print prompt: free aptamer or aptamer + target ligand
        printRecord('Simulation mode: {}'.format(self.params['mode']))
        if (self.actionDict['do docking'] is False):  # no given target ligand, nor do docking
            printRecord('Simulating free aptamer: {}'.format(self.aptamerSeq))
        else:
            printRecord('Simulating {} with {}'.format(self.aptamerSeq, self.targetPDB))

        if (self.params['skip MMB'] is False) and (self.params['device'] == 'local'):
            if (self.params['device platform'] == 'linux') or (self.params['device platform'] == 'WSL'):
                os.environ["LD_LIBRARY_PATH"] = self.params['mmb dir']  # Add MMB library's path to environment variable
                # When running MMB on 'local':
                    # Windows: no action needed. But this pipeline does not support Windows OS, due to NUPACK.
                    # WSL or linux: explicitly add MMB library path.
                    # macos: no action needed here.                                
                # TL warning: do NOT use os.environ["DYLD_LIBRARY_PATH"] for macos! It does not help MMB locate the library AND it confuses OpenMM.
                # TL tip: refer to 'mmb' class in interfaces.py on how to run MMB on macos machine.
                
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
            printRecord('Starting Fresh Run %d' % (prev_max + 1), self.workDir+'/')
        else:
            self.workDir = self.params['workdir'] + '/' + 'run1'
            os.mkdir(self.workDir)
            printRecord('Starting Fresh Run 1', self.workDir+'/')

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
        
        if self.actionDict['do 2d analysis']:  # get secondary structure
            self.pairLists = self.getSecondaryStructure(self.aptamerSeq)
            outputDict['2d analysis'] = self.ssAnalysis
            np.save('opendnaOutput', outputDict)  # save 2d structure results
            printRecord('Running over %d' % len(self.pairLists) + ' possible 2D structures.')
            num_2dSS = len(self.pairLists)
        
        else:  # there are three cases where 2d analysis is not carried out.
            if self.params['pick up from complexChk'] is True:
                printRecord('Resume a sampling of aptamer-ligand complex from .chk file, therefore skip all the steps before it.')
            elif self.params['pick up from freeAptamerChk'] is True:  
                printRecord('Resume a sampling of free aptamer from .chk file, therefore skip all the steps before it.')                
            else:
                printRecord('Starting with an existing folded strcuture, therefore skip 2d analysis and MMB folding.')
            num_2dSS = 1  # quick and dirty

        for self.i in range(num_2dSS):  # loop over all possible secondary structures
            if self.actionDict['do 2d analysis'] is True:  # self.ssAnalysis only exists if we "do 2d analysis"            
                printRecord('2D structure #{} is                              : {}'.format(self.i, self.ssAnalysis['displayed 2d string'][self.i]))
                self.pairList = np.asarray(self.pairLists[self.i])  # be careful!!!: .pairList vs. .pairLists                                                

            if self.actionDict['do MMB']:  # fold 2D into 3D
                self.foldSequence(self.aptamerSeq, self.pairList)
            elif self.params['skip MMB'] is True:   
                # situation here: no MMB + not to resume a MD of either free aptamer or aptamer-ligand complex: ie, (self.params['pick up from freeAptamerChk'] is False) and (self.params['pick up from complexChk'] is False)
                # meaning of the condition: just skip MMb and preceed to next step: MD smoothing or MD sampling or dock (if "coarse dock")
                self.pdbDict['folded aptamer {}'.format(self.i)] = self.params['folded initial structure']
                self.pdbDict['representative aptamer {}'.format(self.i)] = self.params['folded initial structure']  # for "coarse dock" mode.            

            if self.params['skip smoothing'] is True:
                self.actionDict['do smoothing'] = False
                if self.params['mode'] != '2d structure':
                    printRecord('\nNo relaxation (smoothing) of the folded aptamer.')
                self.pdbDict['relaxed aptamer {}'.format(self.i)] = 'foldedAptamer_0.pdb'
            elif self.actionDict['do smoothing']:  # just to double check
                if self.params['skip MMB'] is False:
                    self.MDSmoothing(self.pdbDict['mmb folded aptamer {}'.format(self.i)], relaxationTime=self.params['smoothing time'], implicitSolvent=self.params['implicit solvent'])  # relax for xx nanoseconds
                else:  # if given a folded structure, skip MMB and directly use the provided structure
                    self.MDSmoothing(self.pdbDict['folded aptamer {}'.format(self.i)], relaxationTime=self.params['smoothing time'], implicitSolvent=self.params['implicit solvent'])  # relax for xx nanoseconds

            if self.actionDict['get equil repStructure']:  # definitely did smoothing if want an equil structure                    
                if self.params['pick up from freeAptamerChk'] is False:
                    outputDict['free aptamer results {}'.format(self.i)] = self.freeAptamerDynamics(self.pdbDict['relaxed aptamer {}'.format(self.i)], implicitSolvent=self.params['implicit solvent'])
                else:  # user must make sure they have provided a .chk file
                    outputDict['free aptamer results {}'.format(self.i)] = self.freeAptamerDynamics(self.params['resumed structurePDB'], implicitSolvent=self.params['implicit solvent'])
                    # if using implicit solvent, the .top and .crd files have the same name as .pdb. For ex: relaxed_amb_processed.pdb/top/crd
                np.save('opendnaOutput', outputDict)  # save outputs

            if (self.actionDict['do docking'] is True) and (self.targetPDB != 'False'):  # find docking configuration for the complexed structure
                # coarse dock: no smoothing;
                # smooth dock: smooth + dock;
                # full dock: smooth + dock + equil structure;
                # full binding: smooth + dock + equil structure + sampling dynamics;                
                outputDict['dock scores {}'.format(self.i)] = self.dock(self.pdbDict['representative aptamer {}'.format(self.i)], self.targetPDB)  # eg, "peptide.pdb" which can be created given peptide sequence by buildPeptide in function dock
                # pdbDict['representative aptamer {}' is defined at MMB folding, MD smoothing or freeAptamerDynamics
                # TODO: does lightdock also support Amber implicit solvent model?
                np.save('opendnaOutput', outputDict)  # save outputs
                # N docked structures are specified by user: how many docked structures do we want to investigate
                # num_docked_structure = self.params['N docked structures']
                num_docked_structure = len(self.topDockingScores)
                
                printRecord('Running over %d' % num_docked_structure + ' docked structures.')
            
            if self.params['pick up from complexChk'] is True:  # it means we did not dock but we have a complex structure.
                num_docked_structure = 1
                printRecord('Resuming the sampling of 1 docked aptamer-ligand structure')

            if self.actionDict['do binding'] is True:  # run MD on the aptamer-ligand structure
                for self.j in range(num_docked_structure):  # loop over docking configurations for a given secondary structure                
                    printRecord('Docked structure #{}'.format(self.j))                
                    if self.params['pick up from complexChk'] is False:
                        outputDict['binding results {} {}'.format(self.i, self.j)] = self.complexDynamics(self.pdbDict['binding complex {} {}'.format(self.i, int(self.j))], implicitSolvent=self.params['implicit solvent'])
                        # TODO why need int(self.j)? Why sometimes %d % string, but sometimes {}.format?
                        # np.save('opendnaOutput', outputDict)
                    else:  # user must make sure they have provided a .chk file
                        outputDict['binding results {} {}'.format(self.i, self.j)] = self.complexDynamics(self.params['resumed structurePDB'], implicitSolvent=self.params['implicit solvent'])
                        # if using implicit solvent, the .top and .crd files have the same name as .pdb. For ex: 'complex_1_2_processed.pdb/top/crd'
                    np.save('opendnaOutput', outputDict)
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
        printRecord("Getting Secondary Structure(s)")
        if self.params['secondary structure engine'] == 'seqfold':
            ssString, pairList = getSeqfoldStructure(aptamerSeq, self.params['temperature'])  # seqfold guess
            self.ssAnalysis = {'2d string': ssString, 'pair list': pairList, 'displayed 2d string': ssString}
            # print('seqfold ssAnalysis: ',self.ssAnalysis)
            return pairList
        elif self.params['secondary structure engine'] == 'NUPACK':
            nup = interfaces.nupack(aptamerSeq, self.params['temperature'], self.params['ionicStrength'], self.params['[Mg]'])  # initialize nupack
            self.ssAnalysis = nup.run()  # run nupack analysis of possible 2D structures.

            distances = getSecondaryStructureDistance(self.ssAnalysis['config'])            
            if len(distances) > 1:
                topReps, topProbs, topDists = do2DAgglomerativeClustering(self.ssAnalysis['config'], self.ssAnalysis['state prob'], distances)
                topPairLists = []                
                nLists = min(len(topReps), self.params['N 2D structures'])
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
        printRecord("Folding Aptamer from Sequence. Fold speed={}".format(self.params['fold speed']))        

        mmb = interfaces.mmb(aptamerSeq, pairList, self.params, self.i)
        MDAfoldFidelity, MMBfoldFidelity = mmb.run()
        foldFidelity = MMBfoldFidelity
        printRecord('Initial fold fidelity = %.3f' % foldFidelity)
        printRecord('Initial fold fidelity = %.3f (from MDAnalysis)' % MDAfoldFidelity)
        attempts = 1
        # Skip refold
        # if self.params['fold speed'] != 'quick':  # don't rerun if we're quick folding
        #     intervalLength = 5  # initial interval length
        #     while (foldFidelity <= self.params['foldFidelity']) and (attempts < 5):  # if it didn't fold properly, try again with a longer annealing time - up to XX times
        #         self.params['fold speed'] = 'long'
                
        #         printRecord("Refolding Aptamer from Sequence")
        #         mmb = interfaces.mmb(aptamerSeq, pairList, self.params, self.i, intervalLength)  # extra intervalLength argument                
        #         MDAfoldFidelity, MMBfoldFidelity = mmb.run()
        #         foldFidelity = MMBfoldFidelity
        #         printRecord('Subsequent fold fidelity = %.3f' % foldFidelity)
        #         printRecord('Subsequent MDA fold fidelity = %.3f' % MDAfoldFidelity)
        #         attempts += 1
        #         intervalLength += 5  # on repeated attempts, try harder

        self.pdbDict['mmb folded aptamer {}'.format(self.i)] = mmb.foldedAptamerSeq  # mmb.foldedAptamerSeq = foldedAptamer_{}.pdb: defined in the mmb.run()
        os.system('mv commands.run* ' + mmb.fileDump)  # mmb.fileDump is a directory for intermediate files during running MMB
        printRecord("Folded the aptamer and generated the folded structure: {}".format(self.pdbDict['mmb folded aptamer {}'.format(self.i)]))
        
        # For "coarse dock" mode:
        self.pdbDict['representative aptamer {}'.format(self.i)] = self.pdbDict['mmb folded aptamer {}'.format(self.i)]  # ie, 'relaxedAptamer_{}.pdb'.format(self.i)

    def MDSmoothing(self, structure, relaxationTime=0.01, implicitSolvent=False):  # default relaxation time is 0.01 ns
        """
        Run a short MD sampling to relax the coarse MMB structure
        Relaxation time in nanoseconds, print time in picoseconds
        :param structure:
        :param relaxationTime: params['smoothing time']
        :param implicitSolvent:
        :return:
        """
        structureName = structure.split('.')[0]
        printRecord('\nRunning relaxation (smoothing)')
        if implicitSolvent is False:
            # set up periodic box and condition: pH and ionic strength => protons, ions and their concentrations
            prepPDB(structure, self.params['box offset'], self.params['pH'], self.params['ionicStrength'], MMBCORRECTION=True, waterBox=True)
            print('Done preparing files with waterbox. Start openmm.')
        else:  # prepare prmtop and crd file using LEap in ambertools
            printRecord('Implicit solvent: running LEap to generate .prmtop and .crd for folded aptamer...')
            #os.system('pdb4amber {}.pdb > {}_amb_processed.pdb 2> {}_pdb4amber_out.log'.format(structureName, structureName, structureName))
            copyfile(structure, structureName + '_amb_processed.pdb')
            copyfile('./leap_template.in', 'leap_MDSmoothing.in')
            replaceText('leap_MDSmoothing.in', 'mySTRUCTURE', structureName)
            os.system('tleap -f leap_MDSmoothing.in > leap_MDSmoothing.out')
            os.system('tail -1 leap_MDSmoothing.out')  # show last line
            # printRecord(readFinalLines('leap_MDSmoothing.out', 1))  # show last line. problematic
            structureName += '_amb'  # after pdb4amber and saveAmberParm, will generate structureName_amb_processed.pdb/top/crd

        processedStructure = structureName + '_processed.pdb'
        processedStructureTrajectory = structureName + '_processed_trajectory.dcd'
        omm = interfaces.omm(structurePDB=processedStructure, params=self.params, simTime=relaxationTime, binding=False, implicitSolvent=implicitSolvent)
        self.ns_per_day = omm.doMD()  # run MD in OpenMM framework

        printRecord('Pre-relaxation simulation speed %.1f' % self.ns_per_day + 'ns/day')  # print out sampling speed

        if implicitSolvent is False:
            cleanTrajectory(processedStructure, processedStructureTrajectory)  # remove water and salt from trajectory
            printRecord("Cleaned the trajectory of relaxing free aptamer. Extracting the final frame.")
        else:  # no water or salt to remove
            copyfile(processedStructure, 'clean_' + processedStructure)
            copyfile(processedStructureTrajectory, 'clean_' + processedStructureTrajectory)
            printRecord("No cleaning in implicit solvent, just copied smoothing traj. Extracting the final frame")

        self.dcdDict['relaxed aptamer {}'.format(self.i)] = 'clean_' + processedStructureTrajectory
        self.pdbDict['relaxed aptamer {}'.format(self.i)] = 'relaxedAptamer_{}.pdb'.format(self.i)  # specify the file name for the extracted final frame, ie, relaxed aptamer
        # extractFrame(processedStructure, processedStructureTrajectory, -1, self.pdbDict['relaxed aptamer {}'.format(self.i)])  # extract final frame        
        # # Current warning from MDA: UserWarning: Unit cell dimensions not found. CRYST1 record set to unitary values. I think: MDA cannot find unit cell dimension info. Can we add info when using implicit solvent?            
        # Another possible way to extract the last frame is to use OpenMM:            
        omm.extractLastFrame(self.pdbDict['relaxed aptamer {}'.format(self.i)])
        
        # in "smooth dock" mode:
        self.pdbDict['representative aptamer {}'.format(self.i)] = self.pdbDict['relaxed aptamer {}'.format(self.i)]  # ie, 'relaxedAptamer_{}.pdb'.format(self.i)

    def freeAptamerDynamics(self, aptamerPDB, implicitSolvent=False):
        """
        Run MD sampling for free aptamer and relevant analysis
        :param implicitSolvent:
        :param aptamerPDB: PDB file name of the aptamer to be simulated
        :return:
        """
        structureName = aptamerPDB.split('.')[0]  # aptamerPDB: eg "relaxedAptamer_0.pdb" or "relaxedAptamer_0_processed.pdb" (if resuming sampling)

        if self.params['pick up from freeAptamerChk'] is True:
            printRecord('\nResuming a previous free aptamer dynamics')
            self.autoMD(structurePDB=aptamerPDB, binding=False, implicitSolvent=implicitSolvent)  # In this case, relaxedAptamer_0_processed.pdb (due to resuming)
            processedAptamer = aptamerPDB
            processedAptamerTrajectory = structureName + '_complete_trajectory.dcd'  # this is output file of autoMD
        else:
            printRecord('\nRunning a fresh free aptamer dynamics')
            if implicitSolvent is False:
                # set up periodic box and condition: pH and ionic strength => protons, ions and their concentrations
                prepPDB(aptamerPDB, self.params['box offset'], self.params['pH'], self.params['ionicStrength'], MMBCORRECTION=True, waterBox=True)
            else:  # prepare prmtop and crd file using LEap in ambertools
                printRecord('Implicit solvent: running LEap to generate .prmtop and .crd for relaxed aptamer...')
                # os.system('pdb4amber {}.pdb > {}_amb_processed.pdb 2> {}_pdb4amber_out.log'.format(structureName, structureName, structureName))
                copyfile(aptamerPDB, structureName + '_amb_processed.pdb')
                copyfile('./leap_template.in', 'leap_freeAptamerMD.in')
                replaceText('leap_freeAptamerMD.in', 'mySTRUCTURE', structureName)
                os.system('tleap -f leap_freeAptamerMD.in > leap_freeAptamerMD.out')
                os.system('tail -1 leap_freeAptamerMD.out')  # show last line
                # printRecord(readFinalLines('leap_freeAptamerMD.out', 1))  # show last line. problematic            
                structureName += '_amb'  # after pdb4amber and saveAmberParm, will generate structureName_amb_processed.pdb/top/crd

            processedAptamer = structureName + '_processed.pdb'
            self.autoMD(structurePDB=processedAptamer, binding=False, implicitSolvent=implicitSolvent)  # run MD sampling either for a set amount of time or till converged to equilibrium sampling of RC's
            processedAptamerTrajectory = structureName + '_processed_complete_trajectory.dcd'  # this is output file of autoMD

        printRecord('Free aptamer simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print out sampling speed
        self.checkRuntime()
        print("\nChecked time. Try cleaning")
        if implicitSolvent is False:
            cleanTrajectory(processedAptamer, processedAptamerTrajectory)  # clean up trajectory for later use. by doing what?
            printRecord("Cleaned free aptamer traj. Looking for free aptamer's representative configuration")
        else:  # no water or salt to remove
            copyfile(processedAptamer, 'clean_' + processedAptamer)  # TODO no cleaning for now
            copyfile(processedAptamerTrajectory, 'clean_' + processedAptamerTrajectory)  # TODO no cleaning for now
            printRecord("No cleaning in implicit solvent, just copied free aptamer traj. Looking for free aptamer's representative configuration")
            
        self.dcdDict['sampled aptamer {}'.format(self.i)] = 'clean_' + processedAptamerTrajectory
        self.pdbDict['sampled aptamer {}'.format(self.i)] = 'clean_' + processedAptamer
        aptamerDict = self.analyzeTrajectory(self.pdbDict['sampled aptamer {}'.format(self.i)], self.dcdDict['sampled aptamer {}'.format(self.i)])
        # Within analyzeTrajectory, the last step is also to save an representative frame. We can also replace it using OpenMM??
            # self.omm.extractLastFrame('repStructure_%d' % self.i + '.pdb', representativeIndex) # need more scripting to complete it
        # aptamerDict = {}
        # TODO: analyzeTraj --> getNucDAtraj --> Dihedral: raise ValueError("All AtomGroups must contain 4 atoms")        
        self.pdbDict['representative aptamer {}'.format(self.i)] = 'repStructure_{}.pdb'.format(self.i)        
        printRecord('Generated representative structure of free aptamer: {}.'.format(self.pdbDict['representative aptamer {}'.format(self.i)]))
        printRecord('Free aptamer sampling complete.')

        return aptamerDict

    def dock(self, aptamerPDB, targetPDB):
        """
        Use LightDock to run docking and isolate good complex structures
        :param aptamer: pdb file of aptamer after MD simulation
        :param targetPDB: pdb file of the target ligand
        :return:
        """
        printRecord('\nDocking')
        if bool(self.params['peptide backbone constraint constant']):  # it constant != 0, bool=True
            printRecord('Peptide will be constrained on their dihidral angles.')
        
        # build a peptide from its sequence: no longer needed. Either user provides pdb of target (if it is peptide), or use the example pdb we provide
        # if self.params['build a peptide as target ligand'] is True: buildPeptide(self.targetSeq, self.targetPDB, customAngles=bool(self.params['peptide backbone constraint constant']))
        # No longer use buildPeptide in our pipeline, user must provide a 3D structure of the target ligand. We don't take responsibility of creating a target for them.
        # Now we have a pdb file named self.targetPDB which represents our target, provided by the user

        ld = interfaces.ld(aptamerPDB, targetPDB, self.params, self.i)  # ld is a new class, therefore need to pass in this class's params: self.params
        ld.run()
        topScores = ld.topScores

        for i in range(len(topScores)):
            copyfile('top_%d' % self.i + '/top_%d' % int(i+1) + '.pdb', 'complex_{}_{}.pdb'.format(self.i, int(i)))                
            self.pdbDict['binding complex {} {}'.format(self.i, int(i))] = 'complex_{}_{}.pdb'.format(self.i, int(i))
                
        self.topDockingScores = topScores
        if len(topScores) == 0:
            printRecord('Docking complete, no docking scores are obtained.')
            return 'Not found'
        elif len(topScores) < self.params['N docked structures']:
            printRecord('Number of identified docked structures is less than what you are looking for!')
            printRecord('Docking complete, best score(s): {}.'.format(topScores))
            return topScores
        else:
            printRecord('Docking complete, best score(s): {}.'.format(topScores))
            return topScores

    def complexDynamics(self, complexPDB, implicitSolvent=False):
        """
        Run MD sampling for a complex structure and relevant analysis
        :param implicitSolvent:
        :param complexPDB: complex structure's pdb like 'complex_1_2.pdb'
        :return:
        """    
        structureName = complexPDB.split('.')[0]  # complexPDB: eg "complex_1_2.pdb" or "complex_1_2_processed.pdb" (if resuming sampling)

        if self.params['pick up from complexChk'] is True:
            printRecord('\nResuming a previous aptamer-ligand dynamics')
            self.autoMD(structurePDB=complexPDB, binding=True, implicitSolvent=implicitSolvent)  # In this case, complex_1_2_processed.pdb (due to resuming)
            processedComplex = complexPDB
            processedComplexTrajectory = structureName + '_complete_trajectory.dcd'  # this is output file of autoMD
        else:
            printRecord('\nRunning a fresh aptamer-ligand dynamics')
            if implicitSolvent is False:
                # set up periodic box and condition: pH and ionic strength => protons, ions and their concentrations
                prepPDB(complexPDB, self.params['box offset'], self.params['pH'], self.params['ionicStrength'], MMBCORRECTION=True, waterBox=True)
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
                os.system('tail -1 leap_complexMD.out')  # show last line
                # printRecord(readFinalLines('leap_complexMD.out', 1))  # show last line. problematic            
                structureName += '_amb'  # after pdb4amber and saveAmberParm, will generate structureName_amb_processed.pdb/top/crd

            processedComplex = structureName + '_processed.pdb'          
            self.autoMD(structurePDB=processedComplex, binding=True, implicitSolvent=implicitSolvent)  # run MD sampling either for a set amount of time or till converged to equilibrium sampling of RC's
            processedComplexTrajectory = structureName + '_processed_complete_trajectory.dcd'  # this is output file of autoMD          

        printRecord('Complex simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print out sampling speed
        self.checkRuntime()
        print("\nChecked time. Try cleaning")
        if implicitSolvent is False:
            cleanTrajectory(processedComplex, processedComplexTrajectory)  # clean up trajectory for later use. by doing what?
            printRecord("Cleaned complex traj.")    
        else:  # no water or salt to remove
            copyfile(processedComplex, 'clean_' + processedComplex)  # TODO no cleaning for now
            copyfile(processedComplexTrajectory, 'clean_' + processedComplexTrajectory)  # TODO no cleaning for now
            printRecord("No cleaning in implicit solvent, just copied complex traj.")

        self.dcdDict['sampled complex {} {}'.format(self.i, self.j)] = 'clean_' + processedComplexTrajectory
        self.pdbDict['sampled complex {} {}'.format(self.i, self.j)] = 'clean_' + processedComplex

        # # Not carry out any analysis for current simulator. 
        # bindingDict = self.analyzeBinding(self.pdbDict['sampled complex {} {}'.format(self.i, self.j)],
        #                                   self.dcdDict['sampled complex {} {}'.format(self.i, self.j)],
        #                                   self.pdbDict['sampled aptamer {}'.format(self.i)],
        #                                   self.dcdDict['sampled aptamer {}'.format(self.i)])
        # # print findings
        # printRecord('Binding Results: Contact Persistence = {:.2f}, Contact Score = {:.2f}, Conformation Change = {:.2f}'.format(bindingDict['close contact ratio'], bindingDict['contact score'], bindingDict['conformation change']))
        # # TODO: why {:.2f} instead of {.2f}?
        printRecord('Complex sampling complete. Success! Exiting pipeline.')
        bindingDict = {}
        return bindingDict

    def autoMD(self, structurePDB, binding=False, implicitSolvent=False):
        """
        Run MD either for a set amount of time or till the dynamics 'converge'
        :param structurePDB: pdb file name of the structure; Could be free aptamer or aptamer-ligand complex
        :param binding: ``True`` - complex; ``False``: free aptamer
        :return:
        """
        if binding is True:  # whether it's complex or free aptamer
            maxIter = self.params['max complex sampling iterations']            
        else:
            maxIter = self.params['max aptamer sampling iterations']
        cutoff = self.params['autoMD convergence cutoff']

        structureName = structurePDB.split('.')[0]  # e.g., free aptamer: relaxedAptamer_0_amb_processed.pdb (implicit solvent) or relaxedAptamer_0_processed.pdb (explicit solvent)
        if self.params['auto sampling'] is False:  # just run MD for the given sampling time
            self.targetUnbound = False
            omm = interfaces.omm(structurePDB=structurePDB, params=self.params, binding=binding, implicitSolvent=implicitSolvent)
            self.ns_per_day = omm.doMD()  # run MD in OpenMM framework
            # print('Generated:', structureName + '_trajectory.dcd')
            os.replace(structureName + '_trajectory.dcd', structureName + "_complete_trajectory.dcd")
            # print('Replaced ^ with:', structureName + '_complete_trajectory.dcd')
            printRecord('Generated: ' + structureName + '_complete_trajectory.dcd (autosampling was OFF)')

        elif self.params['auto sampling'] is True:  # run MD till convergence (equilibrium)
            converged = False
            iter = 0
            self.targetUnbound = False

            while (converged is False) and (iter < maxIter):
                iter += 1
                omm = interfaces.omm(structurePDB=structurePDB, params=self.params, binding=binding, implicitSolvent=implicitSolvent)
                self.ns_per_day = omm.doMD()

                if iter > 1:  # if we have multiple trajectory segments, combine them
                    appendTrajectory(structurePDB, structureName + '_trajectory-1.dcd', structureName + '_trajectory.dcd')
                    os.replace('combinedTraj.dcd', structureName + '_trajectory-1.dcd')
                else:
                    os.replace(structureName + '_trajectory.dcd', structureName + '_trajectory-1.dcd')  # in case we need to combine two trajectories

                combinedSlope = checkTrajPCASlope(structurePDB, structureName + '_trajectory-1.dcd', self.params['print step'])
                # TODO what is the slope and what it for?
                if binding is True:
                    self.targetUnbound = checkMidTrajectoryBinding(structurePDB, structureName + '_trajectory-1.dcd', self.targetSeq, self.aptamerSeq, self.params, cutoffTime=1)
                    if self.targetUnbound:
                        printRecord('Target came unbound!')

                if (combinedSlope < cutoff) or (self.targetUnbound is True):
                    converged = True

            os.replace(structureName + '_trajectory-1.dcd', structureName + '_complete_trajectory.dcd')
            printRecord('Generated: ' + structureName + '_complete_trajectory.dcd (autosampling was ON)')

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
        if self.actionDict['do 2d analysis'] is True:
            printRecord('Predicted 2D structure :' + self.ssAnalysis['displayed 2d string'][self.i])  # if we did not fold the structure, we probably do not know the secondary structure.
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

    def analyzeBinding(self, bindStructure, bindTrajectory, freeStructure, freeTrajectory):
        """
        Analyze trajectory of both free aptamer and aptamer-target complex to find out about binding information
        :param bindStructure: complex structure's PDB
        :param bindTrajectory: complex structure's MD trajectory
        :param freeStructure: free aptamer structure's PDB
        :param freeTrajectory: free aptamer structure's MD trajectory
        :return:
        """
        bindu = mda.Universe(bindStructure, bindTrajectory)
        freeu = mda.Universe(freeStructure, freeTrajectory)
        bindingDict = bindingAnalysis(bindu, freeu, self.targetSeq, self.aptamerSeq)  # look for contacts between target ligand and aptamer.
        if self.targetUnbound:
            bindingDict['target came unbound'] = True
        else:
            bindingDict['target came Unbound'] = False

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
        buffer = 0.9  # ideally would like a buffer of 20% of the walltime
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
                printRecord('Binding sampling time chunk reduced to {} ns'.format(self.params['sampling time']))

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
