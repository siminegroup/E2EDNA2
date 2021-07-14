# import statements
import re
import glob
import os
import numpy as np
from shutil import copyfile, copytree
from simtk.openmm import *
import interfaces
from utils import *
from analysisTools import *


# noinspection PyPep8Naming
class opendna():
    def __init__(self, params):
        self.params = params
        self.sequence = self.params['sequence']
        self.peptide = self.params['peptide']

        self.getActionDict()
        if self.actionDict['make workdir']:
            self.setup()  # if we dont need a workdir & MMB files, don't make one
        self.getCheckpoint()

    def getActionDict(self):
        """
        generate a binary sequence of 'to-do's' given user input
        """
        self.actionDict = {}
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
        elif self.params['mode'] == 'coarse dock':
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = False
            self.actionDict['get equil repStructure'] = False
            self.actionDict['do docking'] = True
            self.actionDict['do binding'] = False
        elif self.params['mode'] == 'smooth dock':
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
            self.actionDict['do smoothing'] = False
            self.actionDict['get equil repStructure'] = True
            self.actionDict['do docking'] = False
            self.actionDict['do binding'] = False
        elif self.params['mode'] == 'full docking':
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = False
            self.actionDict['get equil repStructure'] = True
            self.actionDict['do docking'] = True
            self.actionDict['do binding'] = False
        elif self.params['mode'] == 'full binding':
            self.actionDict['make workdir'] = True
            self.actionDict['do 2d analysis'] = True
            self.actionDict['do MMB'] = True
            self.actionDict['do smoothing'] = False
            self.actionDict['get equil repStructure'] = True
            self.actionDict['do docking'] = True
            self.actionDict['do binding'] = True

    def setup(self):
        """
        setup working directory
        copy in relevant xyz and keyfiles
        move to relevant directory
        :return:
        """

        if (self.params['explicit run enumeration'] == True) or (self.params['run num'] == 0):
            if self.params['run num'] == 0:
                self.makeNewWorkingDirectory()
            else:
                self.workDir = self.params['workdir'] + '/run%d' % self.params['run num']  # make a new workdir with the given index
                os.mkdir(self.workDir)

            os.mkdir(self.workDir + '/outfiles')
            # copy structure files
            copyfile(self.params['analyte pdb'], self.workDir + '/analyte.pdb')

            # copy MMB params and command script
            copyfile(self.params['mmb params'], self.workDir + '/parameters.csv')
            copyfile(self.params['mmb template'], self.workDir + '/commands.template.dat')

            # copy lightdock scripts
            copytree('lib/lightdock', self.workDir + '/ld_scripts')

        else:
            self.workDir = self.params['workdir'] + '/' + 'run%d' % self.params['run num']

        # move to working dr
        os.chdir(self.workDir)

        if self.params['device'] == 'local':
            os.environ["LD_LIBRARY_PATH"] = self.params['mmb dir']  # export the path to the MMB library - only necessary in WSL environments

    def makeNewWorkingDirectory(self):  # make working directory
        """
        make a new working directory
        non-overlapping previous entries
        :return:
        """
        workdirs = glob.glob(self.params['workdir'] + '/' + 'run*')  # check for prior working directories
        if len(workdirs) > 0:
            prev_runs = []
            for i in range(len(workdirs)):
                prev_runs.append(int(workdirs[i].split('run')[-1]))

            prev_max = max(prev_runs)
            self.workDir = self.params['workdir'] + '/' + 'run%d' % (prev_max + 1)
            os.mkdir(self.workDir)
            print('Starting Fresh Run %d' % (prev_max + 1))
        else:
            self.workDir = self.params['workdir'] + '/' + 'run1'
            print('Starting Fresh Run 1')
            os.mkdir(self.workDir)

    def getCheckpoint(self):
        """
        identify if any work has been done in this directory, and if so, where to pick up
        :return:
        """
        try:
            f = open('checkpoint.txt', 'r')  # if it does, see how long it is
            text = f.read()
            f.close()
            self.checkpoints = len([m.start() for m in re.finditer('\n', text)])

            printRecord("Resuming Run #%d" % int(self.params['run num']))
        except:
            f = open('checkpoint.txt', 'w')
            f.write('Started!\n')
            f.close()
            self.checkpoints = 1
            self.checkpoint = 'initialized'

    def run(self):
        """
        run the binding simulation end-to-end
        consult checkpoints to not repeat prior steps
        :return:
        """
        outputDict = {}
        outputDict['params'] = self.params
        np.save('opendnaOutput', outputDict)
        if self.actionDict['do 2d analysis']:  # get secondary structure
            self.pairLists = self.getSecondaryStructure(self.sequence)
            outputDict['2d analysis'] = self.ssAnalysis
            np.save('opendnaOutput', outputDict)  # save outputs

        printRecord('Running over %d' % len(self.pairLists) + ' possible 2D structures')

        for self.i in range(len(self.pairLists)):  # loop over possible secondary structures
            self.pairList = np.asarray(self.pairLists[self.i])

            if self.actionDict['do MMB']:  # fold it
                self.foldSequence(self.sequence, self.pairList)
                os.rename('sequence.pdb','sequence_%d'%self.i + '.pdb')

            if self.actionDict['do smoothing']:
                self.MDSmoothing('sequence_%d'%self.i + '.pdb', relaxationTime=0.01)  # relax for xx nanoseconds

            if self.actionDict['get equil repStructure']:
                outputDict['free aptamer results'] = self.runFreeAptamer('sequence_%d'%self.i + '.pdb')
                np.save('opendnaOutput_%d' % self.i, outputDict)  # save outputs

            if self.actionDict['do docking'] and (self.peptide is not False):  # find docking configurations for the complexed structure
                outputDict['dock scores'] = self.dock('repStructure_%d'%self.i + '.pdb', 'peptide.pdb')
                np.save('opendnaOutput_%d' % self.i, outputDict)  # save outputs

            printRecord('Running over %d' % self.params['N docked structures'] + ' docked structures')
            for self.j in range(self.params['N docked structures']):
                if self.actionDict['do binding'] and (self.peptide is not False):  # run MD on the complexed structure
                    outputDict['binding results'] = self.bindingDynamics('complex_{}_{}'.format(self.i,int(self.j + 1)) + '.pdb', 'clean_finished_sequence_%d'%self.i + '.pdb')

                    np.save('opendnaOutput_%d' % self.i + '_%d' % int(self.j + 1), outputDict)  # save outputs

        return outputDict

    # =======================================================
    # =======================================================
    # =======================================================

    def getSecondaryStructure(self, sequence):
        """
        get the secondary structure for a given sequence
        using seqfold here - identical features are available using nupack, though results are sometimes different
        :param sequence:
        :return: a dot-bracket string and list of paired bases (assuming single-strand DNA aptamer)
        """
        printRecord("Get Secondary Structure")
        if self.params['secondary structure engine'] == 'seqfold':
            ssString, pairList = getSeqfoldStructure(sequence, self.params['temperature'])  # seqfold guess
            self.ssAnalysis = [ssString, pairList]
            return pairList
        elif self.params['secondary structure engine'] == 'NUPACK':
            nup = interfaces.nupack(sequence, self.params['temperature'], self.params['ionic strength'])  # initialize nupack
            self.ssAnalysis = nup.run()  # run nupack analysis of possible secondary structures

            distances = getSecondaryStructureDistance(self.ssAnalysis['config'])
            if len(distances) > 1:
                topReps, topProbs, topDists = do2DAgglomerativeClustering(self.ssAnalysis['config'], self.ssAnalysis['state prob'], distances)
                topPairLists = []
                nLists = min(len(topReps),self.params['N 2D structures']) # there may be less clusters than we ask for
                for i in range(nLists):
                    if (i == 0) or (np.amin(topDists[i, :i]) > 0.05):  # only add a config to the list if it's at least 5% different from a previously accepted config
                        topPairLists.append(ssToList(configToString(topReps[i].astype(int))))
                return topPairLists
            else:  # if for some reason there's only a single plausible structure (unlikely)
                return self.ssAnalysis['pair list']

    def foldSequence(self, sequence, pairList):
        """
        generate a coarsely folded structure for the aptamer
        :param sequence: the ssDNA sequence to be folded
        :param pairList: list of binding base pairs
        :return:
        """
        # write pair list as forces to the MMB command file
        printRecord("Folding Sequence")
        mmb = interfaces.mmb(sequence, pairList, self.params)
        mmb.run()
        if (self.params['mode'] == '3d coarse') or (self.params['mode'] == 'coarse dock'):
            self.prepCoarsePDB()

        printRecord("Folded Sequence")

    def prepCoarsePDB(self):
        prepPDB('sequence.pdb', self.params['box offset'], self.params['pH'], self.params['ionic strength'], MMBCORRECTION=True, waterBox=False)
        copyfile('sequence_processed.pdb', 'repStructure_%d'%self.i + '.pdb')  # final structure output

    def MDSmoothing(self, structure, relaxationTime=0.01):
        """
        do a short MD run in water to relax the coarse MMB structure
        relaxation time in nanoseconds, print time in picoseconds
        """
        structureName = structure.split('.')[0]
        printRecord('Running quick relaxation')
        prepPDB(structure, self.params['box offset'], self.params['pH'], self.params['ionic strength'], MMBCORRECTION=True, waterBox=True)
        processedStructure = structureName + '_processed.pdb'
        omm = interfaces.omm(structure, self.params)
        self.ns_per_day = omm.doMD()
        extractFrame(processedStructure, processedStructure.split('.')[0] + '_trajectory.dcd', -1, 'repStructure_%d'%self.i + '.pdb')  # pull the last frame of the relaxation

    def runFreeAptamer(self, aptamer):
        """
        run molecular dynamics for free aptamer only
        do relevant analysis
        :param aptamer:
        :return:
        """
        structureName = aptamer.split('.')[0]
        printRecord('Running free aptamer dynamics')
        prepPDB(aptamer, self.params['box offset'], self.params['pH'], self.params['ionic strength'], MMBCORRECTION=True, waterBox=True)  # add periodic box and appropriate protons
        processedAptamer = structureName + '_processed.pdb'
        self.autoMD(processedAptamer)  # do MD - automatically converge to equilibrium sampling

        printRecord('Free aptamer simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print speed

        finishedPDB = 'finished_' + aptamer
        finishedDCD = 'finished_' + structureName + '.dcd'
        os.replace(structureName + '_processed_complete_trajectory.dcd', finishedDCD)
        os.replace(processedAptamer, finishedPDB)

        cleanTrajectory(finishedPDB, finishedDCD)
        aptamerDict = self.analyzeTrajectory('clean_' + finishedPDB, 'clean_' + finishedDCD, False)

        printRecord('Free aptamer sampling complete')

        return aptamerDict

    def dock(self, aptamer, peptide):
        """
        use LightDock to run docking and isolate good structures
        """
        printRecord('Docking')
        buildPeptide(self.peptide)
        ld = interfaces.ld(aptamer, peptide, self.params, self.i)
        ld.run()
        topScores = ld.topScores

        for i in range(self.params['N docked structures']):
            copyfile('top_%d'%self.i + '/top_%d' % int(i + 1) + '.pdb', 'complex_{}_{}'.format(self.i, int(i + 1)) + '.pdb')

        printRecord('Docking complete')

        return topScores

    def bindingDynamics(self, complex, freeSequence):
        """
        run molecular dynamics
        do relevant analysis
        :param complex:
        :return:
        """
        structureName = complex.split('.')[0]
        printRecord('Running Binding Simulation')
        prepPDB(complex, self.params['box offset'], self.params['pH'], self.params['ionic strength'], MMBCORRECTION=False, waterBox=True)
        processedComplex = complex.split('.')[0] + '_processed.pdb'

        self.autoMD(processedComplex, binding=True)  # do MD - automatically converge to equilibrium sampling

        printRecord('Complex simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print speed

        finishedPDB = 'finished_' + complex
        finishedDCD = 'finished_' + structureName + '.dcd'
        os.replace(structureName + '_processed_complete_trajectory.dcd', finishedDCD)
        os.replace(processedComplex, finishedPDB)

        cleanTrajectory(finishedPDB, finishedDCD)
        bindingDict = self.analyzeBinding('clean_' + finishedPDB, 'clean_' + finishedDCD, freeSequence, freeSequence.split('.')[0] + '.dcd')

        # print findings
        printRecord('Binding Results: Contact Persistence = {:.2f}, Contact Score = {:.2f}, Conformation Change = {:.2f}'.format(bindingDict['close contact ratio'], bindingDict['contact score'], bindingDict['conformation change']))
        printRecord('Complex sampling complete')

        return bindingDict

    def autoMD(self, structure, binding=False):
        """
        run MD until either for a set amount of time or until the dynamics 'converge'
        optionally, sample after we reach convergence ("equilibration") # not implemented yet
        :param structure:
        :return:
        """
        maxIter = self.params['max autoMD iterations']  # some maximum number of allowable iterations
        cutoff = self.params['autoMD convergence cutoff']

        structureName = structure.split('.')[0]
        if self.params['auto sampling'] == False:  # just run MD once
            self.analyteUnbound = False
            omm = interfaces.omm(structure, self.params)
            self.ns_per_day = omm.doMD()
            os.replace(structureName + '_trajectory.dcd', structureName + "_complete_trajectory.dcd")

        elif self.params['auto sampling'] == True:  # run until we detect equilibration
            converged = False
            iter = 0
            self.analyteUnbound = False

            while (converged == False) and (iter < maxIter):
                iter += 1
                omm = interfaces.omm(structure, self.params)
                self.ns_per_day = omm.doMD()

                if iter > 1:  # if we have multiple trajectory segments, combine them
                    appendTrajectory(structure, structureName + '_trajectory-1.dcd', structureName + '_trajectory.dcd')
                    os.replace('combinedTraj.dcd', structureName + '_trajectory-1.dcd')
                else:
                    os.replace(structureName + '_trajectory.dcd', structureName + '_trajectory-1.dcd')

                combinedSlope = checkTrajPCASlope(structure, structureName + '_trajectory-1.dcd')
                if binding:  # if this is a binding simulation, also check to see if the analyte has come unbound from the analyte, and if so, cutoff the simulation
                    self.analyteUnbound = checkMidTrajectoryBinding(structure, structureName + '_trajectory-1.dcd', self.peptide, self.sequence, self.params, cutoffTime=1)
                    if self.analyteUnbound:
                        printRecord('Analyte came unbound!')

                if (combinedSlope < cutoff) or (self.analyteUnbound == True):  # the average magnitude of sloped should be below some cutoff`
                    converged = True

            os.replace(structureName + '_trajectory-1.dcd', structureName + '_complete_trajectory.dcd')  # we'll consider the full trajectory as 'sampling'

    def analyzeTrajectory(self, structure, trajectory, analyte):
        """
        analyze trajectory for aptamer fold + analyte binding information
        :param trajectory:
        :return:
        """
        u = mda.Universe(structure, trajectory)

        wcTraj = getWCDistTraj(u)  # watson-crick base pairing distances (H-bonding) FAST but some errors
        baseDistTraj = getBaseBaseDistTraj(u)  # FAST, base-base center-of-geometry distances
        nucleicAnglesTraj = getNucDATraj(u)  # FAST, new, omits 'chi' angle between ribose and base

        # 2D structure analysis
        pairTraj = getPairTraj(wcTraj)
        secondaryStructure = analyzeSecondaryStructure(pairTraj)  # find equilibrium secondary structure

        # 3D structure analysis
        mixedTrajectory = np.concatenate(
            (baseDistTraj.reshape(len(baseDistTraj), int(baseDistTraj.shape[-2] * baseDistTraj.shape[-1])), nucleicAnglesTraj.reshape(len(nucleicAnglesTraj), int(nucleicAnglesTraj.shape[-2] * nucleicAnglesTraj.shape[-1]))),
            axis=1)  # mix up all our info
        representativeIndex, pcTrajectory, eigenvalues = isolateRepresentativeStructure(mixedTrajectory)

        # save this structure as a separate file
        extractFrame(structure, trajectory, representativeIndex, 'repStructure_%d'%self.i + '.pdb')

        # we also need a function which processes the energies spat out by the trajectory logs

        analysisDict = {}  # compile our results
        analysisDict['2D structure'] = secondaryStructure
        analysisDict['RC trajectories'] = pcTrajectory
        analysisDict['representative structure index'] = representativeIndex

        return analysisDict

    def analyzeBinding(self, bindStructure, bindTrajectory, freeStructure, freeTrajectory):
        """
        analyze trajectory for aptamer fold + analyte binding information
        :param trajectory:
        :return:
        """
        bindu = mda.Universe(bindStructure, bindTrajectory)
        freeu = mda.Universe(freeStructure, freeTrajectory)
        bindingDict = bindingAnalysis(bindu, freeu, self.peptide, self.sequence)  # look for contacts between analyte and aptamer
        if self.analyteUnbound:
            bindingDict['analyte came unbound'] = True
        else:
            bindingDict['analyte came unbound'] = False

        return bindingDict