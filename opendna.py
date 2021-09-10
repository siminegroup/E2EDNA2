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
        self.pdbDict = {}
        self.dcdDict = {}

        self.getActionDict()
        if self.actionDict['make workdir']:
            self.setup()  # if we dont need a workdir & MMB files, don't make one


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
            copyfile(self.params['mmb normal template'], self.workDir + '/commands.template.dat')
            copyfile(self.params['mmb quick template'], self.workDir + '/commands.template_quick.dat')
            copyfile(self.params['mmb long template'], self.workDir + '/commands.template_long.dat')

            # copy lightdock scripts
            copytree('lib/lightdock', self.workDir + '/ld_scripts')

        else:
            self.workDir = self.params['workdir'] + '/' + 'run%d' % self.params['run num']

        # move to working dr
        os.chdir(self.workDir)
        printRecord('Simulating {} with {}'.format(self.sequence, self.peptide))

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

    def computeRuntime(self):
        '''
        estimate the maximum runtime for this code, given a ns/day calculation from existing work done
        omit MMB and docking as being relatively cheap
        '''
        dt = self.params['sampling time']
        ssStructures = self.params['N 2D structures']
        dockingStructures = self.params['N docked structures']
        aptamerSteps = self.params['max aptamer sampling iterations']
        complexSteps = self.params['max complex sampling iterations']
        maxNS = ssStructures * aptamerSteps * dt + ssStructures * dockingStructures * complexSteps * dt # maximum possible sampling time in nanoseconds
        projectedMaxRuntime = maxNS / self.ns_per_day * 24 # number of hours to complete our maximum sampling time

        return projectedMaxRuntime

    def checkRuntime(self):
        '''
        estimate the runtime at the current conditions, and adjust parameters if necessary to fit a given runtime
        '''
        projectedMaxRuntime = self.computeRuntime()
        buffer = 0.8 # we ideally would like a buffer of 20% of the walltime
        doReduction = 0
        while projectedMaxRuntime > (buffer * self.params['max walltime']): # if the code is going to take too long, shorten the sampling time until it fits
            doReduction = 1
            self.params['sampling time'] *= 0.95
            projectedMaxRuntime = self.computeRuntime()

        if doReduction == 1:
            if self.params['sampling time'] < 0.1:
                printRecord('Sampling time reduced to 100 ps - this run is too expensive, and will now terminate')
                self.terminateRun()
            else:
                printRecord('Sampling time chunk reduced to {} ns'.format(self.params['sampling time']))

    def terminateRun(self):
        '''
        for some reason, the run will suddenly end
        '''
        sys.exit()


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
            printRecord('2D structure #{} is {}'.format(self.i, self.ssAnalysis['2d string'][self.i]))
            self.pairList = np.asarray(self.pairLists[self.i])

            if self.actionDict['do MMB']:  # fold it
                self.foldSequence(self.sequence, self.pairList)

            if self.actionDict['do smoothing']:
                self.MDSmoothing(self.pdbDict['mmb folded sequence {}'.format(self.i)], relaxationTime=self.params['smoothing time'])  # relax for xx nanoseconds

            if self.actionDict['get equil repStructure']:
                outputDict['free aptamer results {}'.format(self.i)] = self.runFreeAptamer(self.pdbDict['relaxed sequence {}'.format(self.i)])
                np.save('opendnaOutput', outputDict)  # save outputs

            if self.actionDict['do docking'] and (self.peptide is not False):  # find docking configurations for the complexed structure
                outputDict['dock scores {}'.format(self.i)] = self.dock(self.pdbDict['representative aptamer {}'.format(self.i)], 'peptide.pdb')
                np.save('opendnaOutput', outputDict)  # save outputs

            printRecord('Running over %d' % self.params['N docked structures'] + ' docked structures')
            for self.j in range(self.params['N docked structures']): # loop over docking configurations for a given secondary structure
                printRecord('Docked structure #{}'.format(self.j))
                if self.actionDict['do binding']:  # run MD on the complexed structure
                    outputDict['binding results {} {}'.format(self.i, self.j)] = self.bindingDynamics(self.pdbDict['binding complex {} {}'.format(self.i,int(self.j))])

                    np.save('opendnaOutput', outputDict)  # save outputs

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
            nup = interfaces.nupack(sequence, self.params['temperature'], self.params['ionic strength'], self.params['[Mg]'])  # initialize nupack
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
        mmb = interfaces.mmb(sequence, pairList, self.params, self.i)
        foldFidelity = mmb.run()
        printRecord('Initial fold fidelity = %.3f'%foldFidelity)
        attempts = 1
        if self.params['fold speed'] != 'quick': # don't rerun if we're quick folding
            intervalLength = 5  # initial interval length
            while (foldFidelity <= 0.9) and (attempts < 5): # if it didn't fold properly, try again with a longer annealing time - up to XX times
                self.params['fold speed'] = 'long'
                printRecord("Refolding Sequence")
                mmb = interfaces.mmb(sequence, pairList, self.params, self.i, intervalLength)
                foldFidelity = mmb.run()
                printRecord('Subsequent fold fidelity = %.3f' % foldFidelity)
                attempts += 1
                intervalLength += 5 # on repeated attempts, try harder

        self.pdbDict['mmb folded sequence {}'.format(self.i)] = mmb.foldedSequence
        os.system('mv commands.run* ' + mmb.fileDump)
        printRecord("Folded Sequence")

    def MDSmoothing(self, structure, relaxationTime=0.01):
        """
        do a short MD run in water to relax the coarse MMB structure
        relaxation time in nanoseconds, print time in picoseconds
        """
        structureName = structure.split('.')[0]
        printRecord('Running relaxation')
        prepPDB(structure, self.params['box offset'], self.params['pH'], self.params['ionic strength'], MMBCORRECTION=True, waterBox=True)
        processedStructure = structureName + '_processed.pdb'
        processedStructureTrajectory = processedStructure.split('.')[0] + '_trajectory.dcd'
        omm = interfaces.omm(processedStructure, self.params)
        self.ns_per_day = omm.doMD()

        printRecord('Pre-relaxation simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print speed

        cleanTrajectory(processedStructure, processedStructureTrajectory) # clean up trajectory for later use
        self.dcdDict['relaxed sequence {}'.format(self.i)] = 'clean_' + processedStructureTrajectory
        self.pdbDict['relaxed sequence {}'.format(self.i)] = 'relaxedSequence_{}.pdb'.format(self.i)
        extractFrame(processedStructure,processedStructureTrajectory,-1,self.pdbDict['relaxed sequence {}'.format(self.i)]) # final frame


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
        processedAptamerTrajectory = processedAptamer.split('.')[0] + '_complete_trajectory.dcd' # this is what autoMD will output

        self.autoMD(processedAptamer)  # do MD - automatically converge to equilibrium sampling

        printRecord('Free aptamer simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print speed
        self.checkRuntime()

        cleanTrajectory(processedAptamer, processedAptamerTrajectory)
        self.dcdDict['sampled aptamer {}'.format(self.i)] = 'clean_' + processedAptamerTrajectory
        self.pdbDict['sampled aptamer {}'.format(self.i)] = 'clean_' + processedAptamer
        aptamerDict = self.analyzeTrajectory(self.pdbDict['sampled aptamer {}'.format(self.i)],self.dcdDict['sampled aptamer {}'.format(self.i)])
        self.pdbDict['representative aptamer {}'.format(self.i)] = 'repStructure_{}.pdb'.format(self.i)

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
            copyfile('top_%d'%self.i + '/top_%d' % int(i + 1) + '.pdb', 'complex_{}_{}.pdb'.format(self.i, int(i)))
            self.pdbDict['binding complex {} {}'.format(self.i, int(i))] = 'complex_{}_{}.pdb'.format(self.i,int(i))

        printRecord('Docking complete, best scores are {}'.format(topScores))
        return topScores

    def bindingDynamics(self, complex):
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
        processedComplexTrajectory = processedComplex.split('.')[0] + '_complete_trajectory.dcd' # this is what autoMD will output

        self.autoMD(processedComplex, binding=True)  # do MD - automatically converge to equilibrium sampling

        printRecord('Complex simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print speed
        self.checkRuntime()

        cleanTrajectory(processedComplex, processedComplexTrajectory)

        self.dcdDict['sampled complex {} {}'.format(self.i, self.j)] = 'clean_' + processedComplexTrajectory
        self.pdbDict['sampled complex {} {}'.format(self.i, self.j)] = 'clean_' + processedComplex

        bindingDict = self.analyzeBinding(self.pdbDict['sampled complex {} {}'.format(self.i, self.j)], self.dcdDict['sampled complex {} {}'.format(self.i,self.j)], self.pdbDict['sampled aptamer {}'.format(self.i)], self.dcdDict['sampled aptamer {}'.format(self.i)])

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
        if binding:
            maxIter = self.params['max complex sampling iterations']
        else:
            maxIter = self.params['max aptamer sampling iterations']  # some maximum number of allowable iterations
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

                combinedSlope = checkTrajPCASlope(structure, structureName + '_trajectory-1.dcd', self.params['print step'])
                if binding:  # if this is a binding simulation, also check to see if the analyte has come unbound from the analyte, and if so, cutoff the simulation
                    self.analyteUnbound = checkMidTrajectoryBinding(structure, structureName + '_trajectory-1.dcd', self.peptide, self.sequence, self.params, cutoffTime=1)
                    if self.analyteUnbound:
                        printRecord('Analyte came unbound!')

                if (combinedSlope < cutoff) or (self.analyteUnbound == True):  # the average magnitude of sloped should be below some cutoff`
                    converged = True

            os.replace(structureName + '_trajectory-1.dcd', structureName + '_complete_trajectory.dcd')  # we'll consider the full trajectory as 'sampling'

    def analyzeTrajectory(self, structure, trajectory):
        """
        analyze trajectory for aptamer fold + analyte binding information
        :param trajectory:
        :return:
        """
        u = mda.Universe(structure, trajectory)

        wcTraj = getWCDistTraj(u)  # watson-crick base pairing distances (H-bonding) s
        baseDistTraj = getBaseBaseDistTraj(u)  # FAST, base-base center-of-geometry distances
        nucleicAnglesTraj = getNucDATraj(u)  # FAST, new, omits 'chi' angle between ribose and base

        # 2D structure analysis
        pairTraj = getPairTraj(wcTraj)
        secondaryStructure = analyzeSecondaryStructure(pairTraj)  # find equilibrium secondary structure
        print('Predicted 2D structure :' + self.ssAnalysis['2d string'][self.i])
        print('Actual 2D structure    :' + configToString(secondaryStructure))

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
