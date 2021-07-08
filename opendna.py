# import statements
import re
from utils import *
from simtk.openmm import *
import simtk.unit as unit
from simtk.openmm.app import *
from pdbfixersource import PDBFixer


# noinspection PyPep8Naming
class opendna():
    def __init__(self, sequence, peptide, params):
        self.params = params
        self.sequence = sequence
        self.peptide = peptide
        self.params['sequence'] = sequence
        self.params['peptide'] = peptide

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
            copytree('lib/lightdock', self.workDir + '/lightdock')

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

            print("Resuming Run #%d" % int(self.params['run num']))
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
            self.ssString, self.pairList = self.getSecondaryStructure(self.sequence)
            outputDict['2d string'] = self.ssString
            outputDict['pair list'] = self.pairList
            np.save('opendnaOutput', outputDict)  # save outputs

        if self.actionDict['do MMB']:  # fold it
            self.foldSequence(self.sequence, self.pairList)

        if self.actionDict['do smoothing']:
            self.MDSmoothing('sequence.pdb', relaxationTime=0.01)  # relax for xx nanoseconds

        if self.actionDict['get equil repStructure']:
            outputDict['free aptamer results'] = self.freeAptamerDynamics('sequence.pdb')
            np.save('opendnaOutput', outputDict)  # save outputs

        if self.actionDict['do docking'] and (self.peptide is not False):  # find docking configurations for the complexed structure
            outputDict['dock scores'] = self.runDocking('repStructure.pdb', 'peptide.pdb')
            np.save('opendnaOutput', outputDict)  # save outputs

        if self.actionDict['do binding'] and (self.peptide is not False):  # run MD on the complexed structure
            outputDict['binding results'] = self.bindingDynamics('complex_0.pdb')
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
        print("Get Secondary Structure")
        if os.path.exists("pre_fold.pdb"):  # if we have a pre-folded structure, do nothing
            pass
        else:
            if self.params['secondary structure engine'] == 'seqfold':
                ssString, pairList = getSeqfoldStructure(sequence, self.params['temperature'])  # seqfold guess
            elif self.params['secondary structure engine'] == 'NUPACK':
                ssDict, subopts = doNupackAnalysis(sequence, self.params['temperature'], self.params['ionic strength'])  # NUPACK guess
                ssString = ssDict['2d string']
                pairList = ssDict['pair list']
                self.ssInfo = [ssDict, subopts]  # more comprehensive analysis

            return ssString, np.asarray(pairList)

    def foldSequence(self, sequence, pairList):
        """
        generate a coarsely folded structure for the aptamer
        :param sequence: the ssDNA sequence to be folded
        :param pairList: list of binding base pairs
        :return:
        """
        # write pair list as forces to the MMB command file
        print("Folding Sequence")
        comFile = 'commands.fold.dat'  # name of command file
        copyfile('commands.template.dat', comFile)  # make command file
        replaceText(comFile, 'SEQUENCE', sequence)
        replaceText(comFile, 'TEMPERATURE', str(self.params['temperature'] - 273)) # probably not important, but we can add the temperature in C

        baseString = '#baseInteraction A IND WatsonCrick A IND2 WatsonCrick Cis'
        lineNum = findLine(comFile, baseString)  # find the line number to start enumerating base pairs

        for i in range(len(pairList)):
            filledString = 'baseInteraction A {} WatsonCrick A {} WatsonCrick Cis'.format(pairList[i, 0], pairList[i, 1])
            addLine(comFile, filledString, lineNum + 1)

        # run fold - sometimes this errors for no known reason - keep trying till it works
        Result = None
        attempts = 0
        while (Result is None) and (attempts < 100):
            try:
                attempts += 1
                os.system(self.params['mmb'] + ' -c ' + comFile + ' > outfiles/fold.out')
                os.replace('frame.pdb', 'sequence.pdb')
                if (self.params['mode'] == '3d coarse') or (self.params['mode'] == 'coarse dock'):
                    self.prepPDB('sequence.pdb', MMBCORRECTION=True, waterBox=False)
                    copyfile('sequence_processed.pdb', 'repStructure.pdb')  # final structure output
                Result = 1
                # cleanup
                os.mkdir('mmbFiles')
                os.system('mv last* mmbFiles')
                os.system('mv trajectory.* mmbFiles')
                os.system('mv match.* mmbFiles')
                os.system('mv commands* mmbFiles')

            except:
                pass

        writeCheckpoint("Folded Sequence")

    def prepPDB(self, file, MMBCORRECTION=False, waterBox=True):
        """
        soak pdb file in water box
        :MMBCORRECTION: if the input pdb file is an MMB output, we need to apply a correction, since MMB is a little weird formatting-wise https://simtk.org/plugins/phpBB/viewtopicPhpbb.php?f=359&t=13397&p=0&start=0&view=&sid=bc6c1b9005122914ec7d572999ba945b
        :param file:
        :return:
        """
        if MMBCORRECTION:
            replaceText(file, '*', "'")  # due to a bug in this version of MMB - structures are encoded improperly - this fixes it

        fixer = PDBFixer(filename=file)
        padding, boxSize, boxVectors = None, None, None
        geompadding = float(self.params['box offset']) * unit.nanometer
        padding = geompadding

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()  # may need to optimize bonding here

        fixer.addMissingHydrogens(pH=self.params['pH'])  # add missing hydrogens

        if waterBox == True:
            ionicStrength = float(self.params['ionic strength']) * unit.molar  # not implemented
            positiveIon = 'Na+'  # params['positiveion']+'+'
            negativeIon = 'Cl-'  # params['negativeion']+'-'
            fixer.addSolvent(boxSize, padding, boxVectors, positiveIon, negativeIon, ionicStrength)

        PDBFile.writeFile(fixer.topology, fixer.positions, open(file.split('.pdb')[0] + '_processed.pdb', 'w'))

    def MDSmoothing(self, structure, relaxationTime=0.01):
        """
        do a short MD run in water to relax the coarse MMB structure
        relaxation time in nanoseconds, print time in picoseconds
        """
        structureName = structure.split('.')[0]
        print('Running quick relaxation')
        self.prepPDB(structure, MMBCORRECTION=True, waterBox=True)
        processedStructure = structureName + '_processed.pdb'
        self.openmmDynamics(processedStructure, simTime=relaxationTime)
        extractFrame(processedStructure, processedStructure.split('.')[0] + '_trajectory.dcd', -1, 'repStructure.pdb')  # pull the last frame of the relaxation

    def freeAptamerDynamics(self, aptamer):
        """
        run molecular dynamics
        do relevant analysis
        :param aptamer:
        :return:
        """
        structureName = aptamer.split('.')[0]
        print('Running free aptamer dynamics')
        self.prepPDB(aptamer, MMBCORRECTION=True, waterBox=True)  # add periodic box and appropriate protons
        processedAptamer = structureName + '_processed.pdb'
        self.autoMD(processedAptamer)  # do MD - automatically converge to equilibrium sampling

        print('Free aptamer simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print speed

        os.replace(structureName + '_processed_complete_trajectory.dcd', 'aptamer.dcd')  # rename
        os.replace(processedAptamer, 'aptamer.pdb')

        cleanTrajectory('aptamer.pdb', 'aptamer.dcd')  # make a trajectory without waters and ions and stuff
        aptamerDict = self.analyzeTrajectory('cleanaptamer.pdb', 'cleanaptamer.dcd', False)

        writeCheckpoint('Free aptamer sampling complete')

        return aptamerDict

    def runDocking(self, aptamer, peptide):
        """
        use LightDock to run docking and isolate good structures
        """
        print('Docking')
        buildPeptide(self.peptide)

        topScores = lightDock(aptamer, peptide, self.params)  # get optimal docked structure
        copyfile('top/top_1.pdb', 'complex_0.pdb')

        writeCheckpoint('Docking complete')

        return topScores

    def bindingDynamics(self, complex):
        """
        run molecular dynamics
        do relevant analysis
        :param complex:
        :return:
        """
        structureName = complex.split('.')[0]
        print('Running Binding Simulation')
        self.prepPDB(complex, MMBCORRECTION=False, waterBox=True)
        processedComplex = complex.split('.')[0] + '_processed.pdb'

        self.autoMD(processedComplex, binding=True)  # do MD - automatically converge to equilibrium sampling

        print('Complex simulation speed %.1f' % self.ns_per_day + ' ns/day')  # print speed

        os.replace(structureName + '_processed_complete_trajectory.dcd', 'complex.dcd')
        os.replace(processedComplex, 'complex.pdb')

        cleanTrajectory('complex.pdb', 'complex.dcd')
        bindingDict = self.analyzeTrajectory('cleancomplex.pdb', 'cleancomplex.dcd', True)

        # print findings
        print('Binding Results: Contact Persistence = {:.2f}, Contact Score = {:.2f}, Conformation Change = {:.2f}'.format(bindingDict['close contact ratio'], bindingDict['contact score'], bindingDict['conformation change']))
        writeCheckpoint('Complex sampling complete')

        return bindingDict


    def openmmDynamics(self, structure, simTime=False):
        """
        run OpenMM dynamics
        minimize, equilibrate, and sample
        :param structure:
        :return:
        """

        pdb = PDBFile(structure)
        structureName = structure.split('.')[0]
        waterModel = self.params['water model']
        forcefield = ForceField('amber14-all.xml', 'amber14/' + waterModel + '.xml')

        # System Configuration
        nonbondedMethod = self.params['nonbonded method']
        nonbondedCutoff = self.params['nonbonded cutoff'] * unit.nanometer
        ewaldErrorTolerance = self.params['ewald error tolerance']
        constraints = self.params['constraints']
        rigidWater = self.params['rigid water']
        constraintTolerance = self.params['constraint tolerance']
        hydrogenMass = self.params['hydrogen mass'] * unit.amu

        # Integration Options
        dt = self.params['time step'] / 1000 * unit.picoseconds
        temperature = self.params['temperature'] * unit.kelvin
        friction = self.params['friction'] / unit.picosecond

        # Simulation Options
        if simTime != False:  # we can override the simulation time here
            steps = int(simTime / 2 * 1e6 // self.params['time step'])  # number of steps
            equilibrationSteps = int(simTime / 2 * 1e6 // self.params['time step'])
        else:
            steps = int(self.params['sampling time'] * 1e6 // self.params['time step'])  # number of steps
            equilibrationSteps = int(self.params['equilibration time'] * 1e6 // self.params['time step'])

        if self.params['platform'] == 'CUDA':  # 'CUDA' or 'cpu'
            platform = Platform.getPlatformByName('CUDA')
            if self.params['platform precision'] == 'single':
                platformProperties = {'Precision': 'single'}
            elif self.params['platform precision'] == 'double':
                platformProperties = {'Precision': 'double'}
        else:
            platform = Platform.getPlatformByName('CPU')

        reportSteps = int(self.params['print step'] * 1000 / self.params['time step'])  # report step in ps, time step in fs
        dcdReporter = DCDReporter(structureName + '_trajectory.dcd', reportSteps)
        # pdbReporter = PDBReporter('trajectory.pdb', reportSteps) # much less efficient trajectory format
        dataReporter = StateDataReporter('log.txt', reportSteps, totalSteps=steps, step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, separator='\t')
        checkpointReporter = CheckpointReporter(structure.split('.')[0] + '_state.chk', 10000)

        # Prepare the Simulation
        #print('Building system...') # waste of space
        topology = pdb.topology
        positions = pdb.positions
        system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff, constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance, hydrogenMass=hydrogenMass)
        integrator = LangevinMiddleIntegrator(temperature, friction, dt)
        integrator.setConstraintTolerance(constraintTolerance)
        if self.params['platform'] == 'CUDA':
            simulation = Simulation(topology, system, integrator, platform, platformProperties)
        elif self.params['platform'] == 'CPU':
            simulation = Simulation(topology, system, integrator, platform)
        simulation.context.setPositions(positions)

        if not os.path.exists(structure.split('.')[0] + '_state.chk'):
            # Minimize and Equilibrate
            print('Performing energy minimization...')
            simulation.minimizeEnergy()
            print('Equilibrating...')
            simulation.context.setVelocitiesToTemperature(temperature)
            simulation.step(equilibrationSteps)
        else:
            simulation.loadCheckpoint(structure.split('.')[0] + '_state.chk')

        # Simulate
        print('Simulating...')
        simulation.reporters.append(dcdReporter)
        # simulation.reporters.append(pdbReporter)
        simulation.reporters.append(dataReporter)
        simulation.reporters.append(checkpointReporter)
        simulation.currentStep = 0
        with Timer() as md_time:
            simulation.step(steps)  # run the dynamics

        simulation.saveCheckpoint(structure.split('.')[0] + '_state.chk')

        self.ns_per_day = (steps * dt) / (md_time.interval * unit.seconds) / (unit.nanoseconds / unit.day)

    def autoMD(self, structure, binding=False):
        """
        run MD until either for a set amount of time or until the dynamics 'converge'
        optionally, sample after we reach convergence ("equilibration")
        :param structure:
        :return:
        """
        maxIter = self.params['max autoMD iterations']  # some maximum number of allowable iterations
        cutoff = self.params['autoMD convergence cutoff']
        structureName = structure.split('.')[0]
        if self.params['auto sampling'] == False:  # just run MD once
            self.openmmDynamics(structure)
            os.replace(structureName + '_trajectory.dcd', structureName + "_complete_trajectory.dcd")
        elif self.params['auto sampling'] == True:  # run until we detect equilibration
            converged = False
            iter = 0
            self.analyteUnbound = False
            while (converged == False) and (iter < maxIter):
                iter += 1
                self.openmmDynamics(structure)
                time.sleep(20)  # let .dcd file settle down

                if iter > 1:  # if we have multiple trajectory segments, combine them
                    appendTrajectory(structure, structureName + '_trajectory-1.dcd', structureName + '_trajectory.dcd')
                    os.replace('combinedTraj.dcd', structureName + '_trajectory-1.dcd')
                else:
                    os.replace(structureName + '_trajectory.dcd', structureName + '_trajectory-1.dcd')

                combinedSlope = checkTrajPCASlope(structure, structureName + '_trajectory-1.dcd')
                if binding: # if this is a binding simulation, also check to see if the analyte has come unbound from the analyte, and if so, cutoff the simulation
                    self.analyteUnbound = checkMidTrajectoryBinding(structure, structureName + '_trajectory-1.dcd', self.peptide, self.sequence, self.params, cutoffTime=1)
                    print('Analyte came unbound!')

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
        '''
        types of analysis:
        1) folding info (base-base distances)
        2) binding info (analyte-base distances)
        '''
        # wcTraj1 = dnaBasePairDist(u,sequence) # watson-crick base pairing distances (H-bonding) SLOW
        # nucleicAnglesTraj = dnaBaseDihedrals(u,sequence) # base-base backbone dihedrals # slow, old

        wcTraj = wcTrajAnalysis(u)  # watson-crick base pairing distances (H-bonding) FAST but some errors
        baseDistTraj = dnaBaseCOGDist(u)  # FAST, base-base center-of-geometry distances
        nucleicAnglesTraj = nucleicDihedrals(u)  # FAST, new, omits 'chi' angle between ribose and base

        # 2D structure analysis
        pairTraj = getPairs(wcTraj)
        secondaryStructure = analyzeSecondaryStructure(pairTraj)  # find equilibrium secondary structure
        predictedConfig = np.zeros(len(self.sequence))
        for i in range(1, len(predictedConfig) + 1):
            if i in self.pairList:
                ind1, ind2 = np.where(self.pairList == i)
                if ind2 == 0:
                    predictedConfig[self.pairList[ind1][0][0] - 1] = self.pairList[ind1][0][1]
                    predictedConfig[self.pairList[ind1][0][1] - 1] = self.pairList[ind1][0][0]


        predictionError = getSecondaryStructureDistance([np.asarray(secondaryStructure), predictedConfig])[0, 1]  # compare to predicted structure
        print('Secondary Structure Prediction error = %.2f' % predictionError)

        # 3D structure analysis
        representativeIndex, pcTrajectory, eigenvalues = isolateRepresentativeStructure(nucleicAnglesTraj)

        # save this structure as a separate file
        extractFrame(structure, trajectory, representativeIndex, 'repStructure.pdb')

        # we also need a function which processes the energies spat out by the trajectory logs

        analysisDict = {}  # compile our results
        analysisDict['2nd structure error'] = predictionError
        analysisDict['predicted 2nd structure'] = predictedConfig
        analysisDict['actual 2d structure'] = secondaryStructure
        analysisDict['RC trajectories'] = pcTrajectory
        analysisDict['representative structure index'] = representativeIndex
        if analyte != False:
            analysisDict['aptamer-analyte binding'] = bindingAnalysis(u, self.peptide, self.sequence)  # look for contacts between analyte and aptamer
            if self.analyteUnbound:
                analysisDict['aptamer-analyte binding']['analyte came unbound'] = True
            else:
                analysisDict['aptamer-analyte binding']['analyte came unbound'] = False


        return analysisDict
