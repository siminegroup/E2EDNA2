# import statements
from shutil import copyfile
import re
from seqfold import dg, fold
from secondaryStructureScripts import *
from utils import *
from simtk.openmm import *
import simtk.unit as unit
from simtk.openmm.app import *
from pdbfixersource import PDBFixer


class opendna():
    def __init__(self,sequence,peptide,params):
        self.params = params
        self.sequence = sequence
        self.peptide = peptide

        self.setup()
        self.getCheckpoint()


    def setup(self):
        '''
        setup working directory
        copy in relevant xyz and keyfiles
        move to relevant directory
        :return:
        '''
        if self.params['run num'] == 0:
            self.makeNewWorkingDirectory()
            os.mkdir(self.workDir + '/outfiles')
            # copy structure files
            copyfile(self.params['analyte pdb'], self.workDir + '/analyte.pdb')

            # copy MMB params
            copyfile(self.params['mmb params'], self.workDir + '/parameters.csv')
            copyfile(self.params['mmb template'], self.workDir + '/commands.template.dat')

        else:
            self.workDir = self.params['workdir'] + '/' + 'run%d' %self.params['run num']

        # move to working dr
        os.chdir(self.workDir)

        os.environ["LD_LIBRARY_PATH"] = self.params['mmb dir'] # export the path to the MMB library


    def makeNewWorkingDirectory(self):    # make working directory
        '''
        make a new working directory
        non-overlapping previous entries
        :return:
        '''
        workdirs = glob.glob(self.params['workdir'] + '/' + 'run*') # check for prior working directories
        if len(workdirs) > 0:
            prev_runs = []
            for i in range(len(workdirs)):
                prev_runs.append(int(workdirs[i].split('run')[-1]))

            prev_max = max(prev_runs)
            self.workDir = self.params['workdir'] + '/' + 'run%d' %(prev_max + 1)
            os.mkdir(self.workDir)
            print('Starting Fresh Run %d' %(prev_max + 1))
        else:
            self.workDir = self.params['workdir'] + '/' + 'run1'
            os.mkdir(self.workDir)


    def getCheckpoint(self):
        '''
        identify if any work has been done in this directory, and if so, where to pick up
        :return:
        '''
        try:
            f = open('checkpoint.txt','r')#if it does, see how long it is
            text = f.read()
            f.close()
            self.checkpoints = len([m.start() for m in re.finditer('\n', text)])

            if self.checkpoints == 1:
                self.checkpoint = 'initialized'
            elif self.checkpoints == 2:
                self.checkpoint = 'something else'
                '''
                and so on
                '''

            print("Resuming Run #%d"%int(self.params['run num']))


        except:
            f = open('checkpoint.txt','w')
            f.write('Started!\n')
            f.close()
            self.checkpoints = 1
            self.checkpoint = 'initialized'


    def run(self):
        '''
        run the binding simulation end-to-end
        consult checkpoints to not repeat prior steps
        :return:
        '''
        if True:#self.checkpoints < 2:
            self.ssString, self.pairList = self.getSecondaryStructure(self.sequence)


        if self.checkpoints < 2:  # if we have the secondary structure
            self.foldSequence(self.sequence, self.pairList)
            writeCheckpoint("Folded Sequence")


        if self.checkpoints < 3: # run MD
            print('Running free aptamer dynamics')

            if self.params['water model'] != 'implicit':
                self.prepPDB('sequence.pdb',MMBCORRECTION=True,waterBox=True) # add periodic box and appropriate protons
                #dict = self.autoMD('sequence_processed.pdb') #auto convergence - in progress
                self.runMD('sequence_processed.pdb')
            else:
                self.prepPDB('sequence.pdb',MMBCORRECTION=True,waterBox=False) # add periodic box and appropriate protons
                self.runMD('sequence.pdb')

            print('Free aptamer simulation speed %.1f'%self.ns_per_day+' ns/day')
            os.rename('trajectory.dcd','Aptamer.dcd')
            os.rename('sequence_processed.pdb','Aptamer.pdb')
            cleanTrajectory('aptamer.pdb','aptamer.dcd') # make a trajectory without waters and ions and stuff
            aptamerDict = self.analyzeTrajectory('cleanaptamer.pdb','cleanaptamer.dcd',False)
            writeCheckpoint('Free Aptamer Sampled')


        if (self.checkpoints < 4) and (self.peptide != False): # run MD on the complexed structure
            print('Docking')
            buildPeptide(self.peptide)

            lightDock('repStructure.pdb', 'peptide.pdb',self.params) # get optimal docked structure
            copyfile('top/top_1.pdb','complex_0.pdb')

            print('Running Binding Simulation')
            if self.params['water model'] != 'implicit':
                self.prepPDB('complex_0.pdb',MMBCORRECTION=False,waterBox=True)
            else:
                self.prepPDB('complex_0.pdb',MMBCORRECTION=False,waterBox=False)

            self.runMD('complex_0_processed.pdb') # would be nice here to have an option to change the paramters for the binding run - or make it adaptive based on feedback from the run
            print('Binding complex simulation speed %.1f' % self.ns_per_day + ' ns/day')
            os.rename('trajectory.dcd','complex.dcd')
            os.rename('complex_0_processed.pdb','complextraj.pdb')
            cleanTrajectory('complex.pdb','complextraj.dcd')
            bindingDict = self.analyzeTrajectory('cleancomplex.pdb','cleancomplextraj.dcd',True)

        if self.peptide != False:
            return [aptamerDict,bindingDict]
        else:
            return aptamerDict

#=======================================================
#=======================================================
#=======================================================


    def getSecondaryStructure(self, sequence):
        '''
        get the secondary structure for a given sequence
        using seqfold here - identical features are available using nupack, though results are sometimes different
        :param sequence:
        :return: a dot-bracket string and list of paired bases (assuming single-strand DNA aptamer)
        '''
        print("Get Secondary Structure")
        if os.path.exists("pre_fold.pdb"):  # if we have a pre-folded structure, do nothing
            pass
        else:
            if self.params['secondary structure engine'] == 'seqfold':
                ssString, pairList = getSeqfoldStructure(sequence, self.params['temperature']) # seqfold guess
            elif self.params['secondary structure engine'] == 'NUPACK':
                ssString, pairList = getNupackStructure(sequence,self.params['temperature'],self.params['ionic strength']) # NUPACK guess

            return ssString, pairList

        writeCheckpoint("Got Secondary Structure")


    def foldSequence(self, sequence, pairList):
        '''
        generate a coarsely folded structure for the aptamer
        :param sequence: the ssDNA sequence to be folded
        :param pairList: list of binding base pairs
        :return:
        '''
        # write pair list as forces to the MMB command file
        print("Folding Sequence")
        comFile = 'commands.fold.dat'  # name of command file
        copyfile('commands.template.dat', comFile)  # make command file
        replaceText(comFile, 'SEQUENCE', sequence)

        baseString = '#baseInteraction A IND WatsonCrick A IND2 WatsonCrick Cis'
        lineNum = findLine(comFile, baseString)  # find the line number to start enumerating base pairs

        for i in range(len(pairList)):
            filledString = 'baseInteraction A {} WatsonCrick A {} WatsonCrick Cis'.format(pairList[i,0], pairList[i,1])
            addLine(comFile, filledString, lineNum + 1)

        # run fold - sometimes this errors for no known reason - keep trying till it works
        Result = None
        attempts = 0
        while (Result == None) and (attempts < 100):
            try:
                attempts += 1
                os.system(self.params['mmb'] + ' -c ' + comFile + ' > outfiles/fold.out')
                os.rename('frame.pdb','sequence.pdb')
                Result = 1
            except:
                pass


    def prepPDB(self, file, MMBCORRECTION=False, waterBox=True):
        '''
        soak pdb file in water box
        :param file:
        :return:
        '''
        if MMBCORRECTION:
            replaceText(file, '*', "'") # due to a bug in this version of MMB - structures are encoded improperly - this fixes it

        fixer = PDBFixer(filename=file)
        padding, boxSize, boxVectors = None, None, None
        geompadding = float(self.params['box offset']) * unit.nanometer
        padding = geompadding

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms() # may need to optimize bonding here

        fixer.addMissingHydrogens(pH=self.params['pH']) # add missing hydrogens


        if waterBox == True:
            ionicStrength = float(self.params['ionic strength']) * unit.molar  # not implemented
            positiveIon = 'Na+'  # params['positiveion']+'+'
            negativeIon = 'Cl-'  # params['negativeion']+'-'
            fixer.addSolvent(boxSize, padding, boxVectors, positiveIon, negativeIon, ionicStrength)

        PDBFile.writeFile(fixer.topology, fixer.positions, open(file.split('.pdb')[0] + '_processed.pdb', 'w'))


    def runMD(self,structure):
        '''
        run OpenMM dynamics
        minimize, equilibrate, and sample
        :param structure:
        :return:
        '''

        pdb = PDBFile(structure)
        waterModel = self.params['water model']
        if waterModel != 'implicit':
            forcefield = ForceField('amber14-all.xml', 'amber14/' + waterModel + '.xml')
        else:
            forcefield = ForceField('amber10.xml','amber10_obc.xml')

        # System Configuration
        if self.params['water model'] != 'implicit':
            nonbondedMethod = self.params['nonbonded method']
        else:
            nonbondedMethod = CutoffNonPeriodic
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
        steps = int(self.params['sampling time'] * 1e6 // self.params['time step']) # number of steps
        equilibrationSteps = int(self.params['equilibration time'] * 1e6 // self.params['time step'])
        if self.params['platform'] == 'CUDA': # 'CUDA' or 'cpu'
            platform = Platform.getPlatformByName('CUDA')
            if self.params['platform precision'] == 'single':
                platformProperties = {'Precision': 'single'}
            elif self.params['platform precision'] == 'double':
                platformProperties = {'Precision': 'double'}
        else:
            platform = Platform.getPlatformByName('CPU')

        reportSteps = int(self.params['print step'] * 1000 / self.params['time step']) # report step in ps, time step in fs
        dcdReporter = DCDReporter('trajectory.dcd', reportSteps)
        #pdbReporter = PDBReporter('trajectory.pdb', reportSteps)
        dataReporter = StateDataReporter('log.txt', reportSteps, totalSteps=steps, step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, separator='\t')
        checkpointReporter = CheckpointReporter('checkpoint.chk', 10000)

        # Prepare the Simulation
        print('Building system...')
        topology = pdb.topology
        positions = pdb.positions
        if self.params['water model'] != 'implicit':
            system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff, constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance, hydrogenMass=hydrogenMass)
        else:
            system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance, hydrogenMass=hydrogenMass)
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
        #simulation.reporters.append(pdbReporter)
        simulation.reporters.append(dataReporter)
        simulation.reporters.append(checkpointReporter)
        simulation.currentStep = 0
        with Timer() as md_time:
            simulation.step(steps) # run the dynamics

        simulation.saveCheckpoint(structure.split('.')[0] + '_state.chk')

        self.ns_per_day = (steps * dt) / (md_time.interval * unit.seconds) / (unit.nanoseconds/unit.day)


    def autoMD(self,structure):
        '''
        run MD until either for a set amount of time or until the dynamics 'converge'
        optionally, sample after we reach convergence ("equilibration")
        :param structure:
        :return:
        '''
        if self.params['auto sampling'] == False: # just run MD once
            self.runMD(structure)
        elif self.params['auto sampling'] == True: # run until we detect equilibration, then sample
            converged = 0
            iter = 0
            while converged == 0:
                iter += 1
                self.runMD(structure)

                if iter > 1: # if we have multiple trajectory segments, combine them
                    appendTrajectory(structure,'trajectory-1.dcd','trajectory.dcd')

                converged = checkTrajConvergence(structure,'trajectory-1.dcd')


    def analyzeTrajectory(self,structure,trajectory,analyte):
        '''
        analyze trajectory for aptamer fold + analyte binding information
        :param trajectory:
        :return:
        '''
        u = mda.Universe(structure, trajectory)
        '''
        types of analysis:
        1) folding info (base-base distances)
        2) binding info (analyte-base distances)
        '''
        #baseWC1 = dnaBasePairDist(u,sequence) # watson-crick base pairing distances (H-bonding) SLOW
        #baseAngles = dnaBaseDihedrals(u,sequence) # base-base backbone dihedrals # slow, old

        baseWC = wcTrajAnalysis(u) # watson-crick base pairing distances (H-bonding) FAST but some errors
        baseDists = dnaBaseCOGDist(u) # FAST, base-base center-of-geometry distances
        baseAngles = nucleicDihedrals(u) # FAST, new, omits 'chi' angle between ribose and base

        # 2D structure analysis
        pairingTrajectory = getPairs(baseWC)
        secondaryStructure = analyzeSecondaryStructure(pairingTrajectory)  # find equilibrium secondary structure
        predictedConfig = np.zeros(len(self.sequence))
        for i in range(1, len(predictedConfig) + 1):
            if i in self.pairList:
                ind1, ind2 = np.where(self.pairList == i)
                if ind2 == 0:
                    predictedConfig[self.pairList[ind1][0][0]-1] = self.pairList[ind1][0][1]
                    predictedConfig[self.pairList[ind1][0][1]-1] = self.pairList[ind1][0][0]

        predictionError = getSecondaryStructureDistance([np.asarray(secondaryStructure), predictedConfig])[0,1]  # compare to predicted structure
        print('Secondary Structure Prediction error = %.2f'%predictionError)

        # 3D structure analysis
        representativeIndex, pcTrajectory = isolateRepresentativeStructure(baseAngles)

        # save this structure as a separate file
        extractFrame(structure, trajectory, representativeIndex)

        if analyte != False: # if we have an analyte
            bindingInfo = bindingAnalysis(u, self.peptide, self.sequence) # look for contacts between analyte and aptamer


        # we also need a function which processes the energies spat out by the trajectory thingy


        analysisDict = {} # compile our results
        analysisDict['2nd structure error'] = predictionError
        analysisDict['predicted 2nd structure'] = predictedConfig
        analysisDict['actual 2nd structure'] = secondaryStructure
        analysisDict['RC trajectories'] = pcTrajectory
        analysisDict['representative structure index'] = representativeIndex
        if analyte != False:
            analysisDict['aptamer-analyte binding'] = bindingInfo

        return analysisDict

