# import statements
from shutil import copyfile
import re
from seqfold import dg, fold
from utils import *
from simtk.openmm import *
import simtk.unit as unit
from simtk.openmm.app import *
from pdbfixersource import PDBFixer


'''
script which accepts a DNA sequence and an analyte, and returns the binding affinity

To-Do:
==> testing
==> docking
==> add back PCA to 2nd structure, instead of clustering
==> auto-equilibration
    ==> need to implement checkpointing
==> multi-state comparision in 2d and 3d
==:>> analysis
    ==> print summary output
        
==> WC pairing tools don't agree but I can't find the error - in any case at short range they seem to be pretty close

Notes
==> implicit solvent - amber10 codes don't agree
'''


params = {}
params['device'] = 'local' # 'local' or 'cluster'
params['platform'] = 'CUDA' # no alternative at the moment
params['platform precision'] = 'single'

if params['device'] == 'cluster':
    params['run num'] = get_input()
elif params['device'] == 'local':
    params['run num'] = 0 # manual setting, for 0, do a fresh run, for != 0, pickup on a previous run.

# Simulation parameters
params['force field'] = 'AMBER' # this does nothing
params['water model'] = 'tip3p' # 'tip3p' (runs on amber 14), 'implicit' (runs on amber 10 - not working)
params['equilibration time'] = 0.01 # equilibration time in nanoseconds
params['sampling time'] = 0.05 # sampling time in nanoseconds
params['auto sampling'] = True # 'True' run sampling until RC's equilibrate + 'sampling time', 'False' just run sampling for 'sampling time'
params['time step'] = 3.0 # in fs
params['print step'] = 1 # printout step in ps - want there to be more than 2 and less than 100 total frames in any trajectory

params['box offset'] = 1.0 # nanometers
params['barostat interval'] = 25
params['friction'] = 1.0 # 1/picosecond
params['nonbonded method'] = PME
params['nonbonded cutoff'] = 1.0 # nanometers
params['ewald error tolerance'] = 5e-4
params['constraints'] = HBonds
params['rigid water'] = True
params['constraint tolerance'] = 1e-6
params['hydrogen mass'] = 4.0 # in amu

# physical params
params['pressure'] = 1 # atmospheres
params['temperature'] = 310 # Kelvin
params['ionic strength'] = .163 # mmol
params['pH'] = 7.4 # simulation will automatically protonate the peptide up to this pH

# paths
if params['device'] == 'local':
    params['workdir'] = 'C:/Users\mikem\Desktop/mmruns'
    params['mmb'] = 'C:/Users/mikem/Desktop/Installer.2_14.Windows/MMB.2_14.exe'
    params['mmb params'] = 'lib/parameters.csv'
    params['mmb template'] = 'lib/commands.template.dat'
elif params['device'] == 'cluster':
    params['workdir'] = '/home/kilgourm/scratch/mmruns' # specify your working directory here
    params['mmb'] = '~/projects/def-simine/programs/MMB/Installer.2_14.Linux64/MMB.2_14.Linux64'
    params['mmb params'] = 'lib/parameters.csv'
    params['mmb template'] = 'lib/commands.template.dat'

# structure files
params['analyte pdb'] = 'lib/peptide/peptide.pdb' # optional - currently not used


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
            copyfile(params['mmb template'], self.workDir + '/commands.template.dat')

        else:
            self.workDir = self.params['workdir'] + '/' + 'run%d' %self.params['run num']

        # move to working dr
        os.chdir(self.workDir)


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
        if self.checkpoints < 2:
            self.ssString, self.pairList = self.getSecondaryStructure(sequence)


        if self.checkpoints < 3:  # if we have the secondary structure
            self.foldSequence(sequence, self.pairList)


        if self.checkpoints < 4: # run MD
            print('Running free aptamer dynamics')

            if self.params['water model'] != 'implicit':
                self.prepPDB('sequence.pdb',MMBCORRECTION=True) # add periodic box and appropriate protons
                dict = self.autoMD('sequence_processed.pdb')
                aa = 1
                #self.runMD('sequence_processed.pdb')
            else:
                self.runMD('sequence.pdb')

            print('Free aptamer simulation speed %.1f'%self.ns_per_day+' ns/day')
            os.rename('trajectory.dcd','Aptamer.dcd')
            os.rename('sequence_processed.pdb','Aptamer.pdb')
            cleanTrajectory('aptamer.pdb','aptamer.dcd') # make a trajectory without waters and ions and stuff
            aptamerDict = self.analyzeTrajectory('cleanaptamer.pdb','cleanaptamer.dcd',False)
            writeCheckpoint('Free Aptamer Sampled')


        if (self.checkpoints < 5) and (self.peptide != False): # run MD on the complexed structure
            print('Running Binding Simulation')
            buildPeptide(self.peptide)
            combinePDB('repStructure.pdb','peptide.pdb') # this is where we would implement a docker
            if self.params['water model'] != 'implicit':
                self.prepPDB('combined.pdb',MMBCORRECTION=False) # doesn't count the peptide!!
                self.runMD('combined_processed.pdb') # would be nice here to have an option to change the paramters for the binding run - or make it adaptive based on feedback from the run
            else:
                self.runMD('combined.pdb')
            print('Binding complex simulation speed %.1f' % self.ns_per_day + ' ns/day')
            os.rename('trajectory.dcd','complex.dcd')
            os.rename('combined_processed.pdb','complex.pdb')
            cleanTrajectory('complex.pdb','complex.dcd')
            bindingDict = self.analyzeTrajectory('cleancomplex.pdb','cleancomplex.dcd',True)


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
            temperature = 37.0
            dg(sequence, temp=temperature) # get energy of the structure
            #print(round(sum(s.e for s in structs), 2)) # predicted energy of the final structure

            structs = fold(sequence) # identify structural features
            desc = ["."] * len(sequence)
            pairList = []
            for s in structs:
                pairList.append(s.ij[0])
                pairList[-1]# list of bound pairs indexed from 1
                if len(s.ij) == 1:
                    i, j = s.ij[0]
                    desc[i] = "("
                    desc[j] = ")"

            ssString = "".join(desc)
            pairList = np.asarray(pairList) + 1

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

        # run fold
        os.system(self.params['mmb'] + ' -c ' + comFile + ' > outfiles/fold.out')
        os.rename('frame.pdb','sequence.pdb')
        writeCheckpoint("Folded Sequence")


    def prepPDB(self, file, MMBCORRECTION=False):
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
        ionicStrength = float(self.params['ionic strength']) * unit.molar  # not implemented
        positiveIon = 'Na+'  # params['positiveion']+'+'
        negativeIon = 'Cl-'  # params['negativeion']+'-'

        fixer.addMissingHydrogens(pH=self.params['pH'])
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
        steps = int(self.params['sampling time'] * 1e6 // self.params['time step']) # number of steps
        equilibrationSteps = int(self.params['equilibration time'] * 1e6 // self.params['time step'])
        if self.params['platform'] == 'CUDA': # no alternative at the moment
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
        system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff, constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance, hydrogenMass=hydrogenMass)
        integrator = LangevinMiddleIntegrator(temperature, friction, dt)
        integrator.setConstraintTolerance(constraintTolerance)
        simulation = Simulation(topology, system, integrator, platform, platformProperties)
        simulation.context.setPositions(positions)


        if not os.path.exists('state.chk'):
            # Minimize and Equilibrate
            print('Performing energy minimization...')
            simulation.minimizeEnergy()
            print('Equilibrating...')
            simulation.context.setVelocitiesToTemperature(temperature)
            simulation.step(equilibrationSteps)
        else:
            simulation.loadCheckpoint('state.chk')

        # Simulate
        print('Simulating...')
        simulation.reporters.append(dcdReporter)
        #simulation.reporters.append(pdbReporter)
        simulation.reporters.append(dataReporter)
        simulation.reporters.append(checkpointReporter)
        simulation.currentStep = 0
        with Timer() as md_time:
            simulation.step(steps)

        simulation.saveCheckpoint('state.chk')


        self.ns_per_day = (steps * dt) / (md_time.interval * unit.seconds) / (unit.nanoseconds/unit.day)


    def autoMD(self,structure):
        '''
        run MD until 'equilibration', then sample
        :param structure:
        :return:
        '''
        self.runMD(structure)
        self.runMD(structure)
        cleanTrajectory(structure,'trajectory.dcd')
        runDict = self.analyzeTrajectory('clean'+structure,'cleantrajectory.dcd',False)
        pcTrajectory = runDict['RC trajectories']
        # see if we're equilibrated in the PC basis - see if there's any net slope - if there is a big slope we're not equilibrated AND/OR not sampling enough
        slopes = np.zeros(pcTrajectory.shape[-1])
        for i in range(len(slopes)):
            slopes[i] = np.abs(np.polyfit(np.arange(len(pcTrajectory)),pcTrajectory[:,i],1)[0])

        return runDict


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
        baseDists = dnaBaseCOGDist(u,sequence) # FAST, base-base center-of-geometry distances
        baseAngles = nucleicDihedrals(u) # FAST, new, omits 'chi' angle between ribose and base

        # 2D structure analysis
        pairingTrajectory = getPairs(baseWC)
        secondaryStructure = analyzeSecondaryStructure(pairingTrajectory)  # find equilibrium secondary structure
        predictedConfig = np.zeros(len(sequence))
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
            bindingInfo = bindingAnalysis(u,self.peptide, self.sequence) # look for contacts between analyte and aptamer


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

'''
==============================================================
'''

if __name__ == '__main__':
    sequence = 'CGCTTTGCG' #'ACCTGGGGGAGTATTGCGGAGGAAGGT' #ATP binding aptamer
    peptide = False #'YQT'#'YQTQTNSPRRAR'
    opendna = opendna(sequence,peptide, params)
    opendnaOutput = opendna.run() # retrieve binding score and center-of-mass time-series

    #os.chdir('C:/Users\mikem\Desktop/tinkerruns\clusterTests/fullRuns/run36')
    #comProfile, bindingScore = evaluateBinding('complex_sampling.arc') # normally done inside the run loop
    '''
    timeEq, potEq, kinEq = getEnergy()
    timeSa, potSa, kinSa = getEnergy()

    outputs = {}
    outputs['time equil'] = timeEq
    outputs['pot equil'] = potEq
    outputs['kin equil'] = kinEq
    outputs['time sampling'] = timeSa
    outputs['pot sampling'] = potSa
    outputs['kin sampling'] = kinSa
    outputs['binding score'] = bindingScore
    outputs['com profile'] = comProfile
    outputs['params'] = params
    np.save('bindingOutputs', outputs) # unpack with load then outputs.item()
    '''

