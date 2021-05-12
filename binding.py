# import statements
from shutil import copyfile
import re
from seqfold import dg, fold
from utils import *
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
import simtk.unit as unit
from MDAnalysis.analysis import distances



'''
script which accepts a DNA sequence and an analyte, and returns the binding affinity

To-Do:
==> testing
==> upgrade peptide source
    > incorporate protonation
==> check peptide fold
==> analysis
    -> fold
    -> binding
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
params['equilibration time'] = 0.01 # equilibration time in nanoseconds
params['sampling time'] = 0.1 # sampling time in nanoseconds
params['time step'] = 2.0 # in fs
params['print step'] = 1 # printout step in ps - want there to be more than 2 and less than 100 total frames in any trajectory

params['box offset'] = 1.0 # nm
params['barostat interval'] = 25
params['friction'] = 1.0 # /picosecond
params['nonbonded method'] = PME
params['nonbonded cutoff'] = 1.0 # nanometers
params['ewald error tolerance'] = 5e-4
params['constraints'] = HBonds
params['rigid water'] = True
params['constraint tolerance'] = 1e-6
params['hydrogen mass'] = 1 # in amu

# physical params
params['num charges'] = 0 # number of positive charges (Na+) to add to simulation box
params['pressure'] = 1 # atmospheres
params['temperature'] = 300 # Kelvin
params['ionic strength'] = 0 # mmol
params['pH'] = 7.0

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
params['analyte pdb'] = 'lib/peptide/peptide.pdb'


class binder():
    def __init__(self,sequence,params):
        self.params = params
        self.sequence = sequence

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
        consult checkpoints to not repeat prior outputs
        :return:
        '''
        # get secondary structure
        print("Get Secondary Structure")
        if os.path.exists("pre_fold.pdb"):  # if we have a pre-folded structure, do nothing
            pass
        else:
            self.ssString, self.pairList = self.getSecondaryStructure(sequence)
        writeCheckpoint("Got Secondary Structure")


        if self.checkpoints < 3:  # if we have the secondary structure
            # fold it!
            print("Folding Sequence")
            self.foldSequence(sequence, self.pairList)
            writeCheckpoint("Folded Sequence")


        if self.checkpoints < 4:  # if we have the secondary structure
            # fold it!
            print("Preparing Box")
            self.PrepPDB('sequence.pdb')
            writeCheckpoint("Box Solvated")


        if self.checkpoints < 5: # run MD
            print('Running Dynamics')
            self.runMD('sequence_processed.pdb')
            writeCheckpoint('Complex Sampled')

            baseDists = self.analyzeTrajectory('sequence_processed.pdb','trajectory.dcd')

            return baseDists


#=======================================================
# =======================================================
# =======================================================


    def getSecondaryStructure(self, sequence):
        '''
        get the secondary structure for a given sequence
        using seqfold here - identical features are available using nupack, though results are sometimes different
        :param sequence:
        :return: a dot-bracket string and list of paired bases (assuming single-strand DNA aptamer)
        '''
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


    def foldSequence(self, sequence, pairList):
        '''
        generate a coarsely folded structure for the aptamer
        :param sequence: the ssDNA sequence to be folded
        :param pairList: list of binding base pairs
        :return:
        '''
        # write pair list as forces to the MMB command file
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


    def PrepPDB(self, file):
        '''
        soak pdb file in water box
        :param file:
        :return:
        '''
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
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

        # System Configuration
        nonbondedMethod = self.params['nonbonded method']
        nonbondedCutoff = self.params['nonbonded cutoff'] * nanometer
        ewaldErrorTolerance = self.params['ewald error tolerance']
        constraints = self.params['constraints']
        rigidWater = self.params['rigid water']
        constraintTolerance = self.params['constraint tolerance']
        hydrogenMass = self.params['hydrogen mass'] * amu

        # Integration Options
        dt = self.params['time step'] / 1000 * picoseconds
        temperature = self.params['temperature'] * kelvin
        friction = self.params['friction'] / picosecond

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

        # Minimize and Equilibrate
        print('Performing energy minimization...')
        simulation.minimizeEnergy()
        print('Equilibrating...')
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation.step(equilibrationSteps)

        # Simulate
        print('Simulating...')
        simulation.reporters.append(dcdReporter)
        #simulation.reporters.append(pdbReporter)
        simulation.reporters.append(dataReporter)
        simulation.reporters.append(checkpointReporter)
        simulation.currentStep = 0
        simulation.step(steps)


    def analyzeTrajectory(self,structure,trajectory):
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

        #base-base distances
        # identify DNA residues

        # compute and collate all their center-of-geometry distances
        baseDists = np.zeros((len(u.trajectory),len(self.sequence),len(self.sequence)))
        tt = 0
        for ts in u.trajectory:
            dnaResidues = u.segments[0]
            posMat = np.zeros((len(self.sequence),3))
            for i in range(len(self.sequence)):
                posMat[i] = dnaResidues.residues[i].atoms.center_of_geometry()

            baseDists[tt,:,:] = distances.distance_array(posMat,posMat,box=u.dimensions)
            tt += 1

        # we could then go on to analyze e.g., the closest bases, and see if their identites or stable
        # or look at global reorganization - though in 3D this may be substantial

        # when we get the analyte - we can do something very similar here

        # we also need a function which processes the energies spat out by the trajectory thingy

        return baseDists

'''
==============================================================
'''

if __name__ == '__main__':
    sequence = 'ATGCTAGCGA'
    peptide = 'AAA'
    binder = binder(sequence, params)
    baseDists = binder.run() # retrieve binding score and center-of-mass time-series

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

