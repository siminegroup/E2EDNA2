'''
interfaces for other major software used by OpenDNA
Macromolecule builder (mmb)
OpenMM (omm)
LightDock (ld)
'''
from shutil import copyfile
from utils import *
from analysisTools import *
from simtk.openmm import *
import simtk.unit as unit
from simtk.openmm.app import *
from numpy import pi

'''
=> write docstrings
=> clean up functions & formatting
=>  test
'''


class mmb(): # macromolecule builder
    def __init__(self, sequence, pairList, params, ind1, intervalLength=None):
        if params['fold speed'] == 'quick':
            self.template = 'commands.template_quick.dat'  # only for debugging runs - very short
        elif params['fold speed'] == 'normal':
            self.template = 'commands.template.dat'  # default folding algorithm
        elif params['fold speed'] == 'long':
            self.template = 'commands.template_long.dat'  # extended annealing - for difficult sequences
        self.comFile = 'commands.run_fold.dat'
        self.foldSpeed = params['fold speed']
        self.sequence = sequence
        self.temperature = params['temperature']
        self.pairList = pairList
        self.mmbPath = params['mmb']

        self.ind1 = ind1

        self.foldedSequence = 'foldedSequence_{}.pdb'.format(ind1)

        self.intervalLength = intervalLength
        
        self.fileDump = 'mmbFiles_%d'%self.ind1

    def generateCommandFile(self):
        copyfile(self.template, self.comFile)  # make command file
        replaceText(self.comFile, 'SEQUENCE', self.sequence)
        replaceText(self.comFile, 'TEMPERATURE', str(self.temperature - 273)) # probably not important, but we can add the temperature in C
        if self.foldSpeed == 'long':
            replaceText(self.comFile, 'INTERVAL', str(self.intervalLength))

        baseString = '#baseInteraction A IND WatsonCrick A IND2 WatsonCrick Cis'
        lineNum = findLine(self.comFile, baseString)  # find the line number to start enumerating base pairs

        for i in range(len(self.pairList)):
            filledString = 'baseInteraction A {} WatsonCrick A {} WatsonCrick Cis'.format(self.pairList[i, 0], self.pairList[i, 1])
            addLine(self.comFile, filledString, lineNum + 1)

    def fold(self):
        # run fold - sometimes this errors for no known reason - keep trying till it works
        Result = None
        attempts = 0
        while (Result is None) and (attempts < 100):
            try:
                attempts += 1
                os.system(self.mmbPath + ' -c ' + self.comFile + ' > outfiles/fold.out')
                os.replace('frame.pdb', self.foldedSequence)

                Result = 1
                # cleanup
                self.fileDump = 'mmbFiles_%d'%self.ind1
                os.mkdir(self.fileDump)
                os.system('mv last* ' + self.fileDump)
                os.system('mv trajectory.* ' + self.fileDump)
                os.system('mv match.* ' + self.fileDump)

            except:
                pass

    def check2DAgreement(self):
        '''
        check agreement between prescribed 2D structure and the actual fold
        '''
        ''' do it with mda
        u = mda.Universe(self.foldedSequence)
        wcTraj = getWCDistTraj(u)  # watson-crick base pairing distances (H-bonding)
        pairTraj = getPairTraj(wcTraj)
        trueConfig = pairListToConfig(self.pairList, len(self.sequence))
        foldDiscrepancy = getSecondaryStructureDistance([pairTraj[0],trueConfig])[0,1]
        self.foldFidelity = 1-foldDiscrepancy # return the fraction of correctly paired bases
        '''
        # do it with MMB
        f = open('outfiles/fold.out')
        text = f.read()
        f.close()
        text = text.split('\n')
        scores = []
        for i in range(len(text)):
            if "Satisfied baseInteraction's" in text[i]:
                line = text[i].split(':')
                numerator = float(line[-2].split(' ')[1])
                denominator = float(line[-1])
                if denominator == 0: # if there are no pairs, the denominator will be zero (no constraints)
                    scores.append(1)
                else:
                    scores.append(numerator/denominator)

        scores = np.asarray(scores)
        self.foldFidelity = np.amax(scores)

    def run(self):
        self.generateCommandFile()
        self.fold()
        self.check2DAgreement()

        return self.foldFidelity


class ld(): # lightdock
    def __init__(self, aptamer, peptide, params, ind1):
        self.setupPath = params['ld setup path']
        self.runPath = params['ld run path']
        self.genPath = params['lgd generate path']
        self.clusterPath = params['lgd cluster path']
        self.rankPath = params['lgd rank path']
        self.topPath = params['lgd top path']

        self.aptamerPDB = aptamer
        self.peptidePDB = peptide

        self.glowWorms = 300 # params['glowworms']
        self.dockingSteps = params['docking steps']
        self.numTopStructures = params['N docked structures']

        self.pH = params['pH']

        self.ind1 = ind1

    def getSwarmCount(self):
        # identify aptamer size
        dimensions = getMoleculeSize(self.aptamerPDB)
        # compute surface area of rectangular prism with no offset
        surfaceArea = 2 * (dimensions[0] * dimensions[1] + dimensions[1] * dimensions[2] + dimensions[2] * dimensions[1])  # in angstroms
        # compute number of swarms which cover this surface area - swarms are spheres 2 nm in diameter - this is only approximately the right number of swarms, since have no offset, and are using a rectangle
        nSwarms = np.ceil(surfaceArea / (np.pi * 10 ** 2))

        self.swarms = int(nSwarms)  # number of glowworm swarms

    def prepPDBs(self):
        killH(self.aptamerPDB)  # DNA needs to be deprotonated
        addH(self.peptidePDB, self.pH)  # peptide needs to be hydrogenated
        self.aptamerPDB2 = self.aptamerPDB.split('.')[0] + "_noH.pdb"
        self.peptidePDB2 = self.peptidePDB.split('.')[0] + "_H.pdb"
        changeSegment(self.peptidePDB2, 'A', 'B')

    def runLightDock(self):
        # run setup
        os.system(self.setupPath + ' ' + self.aptamerPDB2 + ' ' + self.peptidePDB2 + ' -s ' + str(self.swarms) + ' -g ' + str(self.glowWorms) + ' >> outfiles/lightdockSetup.out')

        # run docking
        os.system(self.runPath + ' setup.json ' + str(self.dockingSteps) + ' -s dna >> outfiles/lightdockRun.out')

    def generateAndCluster(self):
        # generate docked structures and cluster them
        for i in range(self.swarms):
            os.chdir('swarm_%d' % i)
            os.system(self.genPath + ' ../' + self.aptamerPDB2 + ' ../' + self.peptidePDB2 + ' gso_%d' % self.dockingSteps + '.out' + ' %d' % self.glowWorms + ' > /dev/null 2> /dev/null; >> generate_lightdock.list')  # generate configurations
            os.system(self.clusterPath + ' gso_%d' % self.dockingSteps + '.out >> cluster_lightdock.list')  # cluster glowworms
            os.chdir('../')

    def rank(self):
        os.system(self.rankPath + ' %d' % self.swarms + ' %d' % self.dockingSteps + ' >> outfiles/ld_rank.out')  # rank the clustered docking setups

    def extractTopStructures(self):
        # generate top structures
        os.system(self.topPath + ' ' + self.aptamerPDB2 + ' ' + self.peptidePDB2 + ' rank_by_scoring.list %d' % self.numTopStructures + ' >> outfiles/ld_top.out')
        os.mkdir('top_%d'%self.ind1)
        os.system('mv top*.pdb top_%d'%self.ind1 + '/')  # collect top structures (clustering currently dubiously working)

    def extractTopScores(self):
        self.topScores = readInitialLines('rank_by_scoring.list', self.numTopStructures + 1)[1:]
        for i in range(len(self.topScores)):
            self.topScores[i] = float(self.topScores[i].split(' ')[-1])  # get the last number, which is the score

        return self.topScores

    def run(self):
        self.getSwarmCount()
        self.prepPDBs()
        self.runLightDock()
        self.generateAndCluster()
        self.rank()
        self.extractTopStructures()
        self.extractTopScores()
        # cleanup
        dir = 'lightdockOutputs_{}'.format(self.ind1)
        os.mkdir(dir)
        os.system('mv swarm* ' + dir)
        os.system('mv lightdock_* ' + dir)
        os.system('mv setup.json ' + dir)
        os.system('mv *.list ' + dir)
        os.system('mv lightdock.info ' + dir)
        os.system('mv init ' + dir)

class omm(): # openmm
    def __init__(self, structure, params, simTime=None):
        self.pdb = PDBFile(structure)
        self.structureName = structure.split('.')[0]
        self.waterModel = params['water model']
        self.forcefield = ForceField('amber14-all.xml', 'amber14/' + self.waterModel + '.xml')
        self.peptide = params['peptide']

        # System Configuration
        self.nonbondedMethod = params['nonbonded method']
        self.nonbondedCutoff = params['nonbonded cutoff'] * unit.nanometer
        self.ewaldErrorTolerance = params['ewald error tolerance']
        self.constraints = params['constraints']
        self.rigidWater = params['rigid water']
        self.constraintTolerance = params['constraint tolerance']
        self.hydrogenMass = params['hydrogen mass'] * unit.amu

        # Integration Options
        self.dt = params['time step'] / 1000 * unit.picoseconds
        self.temperature = params['temperature'] * unit.kelvin
        self.friction = params['friction'] / unit.picosecond

        self.quickSim = params['test mode']

        # Simulation Options
        if simTime != None:  # we can override the simulation time here
            self.steps = int(simTime * 1e6 // params['time step'])  #
            self.equilibrationSteps = int(params['equilibration time'] * 1e6 // params['time step'])

        else:
            self.steps = int(params['sampling time'] * 1e6 // params['time step'])  # number of steps
            self.equilibrationSteps = int(params['equilibration time'] * 1e6 // params['time step'])      
        
        #platform
        if params['platform'] == 'CUDA':  # 'CUDA' or 'cpu'
            self.platform = Platform.getPlatformByName('CUDA')
            if params['platform precision'] == 'single':
                self.platformProperties = {'Precision': 'single'}
            elif params['platform precision'] == 'double':
                self.platformProperties = {'Precision': 'double'}
        else:
            self.platform = Platform.getPlatformByName('CPU')

        self.reportSteps = int(params['print step'] * 1000 / params['time step'])  # report step in ps, time step in fs
        self.dcdReporter = DCDReporter(self.structureName + '_trajectory.dcd', self.reportSteps)
        #self.pdbReporter = PDBReporter(self.structureName + '_trajectory.pdb', self.reportSteps) # huge files
        self.dataReporter = StateDataReporter('log.txt', self.reportSteps, totalSteps=self.steps, step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, separator='\t')
        self.checkpointReporter = CheckpointReporter(structure.split('.')[0] + '_state.chk', 10000)

        # Prepare the Simulation
        #printRecord('Building system...') # waste of space
        self.topology = self.pdb.topology
        self.positions = self.pdb.positions
        self.system = self.forcefield.createSystem(self.topology, nonbondedMethod=self.nonbondedMethod, nonbondedCutoff=self.nonbondedCutoff, constraints=self.constraints, rigidWater=self.rigidWater, ewaldErrorTolerance=self.ewaldErrorTolerance, hydrogenMass=self.hydrogenMass)
        
        if params['peptide backbone constraint constant'] != 0:
            self.force = CustomTorsionForce('0.5*K*dtheta^2; dtheta = min(diff, 2*' + str(round(pi, 3)) + '-diff); diff = abs(theta - theta0)')
            self.force.addGlobalParameter('K', params['peptide backbone constraint constant'])
            self.force.addPerTorsionParameter('theta0')
            
            self.angles_to_constrain = findAngles()
            
#             self.da_atoms = [atom for atom in self.topology.atoms() if atom.residue.chain.index == 1 and atom.name in {'N', 'CA', 'C'}] #  assumes the chain id of the peptide is 1 - for future releases, will need to be more dynamic       
            
#             self.atom_ids = [atom.index for atom in self.da_atoms]
            
            self.phi_tup = ('C', 'N', 'CA', 'C')
            self.psi_tup = ('N', 'CA', 'C', 'N')
            self.angle_tups = self.phi_tup, self.psi_tup
            radians = unit.radians
            
            rad_conv = pi / 180
            self.nchains = len(self.topology._chains)
            
            for chain in range(self.nchains):                
                for row in self.angles_to_constrain:
                    aa_id, phi, psi, chain_id = row[0], row[1], row[2], row[3]
                    aa_id, phi, psi, chain_id = int(aa_id), float(phi), float(psi), int(chain_id)

                    self.da_atoms = [atom for atom in self.topology.atoms() if atom.residue.chain.index == chain 
                                     and atom.name in {'N', 'CA', 'C'} 
                                     and atom.residue.index in {aa_id, aa_id - 1}]
                    
                    # aa_id - 1 is included to account for the atoms in the previous residue being part of the current residue's dihedrals

                    for i in range(len(self.da_atoms) - 3):
                        self.tup = tuple([atom.name for atom in da_atoms[i:i + 4]])
                        self.tupIndex = tuple([atom.index for atom in da_atoms[i:i + 4]])

                        if self.da_atoms[i + 3].residue.index == aa_id:
                            if self.tup == self.phi_tup:
                                self.force.addTorsion(self.tupIndex[0], 
                                                      self.tupIndex[1], 
                                                      self.tupIndex[2], 
                                                      self.tupIndex[3], (phi * rad_conv,) * radians)

                            elif self.tup == self.psi_tup:
                                self.force.addTorsion(self.tupIndex[0], 
                                                      self.tupIndex[1], 
                                                      self.tupIndex[2], 
                                                      self.tupIndex[3], (psi * rad_conv,) * radians)
            self.system.addForce(self.force)
            
        self.integrator = LangevinMiddleIntegrator(self.temperature, self.friction, self.dt)
        self.integrator.setConstraintTolerance(self.constraintTolerance)
        if params['platform'] == 'CUDA':
            self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform, self.platformProperties)
        elif params['platform'] == 'CPU':
            self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform)
        self.simulation.context.setPositions(self.positions)
            
        

    def doMD(self):
        if not os.path.exists(self.structureName.split('.')[0] + '_state.chk'):
            # Minimize and Equilibrate
            printRecord('Performing energy minimization...')
            if self.quickSim: # if want to do it fast, loosen the tolerances
                self.simulation.minimizeEnergy(tolerance = 20, maxIterations = 100) # default is 10 kJ/mol - also set a max number of iterations
            else:
                self.simulation.minimizeEnergy() # (tolerance = 1 * unit.kilojoules / unit.mole)
            printRecord('Equilibrating...')
            self.simulation.context.setVelocitiesToTemperature(self.temperature)
            self.simulation.step(self.equilibrationSteps)
        else:
            self.simulation.loadCheckpoint(self.structureName.split('.')[0] + '_state.chk')

        # Simulate
        printRecord('Simulating...')
        self.simulation.reporters.append(self.dcdReporter)
        #self.simulation.reporters.append(self.pdbReporter)
        self.simulation.reporters.append(self.dataReporter)
        self.simulation.reporters.append(self.checkpointReporter)
        self.simulation.currentStep = 0
        with Timer() as md_time:
            self.simulation.step(self.steps)  # run the dynamics

        self.simulation.saveCheckpoint(self.structureName.split('.')[0] + '_state.chk')

        self.ns_per_day = (self.steps * self.dt) / (md_time.interval * unit.seconds) / (unit.nanoseconds / unit.day)

        return self.ns_per_day


class nupack():
    def __init__(self, sequence, temperature, ionicStrength):
        self.sequence = sequence
        self.temperature = temperature
        self.ionicStrength = ionicStrength
        self.R = 0.0019872  # ideal gas constant in kcal/mol/K


    def energyCalc(self):
        '''
        sequence is DNA FASTA
        temperature in C or K
        ionicStrength in M
        return a bunch of analysis
        '''
        if self.temperature > 273:  # auto-detect Kelvins
            gap = 2 * self.R * self.temperature
            CelsiusTemperature = self.temperature - 273
        else:
            gap = 2 * self.R * (self.temperature + 273)  # convert to Kelvin for kT
            CelsiusTemperature = self.temperature

        A = Strand(self.sequence, name='A')
        comp = Complex([A], name='AA')
        set1 = ComplexSet(strands=[A], complexes=SetSpec(max_size=1, include=[comp]))
        model1 = Model(material='dna', celsius=CelsiusTemperature, sodium=self.ionicStrength)
        results = complex_analysis(set1, model=model1, compute=['pfunc', 'mfe', 'subopt', 'pairs'], options={'energy_gap': gap})
        self.output = results[comp]


    def structureAnalysis(self):
        '''
        input nupack structure and return a bunch of analysis
        :param nupackStructure:
        :return:
        '''
        nStructures = len(self.output.subopt)
        self.ssDict = {}  # best structure dictionary
        self.ssDict['2d string'] = []
        self.ssDict['pair list'] = []
        self.ssDict['config'] = np.zeros((nStructures,len(self.sequence)))
        self.ssDict['num pairs'] = np.zeros(nStructures).astype(int)
        self.ssDict['pair frac'] = np.zeros(nStructures)
        self.ssDict['num pins'] = np.zeros(nStructures).astype(int)
        self.ssDict['state prob'] = np.zeros(nStructures)

        for i in range(nStructures):
            ssString = str(self.output.subopt[i].structure)
            self.ssDict['2d string'].append(ssString)
            self.ssDict['pair list'].append(ssToList(ssString))
            self.ssDict['config'][i] = pairListToConfig(self.ssDict['pair list'][-1],len(self.sequence))
            self.ssDict['state prob'][i] = np.exp(-self.output.subopt[i].energy / self.R / self.temperature) / float(self.output.pfunc)
            self.ssDict['num pairs'][i] = ssString.count('(')  # number of pairs
            self.ssDict['pair frac'][i] = 2 * self.ssDict['num pairs'][i] / len(self.sequence)  # fraction of paired bases

            nPins = 0  # number of distinct hairpins
            indA = 0
            for j in range(len(self.sequence)):
                if ssString[j] == '(':
                    indA += 1
                elif ssString[j] == ')':
                    indA -= 1
                    if indA == 0:  # if we come to the end of a distinct hairpin
                        nPins += 1

            self.ssDict['num pins'][i] = nPins

    def run(self):
        self.energyCalc()
        self.structureAnalysis()
        return self.ssDict

