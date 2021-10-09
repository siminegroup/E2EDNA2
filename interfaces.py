"""
interfaces for other major software used by OpenDNA
MacroMolecule Builder (MMB): installed PeptideBuilder on pip. Is it enough?
OpenMM (omm)
LightDock (ld)
NUPACK
"""

from shutil import copyfile
from numpy import pi
from nupack import *

from simtk.openmm import *  # In fact what happens under the hood: from openmm import *
from simtk.openmm.app import *
import simtk.unit as unit

# from openmm import *
# from openmm.app import *
# import openmm.unit as unit

from utils import *
from analysisTools import *

# installed openmm: http://docs.openmm.org/7.5.0/userguide/application.html#installing-openmm
# openmm is installed as an individual package;
# package simtk is created with packages "openmm" and "unit" inside but literally import the stand-alone openmm.


class nupack:
    def __init__(self, sequence, temperature, ionicStrength, mgConc=0):
        self.sequence = sequence
        self.temperature = temperature
        self.ionicStrength = ionicStrength
        self.naConc = ionicStrength
        self.mgConc = mgConc
        self.R = 0.0019872  # ideal gas constant in kcal/mol/K

    def run(self):
        self.energyCalc()
        self.structureAnalysis()
        return self.ssDict

    def energyCalc(self):
        """
        sequence is DNA FASTA format
        temperature in C or K
        ionicStrength in Molar
        Ouput a lot of analysis results in self.output
        :return:
        """
        if self.temperature > 273:  # auto-detect Kelvins
            gap = 2 * self.R * self.temperature
            CelsiusTemprature = self.temperature - 273
        else:
            gap = 2 * self.R * (self.temperature + 273)  # convert to Kelvin fir kT
            CelsiusTemprature = self.temperature

        A = Strand(self.sequence, name='A')
        comp = Complex([A], name='AA')
        set1 = ComplexSet(strands=[A], complexes=SetSpec(max_size=1, include=[comp]))
        model1 = Model(material='dna', celsius=CelsiusTemprature, sodium=self.naConc, magnesium=self.mgConc)
        results = complex_analysis(set1, model=model1, compute=['pfunc', 'mfe', 'subopt', 'pairs'], options={'energy_gap': gap})
        self.output = results[comp]
    # TODO: what do these nupack functions do?

    def structureAnalysis(self):
        """
        Input nupack structure. return a bunch of analysis results
        :return:
        """
        nStructures = len(self.output.subopt)
        self.ssDict = {}
        self.ssDict['2d string'] = []
        self.ssDict['pair list'] = []
        self.ssDict['config'] = np.zeros((nStructures, len(self.sequence)))
        self.ssDict['num pairs'] = np.zeros(nStructures).astype(int)
        self.ssDict['pair frac'] = np.zeros(nStructures)
        self.ssDict['num pins'] = np.zeros(nStructures).astype(int)
        self.ssDict['state prob'] = np.zeros(nStructures)

        for i in range(nStructures):
            ssString = str(self.output.subopt[i].structure)
            self.ssDict['2d string'].append(ssString)
            self.ssDict['pair list'].append(ssToList(ssString))
            self.ssDict['config'][i] = pairListToConfig(self.ssDict['pair list'][-1], len(self.sequence))
            self.ssDict['state prob'][i] = np.exp(-self.output.subopt[i].energy / self.R / self.temperature / float(self.output.pfunc))
            self.ssDict['num pairs'][i] = ssString.count('(')  # number of pairs
            self.ssDict['pair frac'][i] = 2 * self.ssDict['num pairs'][i] / len(self.sequence)

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


class mmb:  # MacroMolecule Builder (MMB)
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
        self.foldedSequence = 'foldedSequence_{}.pdb'.format(ind1)  # output structures of MMB
        self.intervalLength = intervalLength
        self.fileDump = 'mmbFiles_%d' % self.ind1  # directory to save the mmb run files

        self.foldFidelity = 0

    def run(self):
        self.generateCommandFile()
        self.fold()
        self.check2DAgreement()

        return self.foldFidelity

    # # myOwn
    def generateCommandFile(self):
        copyfile(self.template, self.comFile)  # make a command file
        replaceText(self.comFile, 'SEQUENCE', self.sequence)
        replaceText(self.comFile, 'TEMPERATURE', str(self.temperature - 273))  # probably not important, but we can add the temperature in C

        if self.foldSpeed == 'long':
            replaceText(self.comFile, 'INTERVAL', str(self.intervalLength))

        baseString = '#baseInteraction A IND WatsonCrick A IND2 WatsonCrick Cis'                     
        lineNum = findLine(self.comFile, baseString)  # this line defines the attractive forces in MMB

        for i in range(len(self.pairList)):
            filledString = 'baseInteraction A {} WatsonCrick A {} WatsonCrick Cis'.format(self.pairList[i, 0], self.pairList[i, 1])
            addLine(self.comFile, filledString, lineNum + 1)
    # def generateCommandFile(self):
    #     copyfile(self.template, self.comFile)  # make command file
    #     replaceText(self.comFile, 'SEQUENCE', self.sequence)
    #     replaceText(self.comFile, 'TEMPERATURE', str(self.temperature - 273)) # probably not important, but we can add the temperature in C
    #     if self.foldSpeed == 'long':
    #         replaceText(self.comFile, 'INTERVAL', str(self.intervalLength))

    #     baseString = '#baseInteraction A IND WatsonCrick A IND2 WatsonCrick Cis'
    #     lineNum = findLine(self.comFile, baseString)  # find the line number to start enumerating base pairs

    #     for i in range(len(self.pairList)):
    #         filledString = 'baseInteraction A {} WatsonCrick A {} WatsonCrick Cis'.format(self.pairList[i, 0], self.pairList[i, 1])
    #         addLine(self.comFile, filledString, lineNum + 1)

    def fold(self):
        # run fold - sometimes this errors for no known reasons - keep trying till it works
        result = None
        attempts = 0
        while (result is None) and (attempts < 100):
            try:
                attempts += 1
                os.system(self.mmbPath + ' -c ' + self.comFile + ' > outfiles/fold.out')
                os.replace('frame.pdb', self.foldedSequence)

                result = 1
                # clean up
                self.fileDump = 'mmbFiles_%d' % self.ind1
                os.mkdir(self.fileDump)
                os.system('mv last* ' + self.fileDump)
                os.system('mv trajectory.* ' + self.fileDump)
                os.system('mv watch.* ' + self.fileDump)
            except:  # TODO look up this warning: do not use bare except; Too broad except clause
                pass

    def check2DAgreement(self):
        """
        check agreement between prescribed 2D structure and the actual fold
        :return:
        """
        # # do it with MDA: MD_Analysis
        # u = mda.Universe(self.foldedSequence)
        # wcTraj = getWCDistTraj(u)
        # pairTraj = getPairTraj(wcTraj)
        # trueConfig = pairListToConfig(self.pairList, len(self.sequence))
        # foldDiscrepancy = getSecondaryStructureDistance([pairTraj[0], trueConfig])[0,1]
        # self.foldFidelity = 1-foldDiscrepancy

        # do it with MMB
        f = open('outfiles/fold.out')
        text = f.read()
        f.close()
        text = text.split('\n')
        scores = []
        for i in range(len(text)):
            if "Satisfied baseInteraction's" in text[i]:  # eg, Satisfied baseInteraction's:   : 8 out of : 13
                line = text[i].split(':')
                numerator = float(line[-2].split(' ')[1])
                denominator = float(line[-1])
                if denominator == 0:  # if there are no pairs, there are no constraints
                    scores.append(1)  # then ofc, any folded structure satisfies no constraints.
                else:
                    scores.append(numerator/denominator)

        scores = np.asarray(scores)
        self.foldFidelity = np.amax(scores)


# openmm
class omm:
    def __init__(self, structure, params, simTime=None, implicitSolvent=False):
        """
        pass on the pre-set and user-defined params to openmm engine
        """
        self.structureName = structure.split('.')[0]  # e.g., structure: relaxedSequence_0_amb_processed.pdb
        self.peptide = params['peptide']

        if implicitSolvent is False:
            self.pdb = PDBFile(structure)
            self.waterModel = params['water model']
            self.forcefield = ForceField('amber14-all.xml', 'amber14/' + self.waterModel + '.xml')

        # System configuration
        self.nonbondedMethod = params['nonbonded method']
        self.nonbondedCutoff = params['nonbonded cutoff'] * unit.nanometer
        self.ewaldErrorTolerance = params['ewald error tolerance']
        self.constraints = params['constraints']
        self.rigidWater = params['rigid water']
        self.constraintTolerance = params['constraint tolerance']
        self.hydrogenMass = params['hydrogen mass'] * unit.amu

        # Integration options
        self.dt = params['time step'] / 1000 * unit.picosecond

        self.temperature = params['temperature'] * unit.kelvin
        self.friction = params['friction'] / unit.picosecond
        self.quickSim = params['test mode']

        # Simulation options
        if simTime is not None:  # we can override the simulation time here
            self.steps = int(simTime * 1e6 // params['time step'])
            self.equilibrationSteps = int(params['equilibration time'] * 1e6 // params['time step'])
        else:
            self.steps = int(params['sampling time'] * 1e6 // params['time step'])  # number of steps
            self.equilibrationSteps = int(params['equilibration time'] * 1e6 // params['time step'])

        # Platform
        if params['platform'] == 'CUDA':  # 'CUDA' or 'cpu'
            self.platform = Platform.getPlatformByName('CUDA')
            if params['platform precision'] == 'single':
                self.platformProperties = {'Precision': 'single'}
            elif params['platform precision'] == 'double':
                self.platformProperties = {'Precision': 'double'}
        else:
            self.platform = Platform.getPlatformByName('CPU')

        self.reportSteps = int(params['print step'] * 1000 / params['time step'])  # report steps in ps, time step in fs
        self.dcdReporter = DCDReporter(self.structureName + '_trajectory.dcd', self.reportSteps)
        # self.pdbReporter = PDBReporter(self.structureName + '_trajectory.pdb', self.reportSteps)  # huge files
        self.dataReporter = StateDataReporter('log.txt', self.reportSteps, totalSteps=self.steps, step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, separator='\t')        
        self.checkpointReporter = CheckpointReporter(self.structureName + '_state.chk', 10000)

        # Prepare the simulation
        if implicitSolvent is False:
            # print_record('Building system...') # waste of space
            self.topology = self.pdb.topology
            self.positions = self.pdb.positions
            self.system = self.forcefield.createSystem(self.topology, nonbondedMethod=self.nonbondedMethod, nonbondedCutoff=self.nonbondedCutoff, constraints=self.constraints, rigidWater=self.rigidWater, ewaldErrorTolerance=self.ewaldErrorTolerance, hydrogenMass=self.hydrogenMass)
            for atom in self.topology.atoms():
                if atom.residue.name == 'Y' or atom.residue.name == 'TYR':
                    print_record("The first amino acid of the peptide (TYR) belongs to chain ID = " + str(atom.residue.chain.index))
            # TODO why looking for the TYR? covid peptide residue?

        else:  # create a system using prmtop file
            self.solventModel = params['implicit solvent model']
            print('\nCreating a simulation system under implicit solvent model of {}'.format(self.solventModel))
            self.prmtop = AmberPrmtopFile(self.structureName + '.top')  # e.g. foldedSequence_amb_processed.top or relaxedSequence_0_amb_processed.top
            self.inpcrd = AmberInpcrdFile(self.structureName + '.crd')
            self.topology = self.prmtop.topology
            self.positions = self.inpcrd.positions

            self.system = self.prmtop.createSystem(implicitSolvent=self.solventModel, implicitSolventSaltConc=0.0 * (unit.moles / unit.liter), nonbondedCutoff=1 * unit.nanometer, constraints=HBonds)
            # self.simulation.context.setPositions(inpcrd.positions)
            print('Initial positions set.')

            if self.inpcrd.boxVectors is not None:
                self.simulation.context.setPeriodicBoxVectors(*self.inpcrd.boxVectors)

        # Apply constraint if specified so
        if params['peptide backbone constraint constant'] != 0:
            self.force = CustomTorsionForce(
                '0.5*K*dtheta^2; dtheta = min(diff, 2*' + str(round(pi, 3)) + '-diff); diff = abs(theta - theta0)')
            self.force.addGlobalParameter('K', params['peptide backbone constraint constant'])
            self.force.addPerTorsionParameter('theta0')

            self.angles_to_constrain = findAngles()

            self.phi_tup = ('C', 'N', 'CA', 'C')
            self.psi_tup = ('N', 'CA', 'C', 'N')
            self.angle_tups = self.phi_tup, self.psi_tup

            radians = unit.radians

            rad_conv = pi / 180
            self.nchains = len(self.topology._chains)

            print_record("Number of chains = " + str(self.nchains))
            print_record("Beginning to iterate through chains...\n")

            for row in self.angles_to_constrain:
                aa_id, phi, psi, chain_id = row[0], row[1], row[2], row[3]
                aa_id, phi, psi, chain_id = int(aa_id), float(phi), float(psi), int(chain_id)

                print_record(f"aa_id = {aa_id}, phi = {phi}, psi = {psi}, chain_id = {chain_id}")
                print_record("Printing first 5 atoms in topology.atoms()...")

                # first5 = [atom for atom in self.topology.atoms()]
                # first5 = first5[:5]

                # for atom in first5:
                # print_record(f"Atom name={atom.name}, Atom residue chain index = {atom.residue.chain.index}, Atom residue index = {atom.residue.index}")

                self.da_atoms = [atom for atom in self.topology.atoms() if
                                 atom.residue.chain.index == chain_id and atom.name in {'N', 'CA',
                                                                                        'C'} and atom.residue.index in {
                                     aa_id, aa_id - 1}]

                print_record("Identified da_atoms.\n")

                # aa_id - 1 is included to account for the atoms in the previous residue being part of the current residue's dihedrals

                print_record("da_atoms length = " + str(len(self.da_atoms)))  # returns 0 for some reason

                for i in range(len(self.da_atoms) - 3):
                    self.tup = tuple([atom.name for atom in self.da_atoms[i:i + 4]])
                    self.tupIndex = tuple([atom.index for atom in self.da_atoms[i:i + 4]])
                    print_record("Found tup and tupIndex.\n")

                    if self.da_atoms[i + 3].residue.index == aa_id:
                        print_record("Beginning to add torsions...\n")

                        if self.tup == self.phi_tup:
                            self.force.addTorsion(self.tupIndex[0],
                                                  self.tupIndex[1],
                                                  self.tupIndex[2],
                                                  self.tupIndex[3], (phi * rad_conv,) * radians)
                            print_record("Successfully added a phi torsion restraint.\n")

                        elif self.tup == self.psi_tup:
                            self.force.addTorsion(self.tupIndex[0],
                                                  self.tupIndex[1],
                                                  self.tupIndex[2],
                                                  self.tupIndex[3], (psi * rad_conv,) * radians)
                            print_record("Successfully added a phi torsion restraint.\n")

            self.system.addForce(self.force)
            print_record("Successfully added the force.\n")
        # done with constraint

        self.integrator = LangevinMiddleIntegrator(self.temperature, self.friction, self.dt)  # another object
        self.integrator.setConstraintTolerance(self.constraintTolerance)  # What is this tolerance for? For constraint?
        if params['platform'] == 'CUDA':
            self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform, self.platformProperties)
        elif params['platform'] == 'CPU':
            self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform)
        self.simulation.context.setPositions(self.positions)

        print_record("Positions set.\n")

    def doMD(self):  # no need to be aware of the implicitSolvent
        if not os.path.exists(self.structureName + '_state.chk'):
            # Minimize and Equilibrate
            print_record('Performing energy minimization...')
            if self.quickSim:  # if want to do it fast, loosen the tolerances
                self.simulation.minimizeEnergy(tolerance=20, maxIterations=100)  # default is 10 kJ/mol - also set a max number of iterations
            else:
                self.simulation.minimizeEnergy()  # (tolerance = 1 * unit.kilojoules / unit.mole)
            print_record('Equilibrating...')
            self.simulation.context.setVelocitiesToTemperature(self.temperature)
            self.simulation.step(self.equilibrationSteps)
        else:  # no checkpoint file
            self.simulation.loadCheckpoint(self.structureName + '_state.chk')

        # Simulation
        print_record('Simulating...')
        self.simulation.reporters.append(self.dcdReporter)
        # self.simulation.reporters.append(self.pdbReporter)
        self.simulation.reporters.append(self.dataReporter)
        self.simulation.reporters.append(self.checkpointReporter)
        self.simulation.currentStep = 0
        with Timer() as md_time:
            self.simulation.step(self.steps)  # run the dynamics

        self.simulation.saveCheckpoint(self.structureName + '_state.chk')

        self.ns_per_day = (self.steps * self.dt) / (md_time.interval * unit.seconds) / (unit.nanoseconds / unit.day)

        return self.ns_per_day


class ld:  # lightdock
    def __init__(self, aptamer, peptide, params, ind1):
        self.setupPath = params['ld setup path']  # ld_scripts/... there should be a directory ld/scripts in workdir
        self.runPath = params['ld run path']
        self.genPath = params['lgd generate path']
        self.clusterPath = params['lgd cluster path']
        self.rankPath = params['lgd rank path']
        self.topPath = params['lgd top path']

        self.aptamerPDB = aptamer
        self.peptidePDB = peptide

        self.glowWorms = 300  # params['glowworms']
        self.dockingSteps = params['docking steps']
        self.numTopStructures = params['N docked structures']

        self.pH = params['pH']

        self.ind1 = ind1

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

    def getSwarmCount(self):
        # identify aptamer size
        dimensions = getMoleculeSize(self.aptamerPDB)
        # compute surface area of rectangular prism with no offset
        surfaceArea = 2 * (dimensions[0] * dimensions[1] + dimensions[1] * dimensions[2] + dimensions[2] * dimensions[1])  # in angstroms
        # compute number of swarms which cover this surface area - swarms are spheres 2 nm in diameter - this is only approximately the right number of swarms, since have no offset, and are using a rectangle
        nSwarms = np.ceil(surfaceArea / (np.pi * 10 ** 2))

        self.swarms = int(nSwarms)  # number of glowworm swarms

    def prepPDBs(self):  # different from prepPDB in utils.py
        killH(self.aptamerPDB)  # DNA needs to be deprotonated on the phosphate groups
        addH(self.peptidePDB, self.pH)  # peptide needs to be hydrogenated: side chains or terminal amino and carbonate groups?
        self.aptamerPDB2 = self.aptamerPDB.split('.')[0] + "_noH.pdb"
        self.peptidePDB2 = self.peptidePDB.split('.')[0] + "_H.pdb"
        changeSegment(self.peptidePDB2, 'A', 'B')

    def runLightDock(self):
        # Run setup
        os.system(self.setupPath + ' ' + self.aptamerPDB2 + ' ' + self.peptidePDB2 + ' -s ' + str(self.swarms) + ' -g ' + str(self.glowWorms) + ' >> outfiles/lightdockSetup.out')

        # Run docking
        os.system(self.runPath + ' setup.json ' + str(self.dockingSteps) + ' -s dna >> outfiles/lightdockRun.out')

    def generateAndCluster(self):
        # Generate docked structures and cluster them
        for i in range(self.swarms):
            os.chdir('swarm_%d' % i)
            os.system(self.genPath + ' ../' + self.aptamerPDB2 + ' ../' + self.peptidePDB2 + ' gso_%d' % self.dockingSteps + '.out' + ' %d' % self.glowWorms + ' > /dev/null 2> /dev/null; >> generate_lightdock.list')  # generate configurations
            os.system(self.clusterPath + ' gso_%d' % self.dockingSteps + '.out >> cluster_lightdock.list')  # cluster glowworms
            os.chdir('../')

    def rank(self):
        # Rank the clustered docking setups
        os.system(self.rankPath + ' %d' % self.swarms + ' %d' % self.dockingSteps + ' >> outfiles/ld_rank.out')

    def extractTopStructures(self):
        # Generate top structures
        os.system(self.topPath + ' ' + self.aptamerPDB2 + ' ' + self.peptidePDB2 + ' rank_by_scoring.list %d' % self.numTopStructures + ' >> outfiles/ld_top.out')
        os.mkdir('top_%d'%self.ind1)
        os.system('mv top*.pdb top_%d'%self.ind1 + '/')  # collect top structures (clustering currently dubiously working)

    def extractTopScores(self):
        self.topScores = readInitialLines('rank_by_scoring.list', self.numTopStructures + 1)[1:]
        for i in range(len(self.topScores)):
            self.topScores[i] = float(self.topScores[i].split(' ')[-1])  # get the last number, which is the score

        return self.topScores
