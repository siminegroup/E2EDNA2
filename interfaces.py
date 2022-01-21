"""
interfaces for other major software used by OpenDNA
MacroMolecule Builder (MMB)
PeptideBuilder: installed by pip
OpenMM (omm)
LightDock (ld)
NUPACK
"""

from shutil import copyfile
from numpy import pi
from nupack import *

from openmm import *
from openmm.app import *
import openmm.unit as unit

from utils import *
from analysisTools import *

# installed openmm: http://docs.openmm.org/7.5.0/userguide/application.html#installing-openmm
# openmm is installed as an individual package;
# package simtk is created with packages "openmm" and "unit" inside but literally import the stand-alone openmm.


class nupack:
    def __init__(self, aptamerSeq, temperature, ionicStrength, mgConc=0):
        self.aptamerSeq = aptamerSeq
        self.temperature = temperature
        self.naConc = ionicStrength
        self.mgConc = mgConc
        self.R = 0.0019872  # ideal gas constant in kcal/mol/K

    def run(self):
        self.energyCalc()
        self.structureAnalysis()
        return self.ssDict

    def energyCalc(self):
        """
        aptamerSeq is in FASTA format
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

        A = Strand(self.aptamerSeq, name='A')
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
        self.ssDict['config'] = np.zeros((nStructures, len(self.aptamerSeq)))
        self.ssDict['num pairs'] = np.zeros(nStructures).astype(int)
        self.ssDict['pair frac'] = np.zeros(nStructures)
        self.ssDict['num pins'] = np.zeros(nStructures).astype(int)
        self.ssDict['state prob'] = np.zeros(nStructures)

        for i in range(nStructures):
            ssString = str(self.output.subopt[i].structure)
            self.ssDict['2d string'].append(ssString)
            self.ssDict['pair list'].append(ssToList(ssString))
            self.ssDict['config'][i] = pairListToConfig(self.ssDict['pair list'][-1], len(self.aptamerSeq))
            self.ssDict['state prob'][i] = np.exp(-self.output.subopt[i].energy / self.R / self.temperature / float(self.output.pfunc))
            self.ssDict['num pairs'][i] = ssString.count('(')  # number of pairs
            self.ssDict['pair frac'][i] = 2 * self.ssDict['num pairs'][i] / len(self.aptamerSeq)

            nPins = 0  # number of distinct hairpins
            indA = 0
            for j in range(len(self.aptamerSeq)):
                if ssString[j] == '(':
                    indA += 1
                elif ssString[j] == ')':
                    indA -= 1
                    if indA == 0:  # if we come to the end of a distinct hairpin
                        nPins += 1

            self.ssDict['num pins'][i] = nPins


class mmb:  # MacroMolecule Builder (MMB)
    def __init__(self, aptamerSeq, pairList, params, ind1, intervalLength=None):
        if params['fold speed'] == 'quick':
            self.template = 'commands.template_quick.dat'  # can be used for debugging runs - very short
        elif params['fold speed'] == 'normal':
            self.template = 'commands.template.dat'        # default folding algorithm
        elif params['fold speed'] == 'long':
            self.template = 'commands.template_long.dat'   # extended annealing - for difficult sequences
        self.comFile = 'commands.run_fold.dat'
        self.foldSpeed = params['fold speed']
        self.aptamerSeq = aptamerSeq
        self.temperature = params['temperature']
        self.pairList = pairList
        self.mmbPath = params['mmb']
        
        # Conditions used to decide how to run MMB executable
        self.device = params['device']  # 'local' or 'cluster'
        self.devicePlatform = params['device platform']  # 'macos' or 'linux' or 'WSL'
        self.mmbDylibPath = params['mmb dir']
        
        self.ind1 = ind1
        self.foldedAptamerSeq = 'foldedAptamer_{}.pdb'.format(ind1)  # output structures of MMB
        self.intervalLength = intervalLength
        self.fileDump = 'mmbFiles_%d' % self.ind1  # directory to save the mmb run files

        self.foldFidelity = 0

    def run(self):
        self.generateCommandFile()
        self.fold()
        self.check2DAgreement()

        return self.MDAfoldFidelity, self.MMBfoldFidelity

    def generateCommandFile(self):
        copyfile(self.template, self.comFile)  # make a command file
        replaceText(self.comFile, 'SEQUENCE', self.aptamerSeq)
        replaceText(self.comFile, 'TEMPERATURE', str(self.temperature - 273))  # probably not important, but we can add the temperature in C

        if self.foldSpeed == 'long':
            replaceText(self.comFile, 'INTERVAL', str(self.intervalLength))

        baseString = '#baseInteraction A IND WatsonCrick A IND2 WatsonCrick Cis'                     
        lineNum = findLine(self.comFile, baseString)  # this line defines the attractive forces in MMB

        for i in range(len(self.pairList)):
            filledString = 'baseInteraction A {} WatsonCrick A {} WatsonCrick Cis'.format(self.pairList[i, 0], self.pairList[i, 1])
            addLine(self.comFile, filledString, lineNum + 1)

    def fold(self):
        # run fold - sometimes this errors for no known reasons - keep trying till it works
        result = None
        attempts = 0
        while (result is None) and (attempts < 100):
            try:
                attempts += 1                
                if (self.device == 'local') and (self.devicePlatform == 'macos'):  # special care for macos. Assume no cluster is on macos
                    os.system('export DYLD_LIBRARY_PATH=' + self.mmbDylibPath + ';' + self.mmbPath + ' -c ' + self.comFile + ' > outfiles/fold.out')    
                else:  # linux or WSL: "LD_LIBRARY_PATH" is specified in opendna.py ('local') or sub.sh ('cluster')
                    os.system(self.mmbPath + ' -c ' + self.comFile + ' > outfiles/fold.out')  # should we append new outputs to the fold.out?
                os.replace('frame.pdb', self.foldedAptamerSeq)
                result = 1

                # clean up
                self.fileDump = 'mmbFiles_%d' % self.ind1
                if os.path.isdir('./' + self.fileDump):  # this is refolding the aptamer
                    self.fileDump = self.fileDump + '/refolding'
                    if os.path.isdir('./' + self.fileDump):
                        pass
                    else:
                        os.mkdir(self.fileDump)
                else:
                    os.mkdir(self.fileDump)
                os.system('mv last* ' + self.fileDump)
                os.system('mv trajectory.* ' + self.fileDump)                
                # os.system('mv watch.* ' + self.fileDump)
                os.system('mv match.4*.pdb ' + self.fileDump)  # the intermediate files are updated during each stage
            except:  # TODO look up this warning: "do not use bare except; Too broad except clause"
                pass

    def check2DAgreement(self):
        """
        check agreement between prescribed 2D structure and the actual fold
        :return:
        """
        # do it with MDA: MDAnalysis
        u = mda.Universe(self.foldedAptamerSeq)
        # extract distance info through the trajectory
        wcTraj = getWCDistTraj(u)  # watson-crick base pairing distances (H-bonding)
        pairTraj = getPairTraj(wcTraj)        
        
        # 2D structure analysis        
        secondaryStructure = analyzeSecondaryStructure(pairTraj)  # find equilibrium secondary structure
        printRecord('2D structure after MMB folding (from MDAnalysis): ' + configToString(secondaryStructure))

        # MDAnalysis: fidelity score
        trueConfig = pairListToConfig(self.pairList, len(self.aptamerSeq))
        foldDiscrepancy = getSecondaryStructureDistance([pairTraj[0], trueConfig])[0,1]
        self.MDAfoldFidelity = 1-foldDiscrepancy

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
        self.MMBfoldFidelity = np.amax(scores)  # the score would increase with more and more stages


# openmm
class omm:
    def __init__(self, structurePDB, params, simTime=None, binding=False, implicitSolvent=False):
        """
        pass on the pre-set and user-defined params to openmm engine
        """
        self.structureName = structurePDB.split('.')[0]  # e.g., structurePDB: relaxedAptamer_0_amb_processed.pdb or complex_1_2_processed.pdb        
        self.chkFile = params['chk file']  # if not resuming the simulation, it is empty string ""

        if implicitSolvent is False:
            self.pdb = PDBFile(structurePDB)
            self.waterModel = params['water model']
            self.forcefield = ForceField('amber14-all.xml', 'amber14/' + self.waterModel + '.xml')

        # System configuration
        self.nonbondedMethod = params['nonbonded method']  # Currently PME!!! It's used with periodic boundary condition applied
        self.nonbondedCutoff = params['nonbonded cutoff'] * unit.nanometer
        self.ewaldErrorTolerance = params['ewald error tolerance']
        self.constraints = params['constraints']
        self.rigidWater = params['rigid water']
        self.constraintTolerance = params['constraint tolerance']
        self.hydrogenMass = params['hydrogen mass'] * unit.amu
        self.temperature = params['temperature'] * unit.kelvin

        # Integration options
        self.dt = params['time step'] / 1000 * unit.picosecond
        
        self.quickSim = params['test mode']

        # Simulation options
        self.timeStep = params['time step']  # fs
        self.equilibrationSteps = int(params['equilibration time'] * 1e6 // params['time step'])
        if simTime is not None:  # override the default "sampling time": e.g., when running smoothing.
            self.steps = int(simTime * 1e6 // params['time step'])            
            self.simTime = simTime
        else:  # default sampling time = params['sampling time']
            self.steps = int(params['sampling time'] * 1e6 // params['time step'])  # number of steps
            self.simTime = params['sampling time']
        
        # Platform
        if params['platform'] == 'CUDA':  # 'CUDA' or 'cpu'
            self.platform = Platform.getPlatformByName('CUDA')
            if params['platform precision'] == 'single':
                self.platformProperties = {'Precision': 'single'}
            elif params['platform precision'] == 'double':
                self.platformProperties = {'Precision': 'double'}
        else:
            self.platform = Platform.getPlatformByName('CPU')

        # Can resumed run append the old log.txt or .dcd file?
        self.reportSteps = int(params['print step'] * 1000 / params['time step'])  # report steps in ps, time step in fs
        self.dcdReporter = DCDReporter(self.structureName + '_trajectory.dcd', self.reportSteps)
        # self.pdbReporter = PDBReporter(self.structureName + '_trajectory.pdb', self.reportSteps)  # huge files
        if simTime is None:
            if binding is True:
                logFileName = 'log_complex.txt'
            else:
                logFileName = 'log.txt'
        else:  # MD smoothing: simTime is specified as 'smoothing time'
            logFileName = 'log_smoothing.txt'
        self.dataReporter = StateDataReporter(logFileName, self.reportSteps, totalSteps=self.steps, step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True,separator='\t')        
                
        if (params['pick up from freeAptamerChk'] is False) and (params['pick up from complexChk'] is False):  # or if not self.chkFile:
            self.checkpointReporter = CheckpointReporter(self.structureName + '_state.chk', 10000)
        else:  # ie, we are resuming a sampling, chkFile must exist.
            self.checkpointReporter = CheckpointReporter(self.chkFile, 10000)                   
            
        # Prepare the simulation        
        if implicitSolvent is False:
            self.topology = self.pdb.topology
            self.positions = self.pdb.positions
            printRecord('Creating a simulation system under explicit solvent')

            self.system = self.forcefield.createSystem(self.topology, nonbondedMethod=self.nonbondedMethod, nonbondedCutoff=self.nonbondedCutoff, constraints=self.constraints, rigidWater=self.rigidWater, hydrogenMass=self.hydrogenMass, ewaldErrorTolerance=self.ewaldErrorTolerance)
            # ewaldErrorTolerance: as "**args": Arbitrary additional keyword arguments may also be specified. This allows extra parameters to be specified that are specific to particular force fields.

        else:  # create a system using prmtop file and use implicit solvent
            self.prmtop = AmberPrmtopFile(self.structureName + '.top')  # e.g. foldedAptamer_amb_processed.top or relaxedAptamer_0_amb_processed.top or complex_1_2_amb_processed.pdb
            self.inpcrd = AmberInpcrdFile(self.structureName + '.crd')
            self.topology = self.prmtop.topology
            self.positions = self.inpcrd.positions
            printRecord('Creating a simulation system under implicit solvent model of {}'.format(params['implicit solvent model']))
            
            self.implicitSolventModel = params['implicit solvent model']
            self.implicitSolventSaltConc = params['implicit solvent salt conc'] * (unit.moles / unit.liter)
            self.implicitSolventKappa = params['implicit solvent Kappa']  # add this in main_resume/main.py
            self.soluteDielectric = params['soluteDielectric']
            self.solventDielectric = params['solventDielectric']

            self.system = self.prmtop.createSystem(nonbondedMethod=self.nonbondedMethod, nonbondedCutoff=self.nonbondedCutoff, constraints=self.constraints, rigidWater=self.rigidWater,
                                                   implicitSolvent=self.implicitSolventModel, implicitSolventSaltConc=self.implicitSolventSaltConc, implicitSolventKappa=self.implicitSolventKappa, temperature=self.temperature,
                                                   soluteDielectric=self.soluteDielectric, solventDielectric=self.solventDielectric, removeCMMotion=True, hydrogenMass=self.hydrogenMass, ewaldErrorTolerance=self.ewaldErrorTolerance, switchDistance=0.0*unit.nanometer)
            '''
            temperature : 
                Temperture of the system, only used to compute the Debye length from implicitSolventSoltConc
            implicitSolventSaltConc : 
                Concentration of all cations (or anions): because it's to replace ionic strength, which = [cation] = [anion] for any monovalent electrolyte
                The salt concentration for GB calculations (modelled as a debye screening parameter). 
                It's converted to the debye length (kappa) using the provided temperature and solventDielectric
            implicitSolventKappa: 
                float unit of 1/length: 1/angstroms. No need to * unit.xxx; 
                If this value is set, implicitSolventSaltConc will be ignored. 
                If not set, it's calculated as follows in OpenMM:                     
                    implicitSolventKappa = 50.33355 * sqrt(implicitSolventSaltConc / solventDielectric / temperature)
                    # The constant is 1 / sqrt( epsilon_0 * kB / (2 * NA * q^2 * 1000) ),
                    # where NA is avogadro's number, epsilon_0 is the permittivity of free space, q is the elementary charge (this number matches Amber's kappa conversion factor)
                    # The equation above is for Debye screening parameter, Kappa, for a monovalent electrolyte in an electrolyte or a colloidal suspension. Debye length = 1/Kappa
                    # Then Multiply by 0.73 to account for ion exclusions, and multiply by 10 to convert to 1/nm from 1/angstroms
                    implicitSolventKappa *= 7.3
            soluteDielectric :
                The solute dielectric constant to use in the implicit solvent model.
                Default value = 1.0, meaning no screening effect at all.
            solventDielectric :
                The solvent dielectric constant to use in the implicit solvent model.
                Default value = 78.5, which is the value for water.
                This offers a way to model non-aqueous system.
            rigidWater: default = True
                If true, water molecules will be fully rigid regardless of the value passed for the constraints argument.
                Even using implicit solvent models, the simulation system could still contain water molecules, but just not from the solvent.
            switchDistance: 
                The distance at which the potential energy switching function is turned on for Lennard-Jones interactions. 
                If the switchDistance is 0 or evaluates to boolean False, no switching function will be used. 
                Values greater than nonbondedCutoff or less than 0 raise ValueError
            '''
            if self.inpcrd.boxVectors is not None:
                self.simulation.context.setPeriodicBoxVectors(*self.inpcrd.boxVectors)

        # Apply constraint if specified so
        if params['peptide backbone constraint constant'] != 0:
            self.force = CustomTorsionForce('0.5*K*dtheta^2; dtheta = min(diff, 2*' + str(round(pi, 3)) + '-diff); diff = abs(theta - theta0)')
            self.force.addGlobalParameter('K', params['peptide backbone constraint constant'])
            self.force.addPerTorsionParameter('theta0')

            self.angles_to_constrain = findAngles()

            self.phi_tup = ('C', 'N', 'CA', 'C')
            self.psi_tup = ('N', 'CA', 'C', 'N')
            self.angle_tups = self.phi_tup, self.psi_tup

            radians = unit.radians

            rad_conv = pi / 180
            self.nchains = len(self.topology._chains)

            printRecord("Number of chains = " + str(self.nchains))
            printRecord("Beginning to iterate through chains...\n")

            for row in self.angles_to_constrain:
                aa_id, phi, psi, chain_id = row[0], row[1], row[2], row[3]
                aa_id, phi, psi, chain_id = int(aa_id), float(phi), float(psi), int(chain_id)

                printRecord(f"aa_id = {aa_id}, phi = {phi}, psi = {psi}, chain_id = {chain_id}")
                printRecord("Printing first 5 atoms in topology.atoms()...")

                # first5 = [atom for atom in self.topology.atoms()]
                # first5 = first5[:5]

                # for atom in first5:
                # printRecord(f"Atom name={atom.name}, Atom residue chain index = {atom.residue.chain.index}, Atom residue index = {atom.residue.index}")

                self.da_atoms = [atom for atom in self.topology.atoms() if
                                 atom.residue.chain.index == chain_id and atom.name in {'N', 'CA',
                                                                                        'C'} and atom.residue.index in {
                                     aa_id, aa_id - 1}]

                printRecord("Identified da_atoms.\n")

                # aa_id - 1 is included to account for the atoms in the previous residue being part of the current residue's dihedrals

                printRecord("da_atoms length = " + str(len(self.da_atoms)))  # returns 0 for some reason

                for i in range(len(self.da_atoms) - 3):
                    self.tup = tuple([atom.name for atom in self.da_atoms[i:i + 4]])
                    self.tupIndex = tuple([atom.index for atom in self.da_atoms[i:i + 4]])
                    printRecord("Found tup and tupIndex.\n")

                    if self.da_atoms[i + 3].residue.index == aa_id:
                        printRecord("Beginning to add torsions...\n")

                        if self.tup == self.phi_tup:
                            self.force.addTorsion(self.tupIndex[0],
                                                  self.tupIndex[1],
                                                  self.tupIndex[2],
                                                  self.tupIndex[3], (phi * rad_conv,) * radians)
                            printRecord("Successfully added a phi torsion restraint.\n")

                        elif self.tup == self.psi_tup:
                            self.force.addTorsion(self.tupIndex[0],
                                                  self.tupIndex[1],
                                                  self.tupIndex[2],
                                                  self.tupIndex[3], (psi * rad_conv,) * radians)
                            printRecord("Successfully added a phi torsion restraint.\n")

            self.system.addForce(self.force)
            printRecord("Successfully added the force.\n")
        # done with constraint

        # # Need to create a Simulation class, even to load checkpoint file
        # self.integrator = LangevinMiddleIntegrator(self.temperature, self.friction, self.dt)  # another object
        # print("Using LangevinMiddleIntegrator integrator, at T={}.".format(self.temperature))
        self.friction = params['friction'] / unit.picosecond
        self.integrator = LangevinIntegrator(self.temperature, self.friction, self.dt)  # another object
        printRecord("Using LangevinIntegrator integrator, at T={}.".format(self.temperature))
        
        self.integrator.setConstraintTolerance(self.constraintTolerance)  # What is this tolerance for? For constraint?
        
        if params['platform'] == 'CUDA':
            self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform, self.platformProperties)
        elif params['platform'] == 'CPU':
            self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform)
        
        if (params['pick up from freeAptamerChk'] is False) and (params['pick up from complexChk'] is False):
            self.simulation.context.setPositions(self.positions)
            printRecord("Initial positions set.")
        else:
            pass  # if resuming a run, the initial position comes from the chk file.        

    def doMD(self):  # no need to be aware of the implicitSolvent        
        # if not os.path.exists(self.structureName + '_state.chk'):
        if not self.chkFile:  # "self.chkFile is not empty" is equivalent to "either params['pick up from freeAptamerChk'] or params['pick up from complexChk'] is True"
            # User did not specify a .chk file ==> we are doing a fresh sampling, not resuming.
            # Minimize and Equilibrate
            printRecord('Performing energy minimization...')
            if self.quickSim:  # if want to do it fast, loosen the tolerances
                self.simulation.minimizeEnergy(tolerance=20, maxIterations=100)  # default is 10 kJ/mol - also set a max number of iterations
            else:
                self.simulation.minimizeEnergy()  # (tolerance = 1 * unit.kilojoules / unit.mole)
            printRecord('Equilibrating({} steps, time step={} fs)...'.format(self.equilibrationSteps, self.timeStep))            
            self.simulation.context.setVelocitiesToTemperature(self.temperature)
            self.simulation.step(self.equilibrationSteps)
        else:
            # Resume a sampling: no need to minimize and equilibrate
            printRecord("Loading checkpoint file: " + self.chkFile + " to resume sampling")            
            self.simulation.loadCheckpoint(self.chkFile)  # Resume the simulation

        # Simulation
        printRecord('Simulating({} steps, time step={} fs, simTime={} ns)...'.format(self.steps, self.timeStep, self.simTime))
        self.simulation.reporters.append(self.dcdReporter)
        # self.simulation.reporters.append(self.pdbReporter)
        self.simulation.reporters.append(self.dataReporter)
        self.simulation.reporters.append(self.checkpointReporter)
        self.simulation.currentStep = 0
        with Timer() as md_time:
            self.simulation.step(self.steps)  # run the dynamics

        # Update the chk file with info from the final step
        if not self.chkFile:
            # User did not specify a .chk file ==> we are doing a fresh sampling, not resuming.
            self.simulation.saveCheckpoint(self.structureName + '_state.chk')
        else:
            self.simulation.saveCheckpoint(self.chkFile)

        self.ns_per_day = (self.steps * self.dt) / (md_time.interval * unit.seconds) / (unit.nanoseconds / unit.day)
    
        return self.ns_per_day

    def extractLastFrame(self, lastFrameFileName):
        lastpositions = self.simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(self.topology, lastpositions, open(lastFrameFileName, 'w'))
        printRecord('OpenMM: save the last frame into: {}'.format(lastFrameFileName))


class ld:  # lightdock
    def __init__(self, aptamerPDB, targetPDB, params, ind1):
        self.setupPath = params['ld setup path']  # ld_scripts/... there should be a directory ld/scripts in workdir
        self.runPath = params['ld run path']
        self.genPath = params['lgd generate path']
        self.clusterPath = params['lgd cluster path']
        self.rankPath = params['lgd rank path']
        self.topPath = params['lgd top path']

        self.aptamerPDB = aptamerPDB
        self.targetPDB = targetPDB

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
        addH(self.targetPDB, self.pH)  # peptide needs to be hydrogenated: side chains or terminal amino and carbonate groups?
        self.aptamerPDB2 = self.aptamerPDB.split('.')[0] + "_noH.pdb"
        self.targetPDB2 = self.targetPDB.split('.')[0] + "_H.pdb"
        changeSegment(self.targetPDB2, 'A', 'B')

    def runLightDock(self):
        # Run setup
        os.system(self.setupPath + ' ' + self.aptamerPDB2 + ' ' + self.targetPDB2 + ' -s ' + str(self.swarms) + ' -g ' + str(self.glowWorms) + ' >> outfiles/lightdockSetup.out')

        # Run docking
        os.system(self.runPath + ' setup.json ' + str(self.dockingSteps) + ' -s dna >> outfiles/lightdockRun.out')

    def generateAndCluster(self):
        # Generate docked structures and cluster them
        for i in range(self.swarms):
            os.chdir('swarm_%d' % i)
            os.system(self.genPath + ' ../' + self.aptamerPDB2 + ' ../' + self.targetPDB2 + ' gso_%d' % self.dockingSteps + '.out' + ' %d' % self.glowWorms + ' > /dev/null 2> /dev/null; >> generate_lightdock.list')  # generate configurations
            os.system(self.clusterPath + ' gso_%d' % self.dockingSteps + '.out >> cluster_lightdock.list')  # cluster glowworms
            os.chdir('../')

    def rank(self):
        # Rank the clustered docking setups
        os.system(self.rankPath + ' %d' % self.swarms + ' %d' % self.dockingSteps + ' >> outfiles/ld_rank.out')

    def extractTopStructures(self):
        # Generate top structures
        os.system(self.topPath + ' ' + self.aptamerPDB2 + ' ' + self.targetPDB2 + ' rank_by_scoring.list %d' % self.numTopStructures + ' >> outfiles/ld_top.out')
        os.mkdir('top_%d'%self.ind1)
        os.system('mv top*.pdb top_%d'%self.ind1 + '/')  # collect top structures (clustering currently dubiously working)
        # If no top*.pdb at all: there will be an warning message: "mv: rename top*.pdb to top_0/top*.pdb: No such file or directory". Will not stop the pipeline though.

    def extractTopScores(self):
        self.topScores = []
        topScoresFromFile = readInitialLines('rank_by_scoring.list', self.numTopStructures + 1)[1:]  # always read the header but not record to self.topScores
        if len(topScoresFromFile) < self.numTopStructures:
            printRecord('Number of identified docked structures is less than what you are looking for!')
        for i in range(len(topScoresFromFile)):
            score = topScoresFromFile[i].split(' ')[-1]  # get the last number, which is the score
            if type(score) == float:  # there is indeed a score
                self.topScores.append(float(score))
        printRecord('Number of identified docked structures = {}'.format(len(self.topScores)))
        # return self.topScores
