'''
E2EDNA 2.0 - OpenMM Implementation of E2EDNA !

An automated pipeline for simulating DNA aptamers complexed with target ligands (peptide, DNA, RNA or small molecules).
    
Copyright (C) 2022 Michael Kilgour, Tao Liu and Lena Simine

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

"""
interfaces for other major software used by E2EDNA
MacroMolecule Builder (MMB)
PeptideBuilder: installed by pip
OpenMM (omm)
LightDock (ld)
NUPACK
"""
import sys
import shutil
from shutil import copyfile
from numpy import pi

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
        try:
            import nupack
            # from nupack import * #: does not work at this level
        except ImportError as err:
            raise ImportError("Required module 'nupack' not found. Please download from http://www.nupack.org/downloads and install in your Python virtual environment.") from err            
        
        # Run nupack-dependent code
        if self.temperature > 273:  # auto-detect Kelvins
            gap = 2 * self.R * self.temperature
            CelsiusTemprature = self.temperature - 273
        else:
            gap = 2 * self.R * (self.temperature + 273)  # convert to Kelvin fir kT
            CelsiusTemprature = self.temperature

        A = nupack.Strand(self.aptamerSeq, name='A')
        comp = nupack.Complex([A], name='AA')
        set1 = nupack.ComplexSet(strands=[A], complexes=nupack.SetSpec(max_size=1, include=[comp]))
        model1 = nupack.Model(material='dna', celsius=CelsiusTemprature, sodium=self.naConc, magnesium=self.mgConc)
        results = nupack.complex_analysis(set1, model=model1, compute=['pfunc', 'mfe', 'subopt', 'pairs'], options={'energy_gap': gap})
        self.output = results[comp]

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
        if params['fold_speed'] == 'quick':
            self.template = 'commands.template_quick.dat'  # can be used for debugging runs - very short
        elif params['fold_speed'] == 'normal':
            self.template = 'commands.template.dat'        # default folding algorithm
        elif params['fold_speed'] == 'slow':
            self.template = 'commands.template_long.dat'   # extended annealing - for difficult sequences
        self.comFile = 'commands.run_fold.dat'
        self.foldSpeed = params['fold_speed']
        self.aptamerSeq = aptamerSeq
        self.temperature = params['temperature']
        self.pairList = pairList
        self.mmbPath = params['mmb']
        
        # Conditions used to decide how to run MMB executable
        self.device = params['device']  # 'local' or 'cluster'
        self.devicePlatform = params['operating_system']  # 'macos' or 'linux' or 'WSL'
        self.mmbDylibPath = params['mmb_dir']
        
        self.ind1 = ind1
        self.foldedAptamerSeq = 'foldedAptamer_{}.pdb'.format(ind1)  # output structures of MMB
        self.intervalLength = intervalLength
        self.fileDump = 'mmb_temp_files_%d' % self.ind1  # directory to save the mmb run files

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

        if self.foldSpeed == 'slow':
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
                    os.system('export DYLD_LIBRARY_PATH=' + self.mmbDylibPath + ';' + self.mmbPath + ' -c ' + self.comFile + ' > outfiles_of_pipeline_modules/mmb_fold.out')    
                else:  # linux or WSL: "LD_LIBRARY_PATH" is specified in opendna.py ('local') or sub.sh ('cluster')
                    os.system(self.mmbPath + ' -c ' + self.comFile + ' > outfiles_of_pipeline_modules/mmb_fold.out')  # maybe append new outputs to the mmb_fold.out
                os.replace('frame.pdb', self.foldedAptamerSeq)
                result = 1

                # clean up
                self.fileDump = 'mmb_temp_files_%d' % self.ind1
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
        u = mda.Universe(self.foldedAptamerSeq, dt=1.0)  # manually set dt=1.0ps to avoid warning. 
        # MMB trajectory doesn't provide dt. dt is not used in analysis here either.

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
        f = open('outfiles_of_pipeline_modules/mmb_fold.out')
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
    def __init__(self, runfilesDir, structurePDB, params, simTime=None, binding=False, implicitSolvent=False):
        """
        pass on the pre-set and user-defined params to openmm engine
        runfilesDir: a directory to store runfiles generated during MD (aptamer_relaxation, aptamer_sampling or complex sampling) 
                    eg, "md_aptamer_relaxation_runfiles_0", "md_aptamer_sampling_runfiles_0", "md_complex_sampling_runfiles_0"
        """
        self.runfilesDir = runfilesDir
        self.structureName = structurePDB.split('.')[0]  # e.g., structurePDB: relaxedAptamer_0_amb_processed.pdb or complex_1_2_processed.pdb        
        self.chkFile = params['chk_file']  # if not resuming the simulation, it is Python None; otherwise, basename of the .chk file, no directory info

        if implicitSolvent is False:
            self.pdb = PDBFile(structurePDB)
            self.forceFieldFilename = params['force_field']
            self.waterModelFilename = params['water_model']
            self.forcefield = ForceField(self.forceFieldFilename + '.xml', self.waterModelFilename + '.xml')

        # System configuration
        self.nonbondedMethod = params['nonbonded_method']  # Currently PME!!! It's used with periodic boundary condition applied
        self.nonbondedCutoff = params['nonbonded_cutoff'] * unit.nanometer
        self.ewaldErrorTolerance = params['ewald_error_tolerance']
        self.constraints = params['constraints']
        self.rigidWater = params['rigid_water']
        self.constraintTolerance = params['constraint_tolerance']
        self.hydrogenMass = params['hydrogen_mass'] * unit.amu
        self.temperature = params['temperature'] * unit.kelvin

        # Integration options
        self.dt = params['time_step'] / 1000 * unit.picosecond
        
        self.quickSim = params['quick_check_mode']

        # Simulation options
        self.timeStep = params['time_step']  # fs
        self.equilibrationSteps = int(params['equilibration_time'] * 1e6 // params['time_step'])
        
        if simTime is None:
            raise ValueError("The argument of MD sampling time is required: 'simTime'")
        else:
            self.steps = int(simTime * 1e6 // params['time_step'])            
            self.simTime = simTime
        
        # Platform
        if params['platform'] == 'CUDA':  # 'CUDA' or 'cpu'
            self.platform = Platform.getPlatformByName('CUDA')
            if params['CUDA_precision'] == 'single':
                self.platformProperties = {'Precision': 'single'}
            elif params['CUDA_precision'] == 'double':
                self.platformProperties = {'Precision': 'double'}
        else:
            self.platform = Platform.getPlatformByName('CPU')

        # Can resumed run append the old log.txt or .dcd file?
        self.reportSteps = int(params['print_step'] * 1000 / params['time_step'])  # report steps in ps, time step in fs
        # self.dcdReporter = DCDReporter(self.structureName + '_trajectory.dcd', self.reportSteps)
        dcdTrajFileName = os.path.join(self.runfilesDir, self.structureName+'_trajectory.dcd')
        self.dcdReporter = DCDReporter(dcdTrajFileName, self.reportSteps)
        # self.pdbReporter = PDBReporter(os.path.join(self.runfilesDir, self.structureName+'_trajectory.pdb'), self.reportSteps)  # huge file

        if params['smoothing_time'] is None:
            if binding is True:
                mdLogFileName = os.path.join(self.runfilesDir, 'MDlog_complexSampling.txt')
            else:
                mdLogFileName = os.path.join(self.runfilesDir, 'MDlog_freeAptamerSampling.txt')
        else:
            mdLogFileName = os.path.join(self.runfilesDir, 'MDlog_freeAptamerSmoothing.txt')
        self.dataReporter = StateDataReporter(mdLogFileName, self.reportSteps, totalSteps=self.steps, step=True, speed=True, progress=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True,separator='\t')        
                
        if (params['pickup_from_freeAptamerChk'] is False) and (params['pickup_from_complexChk'] is False):
            self.newChkFileName = os.path.join(self.runfilesDir, self.structureName+'_state.chk')
            self.checkpointReporter = CheckpointReporter(self.newChkFileName, 10000)

        else:  # ie, we are resuming a sampling, chkFile must exist in the run{run_num} directory
            if self.chkFile is not None: # just for redundancy
                self.chkFile_in_runfilesDir = os.path.join(self.runfilesDir, self.chkFile)
                os.copyfile(self.chkFile, self.chkFile_in_runfilesDir)
                printRecord('Copying the checkpoint file (was copied to the run folder): {} to subfolder: {}...'.format(self.chkFile, self.runfilesDir))
                self.checkpointReporter = CheckpointReporter(self.chkFile_in_runfilesDir, 10000)
            
        # Prepare the simulation        
        if implicitSolvent is False:
            self.topology = self.pdb.topology
            self.positions = self.pdb.positions
            printRecord('Creating a simulation system using force field bundled with OpenMM: {} and {}'.format(self.forceFieldFilename, self.waterModelFilename))

            self.system = self.forcefield.createSystem(self.topology, nonbondedMethod=self.nonbondedMethod, nonbondedCutoff=self.nonbondedCutoff, constraints=self.constraints, rigidWater=self.rigidWater, hydrogenMass=self.hydrogenMass, ewaldErrorTolerance=self.ewaldErrorTolerance)
            # ewaldErrorTolerance: as "**args": Arbitrary additional keyword arguments may also be specified. This allows extra parameters to be specified that are specific to particular force fields.

        else:  # create a system using prmtop file and use implicit solvent
            self.prmtop = AmberPrmtopFile(self.structureName + '.top')  # e.g. foldedAptamer_amb_processed.top or relaxedAptamer_0_amb_processed.top or complex_1_2_amb_processed.pdb
            self.inpcrd = AmberInpcrdFile(self.structureName + '.crd')
            self.topology = self.prmtop.topology
            self.positions = self.inpcrd.positions
            printRecord('Creating a simulation system using implicit solvent model of {}'.format(params['implicit_solvent_model']))
            
            self.implicitSolventModel = params['implicit_solvent_model']
            self.implicitSolventSaltConc = params['implicit_solvent_salt_conc'] * (unit.moles / unit.liter)
            self.implicitSolventKappa = params['implicit_solvent_Kappa'] # default is None, unless specified by user
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

        # # Need to create a Simulation class, even to load checkpoint file
        # self.integrator = LangevinMiddleIntegrator(self.temperature, self.friction, self.dt)  # another object
        # printRecord("Using LangevinMiddleIntegrator integrator, at T={}.".format(self.temperature))
        self.friction = params['friction'] / unit.picosecond
        self.integrator = LangevinIntegrator(self.temperature, self.friction, self.dt)  # another object
        printRecord("Using LangevinIntegrator integrator, at T={}.".format(self.temperature))
        
        self.integrator.setConstraintTolerance(self.constraintTolerance)
        
        if params['platform'] == 'CUDA':
            self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform, self.platformProperties)
        elif params['platform'] == 'CPU':
            self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform)
        
        # if (params['pickup_from_freeAptamerChk'] is False) and (params['pickup_from_complexChk'] is False):
        if params['pickup'] is False:
            self.simulation.context.setPositions(self.positions)
            printRecord("Initial positions set.")
        else:
            printRecord("When resuming a run: initial positions come from the .chk file.")

    def doMD(self):  # no need to be aware of the implicitSolvent        
        # if not os.path.exists(self.structureName + '_state.chk'):
        if self.chkFile is None:
            # User did not specify a .chk_file ==> we are doing a fresh sampling, not resuming.
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
            printRecord("Loading checkpoint file: " + self.chkFile_in_runfilesDir + " to resume sampling")            
            self.simulation.loadCheckpoint(self.chkFile_in_runfilesDir)  # Resume the simulation

        # Simulation
        printRecord('Simulating({} steps, time step={} fs, simTime={} ns)...'.format(self.steps, self.timeStep, self.simTime))
        self.simulation.reporters.append(self.dcdReporter)
        # self.simulation.reporters.append(self.pdbReporter)
        self.simulation.reporters.append(self.dataReporter)
        self.simulation.reporters.append(self.checkpointReporter)
        self.simulation.currentStep = 0
        with Timer() as md_time:
            self.simulation.step(self.steps)  # run the dynamics

        # Update the chk_file with info from the final step
        if self.chkFile is None:
            # User did not specify a .chk_file ==> we are doing a fresh sampling, not resuming.
            # self.simulation.saveCheckpoint(self.structureName + '_state.chk')
            self.simulation.saveCheckpoint(self.newChkFileName)
        else:
            self.simulation.saveCheckpoint(self.chkFile_in_runfilesDir)

        self.ns_per_day = (self.steps * self.dt) / (md_time.interval * unit.seconds) / (unit.nanoseconds / unit.day)
    
        return self.ns_per_day

    def extractLastFrame(self, lastFrameFileName):
        lastpositions = self.simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(self.topology, lastpositions, open(lastFrameFileName, 'w'))
        # printRecord('OpenMM: saved the last frame into: {}'.format(lastFrameFileName))


class ld:  # lightdock
    def __init__(self, aptamerPDB, targetPDB, params, ind1):
        self.setupScript   = params['ld_setup']
        self.runScript     = params['ld_run']
        self.genScript     = params['lgd_generate']
        self.clusterScript = params['lgd_cluster']
        self.rankScript    = params['lgd_rank']
        self.topScript     = params['lgd_top']

        self.aptamerPDB = aptamerPDB
        self.targetPDB = targetPDB

        self.glowWorms = 300  # params['glowworms']
        self.dockingSteps = params['docking_steps']
        self.numTopStructures = params['N_docked_structures']

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
        dir = 'lightdock_temp_files_{}'.format(self.ind1)
        os.mkdir(dir)
        os.system('mv swarm* ' + dir)
        os.system('mv setup.json ' + dir)
        os.system('mv *.list ' + dir)
        os.system('mv lightdock.info ' + dir)
        os.system('mv init ' + dir)
        os.system('mv lightdock_*.npy ' + dir)
        os.system('mv lightdock_*.pdb ' + dir)
        shutil.move(self.aptamerPDB2, dir)
        shutil.move(self.targetPDB2, dir)

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
        os.system(self.setupScript + ' ' + self.aptamerPDB2 + ' ' + self.targetPDB2 + ' -s ' + str(self.swarms) + ' -g ' + str(self.glowWorms) + ' >> outfiles_of_pipeline_modules/lightdock_setup.out')

        # Run docking
        os.system(self.runScript + ' setup.json ' + str(self.dockingSteps) + ' -s dna >> outfiles_of_pipeline_modules/lightdock_run.out')

    def generateAndCluster(self):
        # Generate docked structures and cluster them
        for i in range(self.swarms):
            os.chdir('swarm_%d' % i)
            os.system(self.genScript + ' ../' + self.aptamerPDB2 + ' ../' + self.targetPDB2 + ' gso_%d' % self.dockingSteps + '.out' + ' %d' % self.glowWorms + ' > /dev/null 2> /dev/null; >> generate_lightdock.list')  # generate configurations
            os.system(self.clusterScript + ' gso_%d' % self.dockingSteps + '.out >> cluster_lightdock.list')  # cluster glowworms
            os.chdir('../')

    def rank(self):
        # Rank the clustered docking setups
        os.system(self.rankScript + ' %d' % self.swarms + ' %d' % self.dockingSteps + ' >> outfiles_of_pipeline_modules/lightdock_rank.out')

    def extractTopStructures(self):
        # Generate top structures
        os.system(self.topScript + ' ' + self.aptamerPDB2 + ' ' + self.targetPDB2 + ' rank_by_scoring.list %d' % self.numTopStructures + ' >> outfiles_of_pipeline_modules/lightdock_top.out')
        os.mkdir('lightdock_topscores_%d'%self.ind1)
        if os.path.exists('./top_1.pdb'):  # if there is any docked structure
            os.system('mv top*.pdb lightdock_topscores_%d'%self.ind1 + '/')  # collect top structures (clustering currently dubiously working)        
        else:
            printRecord('Not found good docked structures with top scores!')
            # sys.exit()

    def extractTopScores(self):
        self.topScores = []
        topScoresFromFile = readInitialLines('rank_by_scoring.list', self.numTopStructures + 1)[1:]  # always read the header but not record to self.topScores
        if len(topScoresFromFile) < self.numTopStructures:
            printRecord('Number of identified docked structures is less than what you are looking for!')
        for i in range(len(topScoresFromFile)):
            scoreStr = topScoresFromFile[i].split(' ')[-1]  # get the last number, which is the score in string format
            if len(scoreStr) != 0:       # if not an empty string, there is indeed a score. Other options: if scoreStr or if os.path.exists('top_0/top_1.pdb')  # if there is 1 docked structure
                score = float(scoreStr)  # convert from str to float
                self.topScores.append(score)
        printRecord('Number of identified docked structures = {}'.format(len(self.topScores)))
        # return self.topScores
