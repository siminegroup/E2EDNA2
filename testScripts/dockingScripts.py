# scripts for aptamer-analyte docking
import os
import simtk.openmm.app as app
from utils import *

'''
=> optimize
=> automate picking by score
    --> add back protons to dna topology after picking
=> add restraints
=> add anm
'''

os.chdir('/mnt/c/Users/mikem/Desktop/mmruns/lightdock')


def killH(structure):
    pdb = app.PDBFile(structure)
    topology = pdb.topology
    positions = pdb.positions
    modeller = app.Modeller(topology,positions)
    modeller.delete(atom for atom in topology.atoms() if "H" in atom.name)
    app.PDBFile.writeFile(modeller.topology,modeller.positions,open(structure.split('.pdb')[0] + '_noH.pdb','w'))


def recenterPDB(structure):
    pdb = app.PDBFile(structure)
    topology = pdb.topology
    positions = pdb.positions
    modeller = app.Modeller(topology, positions)
    app.PDBFile.writeFile(modeller.topology, modeller.positions, open(structure.split('.pdb')[0] + '_recentered.pdb','w'))

def getMoleculeSize(structure):
    '''
    use MDAnalysis to determine the cubic xyz dimensions for a given molecule
    '''
    u = mda.Universe(structure)
    positions = u.atoms.positions
    dimensions = np.zeros(3)
    for i in range(3):
        dimensions[i] = np.ptp(positions[:,i])

    return dimensions


def lightDock(aptamer, analyte):
    '''
    setup lightdock run, do the run, generate output models, hierarchical clustering and final concatenation of results
    '''
    params = {}
    params['setup path'] = '/home/mkilgour/miniconda3/bin/lightdock3_setup.py'
    params['lightdock path'] = '/home/mkilgour/miniconda3/bin/lightdock3.py'
    params['lgd generate path'] = '/home/mkilgour/miniconda3/bin/lgd_generate_conformations.py'
    params['lgd cluster path'] = '/home/mkilgour/miniconda3/bin/lgd_cluster_bsas.py'
    params['anthony py'] = '/home/mkilgour/miniconda3/bin/ant_thony.py'
    params['lgd rank'] = '/home/mkilgour/miniconda3/bin/lgd_rank.py'
    params['lgd top'] = '/home/mkilgour/miniconda3/bin/lgd_top.py'

    # identify aptamer size
    dimensions = getMoleculeSize(aptamer)
    # compute surface area of rectangular prism with no offset
    surfaceArea = 2 * (dimensions[0]*dimensions[1] + dimensions[1]*dimensions[2] + dimensions[2]*dimensions[1]) # in angstroms
    # compute number of swarms which cover this surface area - swarms are spheres 2 nm in diameter - this is only approximately the right number of swarms, since have no offset, and are using a rectangle
    nSwarms = np.ceil(surfaceArea / (np.pi * 10 ** 2))

    params['swarms'] = int(nSwarms) # number of glowworm swarms
    params['glowworms'] = 300 # number of glowworms per swarm
    params['docking steps'] = 200 # number of steps per docking run
    params['N docked structures'] = 5 # number of docked structures to output

    killH(aptamer)
    addH(analyte.split('.')[0] + '_reIndexed.pdb',7.4)
    changeSegment(analyte,'A','B')

    aptamer2 = aptamer.split('.')[0] + "_noH.pdb"
    analyte2 = analyte.split('.')[0] + "_H.pdb" #analyte.split('.pdb')[0] + "_noH.pdb"

    # run setup
    os.system(params['setup path'] + ' ' + aptamer2 + ' ' + analyte2 + ' -s ' + str(params['swarms']) + ' -g ' + str(params['glowworms']))

    # run docking
    os.system(params['lightdock path'] + ' setup.json ' + str(params['docking steps']) + ' -s dna')

    # generate docked structures and cluster them
    for i in range(params['swarms']):
        os.chdir('swarm_%d'%i)
        os.system(params['lgd generate path'] + ' ../' + aptamer2 + ' ../' + analyte2 + ' gso_%d'%params['docking steps'] + '.out' + ' %d'%params['docking steps'] + ' > /dev/null 2> /dev/null; >> generate_lightdock.list') # generate configurations
        os.system(params['lgd cluster path'] + ' gso_%d'%params['docking steps'] + '.out >> cluster_lightdock.list') # cluster glowworms
        # os.system(params['anthony py'] + ' generate_lightdock.list')
        # os.system(params['anthony py'] + ' cluster_lightdock.list')

        os.chdir('../')


    os.system(params['lgd rank'] + ' %d' % params['swarms'] + ' %d' % params['docking steps']) # rank the clustered docking setups

    # generate top structures
    os.system(params['lgd top'] + ' ' + aptamer2 + ' ' + analyte2 + ' rank_by_scoring.list %d'%params['N docked structures'])
    os.mkdir('top')
    os.system('mv top*.pdb top/') # collect top structures



aptamer = 'repStructure.pdb'
analyte = 'peptide.pdb'
bestDockedStructures = lightDock(aptamer, analyte)
for i in range(len(bestDockedStructures)):
    prepDockedStructure(bestDockedStructures[i,0],bestDockedStructures[i,1])

