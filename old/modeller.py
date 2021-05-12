from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *


pdb = PDBFile('orig_half.pdb')
modeller = Modeller(pdb.topology, pdb.positions)
#forcefield = ForceField('AMOEBA-BIO-2018.xml')#amber14-all.xml,'amber14/tip3p.xml', AMOEBA-BIO-2018.xml, amoeba2013.xml
forcefield =  ForceField('amber10.xml')
modeller.addSolvent(forcefield,model='tip3p', padding=1.0*nanometers)