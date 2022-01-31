# Implementing constraints
### Background

Sometimes it is necessary to implement angular constraints for a peptide, in order to conserve its geometric shape. 
This is particularly important for docking simulations, where a certain ligand will only bind to a peptide if its in a certain conformation.

One of the most popular (and simple-to-use) constraints is the harmonic one. It has the following formula:
<img src="https://latex.codecogs.com/png.latex?%5Cbg_white%20V(x)=0.5k(x-x_0)^2 " />

Where <img src="https://latex.codecogs.com/png.latex?%5Cbg_white%20V " />
is the potential energy, 
<img src="https://latex.codecogs.com/png.latex?%5Cbg_white%20k " /> 
is the force constant, 
<img src="https://latex.codecogs.com/png.latex?%5Cbg_white%20x_0 " /> 
is an initial variable (the value to which the constraint will try and fix the variable) and 
<img src="https://latex.codecogs.com/png.latex?%5Cbg_white%20x " />
is the variable in question, which can be anything from the distance between 2 atoms, to the torsion angle of a dihedral. 

OpenDNA now allows you to perform and set custom angular constraints (for proteins/peptides only). Here's how.

### Setting up the angular constraints

The file "backbone_dihedrals.csv" will be used as input in order to figure out which amino acids need to be constrained.
If we take a look at the file as it is currently, using our favorite text editor, we can see that the following is written in it:

```
residue_num,phi,psi,chain_id
```

This is the template format that you should follow when inputting the angles into this .csv file.
Please do not remove this first line. Also, please do not change the name of this .csv file.

The strength of the angular constraint can be set with the ```params['peptide backbone constraint constant'] = YOUR_INTEGER_VALUE``` statement in ```main.py```.
It is set at ```10000``` currently because that was determined to be the best value for constraining a peptide, when all other parameters are held constant.

### The variables
```residue_num``` is the residue number. It starts from 0, to ensure it matches with Python's (and OpenMM's) indexing style. 
In other words, if you need to put a constraint on amino acid 3 for instance (as it appears in the pdb), then you would write ```2``` for the ```residue_num```.
Since ```residue_num``` essentially functions as an ID for the residue, having more than one input rows with the same ```residue_num``` will lead to an error.

```phi``` is the angle value that you want to set for phi (in degrees). Likewise for ```psi```. 
Support for an ```omega``` variable does not exist currently. The assignments for each dihedral angle follows the format laid out by the ```PeptideBuilder``` package - see the ```PeptideBuilder.py``` [file on their GitHub](https://github.com/clauswilke/PeptideBuilder/blob/6d38a167b9992c27adc86f64370f7083303ce877/PeptideBuilder/PeptideBuilder.py) for more info.
Therefore, it is impossible to add a constraint to psi and phi on the first amino acid at the start of any peptide chain. 

```chain_id``` is the ID of the chain to which the atom in question belongs. It can be any finite positive integer, including zero. 
The chain ID for an atom can be found by looking at the chain letter in the pdb file, then converting it to a number. The format for doing so is as follows:

```
Chain Letter A = chain_id 0
Chain Letter B = chain_id 1
Chain Letter C = chain_id 2
.
.
.
etc...
```

### An example of a valid input
The following is an example of a valid input (valid backbone_dihedrals.csv file) for constraining phi and psi values for amino acids 2 and 6 found in Chain A of an arbitrary peptide:
```
residue_num,phi,psi,chain_id
1,-90,-110,0
5,-85,-120,0
```

### Common errors that might occur
*Assuming everything above is followed*, some common errors that might show up are:

1. ```Particle Coordinate is nan```

This means the simulation exploded, due to the time step being too large to properly calculate trajectories using the current
```params['peptide backbone constraint constant']``` value. An easy way to rectify this is to simply make ```params['time step']``` shorter.
It is currently set at 2.0 fs because that was the best time-step for the current ```params['peptide backbone constraint constant']``` value.

2. Dihedral analysis issues

Depending on the type of peptide simulated (and the type of analysis package used to find the dihedrals), there may be a few inconsistencies between
PeptideBuilder and MDAnalysis. The constraints in OpenDNA are written to be as compatible as possible, but there may still need to be some things to account for.
For instance, creating a GLY-GLY dimer and simulating it with OpenDNA while attempting to use ```1,-60,140,0``` to put 2 constraints on the 2nd glycine will actually put 1 constraint (phi) on the 2nd glycine and 1 constraint (psi) on the 1st glycine, according to the dihedral analysis done via MDAnalysis. 
This is because there is no phi for the first residue (```aa_id = 0```) for any peptide, 
according to how MDAnalysis analyzes the dihedrals. Likewise, there is no psi on the last residue for any peptide, for the same reason. 
If one were to print out the dihedrals from the GLY-GLY dimer as a list using the phi_selection() method in MDAnalysis , one would get the following:
```
Phi: [None, <AtomGroup with 4 atoms>]
Psi: [<AtomGroup with 4 atoms>, None]
```
Which shows the problem mentioned earlier in this paragraph.

****This is by no means an exhaustive list. If more errors are found, they will be added here.**
