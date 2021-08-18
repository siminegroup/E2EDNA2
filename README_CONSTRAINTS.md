# Implementing constraints
### Background

Sometimes it is necessary to implement angular constraints for a peptide, in order to conserve its geometric shape. 
This is particularly important for docking simulations, where a certain ligand will only bind to a peptide if its in a certain conformation.

OpenDNA now allows you to perform and set custom angular constraints. Here's how.

### Setting up the angular constraints

The file "backbone_dihedrals.csv" will be used as input in order to figure out which amino acids need to be constrained.
If we take a look at the file as it is currently, using our favorite text editor, we can see that the following is written in it:

```
residue_num,phi,psi
```

This is the template format that you should follow when inputting the angles into this .csv file.
Please do not remove this first line. Also, please do not change the name of the .csv file.

### The variables
```residue_num``` is the residue number. It starts from 0, to ensure it matches with Python's (and OpenMM's) indexing style. 
In other words, if you need to put a constraint on amino acid 3 for instance (as it appears in the pdb), then you would write ```2``` for the ```residue_num```.

```phi``` is the angle value that you want to set for phi. Likewise for ```psi```. 
Support for an ```omega``` variable does not exist currently.

### An example of a valid input
```
residue_num,phi,psi
1,-90,-110
5,-85,-120
```

