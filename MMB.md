# MacroMolecule Builder (MMB)

MMB is a multifunctional software from simTK which allows for rapid directed folding of oligonucleotides and peptides via straightforward inputs on a wide variety of platforms.

In this code, we use pre-built binaries downloaded from simTK.org. We have tested versions as old as 2.14 and as new as 3.4 variously on Windows 10, WSL2 and Linux, and found our scripts to function well with all of them.
This is perhaps due to the fact that our MMB scripts are extremely simple, and use only the basic level of functionality MMB affords. 

In short, MMB initializes a given ssDNA sequence in a single-helix configuration, and folds it according to user-specified base-pairing conditions via simulation with strong ficticious forces which pull the respective bases together. 
We include in this code scripts which automatically take in secondary structure instructions from seqfold or NUPACK in the form of lists of paired bases, and generate MMB command files accordingly. 

In general, it is a good idea to ramp up the strength of ficticious base-base forces over the course of an MMB folding run (think of this as annealing), and the default protocol we use in this code in `lib/mmb/commands.template.dat` strikes a good balance of speed and quality for folding ssDNA sequences at least up to 60-bases in length.
If results are unsatisfactory, one can easily adjust this script with guidance from the MMB reference guide enclosed in every MMB download.

While MMB simulations generally take only ~5-20 minutes on a single core, one could parallelize this step, or take advantage of more advanced MMB folding commands - we leave that to users' discretion for now.

It is worth noting that MMB outputs, while technically valid 3D representations which agree with the prescribed secondary structure given by seqfold or NUPACK, are generally not all that physical. 
One may consider the raw MMB output as a 'coarse' or highly strained version of the 3D structure of the aptamer, which generally benefits from relaxation via MD simulation. 

Finally, MMB output .pdb files are a bit quirky, and require post-processing to run properly in OpenMM. 
See the comments in opendna.prepPDB for details.