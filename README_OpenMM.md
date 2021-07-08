# OpenMM

OpenMM is the molecular dynamics engine which powers OpenDNA. 
It is available as an easy-to-install python package, runs on CPU or GPU, is compatible with a wide variety of force fields and simulation types, and admits the use of many powerful analysis utilities. 

We make heavy use of OpenMM utilities such as pdbfixer (the engine behind OpenMM-setup) for preparing our input files for simulation, in addition to OpenMM itself. 

In the future, we hope to implement implicit solvation in OpenMM, which unfortunately is not possible for nucleotides using the included force fields, but requires that a given simulation is built from e.g., AMBER .prmtop files.

The OpenMM code itself is fairly self-explanatory, and we direct users to the excellent [documentation](http://docs.openmm.org/latest/userguide/index.html).

One final point worth noting is that is is not clear whether the AMOEBA force field functional form implemented in OpenMM (corresponding the AMOEBA13) is fully compatible with e.g., the new nucleic acid parameters in AMOEBABIO18, even if one appropriately converts the TINKER .prm files to OpenMM .xml parameter files.
For this reason, we advise users interested in AMOEBA simulations to look at our TINKER-based implementation of [E2EDNA](https://github.com/InfluenceFunctional/E2EDNA), which while less automated, is tailor-made for usage with AMOEBA.