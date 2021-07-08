# MDAnalysis and MDTraj

MDAnalysis (MDA) and MDTraj are the basic tools for automated analysis in OpenDNA.
While we have tried to thoroughly comment our usage, particularly of MDAnalysis, the [documentation](https://docs.mdanalysis.org/stable/index.html) may still be helpful to those learning it for the first time.

The basic object in MDA is the Universe, usually denoted by "_u_", which is loaded as a combination of a topology (.pdb file) and trajectory (.dcd file). 
We undertake a wide range of manipulations on this object, using mostly straightforward MDA functions.

One point worth mentioning is on the functions for analysis of nucleic acids. 
While there are preexisting functions for analysis e.g., of Watson-Crick hydrogen bond lengths in mda.analysis.nuclinfo, these functions are somewhat old, and prohibitively slow for analysis of lengthy trajectories.
To remedy this, we have written the functions `wcTrajAnalysis` and `nucleicDihedrals` using modern, very fast MDA methods, for Watson-Crick hydrogen bond length and nucleic dihedral angle analysis respectively. 
When we have time, we hope to clean these functions up for incorporation into MDA proper.