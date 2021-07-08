# LightDock

We use the LightDock python package for docking aptamer to analyte molecules.

In practice, we run LightDock on a representative 3D structure of the free aptamer, and attempt to dock a (usually smaller) peptide molecule on to it. 

LightDock uses a glowworm swarm algorithm, which attacks the aptamer with swarms of 'glowworms' representing the analyte molecule, coming from various directions.
We compute the number of swarms required as being proportionate to the approximate surface area of a cube with the dimensions of the target aptamer (more surface area, more swarms). 
For the number of glowworms and docking steps, we use default values which seem to give good results (300 & 100, respectively). 

After docking, we generate configurations for each glowworm in each swarm, cluster similar structures together, and score them. 
It is then straightforward to select the best structures and send them back to OpenMM for further MD simulation.

There is a certain amount of alchemy involved in adding & rmoving hydrogens and terminal groups in order to get LightDock and OpenMM to work together. 
See the relevant scripts for detailed comments. 