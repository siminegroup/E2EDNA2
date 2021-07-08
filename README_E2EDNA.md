# E2EDNA

We here will lay out the algorithmic logic (a.k.a. the E2EDNA protocl) underpinning this code.
Automation is at the heart of OpenDNA, and so the code makes many decisions behind the scenes which users should be aware of.

## Basic Algorithm

* get 2D fold + analysis
* 3d fold
* 3d relaxation
* docking
* binding

## Detailed Decision-Making

For the purposes of this section, we will assume we are running in `full binding` mode - the most detailed and expensive way to run the pipeline.
The relevant subsections still apply in other modes, but only to their specific subroutines.

### Secondary Structure
Currently, we use only the most-probable structure according to NUPACK (or seqfold), and ignore any possible suboptimal structures.
In future, we hope/plan to leverage NUPACK's ability to predict all the structures likely to be observed within kT of the best structure, cluster them into a few representative groups, and analyze all of them separately.

### Dimensionality Reduction
To identify e.g., the most likely 2D or 3D structure for a given DNA aptamer, we must first reduce the dimensionality of the problem to something manageable for multidimensional histogramming.
Dimensionality reduction for multidimensional probability analysis is done via principal component analysis (PCA).
PCA is fast and foolproof to implement, but limited in its ability to capture all the relevant dynamics of the system. 
In future, we hope/plan to replace/support PCA with nonlinear dimensionality reduction tools e.g., autoencoders.

In terms of implementation, since PCA is cheap, we run it repeatedly until we accrue enough principal components to explain 85% of total problem variance. 
Then, we take all of the principal components with eigenvalues higher than the average of all components. 

In practice, to facilitate multidimensional histogramming, we limit the total number of principal components to a maximum of 5.


### Representative Structure Identification
We identify 'representative' 3D structures via dimensionality reduction of aptamer trajectories (usually in the form of dihedral angles), and multidimensional histogramming.
The number of bins for histogramming is iteratively selected as the maximum amount for which np.histogramdd does not crash.
The resulting multidimensional histogram (probability distribution) is smoothed via multidimensional Gaussian kernel, with standard deviation set adaptively according to the bin size.

A representative sample from the trajectory is found by first identifying the maximum value of the multidimensional probability distribution (a.k.a. the global free energy minimum), then using `scipy.spatial.cKDTree` to identify the trajectory index closest to this point.

One could alternatively generate a structure directly on the probability maximum by back-transforming from PCA space to dihderal, and constructing an aptamer with the exactly correct dihedral angles, but that is a lot of unecessary work, given the fact that the probabilty maximum (if well sampled) should be densely surrounded by existing points in the trajectory (otherwise, it would not be the probability maximum).

### Docking

As with 2D structure, we accept the LightDock docking configuration with the best score, and discard the rest (typically hundreds). 

It would be relatively easy to keep e.g., the top 5-10 docking structures and sample them all with detailed MD. 
Likely we will do something like this in a future development.

### Binding Analysis
To determine whether the aptamer and oligopeptide (or 'peptide') are in contact during a simulation, we start by computing the center-of-geometry distances between all amino acids and nucleotides over a full binding trajectory.

We then compute the number of 'contacts' by counting the number of nucleotides and amino acids which are within a certain cutoff range of one another (in practice, we use a range of ranges, e.g., {6, 8, 10, 12} Angstroms).
We say an aptamer and amino acid are 'in contact' if they have at least one contact within the minimum range at a given time-step. 

The relevant observables are:
* close contact ratio: the proportion of time, after first making contact, during which the aptamer and peptide are in close contact.
* the contact score: an average of the number of contacts over a range of ranges (e.g., {6, 8, 10, 12} Angstroms), after first making close contact, normalized against the number of peptides. The idea here is to relate how closely the peptide is ingratiating itself with the aptamer.
* conformation change: the average values of the principal component trajectories after complexation with the peptide (the free aptamer trajectories have mean 0 by construction, so this is equivalent to taking the difference of the means of the trajectories). 
  Note the principal components used are the ones derived from the free aptamer simulation, the binding trajectory is transformed into this basis for the purpose of this analysis. 
  This is a simplistic way of getting at how the peptide may shift the probability maximum for the aptamer conformation as a whole. 
  In future, we may develop additional metrics which get at larger or smaller scale motion in specific regions of the aptamer which may be relevant e.g., for reporting/signal transduction.
  
### Automatic MD Stopping

If `auto sampling` is turned on, molecular dynamics simulations are done as a seires of runs with fixed length. 
OpenDNA will repeat simulations of this length until either a predetermined maximum number of segments has been reached, or a finishing condition is reached.

The two current conditions for calling off an MD run are:
* The system has reached equilibrium, within a target thereshold. 
Equilibration here is defined by the eigenvalue-weighted average slopes of the principal components of the 3D trajectory approaching zero, indicating either static behaviour or a sufficiently well-explored potential energy surface.
* Unbinding: if OpenDNA detects that the analyte has detached from the aptamer and does not make contact for a preset amount of time (default = 1 nanosecond), the binding simulation will cut off, and the analyte-peptide pair are considered to be non-binding.