# NUPACK (and seqfold)

NUPACK 4.x and seqfold are the two options for secondary structure prediction engines in this code.
NUPACK versions > 4.0 are extremely rich in features for quantitative analysis of secondary structure and probability, particularly for longer aptamers.
Indeed, there are a great number of NUPACK features which go beyond the current scope of this code.

Seqfold by comparison is easier to install, simpler, and considerably more bare-bones feature-wise. 
It is further worth noting that these two packages do not always agree on the most probable secondary structure for a given sequence.

A few key NUPACK features we take advantage of are:
* Tune predictions according to temperature and ionic strength
* Output explicit probability of observing the most likely secondary structure for a given sequence.
* Output suboptimal structures and their probabilities for a given sequence.

NUPACK has also a wealth of design and utility functions we encourage users to play with on their own.
