Changelog
=========================

### Biorseo 2.1 (Nov 2021)
This is an official, tested, release of Biorseo 2:
- replacing Nupack's dynamic programming scheme supporting simple pseudoknots by ViennaRNA's window-based scheme, which does not support pseudoknots or long-distance contacts but allows to test much longer sequences,
- supporting RINs with no issues, 
- supporting custom modules in JSON format (to be detected in sequences by regular expressions), thanks to Nathalie Bernard
- not running Jar3d or BayesPairing for you anymore. This simplifies a lot the code management (replacing a pipeline by the C++ tool only). Jar3d is getting older, does not support very complex modules, and is biaised because it takes as input loops (not the whole sequence). Therefore, you have to give biorseo the answer as input ! BayesPairing 2.0 is evolving itself into a module-placement tool in secondary structures taking eneries into account (and now comparative information), it is a non sense to include it *within* Biorseo. Approaches should be compared and benchmarked instead. But, you can still use the ouputs of this tools as input for biorseo if you like.
- introducing the MFE criterion (thanks to Nathalie Bernard),
- introducing the Biokop-mode,
- with a much simpler and lighter installation process.

Biorseo 2.1 is availbale as a docker container and as a git branch called "biorseo2".
It is the last version supported by Louis Becquey.

### Biorseo 2.0
This was an unofficial, unsupported and unpublished version after the internship of Lénaic Durand at IBISC.
This version 
- replaces Nupack 3.2 with ViennaRNA to compute the pairing probabilities, thanks to Lénaic,
- introduces early support for BayesPairing 2.0, which was still unofficial too at the time,
- supports CaRNAval RINs,
- but has issues with the constraints to assert RIN basepairs are respected.

Results from this version are published in [Louis Becquey's thesis](https://tel.archives-ouvertes.fr/tel-03440181).

### Biorseo 1.2 (2019) and Biorseo 1.5 (2020)
These brought some improvements, fixing numerical issues, and other technical improvements.
Biorseo 1.2 is still available as a docker, and the 1.5 is available as a Git branch called 'biorseo1'.

### Biorseo 1.0 (2018)
This was the first version published for the paper [*Becquey et al. 2020*.](https://doi.org/10.1093/bioinformatics/btz962)