This is a bi-objective integer programming algorithm.
It predicts the secondary structure of a RNA sequence with pieces of 3D information (non-canonical contacts) at some places, 
by identifying zones that can fold like known motifs from the RNA 3D Motif Atlas.

1/ How it works
===================================
INPUT:
- An RNA sequence (tested with sequences ~100 bases)

THEN
- Identifies possible 2D folds with RNAsubopt.
- Knowing possible 2D folds, locate every possibly unpaired loop (hairpin loop, internal loop, multiple loop...)
- align each unpaired loop to the catalogue of models of known RNA motifs (The 3D Motif Atlas of the BGSU RNA group)
- retrieve a list of potential motif-insertion-sites in the RNA sequence. Use them to define the constraints for the IP problem.
- Solve a bi-objective IP problem: 
    * Maximize the expected accuracy of the secondary structure,
    * Maximize the number and size of motifs inserted in the structure.

OUTPUT:
- A set of secondary structures from the pareto front,
- The list of known motif inserted in the corresponding structures (and the non-canonical contacts)
- (lower score structures from k-Pareto sets, not implemented yet.)

2/ Installation
==================================
- Download and install RNAsubopt from the ViennaRNA package (https://www.tbi.univie.ac.at/RNA/)
- Download and install IBM ILOG Cplex optimization studio (https://www.ibm.com/analytics/cplex-optimizer), free academic account required
- Download and install Java runtime (Tested with Java 10)
- Download and install the latest JAR3D executable "jar3d_releasedate.jar" and motif models in this folder (http://rna.bgsu.edu/data/jar3d/models/)
  Note that for HL and ILs, only the latest version is required (not all the versions provided in the folders).
- Download and install a C++ compiler and building dependencies and utilities (g++ or clang, automake, libboost)
