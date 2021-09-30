Biorseo (Bi-Objective RNA Structure Efficient Optimizer)
===================================

This tool predicts the secondary structure of a RNA sequence with pieces of 3D information (non-canonical contacts) at some places, 
by identifying zones that can fold like known modules from data like the RNA 3D Motif Atlas or Rna3Dmotifs.

Contact : louis.becquey@univ-evry.fr

1/ How it works
===================================
INPUT:
- An RNA sequence (with 16 GB of RAM you can go up to ~230 bases)

THEN
- **Pattern-matching step** : Find all possible occurrences of known RNAmodules in the query sequence, by finding subsequences of the query that score well with the probabilistic models of the modules (like JAR3D, or BayesPairing)
- **Constraints  definition  step** : Define constraints on the secondary structure imposed by modules if they would be included (in this case, the loop closing base-pairs are mandatory)
- **Solve a bi-objective IP problem** : Find a secondary structure that satisfies as much as possible both the expected accuracy of the structure and a criterion taking into account module inclusions, by solving a bi-objective integer linear programming problem, using the previous constraints defined in the previous step.

OUTPUT:
- A set of secondary structures from the Pareto front,
- The list of known modules inserted inplace in the corresponding structures

2/ The different models
==================================
MODULE SOURCES

Biorseo can be used with two modules datasets (yet):
* Rna3Dmotifs (from the work of *Djelloul & Denise, 2008*)
* The RNA 3D Motif Atlas of BGSU's RNA lab (*Petrov et al, 2013*, see http://rna.bgsu.edu/rna3dhub/motifs/)
* CaRNAval 1.0 (*Reinhartz et al, 2018*)
* RNA-Bricks 2, RNAMC, CaRNAval 2.0, and others could theoretically be used, but are not supported (yet). You might write your own API.

PATTERN MATCHING STEP
- Use **simple pattern matching**. Rna3Dmotifs modules are available with sequence information. We use regular expressions to find those known loops in your query. This is the approach of RNA-MoIP (*Reinharz et al, 2012*), we deal the same way with short components and wildcards.

- Use **JAR3D**. The RNA 3D Motif Atlas modules can be scored against a given loop sequence by an hybrid SCFG/MRF method (*Zirbel et al, 2015*). This first requires to identify potential loops, which is achieved by a run of RNAsubopt first.

- Use **Bayesian networks with BayesPairing**. To accurately model probability dependancies between nucleotides, one can use BayesPairing to build bayesian networks of the modules (the RNA 3D Motif Atlas and Rna3Dmotifs are both supported). Then, sequences are sampled with the Bayesian network of a module, and we use regular expressions to find them in your query.

OBJECTIVE FUNCTIONS FOR THE MODULE INSERTION CRITERIA

* **Function A** : weights a module by its squared number of nucleotides (like RNA-MoIP).
* **Function B** : weights a module by its number of components (strands) and penalizes it by the log^(_2) of its nucleotide size.
* **Function C** : weights a module by its insertion site score (JAR3D or BayesPairing score).
* **Function D** : weights a module by its number of components (strands) and insertion site score (JAR3D or BayesPairing score), and penalizes it by the log^(_2) of its nucleotide size.

3/ Installation
==================================
Check the file [INSTALL.md](INSTALL.md) for installation instructions.

4/ Recommended uses
==================================
- If **you know you have no pseudoknot**:
    * Benchmarks show Biorseo does not perform better than simpler tools like RNAsubopt alone. Please use RNAsubopt (ViennaRNA package) or Fold (RNAstructure package).

- If you **might expect a pseudoknot, or don't know**:
    * The most promising method is the use of direct pattern matching with Rna3Dmotifs and function A. But this method is sometimes subject to combinatorial explosion issues. If you have a long RNA or a large number of loops, don't use it. Example:
    `./biorseo.py -i PDB_00304.fa -O resultsFolder/ --rna3dmotifs --patternmatch --func A`
    
    * The use of the RNA 3D Motif Atlas placed by JAR3D and scored with function A is not subject to combinatorial issues, but performs a bit worse. It also returns less solutions. Example:
    `./biorseo.py -i PDB_00304.fa -O resultsFolder/ --3dmotifatlas --jar3d --func A

5/ List of Options
==================================
```
Usage:  You must provide:
        1) a FASTA input file with -i,
        2) a module type with --rna3dmotifs, --carnaval, --3dmotifatlas or --contacts
        3) one module placement method in { --patternmatch, --jar3d, --bayespairing }
        4) one scoring function with --func A, B, C, D, E ou F

        If you are not using the Docker image: 
        5) --modules-path, --biorseo-dir and (--jar3d-exec or --bypdir)

Options:
-h [ --help ]                   Print this help message
--version                       Print the program version
-i [ --seq=… ]                  FASTA file with the query RNA sequence
--rna3dmotifs                   Use DESC modules from Djelloul & Denise, 2008
--carnaval                      Use RIN modules from Reinharz & al, 2018
--3dmotifatlas                  Use the HL and IL loops from BGSU's 3D Motif Atlas (updated)
--contacts			Use the library of motifs, created from RNA sequences linked to proteins provided by I. Chauvot de Beauchene of LORIA laboratory
-p [ --patternmatch ]           Use regular expressions to place modules in the sequence (requires --rna3dmotifs or --carnaval)
-j [ --jar3d ]                  Use JAR3D to place modules in the sequence (requires --3dmotifatlas)
-b [ --bayespairing ]           Use BayesPairing2 to place modules in the sequence (requires --rna3dmotifs or --3dmotifatlas)
-o [ --output=… ]               File to summarize the results
-O [ --outputf=… ]              Folder where to output result and temp files
-f [ --func=… ]                 (A, B, C or D, default is B) Objective function to score module insertions:
                                  (A) insert big modules (B) insert light, high-order modules
                                  (c) insert modules which score well with the sequence
                                  (D) insert light, high-order modules which score well with the sequence.
                                  C and D require cannot be used with --patternmatch.
-c [ --first-objective=… ]      (default 1) Objective to solve in the mono-objective portions of the algorithm.
                                  (1) is the module objective given by --func, (2) is the expected accuracy of the structure.
-l [ --limit=… ]                (default 500) Number of solutions in the Pareto set from which
                                  we give up the computation, before eliminating secondary structure doublons.
-t [ --theta=… ]                (default 0.001) Pairing-probability threshold to consider or not the possibility of pairing
-n [ --disable-pseudoknots ]    Add constraints to explicitly forbid the formation of pseudoknots
-v [ --verbose ]                Print what is happening to stdout
--modules-path=…                Path to the modules data.
                                  The folder should contain modules in the DESC format as output by Djelloul & Denise's
                                  'catalog' program for use with --rna3dmotifs, or the IL/ and HL/ folders
                                  from a release of the RNA 3D Motif Atlas for use with --3dmotifatlas, or the
                                  data/modules/RIN/Subfiles/ folder for use with --carnaval.
                                  Consider placing these files on a fast I/O device (NVMe SSD, ...)
--jar3d-exec=…                  Path to the jar3d executable.
                                  Default is /jar3d_2014-12-11.jar, you should use this option if you run
                                  BiORSEO from outside the docker image.
--bypdir=…                      Path to the BayesParing src folder.
                                  Default is /byp/src, you should use this option if you run
                                  BiORSEO from outside the docker image.
--biorseo-dir=…                 Path to the BiORSEO root directory.
                                  Default is /biorseo, you should use this option if you run
                                  BiORSEO from outside the docker image. Use the FULL path.

Examples:
biorseo.py -i myRNA.fa -O myResultsFolder/ --rna3dmotifs --patternmatch --func B
biorseo.py -i myRNA.fa -O myResultsFolder/ --3dmotifatlas --jar3d --func B -l 800
biorseo.py -i myRNA.fa -v --3dmotifatlas --bayespairing --func D

The allowed module/placement-method/function combinations are:

                --patternmatch  --bayespairing    --jar3d
--rna3dmotifs     A. B.           A. B. C. D.
--3dmotifatlas                    A. B. C. D.     A. B. C. D.
--carnaval        A. B.
--contacts 	  E. F.

```
