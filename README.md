Biorseo (Bi-Objective RNA Structure Efficient Optimizer)
===================================

This tool predicts the secondary structure of a RNA sequence with pieces of 3D information (non-canonical contacts) at some places, 
by identifying zones that can fold like known modules from data like the RNA 3D Motif Atlas or Rna3Dmotifs.

Contact : louis.becquey@univ-evry.fr

1/ How it works
===================================
INPUT:
- An RNA sequence (tested with sequences ~100 bases)

THEN
- **Pattern-matching step** : Find all possible occurrences of known RNAmodules in the query sequence, by finding subsequences of the querythat score well with the probabilistic models of the modules (like JAR3D, or BayesPairing)
- **Constraints  definition  step** : Define constraints on the secondary structure imposed by modules if they would be included (in this case, some of the canonical base-pairs are forbidden)
- **Solve a bi-objective IP problem** : Find a secondary structure that satisfies as much as possible both the expected accuracy of the structure and a criterion taking into account module inclusions, by solving a bi-objective integer linear programming problem, using the previous constraints defined in the previous step.

OUTPUT:
- A set of secondary structures from the Pareto front,
- The list of known modules inserted inplace in the corresponding structures

2/ The different models
==================================
MODULE SOURCES

Biorseo can be used with two modules datasets (yet):
* Rna3Dmotifs (from the work of *Djelloul & Denise, 2008*), but with the 3D data of 2018
* The RNA 3D Motif Atlas of BGSU's RNA lab (*Petrov et al, 2013*, see http://rna.bgsu.edu/rna3dhub/motifs/)
* RNA-Bricks 2 or CaRNAval could theoretically be used, but are not supported (yet). You might write your own API.

PATTERN MATCHING STEP
- Use **simple pattern matching**. Rna3Dmotifs modules are available with sequence information. We use regular expressions to find those known loops in your query. This is the approach of RNA-MoIP (*Reinharz et al, 2012*), we deal the same way with short components and wildcards.

- Use **JAR3D**. The RNA 3D Motif Atlas modules can be scored against a given loop sequence by an hybrid SCFG/MRF method (*Zirbel et al, 2015*). This first requires to identify potential loops, which is achieved by a run of RNAsubopt first.

- Use **Bayesian networks with BayesPairing**. To accurately model probability dependancies between nucleotides, one can use BayesPairing to build bayesian networks of the modules (the RNA 3D Motif Atlas and Rna3Dmotifs are both supported). Then, sequences are sampled with the Bayesian network of a module, and we use regular expressions to find them in your query.

OBJECTIVE FUNCTIONS FOR THE MODULE INSERTION CRITERIA

* **Function A** : weights a module by its squared number of nucleotides (like RNA-MoIP).
* **Function B** : weights a module by its number of components (strands) and penalizes it by the log^(_2) of its nucleotide size.
* **Function C** : weights a module by its insertion site score (JAR3D or BayesPairing score).
* **Function D** : weights a module by its number of components (strands) and insertion site score (JAR3D or BayesPairing score), and penalizes it by the log^(_2) of its nucleotide size.

3/ Recommended uses
==================================
- If **you know you have no pseudoknot**:
    * Benchmarks show Biorseo does not perform better than simpler tools like RNAsubopt alone. Please use RNAsubopt (ViennaRNA package) or Fold (RNAstructure package).

- If you **might expect a pseudoknot, or don't know**:
    * The most promising method is the use of direct pattern matching with Rna3Dmotifs and function B. But this method is sometimes subject to combinatorial explosion issues. If you have a long RNA or a large number of loops, don't use it. Example:
    `./bin/biorseo -s PDB_00304.fa --descfolder ./data/modules/DESC --type B -o PDB_00304.rawB `
    
    * The use of the RNA 3D Motif Atlas placed by JAR3D and scored with function B is not subject to combinatorial issues, but performs a bit worse. It also returns less solutions. Example:
    `./bin/biorseo -s PDB_00304.fa --jar3dcsv PDB_00304.sites.csv --type B -o PDB_00304.jar3dB`


4/ Installation
==================================
### DEPENDENCIES
- Make sure you have Python 3.5+, Cmake, and a C++ compiler installed on your distribution.
- Install automake and libboost-filesystem.
- Download and install [IBM ILOG Cplex optimization studio](https://www.ibm.com/analytics/cplex-optimizer), an academic account is required. The free version is too limited, you must register as academic. This is also free.
- Download and install Eigen: Get the latest Eigen archive from http://eigen.tuxfamily.org. Unpack it, and install it.
```bash
wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz -O eigen_src.tar.gz
tar -xf eigen_src.tar.gz
cd eigen-eigen-323c052e1731
mkdir build
cd build
cmake ..
sudo make install
```
- Download and install NUPACK: Register on [Nupack's website](http://www.nupack.org/downloads/source), download the source, unpack it, build it, and install it:
```bash
wget http://www.nupack.org/downloads/serve_file/nupack3.2.2.tar.gz
tar -xf nupack3.2.2.tar.gz
cd nupack3.2.2
mkdir build
cd build
cmake ..
make -j4
sudo make install
```

### OPTIONAL DEPENDENCIES FOR USE OF JAR3D
- Download and install RNAsubopt from the [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/).
- Download and install Java runtime (Tested with Java 10)
- Download the latest JAR3D executable "*jar3d_releasedate.jar*", and latest IL and HL models from [here](http://rna.bgsu.edu/data/jar3d/models/).
  Note that only the latest version is required (not all the versions provided in the folders).

### OPTIONAL DEPENDENCIES FOR USE OF BAYESPAIRING
- Download and install RNAfold from the [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/).
- Make sure you have Python 3.5+ with packages networkx, numpy, regex, wrapt and biopython
- Clone the latest BayesPairing Git repo, and install it : 
```
git clone http://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing.git BayesPairing
cd BayesPairing
pip install .
```

### RNA3DMOTIFS DATA

If you use Rna3Dmotifs, you need to get RNA-MoIP's .DESC dataset: download it from [GitHub](https://github.com/McGill-CSB/RNAMoIP/blob/master/CATALOGUE.tgz). Put all the .desc from the `Non_Redundant_DESC` folder into `./data/modules/DESC`. Otherwise, you also can run Rna3Dmotifs' `catalog` program to get your own DESC modules collection from updated 3D data (download [Rna3Dmotifs](https://rna3dmotif.lri.fr/Rna3Dmotif.tgz)). You also need to move the final DESC files into `./data/modules/DESC`.

### THE RNA 3D MOTIF ATLAS DATA

If not done during the installation of JAR3D, get the latest version of the HL and IL module models from the [BGSU website](http://rna.bgsu.edu/data/jar3d/models/) and extract the Zip files. Put the HL and IL folders into `./data/modules/BGSU`.

### BUILDING
* Clone this git repository : `git clone https://github.com/persalteas/biorseo.git` and `cd biorseo`.
* Edit the file `EditMe` to set the paths of the above dependencies and data. Fileds that you will not use can be ignored (ex: bypdir if you do not use BayesPairing). Example of my setup:
    * CPLEXDir="/opt/ibm/ILOG/CPLEX_Studio128_Student"
    * IEIGEN="/usr/local/include/eigen3"
    * INUPACK="/usr/local/include/nupack"
    * jar3dexec="/nhome/siniac/lbecquey/Software/jar3dbin/jar3d_2014-12-11.jar"
    * ILmotifDir="/nhome/siniac/lbecquey/Data/RNA/motifs/BGSU/Matlab_results/IL/3.2/lib"
    * HLmotifDir="/nhome/siniac/lbecquey/Data/RNA/motifs/BGSU/Matlab_results/HL/3.2/lib"
    * descfolder="/nhome/siniac/lbecquey/Data/RNA/motifs/Rna3Dmotifs/No_Redondance_DESC/"
    * bypdir="/nhome/siniac/lbecquey/Software/BayesPairing/bayespairing/src"
    * biorseoDir="/nhome/siniac/lbecquey/Software/biorseo"
* You might want to edit `Makefile` if you are not using clang as compiler. For example, if you use g++, replace clang++ by g++.
* Build it: `make -j4`
* The working executable file is `./bin/biorseo`.

### BAYESPAIRING USERS: PREPARE BAYESIAN NETWORKS
We run an example job for it to build the bayesian networks of our modules.
```
cd rnabayespairing/src
python3 parse_sequences.py -d rna3dmotif -seq ACACGGGGUAAGAGCUGAACGCAUCUAAGCUCGAAACCCACUUGGAAAAGAGACACCGCCGAGGUCCCGCGUACAAGACGCGGUCGAUAGACUCGGGGUGUGCGCGUCGAGGUAACGAGACGUUAAGCCCACGAGCACUAACAGACCAAAGCCAUCAU -ss ".................................................................((...............)xxxx(...................................................)xxx).............."
```
Use `-d rna3dmotif` or `-d 3dmotifatlas` depending on the module source you are planning to use.
This is a quite long step, but the bayesian networks will be ready for all the future uses.
