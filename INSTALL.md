Installation
==================================
### CLONING
* Clone this git repository : `git clone https://github.com/persalteas/biorseo.git` and `cd biorseo`.
* Create folders for the modules you will use: `mkdir -p data/modules/`. If you plan to use several module sources, add subdirectories : `mkdir -p data/modules/DESC` and `mkdir -p data/modules/BGSU`

### RNA3DMOTIFS DATA

If you use Rna3Dmotifs, you need to get RNA-MoIP's .DESC dataset: download it from [GitHub](https://github.com/McGill-CSB/RNAMoIP/blob/master/CATALOGUE.tgz). Put all the .desc from the `Non_Redundant_DESC` folder into `./data/modules/DESC`. Otherwise, you also can run Rna3Dmotifs' `catalog` program to get your own DESC modules collection from updated 3D data (download [Rna3Dmotifs](https://rna3dmotif.lri.fr/Rna3Dmotif.tgz)). You also need to move the final DESC files into `./data/modules/DESC`.

### THE RNA 3D MOTIF ATLAS DATA

If not done during the installation of JAR3D, get the latest version of the HL and IL module models from the [BGSU website](http://rna.bgsu.edu/data/jar3d/models/) and extract the Zip files. Put the HL and IL folders from inside the Zip files into `./data/modules/BGSU`. Note that only the latest Zip is required.


### DEPENDENCIES
- Make sure you have Python 3.5+, Cmake, and a C++ compiler installed on your distribution. Please, it's 2019, use a recent one, we use the 2017 C++ standard. The compilation will not work with Ubuntu 16's GCC 5.4 for example. Tested with libstdc++-dev >= 6.0, so use GCC >=6.0 or Clang >= 6.0.
- Install automake, libboost-program-options and libboost-filesystem.
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
You will notice that the installation process is not complete, some of the headers are not well copied to /usr/local. Solve it manually:
```
sudo cp nupack3.2.2/src/thermo/*.h /usr/local/include/nupack/thermo/
```
### OPTIONAL DEPENDENCIES FOR USE OF JAR3D
- Download and install RNAsubopt from the [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/).
- Download and install Java runtime (Tested with Java 10)
- Download the latest JAR3D executable "*jar3d_releasedate.jar*" from [the BGSU website](http://rna.bgsu.edu/data/jar3d/models/). 
  

### OPTIONAL DEPENDENCIES FOR USE OF BAYESPAIRING
- Download and install RNAfold from the [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/) (if not already done at the previous step).
- Make sure you have Python 3.5+ with packages networkx, numpy, regex, wrapt and biopython. You can install them with pip, you will need the python3-dev package to build them.
- Clone the latest BayesPairing Git repo, and install it : 
```
git clone http://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing.git BayesPairing
cd BayesPairing
pip install .
```

### BUILDING
* Edit the file `EditMe` to set the paths of the above dependencies and data. Fields that you will not use can be ignored (ex: bypdir if you do not use BayesPairing). Example of my setup:
    * CPLEXDir="/opt/ibm/ILOG/CPLEX_Studio128_Student"
    * IEIGEN="/usr/local/include/eigen3"
    * INUPACK="/usr/local/include/nupack"
    * jar3dexec="/nhome/siniac/lbecquey/Software/jar3dbin/jar3d_2014-12-11.jar"
    * bypdir="/nhome/siniac/lbecquey/Software/BayesPairing/bayespairing/src"
    * biorseoDir="/nhome/siniac/lbecquey/Software/biorseo"
* You might want to edit `Makefile` if you are not using clang as compiler. For example, if you use g++, replace clang++ by g++.
* Build it: `make -j4`
* Check if the executable file exists: `./bin/biorseo --version`.

### BAYESPAIRING USERS: PREPARE BAYESIAN NETWORKS
We run an example job for it to build the bayesian networks of our modules.
```
cd rnabayespairing/src
python3 parse_sequences.py -d rna3dmotif -seq ACACGGGGUAAGAGCUGAACGCAUCUAAGCUCGAAACCCACUUGGAAAAGAGACACCGCCGAGGUCCCGCGUACAAGACGCGGUCGAUAGACUCGGGGUGUGCGCGUCGAGGUAACGAGACGUUAAGCCCACGAGCACUAACAGACCAAAGCCAUCAU -ss ".................................................................((...............)xxxx(...................................................)xxx).............."
```
Use `-d rna3dmotif` or `-d 3dmotifatlas` depending on the module source you are planning to use.
This is a quite long step, but the bayesian networks will be ready for all the future uses.
