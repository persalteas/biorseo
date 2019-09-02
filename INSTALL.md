Option 1 : Installation using docker image (Windows, Mac, Linux)
==================================
* Clone this git repository : `git clone https://github.com/persalteas/biorseo.git` 
* Move into the cloned repository ( `cd biorseo` )

### Install Docker:
* See the officiel instructions depending on your OS here : https://docs.docker.com/install/
* Example for Ubuntu Linux:
```
sudo apt-get install apt-transport-https ca-certificates curl gnupg-agent software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
sudo apt update && sudo apt install docker-ce docker-ce-cli containerd.io
```

### Download and install the RNA motifs data files:
* If you use Rna3Dmotifs, you need to get RNA-MoIP's .DESC dataset: download it from [GitHub](https://github.com/McGill-CSB/RNAMoIP/blob/master/CATALOGUE.tgz). Put all the .desc from the `Non_Redundant_DESC` folder into `./data/modules/DESC`. Otherwise, you also can run Rna3Dmotifs' `catalog` program to get your own DESC modules collection from updated 3D data (download [Rna3Dmotifs](https://rna3dmotif.lri.fr/Rna3Dmotif.tgz)). You also need to move the final DESC files into `./data/modules/DESC`.

* Get the latest version of the HL and IL module models from the [BGSU website](http://rna.bgsu.edu/data/jar3d/models/) and extract the Zip files. Put the HL and IL folders from inside the Zip files into `./data/modules/BGSU`. Note that only the latest Zip is required.

### Install and train BayesPairing
To use Bayespairing, you need to install it on the host machine to train it (build bayesian networks for every RNA motif in your database). It is already installed in the Docker container, but not trained, so you need to train it on your data and tell docker to mount it (see the docker run command below).

Make sure you have Python 3.5+ with packages networkx, numpy, regex, wrapt and biopython. You can install them with pip, on Linux you will need the python3-dev package to build them.

On Windows, on Mac : script coming soon

On Linux:
```
$ sudo -H pip3 install --upgrade pip
$ sudo -H pip3 install networkx numpy regex wrapt biopython
$ git clone http://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing.git BayesPairing
$ cd BayesPairing
$ sudo -H pip3 install .

$ cd bayespairing/src
$ python3 parse_sequences.py -d rna3dmotif -seq ACACGGGGUAAGAGCUGAACGCAUCUAAGCUCGAAACCCACUUGGAAAAGAGACACCGCCGAGGUCCCGCGUACAAGACGCGGUCGAUAGACUCGGGGUGUGCGCGUCGAGGUAACGAGACGUUAAGCCCACGAGCACUAACAGACCAAAGCCAUCAU -ss ".................................................................((...............)xxxx(...................................................)xxx).............."
$ python3 parse_sequences.py -d 3dmotifatlas -seq ACACGGGGUAAGAGCUGAACGCAUCUAAGCUCGAAACCCACUUGGAAAAGAGACACCGCCGAGGUCCCGCGUACAAGACGCGGUCGAUAGACUCGGGGUGUGCGCGUCGAGGUAACGAGACGUUAAGCCCACGAGCACUAACAGACCAAAGCCAUCAU -ss ".................................................................((...............)xxxx(...................................................)xxx).............."
```
The training is quite long, but has to be run only once.

### Download the docker image
`docker pull persalteas/biorseo:latest`

### Run the docker image
Use the following command to run the docker image:
```
$ docker run 
-v `pwd`/data/modules:/modules 
-v `pwd`/BayesPairing/bayespairing:/byp 
-v `pwd`/results:/biorseo/results 
biorseo 
yourexamplejobcommandhere
```
You can replace \`pwd\` by the full path of the biorseo/ root folder. Here we launch the biorseo image with 3 volumes : A first to give BiORSEO access to the module files, a second to give it access to your trained BayesPairing, and a third for it to output the result files of your job. An example job command can be ` ./biorseo.py -i /biorseo/data/fasta/applications.fa --rna3dmotifs --patternmatch --func B`, so the full run command would be 
```
$ docker run -v `pwd`/data/modules:/modules -v `pwd`/BayesPairing/bayespairing:/byp -v `pwd`/results:/biorseo/results biorseo ./biorseo.py -i /biorseo/data/fasta/applications.fa --rna3dmotifs --patternmatch --func B
```

Option 2 : Compile and Install from source (without docker, Linux only)
==================================

### CLONING
* Clone this git repository : `git clone https://github.com/persalteas/biorseo.git` and `cd biorseo`.

* Create folders for the modules you will use: `mkdir -p data/modules/`. If you plan to use several module sources, add subdirectories : `mkdir -p data/modules/DESC` and `mkdir -p data/modules/BGSU`

### RNA3DMOTIFS DATA

If you use Rna3Dmotifs, you need to get RNA-MoIP's .DESC dataset: download it from [GitHub](https://github.com/McGill-CSB/RNAMoIP/blob/master/CATALOGUE.tgz). Put all the .desc from the `Non_Redundant_DESC` folder into `./data/modules/DESC`. Otherwise, you also can run Rna3Dmotifs' `catalog` program to get your own DESC modules collection from updated 3D data (download [Rna3Dmotifs](https://rna3dmotif.lri.fr/Rna3Dmotif.tgz)). You also need to move the final DESC files into `./data/modules/DESC`.

### THE RNA 3D MOTIF ATLAS DATA

Get the latest version of the HL and IL module models from the [BGSU website](http://rna.bgsu.edu/data/jar3d/models/) and extract the Zip files. Put the HL and IL folders from inside the Zip files into `./data/modules/BGSU`. Note that only the latest Zip is required.


### DEPENDENCIES
- Make sure you have Python 3.5+, Cmake, and a C++ compiler installed on your distribution. Please, it's 2019, use a recent one, we use the 2017 C++ standard. The compilation will not work with Ubuntu 16's GCC 5.4 for example. Tested with libstdc++-dev >= 6.0, so use GCC >=6.0 or Clang >= 6.0.
- Install automake, libboost-program-options and libboost-filesystem.
- Download and install [IBM ILOG Cplex optimization studio](https://www.ibm.com/analytics/cplex-optimizer), through their academic initiative. The Student and the Community Edition versions are fine, but the free version is too limited. registering as academic is free.
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
* You might want to edit `Makefile`:
    - if you are not using clang as compiler. For example, if you use g++, replace clang++ by g++.
    - if you did not install CPLEX, Eigen or Nupack in the default locations. Please update the top variables $ICONCERT, $ICPLEX, $INUPACK, $IEIGEN, $LCONCERT, and $LCPLEX with the correct locations.
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

### RUN BIORSEO
Now you can run biorseo.py, but, as you are not into the Docker environment, you MUST provide the options to tell it the jar3d or BayesPairing locations, for example:
```
$ ./biorseo.py 
-i /biorseo/data/fasta/applications.fa 
--rna3dmotifs --patternmatch --func B 
--modules-path=./data/modules/DESC 
--jar3dexec=./jar3d_releasedate.jar OR --bypdir=./BayesPairing/bayespairing/src
```