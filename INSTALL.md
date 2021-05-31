Option 1 : Installation using docker image (Windows, Mac, Linux)
==================================
* Clone this git repository : `git clone https://github.com/persalteas/biorseo.git` , or download the .zip archive from a BiORSEO release and extract it.
* Move into the repository ( `cd biorseo` )

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

### Download the docker image from Docker Hub
`docker pull persalteas/biorseo:latest`

### Run the docker image
Use the following command to run the docker image:
```
$ docker run 
-v `pwd`/data/modules:/modules 
-v `pwd`/data/fasta:/biorseo/data/fasta
-v `pwd`/BayesPairing/bayespairing:/byp 
-v `pwd`/results:/biorseo/results 
persalteas/biorseo 
yourexamplejobcommandhere
```
You can replace \`pwd\` by the full path of the biorseo/ root folder. Here we launch the biorseo image with 4 volumes : A first to give BiORSEO access to the module files, a second to give it access to your input file(s), a third for your trained BayesPairing, and a last for it to output the result files of your job. Considering you place your input file 'MyFastaFile.fa' into the `data/fasta` folder, an example job command can be ` ./biorseo.py -i /biorseo/data/fasta/myFastaFile.fa  --rna3dmotifs --patternmatch --func B`, so the full run command would be 
```
$ docker run -v `pwd`/data/modules:/modules -v `pwd`/data/fasta:/biorseo/data/fasta -v `pwd`/BayesPairing/bayespairing:/byp -v `pwd`/results:/biorseo/results persalteas/biorseo ./biorseo.py -i /biorseo/data/fasta/applications.fa --rna3dmotifs --patternmatch --func B
```

Note that the paths to the input and output files are paths *inside the Docker container*, and those paths are mounted to folders of the host machine with -v options.

Option 2 : Compile and Install from source (without docker, Linux only)
==================================

### CLONING
* Clone this git repository : `git clone https://forge.ibisc.univ-evry.fr/lbecquey/biorseo.git` (from the IBISC forge) or `git clone https://github.com/persalteas/biorseo.git` (from my personal GitHub, only while i am the current developer !) and `cd biorseo`.

* Create folders for the modules you will use: `mkdir -p data/modules/`. If you plan to use several module sources, add subdirectories :
```bash
mkdir -p data/modules/DESC
mkdir -p data/modules/BGSU
mkdir -p data/modules/RIN
mkdir -p data/modules/ISAURE/Motifs_derniere_version
```

### RNA3DMOTIFS DATA

If you use Rna3Dmotifs, you need to get RNA-MoIP's .DESC dataset: download it from [GitHub](https://github.com/McGill-CSB/RNAMoIP/blob/master/CATALOGUE.tgz). Put all the .desc from the `Non_Redundant_DESC` folder into `./data/modules/DESC`. Otherwise, you also can run Rna3Dmotifs' `catalog` program to get your own DESC modules collection from updated 3D data (download [Rna3Dmotifs](https://rna3dmotif.lri.fr/Rna3Dmotif.tgz)). You also need to move the final DESC files into `./data/modules/DESC`.

### THE RNA 3D MOTIF ATLAS DATA

Get the latest version of the HL and IL module models from the [BGSU website](http://rna.bgsu.edu/data/jar3d/models/) and extract the Zip files. Put the HL and IL folders from inside the Zip files into `./data/modules/BGSU`. Note that only the latest Zip is required.

### CARNAVAL DATA

You first need to have the `unzip` command installed on your machine and the `networkx` package installed for Python 3. Then just run the script `Install_CaRNAval_RINs.py`, this will create files into `./data/modules/RIN/Subfiles` :
```bash
cd scripts
python3 Install_CaRNAval_RINs.py
```
If you do not have the unzip command, download and extract manually the [CaRNAval dataset](http://carnaval.lri.fr/carnaval_dataset.zip) and place the files `RIN.py` and `CaRNAval_1_as_dictionnary.nxpickled` in the folder `data/modules/RIN/`, and run the python script.

### CONTACTS DATA

If you use contacts, you need to put the motifs.json of Isaure in `data/modules/ISAURE`.

### DEPENDENCIES
- Make sure you have Python 3.7+ and a C++ compiler (tested with GCC and clang) installed on your distribution. Use a recent one, we use the 2017 C++ standard. The compilation will not work with Ubuntu 16's GCC 5.4 for example.
- Install automake, libeigen3-dev, libboost-program-options-dev and libboost-filesystem-dev, or equivalent packages in your distribution (Eigen 3 and Boost headers).
- Download and install the [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/). We are interested in the libRNA library to build Biorseo first, and then, Jar3D and BayesPairing may use some of the executables (RNAsubopt, RNAfold). 
    - If you use the pre-complied packages, you need both the "Core" package and the development files. 
    - If you compile it from source, extract the archive (`tar -xvzf ViennaRNA-2.4.15.tar.gz`), go into the folder (`cd ViennaRNA-2.4.15`), configure and build it (`./configure` and then `make -j 4`) and finally install it with root permissions (`sudo make install`). Everything should be fine. This takes ~15 min. Optionnally, you may want to install optional libraries GSL and MPFR before (libgsl-dev and libmpfr-dev on Ubuntu), for better performance.
- Download and install [IBM ILOG Cplex optimization studio](https://www.ibm.com/analytics/cplex-optimizer), through their [academic initiative](https://www.ibm.com/academic/home). The Student and the Community Edition versions are fine, but the free version is too limited. Registering as academic is free. We actually don't use the Studio, but we use the development files to compile Biorseo.


### OPTIONAL DEPENDENCIES FOR USE OF JAR3D
- Download and install Java runtime (Tested with Java 14)
- Download the latest JAR3D executable "*jar3d_releasedate.jar*" from [the BGSU website](http://rna.bgsu.edu/data/jar3d/models/). 
  

### OPTIONAL DEPENDENCIES FOR USE OF BAYESPAIRING
- Make sure you can `import RNA` in your favorite python3.X. Otherwise, there might be a ViennaRNA installation issue.
    - Depending on your distribution, this might require that you add the /usr/local/lib/python3.X/site-packages/ folder, or the libRNA.so library, to your $PYTHONPATH variable. To do so, edit your `~/.bashrc` file and write `export PYTHONPATH="$PYTHONPATH:/usr/local/lib/python3.X/site-packages"` at the end of the file, replacing 3.X by your actual version of Python. Then close and re-open the terminal.
    - On Ubuntu/Debian distros, which have stupid conventions and ignore the site-packages/ folder, you need to move the python package to the dist-packages/ folder. E.g., if you compiled ViennaRNA from source, `sudo mv /usr/local/lib/python3.X/site-packages/RNA /usr/local/lib/python3.X/dist-packages/`.

- Make sure you have Python 3.6+ with pip (packages python3-pip and python3-dev on most distros)
- Clone the latest BayesPairing 2 Git repo, and install it : 
```
git clone http://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing2.git BayesPairing2
cd BayesPairing2
python3 -m pip install .
cd ..
```

### BUILDING
* You might want to edit `Makefile` if you did not install CPLEX or the libRNA in the default location. Please update the top variables $ICONCERT, $ICPLEX, $LCONCERT, and $LCPLEX with the correct locations.
* Build it: `make -j4`
* Check if the executable file exists: `./bin/biorseo --version`.

### BAYESPAIRING USERS: PREPARE BAYESIAN NETWORKS
We run an example job for it to build the bayesian networks of our modules.
```
cd BayesPairing2/bayespairing/src
python3 parse_sequences.py -seq "UUUUUUAAGGAAGAUCUGGCCUUCCCACAAGGGAAGGCCAAAGAAUUUCCUU" -t 4 -d rna3dmotif -ss "......(((((((.((((((((((((....)))))))))..))).)))))))"
```
Use `-d rna3dmotif`, `-d 3dMotifAtlas_ALL` depending on the module source you are planning to use.
This is a quite long step, but the bayesian networks will be ready for all the future uses.

### RUN BIORSEO
Now you can run biorseo.py, but, as you are not into the Docker environment, you MUST provide the options to tell it the jar3d or BayesPairing locations, for example:
```
$ ./biorseo.py 
-i ./data/fasta/applications.fa 
-O ./results/
--rna3dmotifs --patternmatch --func B 
--biorseo-dir /FULL/path/to/the/root/biorseo/dir
--modules-path=./data/modules/DESC 
--jar3d-exec=./jar3d_releasedate.jar OR --bypdir=./BayesPairing2/bayespairing/src
```
