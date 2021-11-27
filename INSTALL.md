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

### Download the docker image from Docker Hub
`docker pull persalteas/biorseo:latest`

### Run the docker image
Use the following command to run the docker image:
```
$ docker run 
-v `pwd`/data/modules:/modules 
-v `pwd`/data/fasta:/biorseo/data/fasta
-v `pwd`/results:/biorseo/results 
persalteas/biorseo 
yourexamplejobcommandhere
```
You can replace \`pwd\` by the full path of the biorseo/ root folder. Here we launch the biorseo image with 4 volumes : A first to give BiORSEO access to the module files, a second to give it access to your input file(s), a third for your trained BayesPairing, and a last for it to output the result files of your job. Considering you place your input file 'MyFastaFile.fa' into the `data/fasta` folder, an example job command can be ` ./biorseo.py -i /biorseo/data/fasta/myFastaFile.fa  --rna3dmotifs --patternmatch --func B`, so the full run command would be 
```
$ docker run -v `pwd`/data/modules:/modules -v `pwd`/data/fasta:/biorseo/data/fasta -v `pwd`/results:/biorseo/results persalteas/biorseo ./biorseo.py -i /biorseo/data/fasta/applications.fa --rna3dmotifs --patternmatch --func B
```

Note that the paths to the input and output files are paths *inside the Docker container*, and those paths are mounted to folders of the host machine with -v options.

Option 2 : Compile and Install from source (without docker, Linux only)
==================================

### CLONING
* Clone this git repository : `git clone https://forge.ibisc.univ-evry.fr/lbecquey/biorseo.git` (from the IBISC forge) or `git clone https://github.com/persalteas/biorseo.git` (from my personal GitHub, only while i am the current developer !) and `cd biorseo`.

* Create folders for the modules you will use: `mkdir -p data/modules/`. If you plan to use several module sources, add subdirectories :
```bash
mkdir -p data/modules/BGSU
mkdir -p data/modules/RIN
mkdir -p data/modules/DESC
mkdir -p data/modules/JSON
```

### THE RNA 3D MOTIF ATLAS DATA

Get the latest version of the HL and IL module models from the [BGSU website](http://rna.bgsu.edu/data/jar3d/models/) and extract the Zip files. Put the HL and IL folders from inside the Zip files into `./data/modules/BGSU`. Note that only the latest Zip is required.

### CARNAVAL DATA

You first need to have the `unzip` command installed on your machine and the `networkx` package installed for Python 3. Then just run the script `Install_CaRNAval_RINs.py`, this will create files into `./data/modules/RIN/Subfiles` :
```bash
cd scripts
python3 Install_CaRNAval_RINs.py
```
If you do not have the unzip command, download and extract manually the [CaRNAval dataset](http://carnaval.lri.fr/carnaval_dataset.zip) and place the files `RIN.py` and `CaRNAval_1_as_dictionnary.nxpickled` in the folder `data/modules/RIN/`, and run the python script.

### RNA3DMOTIFS DATA (DEPRECATED)

If you use Rna3Dmotifs, you need to get RNA-MoIP's .DESC dataset: download it from [GitHub](https://github.com/McGill-CSB/RNAMoIP/blob/master/CATALOGUE.tgz). Put all the .desc from the `Non_Redundant_DESC` folder into `./data/modules/DESC`. Otherwise, you also can run Rna3Dmotifs' `catalog` program to get your own DESC modules collection from updated 3D data (download [Rna3Dmotifs](https://rna3dmotif.lri.fr/Rna3Dmotif.tgz)). You also need to move the final DESC files into `./data/modules/DESC`.

### DEPENDENCIES
- Make sure you have Python 3.7+ and a C++ compiler (tested with GCC and clang) installed on your distribution. Use a recent one, we use the 2017 C++ standard. The compilation will not work with Ubuntu 16's GCC 5.4 for example.
- Install automake, libeigen3-dev, libboost-program-options-dev and libboost-filesystem-dev, or equivalent packages in your distribution (Eigen 3 and Boost headers).
- Download and install the [ViennaRNA package](https://www.tbi.univie.ac.at/RNA/). We use their libRNA C++ API to build Biorseo.
    - If you use the pre-complied packages, you need the "Core" package and the development files.
    - If you compile it from source, extract the archive (`tar -xvzf ViennaRNA-2.5.0.tar.gz`), go into the folder (`cd ViennaRNA-2.5.0`), configure and build it (`./configure` and then `make -j 4`) and finally install it with root permissions (`sudo make install`). Everything should be fine. This takes ~15 min. Optionnally, you may want to install optional libraries GSL and MPFR before (libgsl-dev and libmpfr-dev on Ubuntu), for better performance.
- Download and install [IBM ILOG Cplex optimization studio](https://www.ibm.com/analytics/cplex-optimizer), through their [academic initiative](https://www.ibm.com/academic/home). The Student and the Community Edition versions are fine, but the free version is too limited. Registering as academic is free. We actually don't use the Studio, but we use the development files to compile Biorseo.

### BUILDING
* You might want to edit `Makefile` if you did not install CPLEX or the libRNA in the default location. Please update the top variables $ICONCERT, $ICPLEX, $LCONCERT, and $LCPLEX with the correct locations.
* Build it: `make -j4`
* Check if the executable file exists: `./bin/biorseo --version`.

### RUN BIORSEO
Now you can run biorseo.py, but, as you are not into the Docker environment, you MUST provide the options to tell it the jar3d or BayesPairing locations, for example:
```
$ ./biorseo.py 
-i ./data/fasta/applications.fa 
-O ./results/
--rna3dmotifs --patternmatch --func B 
--biorseo-dir /FULL/path/to/the/root/biorseo/dir
--modules-path=./data/modules/DESC 
```
