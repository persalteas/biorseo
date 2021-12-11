Option 1 : Installation using docker image (Windows, Mac, Linux)
==================================

### Install Docker:
* See the officiel instructions depending on your OS here : https://docs.docker.com/install/
* Example for Ubuntu Linux:
```
sudo apt-get install apt-transport-https ca-certificates curl gnupg-agent software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
sudo apt update && sudo apt install docker-ce docker-ce-cli containerd.io
```

### Create a separated folder
You will need a folder that is mounted in the Docker container (shared between your filesystem and the container's one). For example, create a `data` folder. It should be the place where you place your FASTA input, your modules, and where Biorseo should place the output.
```
mkdir data 
```

### Download and install the RNA motifs data files:
* Move your JSON-formatted or CSV-formatted files containing motifs in a `./data/JSON` or `./data/CSV` folder.
* If you use Rna3Dmotifs (easy, but outdated), you need to get RNA-MoIP's .DESC dataset: download it from [GitHub](https://github.com/McGill-CSB/RNAMoIP/blob/master/CATALOGUE.tgz). Put all the .desc from the `Non_Redundant_DESC` folder into `./data/DESC`. Otherwise, you also can run Rna3Dmotifs' `catalog` program to get your own DESC modules collection from updated 3D data (download [Rna3Dmotifs](https://rna3dmotif.lri.fr/Rna3Dmotif.tgz)). You also need to move the final DESC files into `./data/DESC`.
* If you use CaRNAval (which is supposed to be a long-distance contact module dataset, not a SSE module dataset), use the script [scripts/Install_CaRNAval_RINs.py](scripts/Install_CaRNAval_RINs.py) : 
`python3 Install_CaRNAval_RINs.py`. This will create a `/../data/modules/RIN/` folder (because the script is supposed to be run from the repo's `scripts` subfolder)

### Download the docker image from Docker Hub
```
docker pull persalteas/biorseo:latest
```

### Run the docker image
Use the following command to run the docker image:
```
$ docker run -v `pwd`/data:/workdir/data persalteas/biorseo [optionshere]
```
You can replace \`pwd\`/data by the full path to your data folder. Assuming you place your input file 'MyFastaFile.fa' into the `data/fasta` folder, an example job command can be :
```
$ docker run -v `pwd`/data:/workdir/data persalteas/biorseo -s data/fasta/MyFastaFile.fa --descfolder data/DESC --func B -v -o data/MyOutput.biorseo
```

Note that the paths to the input and output files are paths *inside the Docker container*, and those paths are mounted to the data folder of the host machine with the -v option.

Option 2 : Compile and Install from source (without docker, Linux only)
==================================

### CLONING
* Clone this git repository : `git clone https://forge.ibisc.univ-evry.fr/lbecquey/biorseo.git` (from the IBISC forge) or `git clone https://github.com/persalteas/biorseo.git` (from my personal GitHub, only while i (Louis Becquey) am the current developer !) and `cd biorseo`.

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
Now you can run biorseo, but, as you are not into the Docker environment, you MUST provide the options to tell it the jar3d or BayesPairing locations, for example:
```
$ ./bin/biorseo
-s ./data/fasta/applications.fa 
-o result.bi
--func B 
--descfolder=./data/modules/DESC 
```
