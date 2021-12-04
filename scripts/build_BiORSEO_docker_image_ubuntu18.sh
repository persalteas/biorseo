#!/bin/bash

echo "WARNING: The purpose of this file is to document how the docker image was built.";
echo "You cannot execute it directly, because of licensing reasons. Please get your own:";
echo "- CPLEX academic version: cplex_installer_12.8_Student.bin";
exit 0;

cd ../
THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

####################################################### Dependencies ##############################################################
sudo apt install -y make automake libgsl-dev libmpfr-dev libeigen3-dev libboost-program-options-dev libboost-filesystem-dev

# CPLEX: only to build biorseo
# HERE YOU SHOULD GET YOUR OWN cplex_installer_12.8_Student.bin ! I am not allowed to share mine anymore.
chmod +x cplex_installer_12.8_Student.bin
printf "4\n\n1\n\n\n\n\n" | sudo ./cplex_installer_12.8_Student.bin
rm cplex_installer_12.8_Student.bin

# ViennaRNA (to build Biorseo with libRNA)
wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_5_x/ViennaRNA-2.5.0.tar.gz
tar xzf ViennaRNA-2.5.0.tar.gz
cd ViennaRNA-2.5.0
./configure
make -j 8
sudo make install

######################################################### Build Biorseo ###########################################################
# build here, install later on the docker image (done by the Dockerfile)
mkdir -p results
make -j 8
make clean
rm -rf obj/ figures/

######################################################## Build Docker container ##################################################
# Execute the Dockerfile and build the image
docker build . -t biorseo
