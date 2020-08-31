#!/bin/bash

echo "WARNING: The purpose of this file is to document how the docker image was built.";
echo "You cannot execute it directly, because of licensing reasons. Please get your own:";
echo "- CPLEX academic version: cplex_installer_12.8_Student.bin";
echo "- Nupack header files: nupack_3.2.2.tar.gz";
exit 0;

cd ../
THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

####################################################### Dependencies ##############################################################
sudo apt install -y clang-7 cmake make automake libboost-program-options-dev libboost-filesystem-dev openjdk-11-jre
sudo update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-7 100
sudo update-alternatives --install /usr/bin/clang clang /usr/bin/clang-7 100

# CPLEX: only to build biorseo
# HERE YOU SHOULD GET YOUR OWN cplex_installer_12.8_Student.bin ! I am not allowed to share mine anymore.
chmod +x cplex_installer_12.8_Student.bin
printf "4\n\n1\n\n\n\n\n" | sudo ./cplex_installer_12.8_Student.bin
rm cplex_installer_12.8_Student.bin

# Eigen: only to build biorseo (no need to give it to the docker image)
wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz -O eigen_src.tar.gz
tar -xf eigen_src.tar.gz
cd eigen-eigen-323c052e1731
mkdir build
cd build
cmake ..
sudo make install
cd ../..
rm -rf eigen_src.tar.gz eigen-eigen-323c052e1731

# Nupack: only to build biorseo (no need to give it to the docker image)
#curl -u yourname@yourUni.com:yourPassword http://www.nupack.org/downloads/serve_file/nupack3.2.2.tar.gz --output nupack3.2.2.tar.gz
tar -xf nupack3.2.2.tar.gz
cd nupack3.2.2
mkdir build
cd build
cmake ..
make -j8
sudo make install
cd ../..
sudo cp nupack3.2.2/src/thermo/*.h /usr/local/include/nupack/thermo/
rm -rf nupack3.2.2.tar.gz nupack3.2.2/

# BayesPairing: install on the docker image (done by the Dockerfile)
git clone http://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing.git BayesPairing

######################################################### Build Biorseo ###########################################################
# build here, install later on the docker image (done by the Dockerfile)
mkdir -p results
make -j 8
make clean
rm -rf doc/ obj/

######################################################## Build Docker container ##################################################
# Execute the Dockerfile and build the image
docker build . -t biorseo
