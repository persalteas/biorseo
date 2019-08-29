#!/bin/bash
THISDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

####################################################### Dependencies ##############################################################
sudo apt install -y unzip clang-7 cmake automake python3-dev libboost-program-options-dev libboost-filesystem-dev openjdk-11-jre libgsl23
sudo update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-7 100
sudo update-alternatives --install /usr/bin/clang clang /usr/bin/clang-7 100

# CPLEX: only to build biorseo, or needed at runtime ?
wget --no-check-certificate -O cplex_installer_12.8_Student.bin "https://onedrive.live.com/download?cid=C9F6894F7941BBFB&resid=C9F6894F7941BBFB%21416367&authkey=AC-d-946knh0Yo0"
chmod +x cplex_installer_12.8_Student.bin
printf "4\n\n1\n\n\n\n\n" | sudo ./cplex_installer_12.8_Student.bin

# Eigen: only to build biorseo (no need to give it to the docker image)
wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz -O eigen_src.tar.gz
tar -xf eigen_src.tar.gz
cd eigen-eigen-323c052e1731
mkdir build
cd build
cmake ..
sudo make install
cd ../..
rm -r eigen_src.tar.gz 

# ViennaRNA: needed for runtime. Build here on host, install later in image by the Dockerfile.
wget "https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.13.tar.gz"
tar -zxvf ViennaRNA-2.4.13.tar.gz
cd ViennaRNA-2.4.13
./configure 
make -j 4 
cd ..
rm ViennaRNA-2.4.13.tar.gz

# Nupack: only to build biorseo (no need to give it to the docker image)
curl -u louis.becquey@univ-evry.fr:tQm03G29 http://www.nupack.org/downloads/serve_file/nupack3.2.2.tar.gz --output nupack3.2.2.tar.gz
tar -xf nupack3.2.2.tar.gz
cd nupack3.2.2
mkdir build
cd build
cmake ..
make -j4
sudo make install
cd ../..
rm nupack3.2.2.tar.gz
sudo cp nupack3.2.2/src/thermo/*.h /usr/local/include/nupack/thermo/

# BayesPairing
sudo -H pip3 install --upgrade pip
sudo -H pip3 install networkx numpy regex wrapt biopython
git clone http://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing.git BayesPairing
cd BayesPairing
sudo -H pip3 install .
cd ..


######################################################### Build Biorseo ###########################################################
git clone https://github.com/persalteas/biorseo.git
cd biorseo
mkdir -p results
make -j 4
cd ..

######################################################## RNA modules ##############################################################

# Rna3Dmotifs data
mkdir -p modules/DESC
wget https://github.com/McGill-CSB/RNAMoIP/raw/master/CATALOGUE.tgz
tar -xvzf CATALOGUE.tgz 
mv No_Redondance_DESC/*.desc modules/DESC/
rm -r No_Redondance_VIEW3D No_Redondance_DESC CATALOGUE.tgz

# The RNA 3D Motif Atlas
mkdir -p modules/BGSU
wget http://rna.bgsu.edu/data/jar3d/models/HL/HL_3.2_models.zip
unzip HL_3.2_models.zip
mv HL modules/BGSU
rm HL_3.2_models.zip
wget http://rna.bgsu.edu/data/jar3d/models/IL/IL_3.2_models.zip
unzip IL_3.2_models.zip
mv IL modules/BGSU
rm IL_3.2_models.zip

######################################################## Build Docker container ##################################################
docker build -t biorseo
#  docker run -v `pwd`/modules:/modules -v `pwd`/BayesPairing/bayespairing:/byp -v `pwd`/results:/biorseo/results biorseo ls /byp/models

# Train Bayes Pairing
cd bayespairing/src
python3 parse_sequences.py -d rna3dmotif -seq ACACGGGGUAAGAGCUGAACGCAUCUAAGCUCGAAACCCACUUGGAAAAGAGACACCGCCGAGGUCCCGCGUACAAGACGCGGUCGAUAGACUCGGGGUGUGCGCGUCGAGGUAACGAGACGUUAAGCCCACGAGCACUAACAGACCAAAGCCAUCAU -ss ".................................................................((...............)xxxx(...................................................)xxx).............."
python3 parse_sequences.py -d 3dmotifatlas -seq ACACGGGGUAAGAGCUGAACGCAUCUAAGCUCGAAACCCACUUGGAAAAGAGACACCGCCGAGGUCCCGCGUACAAGACGCGGUCGAUAGACUCGGGGUGUGCGCGUCGAGGUAACGAGACGUUAAGCCCACGAGCACUAACAGACCAAAGCCAUCAU -ss ".................................................................((...............)xxxx(...................................................)xxx).............."
cd ../../..




########################################################### More stuff ###########################################################

# # BiokoP
# mkdir biokop
# cd biokop
# wget --no-check-certificate "https://tanuki.ibisc.univ-evry.fr/media/filer_public/b5/ea/b5eadd8b-ebc6-49f7-821e-b1e30137cc21/biokop.tar"
# tar -xvf biokop.tar
# rm biokop.tar
# cd ..

# # RNA examples
# mkdir data
# echo ">__'CRYSTAL_STRUCTURE_OF_A_TIGHT-BINDING_GLUTAMINE_TRNA_BOUND_TO_GLUTAMINE_AMINOACYL_TRNA_SYNTHETASE_'_(PDB_00376)" > data/tRNA.fa
# echo "GGGGUAUCGCCAAGCGGUAAGGCACCGGAUUCUGAUUCCGGAGGUCGAGGUUCGAAUCCUCGUACCCCAGCCA" >> data/tRNA.fa
# echo ">__'GUANINE_RIBOSWITCH_U22C,_A52G_MUTANT_BOUND_TO_HYPOXANTHINE_'_(PDB_01023)" > data/Griboswitch.fa
# echo "GGACAUACAAUCGCGUGGAUAUGGCACGCAAGUUUCUGCCGGGCACCGUAAAUGUCCGACUAUGUCCA" >> data/Griboswitch.fa
# echo ">__'SOLUTION_STRUCTURE_OF_THE_P2B-P3_PSEUDOKNOT_FROM_HUMAN_TELOMERASE_RNA_'_(PDB_00857)" > data/pseudoknot.fa
# echo "GGGCUGUUUUUCUCGCUGACUUUCAGCCCCAAACAAAAAAGUCAGCA" >> data/pseudoknot.fa
# mv biokop/PKB120.fasta data/

