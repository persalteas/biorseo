
#!/bin/bash
######################################################## RNA modules ##############################################################

cd ../

# Rna3Dmotifs data
mkdir -p data/modules/DESC
wget https://github.com/McGill-CSB/RNAMoIP/raw/master/CATALOGUE.tgz
tar -xvzf CATALOGUE.tgz 
mv No_Redondance_DESC/*.desc data/modules/DESC/
rm -r No_Redondance_VIEW3D No_Redondance_DESC CATALOGUE.tgz

# The RNA 3D Motif Atlas
mkdir -p data/modules/BGSU
wget http://rna.bgsu.edu/data/jar3d/models/HL/HL_3.2_models.zip
unzip HL_3.2_models.zip
mv HL data/modules/BGSU
rm HL_3.2_models.zip
wget http://rna.bgsu.edu/data/jar3d/models/IL/IL_3.2_models.zip
unzip IL_3.2_models.zip
mv IL data/modules/BGSU
rm IL_3.2_models.zip

# Install BayesPairing
sudo -H pip3 install --upgrade pip
sudo -H pip3 install networkx numpy regex wrapt biopython
git clone http://jwgitlab.cs.mcgill.ca/sarrazin/rnabayespairing.git BayesPairing
cd BayesPairing
sudo -H pip3 install .

# Train Bayes Pairing (it has been installed on the image and the source has been deleted, we train the models now, and will remount it as volume at run time)
cd bayespairing/src
python3 parse_sequences.py -d rna3dmotif -seq ACACGGGGUAAGAGCUGAACGCAUCUAAGCUCGAAACCCACUUGGAAAAGAGACACCGCCGAGGUCCCGCGUACAAGACGCGGUCGAUAGACUCGGGGUGUGCGCGUCGAGGUAACGAGACGUUAAGCCCACGAGCACUAACAGACCAAAGCCAUCAU -ss ".................................................................((...............)xxxx(...................................................)xxx).............."
python3 parse_sequences.py -d 3dmotifatlas -seq ACACGGGGUAAGAGCUGAACGCAUCUAAGCUCGAAACCCACUUGGAAAAGAGACACCGCCGAGGUCCCGCGUACAAGACGCGGUCGAUAGACUCGGGGUGUGCGCGUCGAGGUAACGAGACGUUAAGCCCACGAGCACUAACAGACCAAAGCCAUCAU -ss ".................................................................((...............)xxxx(...................................................)xxx).............."
cd ../../..

######################################################## Run it ##############################################################

# docker run -v `pwd`/data/modules:/modules -v `pwd`/BayesPairing/bayespairing:/byp -v `pwd`/results:/biorseo/results biorseo ./biorseo.py -i /biorseo/data/fasta/applications.fa --rna3dmotifs --patternmatch --func B