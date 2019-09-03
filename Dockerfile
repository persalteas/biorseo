FROM ubuntu:bionic

# installing dependencies
RUN apt-get update -yq && \
    apt-get upgrade -y && \
    apt-get install -y python3-dev python3-pip openjdk-11-jre libgsl23 libgslcblas0 libboost-program-options-dev libboost-filesystem-dev && \
    rm -rf /var/lib/apt/lists/*

# compiled biorseo
COPY . /biorseo 
# ViennaRNA installer
ADD "https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_18_04/viennarna_2.4.14-1_amd64.deb" /
# jar3d archive
ADD http://rna.bgsu.edu/data/jar3d/models/jar3d_2014-12-11.jar /

# install codes
RUN dpkg -i /viennarna_2.4.14-1_amd64.deb && \
    apt-get install -f          && \
    \
    pip3 install networkx numpy regex wrapt biopython /biorseo/BayesPairing && \
    \
    cd / && \
    rm -rf /biorseo/BayesPairing /ViennaRNA-2.4.13 /ViennaRNA-2.4.13.tar.gz
WORKDIR /biorseo