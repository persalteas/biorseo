FROM ubuntu:bionic

# installing dependencies
RUN apt-get update -yq && \
    apt-get upgrade -y && \
    apt-get install -y build-essential make python3-dev python3-pip libmpfr-dev libgsl-dev libboost-program-options-dev libboost-filesystem-dev openjdk-11-jre && \
    rm -rf /var/lib/apt/lists/*

# compiled biorseo
COPY . /biorseo 
# uncompiled ViennaRNA
ADD "https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.13.tar.gz" /
# jar3d archive
ADD http://rna.bgsu.edu/data/jar3d/models/jar3d_2014-12-11.jar /

# install codes
RUN tar -xvzf /ViennaRNA-2.4.13.tar.gz && \
    cd /ViennaRNA-2.4.13        && \
    ./configure                 && \
    make -j 15                  && \
    make install                && \
    \
    pip3 install networkx numpy regex wrapt biopython /biorseo/BayesPairing && \
    \
    cd / && \
    rm -rf /biorseo/BayesPairing /ViennaRNA-2.4.13 /ViennaRNA-2.4.13.tar.gz
WORKDIR /biorseo