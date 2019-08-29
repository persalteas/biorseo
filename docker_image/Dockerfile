FROM ubuntu:cosmic

# installing dependencies
RUN apt-get update -yq && \
    apt-get upgrade -y && \
    apt-get install -y make python3-dev python3-pip libmpfr-dev libgsl-dev libboost-program-options-dev libboost-filesystem-dev openjdk-11-jre && \
    rm -rf /var/lib/apt/lists/*


# import compiled code
COPY ./ViennaRNA-2.4.13 /ViennaRNA
COPY ./BayesPairing /BayesPairing
COPY ./biorseo /biorseo
ADD http://rna.bgsu.edu/data/jar3d/models/jar3d_2014-12-11.jar /

# install codes
RUN cd /ViennaRNA               && \
    make install                && \
    pip3 install --upgrade pip  && \
    pip3 install networkx numpy regex wrapt biopython /BayesPairing && \
    cd / && \
    rm -rf /BayesPairing /ViennaRNA 
WORKDIR /biorseo