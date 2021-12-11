# You can pick the Ubuntu version that suits you instead, according to the version of the boost libraries
# that you are using to compile biorseo.
#
# Typically, on the machine where you typed 'make', check :
# ls /usr/lib/libboost_filesystem.so.*
# this will give you the file name of your boost library, including the version number.
# Use the docker basis image of the Ubuntu which has this version of boost in the apt sources.
FROM ubuntu:focal

# compiled biorseo
COPY ./bin /workdir/

# Install runtime dependencies
RUN apt-get update -yq && \
    apt-get upgrade -y && \
    apt-get install -y libboost-program-options-dev libboost-filesystem-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /workdir
ENTRYPOINT ["/workdir/biorseo"]
