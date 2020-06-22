#!/usr/bin/python3
#coding=utf-8

# typical usage : ./benchmark.py data/sec_structs/verified_secondary_structures_database.dbn data/sec_structs/pseudoknots.dbn data/sec_structs/applications.dbn

# the .dbn files should be formatted the following way:
# > header of the sequence (somecode)
# ACGUACGUACGUACGUACGU
# ...(((...((...))))).
# > header of the next sequence (somecode2)
# ACGUACGUACGGGCGUACGU
# ...(((..........))).



# ============================ IMPORTS ====================================
from os import path
from sys import argv
import subprocess


# ================== DEFINITION OF THE PATHS ==============================
biorseoDir = path.realpath(".")
jar3dexec = "/opt/jar3d_2014-12-11.jar"
bypdir = "/opt/BayesPairing/bayespairing/src"
moipdir = "/opt/RNAMoIP/Src/RNAMoIP.py"
biokopdir = "/opt/biokop"
runDir = path.dirname(path.realpath(__file__))

RNAStrandFile = argv[1]
PseudobaseFile = argv[2]
StudyCaseFile = argv[3]

outputDir = biorseoDir + "/benchmark_results/"
HLmotifDir = biorseoDir + "/data/modules/BGSU/HL/3.2/lib"
ILmotifDir = biorseoDir + "/data/modules/BGSU/IL/3.2/lib"
descfolder = biorseoDir + "/data/modules/DESC"

# Create some folders to store the results
subprocess.call(["mkdir", "-p", outputDir])
subprocess.call(["mkdir", "-p", outputDir + "PK/"])
subprocess.call(["mkdir", "-p", outputDir + "noPK/"])