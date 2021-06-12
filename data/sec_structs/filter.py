#!/usr/bin/python3

# This script filters a DBN file to keep only instances
# between 10 and 100 nucleotides, and at least one basepair.

# Pseudobase++ : 398 --> 351
# bpRNA-1m_90 : 28370 --> 13716
import sys

db = open(sys.argv[1], "r")
all = db.readlines()
db.close()
newlines = []

for line in all:
    if line[0] == ">":
        header = line
        c = 0
        keepthis = False
    else:
        c += 1
    if c == 1:
        if len(line) > 11 and len(line) < 101:
            keepthis = True
            seq = line
    if c == 2 and keepthis and '(' in line:
        newlines.append(header)
        newlines.append(seq)
        newlines.append(line)

db = open(sys.argv[1].strip(".dbn")+"_short.dbn", "w")
db.writelines(newlines)
db.close()