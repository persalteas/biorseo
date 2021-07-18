from math import sqrt, ceil
import numpy as np
import matplotlib.pyplot as plt
import re


def get_result(filename):
    ext = "json_pmE"
    file2 = open( "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/results/" + filename + ext, "r")

    name = file2.readline()
    rna = file2.readline()
    twod = file2.readline()
    pred = re.findall(r'\S+', twod)

    score = '-' + pred[len(pred)-1]
    min = float(score)
    contacts = file2.readline()
    while twod:
        twod = file2.readline()
        pred = re.findall(r'\S+', twod)
        if len(pred) > 0:
            score = '-' + pred[len(pred) - 1]
            if float(score) < min:
                min = float(score)
        contacts = file2.readline()
    file2.close()
    return min

file = open("/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/benchmark.dbn", "r")
name = file.readline().strip()
rna = file.readline()
twod = file.readline()
contacts = file.readline()
while name:
    print(name)
    best_mea = get_result(name)
    print(str(best_mea) + "\n")
    name = file.readline().strip()
    rna = file.readline()
    twod = file.readline()
    contacts = file.readline()

file.close()
