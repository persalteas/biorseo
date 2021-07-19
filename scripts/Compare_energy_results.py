from math import sqrt, ceil
import numpy as np
import matplotlib.pyplot as plt
import re
import seaborn as sns
import pandas as pd

def get_result_MEA(filename):
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

fileMFE = open( "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/RNAfold_bm.log", "r")
lineRna = fileMFE.readline()
lineStruct = fileMFE.readline()

fileEval = open( "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/RNAeval_bm.log", "r")
lineRna2 = fileEval.readline()
lineStruct2 = fileEval.readline()

file = open("/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/benchmark.dbn", "r")
name = file.readline().strip()
rna = file.readline()
twod = file.readline()
contacts = file.readline()
list_name = []
print(np)
while name:
    #print(name)
    if lineRna != rna:
        while lineRna != rna:
            lineRna = fileMFE.readline()
            lineStruct = fileMFE.readline()
    MFE = float(lineStruct[len(lineStruct)-8:len(lineStruct)-2])
    list_name.append(name[5:len(name)-1])
    list_MFE.append(MFE)
    #print("MFE:" + str(MFE))
    lineRna = fileMFE.readline()
    lineStruct = fileMFE.readline()

    if lineRna2 != rna:
        while lineRna2 != rna:
            lineRna2 = fileEval.readline()
            lineStruct2 = fileEval.readline()
    eval = float(lineStruct2[len(lineStruct2)-8:len(lineStruct2)-2])
    list_name.append(name[5:len(name) - 1])
    list_eval.append(eval)
    #print("Eval:" + str(eval))
    lineRna2 = fileEval.readline()
    lineStruct2 = fileEval.readline()

    best_mea = get_result_MEA(name)
    #print("MEA: " + str(best_mea) + "\n")
    list_name.append(name[5:len(name) - 1])
    list_MEA.append(best_mea)
    name = file.readline().strip()
    rna = file.readline()
    twod = file.readline()
    contacts = file.readline()

file.close()
fileMFE.close()
fileEval.close()

'''print(list_MFE)
print(list_MEA)
print(list_eval)'''

#np = [["rna", "type_score", "score"]]
df = pd.DataFrame({'rna':list_name})
print(df)

sns.set_theme(style="ticks")
plt.savefig("output.png")
t=sns.load_dataset('tips')
#sns.stripplot(x="rna",y="score",data=t,jitter=True,hue='type_score',palette='Set1')


