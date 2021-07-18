from math import sqrt, ceil
import numpy as np
import matplotlib.pyplot as plt

file = open("/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_version_initiale/benchmark.dbn", "r")
name = file.readline()
rna = file.readline()
twod = file.readline()
contacts = file.readline()
length = len(rna)
nb_ctc = contacts.count('*')
print("--------------------------------------------------------")

ctc_max = nb_ctc
ctc_min = nb_ctc

np_lgt = []
np_lgt.append(length)

np_ctc = []
np_ctc.append(nb_ctc)

np = []
np.append([length, nb_ctc])

while name:
    print(contacts)
    print(length)
    print(nb_ctc)
    print("--------------------------------------------------------")

    name = file.readline()
    rna = file.readline()
    length = len(rna)
    if length != 0 :
        np_lgt.append(length)
    twod = file.readline()
    contacts = file.readline()
    nb_ctc = contacts.count('*')
    if nb_ctc != 0:
        np_ctc.append(nb_ctc)
        np.append([length, nb_ctc])
    if nb_ctc > ctc_max:
        ctc_max = nb_ctc
    if nb_ctc < ctc_min and nb_ctc != 0:
        ctc_min = nb_ctc
file.close()
print(np_lgt)
print(np_ctc)
print(np)

x = np_lgt
y = np_ctc

index = np_ctc.index(ctc_max)
index2 = np_ctc.index(ctc_min)

plt.scatter(x, y, c = 'blue')
plt.annotate("(" + str(np_lgt[index]) + "," + str(ctc_max) + ")", (np_lgt[index], ctc_max),c ='red')
plt.scatter(np_lgt[index], ctc_max,c = 'red')
plt.annotate("(" + str(np_lgt[index2]) + "," + str(ctc_min) + ")", (np_lgt[index2], ctc_min),c ='green')
plt.scatter(np_lgt[index2], ctc_min,c = 'green')

plt.xlabel('longeur de l\'arn')
plt.ylabel('nombre de contacts')
plt.savefig('stats.png')
