import time
import subprocess
import os
import os.path
from math import sqrt, ceil
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


log_path = "test.log"
log = open(log_path, 'a')

def run_test(cmd, log):
    log.write(time.asctime(time.localtime(time.time())) + " : Run process \"" + cmd + "\"\n")
    log.flush()
    process = subprocess.Popen(cmd.split(' ') ,shell=False,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # Poll process.stdout to show stdout live
    while process.poll() is None:
      output = process.stdout.readline()
      if output:
        log.write(output.decode())
        log.flush()
    rc = process.poll()

def create_command_E(name, estimator):
    #cmd = ("python3 /mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/biorseo.py -i " +
    cmd = ("python3 /local/local/BiorseoNath/biorseo.py -i " +
      "/local/local/BiorseoNath/data/fasta/" +
      name + ".fa  " +
      "-O results/ " +
      "--contacts " +
      "--patternmatch " +
      "--func E --" + estimator + " -v " +
      "--biorseo-dir /local/local/BiorseoNath " +
      "--modules-path /local/local/BiorseoNath/data/modules/ISAURE/Motifs_derniere_version ")
    return cmd

def create_command_F(name, estimator):
    #cmd = ("python3 /mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/biorseo.py -i " +
    cmd = ("python3 /local/local/BiorseoNath/biorseo.py -i " +
      "/local/local/BiorseoNath/data/fasta/" +
      name + ".fa  " +
      "-O results/ " +
      "--contacts " +
      "--patternmatch " +
      "--func F --" + estimator + " -v " +
      "--biorseo-dir /local/local/BiorseoNath " +
      "--modules-path /local/local/BiorseoNath/data/modules/ISAURE/Motifs_derniere_version ")
    return cmd

# ================== Code from Louis Beckey Benchark.py ==============================
def dbn_to_basepairs(structure):
    parenthesis = []
    brackets = []
    braces = []
    rafters = []
    basepairs = []
    As = []
    Bs = []
    for i, c in enumerate(structure):
        if c == '(':
            parenthesis.append(i)
        if c == '[':
            brackets.append(i)
        if c == '{':
            braces.append(i)
        if c == '<':
            rafters.append(i)
        if c == 'A':
            As.append(i)
        if c == 'B':
            Bs.append(i)
        if c == '.':
            continue
        if c == ')':
            basepairs.append((i, parenthesis.pop()))
        if c == ']':
            basepairs.append((i, brackets.pop()))
        if c == '}':
            basepairs.append((i, braces.pop()))
        if c == '>':
            basepairs.append((i, rafters.pop()))
        if c == 'a':
            basepairs.append((i, As.pop()))
        if c == 'b':
            basepairs.append((i, Bs.pop()))
    return basepairs

def compare_two_contacts(true_ctc, prediction):
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    for i in range(len(true_ctc)):
        if true_ctc[i] == '*' and prediction[i] == '*':
            tp += 1
        elif true_ctc[i] == '.' and prediction[i] == '.':
            tn += 1
        elif true_ctc[i] == '.' and prediction[i] == '*':
            fp += 1
        elif true_ctc[i] == '*' and prediction[i] == '.':
            fn += 1
    """if ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) == 0):
        print(str(tp) + " " + str(tn) + " " + str(fp) + " " + str(fn) + "\n")
        print(true_ctc)
        print(prediction)"""
    return [tp, tn, fp, fn]

def compare_two_structures(true2d, prediction):
    true_basepairs = dbn_to_basepairs(true2d)
    pred_basepairs = dbn_to_basepairs(prediction)
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    for bp in true_basepairs:
        if bp in pred_basepairs:
            tp += 1
        else:
            fn += 1
    for bp in pred_basepairs:
        if bp not in true_basepairs:
            fp += 1
    tn = len(true2d) * (len(true2d) - 1) * 0.5 - fp - fn - tp
    return [tp, tn, fp, fn]

def mattews_corr_coeff(tp, tn, fp, fn):
    if ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn) == 0):
        print("warning: division by zero! no contact in the prediction")
        #print("tp: " + str(tp) + " fp: " + str(fp) + " tn: " + str(tn) + " fn: " + str(fn))
        return -1
    elif (tp + fp == 0):
        print("We have an issue : no positives detected ! (linear structure)")
    return (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

def f1_score(tp, tn, fp, fn):
    return 2 * tp / (2 * tp + fp + fn)

def specificity(tp, tn, fp, fn):
    return tn / (tn + fp)

# ================== Code from Louis Beckey Benchark.py ==============================

def write_mcc_in_file_E(sequence_id, true_contacts, true_structure, estimator):
    read_prd = open("results/test_" + sequence_id + ".json_pmE_"+ estimator, "r")
    write = open("results/test_" + sequence_id + ".mcc_E_" + estimator, "w")

    max_mcc_str = -1;
    max_mcc_ctc = -1;

    title_exp = ">test_" + sequence_id + ": "
    write.write(title_exp)
    contacts_exp = true_contacts
    structure_exp = true_structure
    write.write("structure 2d attendue:\n" + structure_exp + "\n")
    write.write("contacts attendus:\n" + contacts_exp + "\n" + len(structure_exp) * "-")

    title_prd = read_prd.readline()
    structure_prd = read_prd.readline()
    sequence_prd = structure_prd
    while structure_prd:
        structure_prd = read_prd.readline()
        if (len(structure_prd) != 0):
            write.write("\nstructure 2d predite:\n" + structure_prd[:len(sequence_prd)] + "\n")
            mcc_tab = compare_two_structures(structure_exp, structure_prd[:len(sequence_prd)])
            mcc_str = mattews_corr_coeff(mcc_tab[0], mcc_tab[1], mcc_tab[2], mcc_tab[3])
            if (max_mcc_str < mcc_str):
                max_mcc_str = mcc_str
            write.write("mcc: " + str(mcc_str) + "\n")

            contacts_prd = read_prd.readline()
            write.write("\ncontacts predits:\n" + contacts_prd)
            if (len(contacts_prd) == len(contacts_exp)):
                mcc_tab = compare_two_contacts(contacts_exp, contacts_prd)
                mcc_ctc = mattews_corr_coeff(mcc_tab[0], mcc_tab[1], mcc_tab[2], mcc_tab[3])
                if (max_mcc_ctc < mcc_ctc):
                    max_mcc_ctc = mcc_ctc
                write.write("mcc: " + str(mcc_ctc) + "\n\n")
            else:
                write.write("mcc: no expected contacts sequence or not same length between expected and predicted\n\n")
    write.write("max mcc 2D:" + str(max_mcc_str))
    write.write("max mcc ctc:" + str(max_mcc_ctc))
    read_prd.close()
    write.close()
    return [max_mcc_ctc, max_mcc_str]

def write_mcc_in_file_F(sequence_id, true_contacts, true_structure, estimator):

    read_prd = open("results/test_" + sequence_id + ".json_pmF_" + estimator, "r")
    write = open("results/test_" + sequence_id + ".mcc_F_" + estimator, "w")

    max_mcc_str = -1;
    max_mcc_ctc = -1;

    title_exp = ">test_" + sequence_id + ": "
    write.write(title_exp)
    contacts_exp = true_contacts
    structure_exp = true_structure
    write.write("structure 2d attendue:\n" + structure_exp + "\n")
    write.write("contacts attendus:\n" + contacts_exp + "\n" + len(structure_exp) * "-")

    title_prd = read_prd.readline()
    structure_prd = read_prd.readline()
    sequence_prd = structure_prd
    while structure_prd:
        structure_prd = read_prd.readline()
        if (len(structure_prd) != 0):
            write.write("\nstructure 2d predite:\n" + structure_prd[:len(sequence_prd)] + "\n")
            mcc_tab = compare_two_structures(structure_exp, structure_prd[:len(sequence_prd)])
            mcc_str = mattews_corr_coeff(mcc_tab[0], mcc_tab[1], mcc_tab[2], mcc_tab[3])
            if (max_mcc_str < mcc_str):
                max_mcc_str = mcc_str
            write.write("mcc: " + str(mcc_str) + "\n")

            contacts_prd = read_prd.readline()
            write.write("\ncontacts predits:\n" + contacts_prd)
            if (len(contacts_prd) == len(contacts_exp)):
                mcc_tab = compare_two_contacts(contacts_exp, contacts_prd)
                mcc_ctc = mattews_corr_coeff(mcc_tab[0], mcc_tab[1], mcc_tab[2], mcc_tab[3])
                if (max_mcc_ctc < mcc_ctc):
                    max_mcc_ctc = mcc_ctc
                write.write("mcc: " + str(mcc_ctc) + "\n\n")
            else:
                write.write("mcc: no expected contacts sequence or not same length between expected and predicted\n\n")
    write.write("max mcc 2D:" + str(max_mcc_str))
    write.write("max mcc ctc:" + str(max_mcc_ctc))
    read_prd.close()
    write.close()
    return [max_mcc_ctc, max_mcc_str]

def set_axis_style(ax, labels):
    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')

def visualization_best_mcc(list_struct2d, list_contacts, estimator, function, color, lines_color):

    np_struct2d = np.array(list_struct2d)
    np_contacts = np.array(list_contacts)

    data_to_plot = [np_struct2d, np_contacts]
    median_2d = np.median(np_struct2d)
    median_ctc = np.median(np_contacts)
    print("mediane 2D: " + str(median_2d) + "\n")
    print("mediane ctc: " + str(median_ctc) + "\n")

    fig = plt.figure()

    ax = fig.add_axes([0, 0, 1, 1])

    labels = ['structure 2D', 'contacts']
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlabel(function)
    ax.set_ylabel('MCC')

    violins = ax.violinplot(data_to_plot, showmedians=True)

    for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        vp = violins[partname]
        vp.set_edgecolor(lines_color)
        vp.set_linewidth(1)

    for v in violins['bodies']:
        v.set_facecolor(color)
    plt.savefig('visualisation_16_06_' + estimator + '_' + function + '.png', bbox_inches='tight')

def get_list_structs_contacts(path_benchmark, estimator, function):
    myfile = open(path_benchmark, "r")
    list_name = []

    complete_list_struct2d_F = []
    complete_list_contacts_F = []

    name = myfile.readline()
    contacts = myfile.readline()
    seq = myfile.readline()
    structure2d = myfile.readline()
    count = 0
    while seq:
        name = name[6:].strip()
        count = count + 1
        file_path = "results/test_" + name + ".json_pm" + function +"_" + estimator
        if os.path.isfile(file_path):
            file_result = open(file_path, "r")
            list_struct2d_F = []
            list_contacts_F = []
            list_name.append(name)
            title_prd = file_result.readline()
            structure_prd = file_result.readline()
            sequence = structure_prd
            while structure_prd:
                structure_prd = file_result.readline()
                if (len(structure_prd) != 0):
                    mcc_tab = compare_two_structures(structure2d, structure_prd[:len(sequence)])
                    mcc_str = mattews_corr_coeff(mcc_tab[0], mcc_tab[1], mcc_tab[2], mcc_tab[3])
                    list_struct2d_F.append(mcc_str)

                    contacts_prd = file_result.readline()
                    if (len(contacts_prd) == len(contacts)):
                        mcc_tab = compare_two_contacts(contacts, contacts_prd)
                        mcc_ctc = mattews_corr_coeff(mcc_tab[0], mcc_tab[1], mcc_tab[2], mcc_tab[3])
                        list_contacts_F.append(mcc_ctc)
            complete_list_struct2d_F.append(list_struct2d_F)
            complete_list_contacts_F.append(list_contacts_F)
        name = myfile.readline()
        contacts = myfile.readline()
        seq = myfile.readline()
        structure2d = myfile.readline()
    return [list_name, complete_list_struct2d_F, complete_list_contacts_F]
    myfile.close()

def get_half(list):

    first_half = []
    second_half = []
    if (len(list) % 2 == 0):
        middle = len(list) / 2
    else:
        middle = len(list) / 2 + 0.5

    for i in range (int(middle)):
        first_half.append(list[i])

    for i in range (int(middle)):
        if i + int(middle) < len(list):
         second_half.append(list[i + int(middle)])

    return [first_half, second_half]

def visualization_all_mcc(path_benchmark, estimator, function):

    list_name = get_list_structs_contacts(path_benchmark, estimator, function)[0]
    tab_struct2d = get_list_structs_contacts(path_benchmark, estimator, function)[1]
    tab_contacts = get_list_structs_contacts(path_benchmark, estimator, function)[2]

    np_struct2d = np.array(tab_struct2d)
    size = len(tab_struct2d)
    list_median_str = []
    for i in range(size):
        list_median_str.append(np.median(np_struct2d[i]))

    data = [x for _, x in sorted(zip(list_median_str, tab_struct2d))]
    boxName = [x for _, x in sorted(zip(list_median_str, list_name))]

    if (len(data) % 2 == 0):
        absciss = len(data) / 2
    else:
        absciss = len(data) / 2 + 0.5

    divide_tab_name = get_half(boxName)
    divide_tab_data = get_half(data)

    plt.figure(figsize=(15,4),dpi=200)
    plt.xticks(rotation=90)
    plt.boxplot(divide_tab_data[0], medianprops=dict(color='black'))
    for i in range(int(absciss)):
        y =data[i]
        x = np.random.normal(1 + i, 0.04, size=len(y))
        plt.scatter(x, y)
        plt.xticks(np.arange(1, absciss + 1), divide_tab_name[0])

    plt.xlabel('nom de la séquence')
    plt.ylabel('MCC (appariements)')
    plt.savefig('visualisation_128arn_structure2d_' + estimator + "_" + function + '.png', bbox_inches='tight')

    plt.figure(figsize=(15, 4), dpi=200)
    plt.xticks(rotation=90)
    plt.boxplot(divide_tab_data[1], medianprops=dict(color='black'))
    for i in range(len(data)):
        if i + int(absciss) < len(data):
            y = data[i + int(absciss)]
            x = np.random.normal(1 + i, 0.04, size=len(y))
            plt.scatter(x, y)
            plt.xticks(np.arange(1, absciss + 1), divide_tab_name[1])

    plt.xlabel('nom de la séquence')
    plt.ylabel('MCC')
    plt.savefig('visualisation_128arn_structure2d_' + estimator + "_" + function + '_2.png', bbox_inches='tight')

    np_contacts = np.array(tab_contacts)
    size = len(tab_contacts)
    list_median_ctc = []
    for i in range(size):
        list_median_ctc.append(np.median(np_contacts[i]))

    data = [x for _, x in sorted(zip(list_median_ctc, tab_contacts))]
    boxName = [x for _, x in sorted(zip(list_median_ctc, list_name))]

    if (len(data) % 2 == 0) :
        absciss = len(data)/2
    else :
        absciss = len(data)/2 + 0.5

    divide_tab_name = get_half(boxName)
    divide_tab_data = get_half(data)

    plt.figure(figsize=(15, 4), dpi=200)
    plt.xticks(rotation=90)
    plt.boxplot(divide_tab_data[0], medianprops=dict(color='black'))
    for i in range(int(absciss)):
        y = data[i]
        x = np.random.normal(1 + i, 0.04, size=len(y))
        plt.scatter(x, y)
        plt.xticks(np.arange(1, absciss + 1), divide_tab_name[0])

    plt.xlabel('nom de la séquence')
    plt.ylabel('MCC (contacts)')
    plt.savefig('visualisation_128arn_contacts_' + estimator + "_" + function + '.png', bbox_inches='tight')

    plt.figure(figsize=(15, 4), dpi=200)
    plt.xticks(rotation=90)
    plt.boxplot(divide_tab_data[1], medianprops=dict(color='black'))
    for i in range(len(data)):
        if i + int(absciss) < len(data) :
            y = data[i + int(absciss)]
            x = np.random.normal(1 + i, 0.04, size=len(y))
            plt.scatter(x, y)
            plt.xticks(np.arange(1, absciss + 1), divide_tab_name[1])

    plt.xlabel('nom de la séquence')
    plt.ylabel('MCC')
    plt.savefig('visualisation_128arn_contacts_' + estimator + "_" + function + '_2.png', bbox_inches='tight')

#cmd = ("cppsrc/Scripts/create")
#cmd0 = ("cppsrc/Scripts/addDelimiter")
#cmd1 = ("cppsrc/Scripts/countPattern")
#cmd2 = ("cppsrc/Scripts/deletePdb")

myfile = open("data/modules/ISAURE/Motifs_version_initiale/benchmark.txt", "r")
name = myfile.readline()
contacts = myfile.readline()
seq = myfile.readline()
structure2d = myfile.readline()

list_struct2d_E_MFE = []
list_contacts_E_MFE = []
list_struct2d_F_MFE = []
list_contacts_F_MFE = []

list_struct2d_E_MEA = []
list_contacts_E_MEA = []
list_struct2d_F_MEA = []
list_contacts_F_MEA = []

countE_MFE = 0
countF_MFE = 0

countE_MEA = 0
countF_MEA = 0
"""
while seq:
    name = name[6:].strip()
    print(name)

    cmd2 = ("cppsrc/Scripts/deletePdb " + name)
    
    cmd3 = create_command_E(name, 'MFE')
    os.system(cmd3)

    file_path = "results/test_" + name + ".json_pmE_MFE"
    if os.path.isfile(file_path):
        tabE_MFE = write_mcc_in_file_E(name, contacts, structure2d, 'MFE')
        list_contacts_E_MFE.append(tabE_MFE[0])
        list_struct2d_E_MFE.append(tabE_MFE[1])
        countE_MFE = countE_MFE + 1

    cmd3 = create_command_F(name, 'MFE')
    os.system(cmd3)

    file_path = "results/test_" + name + ".json_pmF_MFE"
    if os.path.isfile(file_path):
        tabF_MFE = write_mcc_in_file_F(name, contacts, structure2d, 'MFE')
        list_contacts_F_MFE.append(tabF_MFE[0])
        list_struct2d_F_MFE.append(tabF_MFE[1])
        countF_MFE = countF_MFE + 1

    cmd3 = create_command_E(name, 'MEA')
    os.system(cmd3)

    file_path = "results/test_" + name + ".json_pmE_MEA"
    if os.path.isfile(file_path):
        tabE_MEA = write_mcc_in_file_E(name, contacts, structure2d, 'MEA')
        list_contacts_E_MEA.append(tabE_MEA[0])
        list_struct2d_E_MEA.append(tabE_MEA[1])
        countE_MEA = countE_MEA + 1

    cmd3 = create_command_F(name, 'MEA')
    os.system(cmd3)

    file_path = "results/test_" + name + ".json_pmF_MEA"
    if os.path.isfile(file_path):
        tabF_MEA = write_mcc_in_file_F(name, contacts, structure2d, 'MEA')
        list_contacts_F_MEA.append(tabF_MEA[0])
        list_struct2d_F_MEA.append(tabF_MEA[1])
        countF_MEA = countF_MEA + 1

    name = myfile.readline()
    contacts = myfile.readline()
    seq = myfile.readline()
    structure2d = myfile.readline()

visualization_best_mcc(list_struct2d_E_MFE, list_contacts_E_MFE, 'MFE', 'E', 'red', '#900C3F')
visualization_best_mcc(list_struct2d_F_MFE, list_contacts_F_MFE, 'MFE', 'F', 'blue', '#0900FF')
visualization_best_mcc(list_struct2d_E_MEA, list_contacts_E_MEA, 'MEA', 'E', 'red', '#900C3F')
visualization_best_mcc(list_struct2d_F_MEA, list_contacts_F_MEA, 'MEA', 'F', 'blue', '#0900FF')

print("countE_MFE: " + str(countE_MFE) + "\n")
print("countF_MFE: " + str(countF_MFE) + "\n")
print("countE_MEA: " + str(countE_MEA) + "\n")
print("countF_MEA: " + str(countF_MEA) + "\n")"""
myfile.close()
path_benchmark = "data/modules/ISAURE/Motifs_version_initiale/benchmark.txt"
visualization_all_mcc(path_benchmark,'MEA', 'F')
visualization_all_mcc(path_benchmark,'MEA', 'E')
visualization_all_mcc(path_benchmark,'MFE', 'E')
visualization_all_mcc(path_benchmark,'MFE', 'F')