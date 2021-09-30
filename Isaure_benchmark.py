import time
import subprocess
import os
import os.path
from math import sqrt, ceil
import numpy as np
import matplotlib.pyplot as plt



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


#create the command line to run BiORSEO with the modules library of CaRNAval
def create_command_rin(path, name, function, estimator):
    cmd = ("python3 " + path + "/biorseo.py -i " +
      path + "/data/fasta/" +
      name + ".fa  " +
      "-O results/ " +
      "--carnaval " +
      "--patternmatch " +
      "--func " + function + " --" + estimator + " -v " +
      " --biorseo-dir " + path + " " +
      "--modules-path " + path + "/data/modules/RIN/Subfiles")
    return cmd

#create the command line to run BiORSEO with the modules library of RNA3dMotifs Atlas
def create_command_bgsu(path, name, function, estimator):
    cmd = ("python3 " + path + "/biorseo.py -i " +
      path + "/data/fasta/" +
      name + ".fa  " +
      "-O results/ " +
      "--3dmotifatlas " +
      "--jar3d " +
      "--func " + function + " --" + estimator + " -v " +
      "--jar3d-exec " + " /local/local/localopt/jar3d_2014-12-11.jar" +
      " --biorseo-dir " + path + " " +
      "--modules-path " + path + "/data/modules/BGSU")
    return cmd

#create the command line to run BiORSEO with the modules library of RNA3dMotifs
def create_command_desc(path, name, function, estimator):
    cmd = ("python3 " + path + "/biorseo.py -i " +
      path + "/data/fasta/" +
      name + ".fa  " +
      "-O results/ " +
      "--rna3dmotifs " +
      "--patternmatch " +
      "--func " + function + " --" + estimator + " -v " +
      " --biorseo-dir " + path + " " +
      "--modules-path " + path + "/data/modules/DESC")
    return cmd

#create the command line to run BiORSEO with the motifs library of Isaure in .json
def create_command_isaure(path, name, function, estimator):
    cmd = ("python3 " + path + "/biorseo.py -i " +
      path + "/data/fasta/" +
      name + ".fa  " +
      "-O results/ " +
      "--contacts " +
      "--patternmatch " +
      "--func " + function + " --" + estimator +
      " --biorseo-dir " + path + " " +
      "--modules-path " + path + "/data/modules/ISAURE/bibliotheque_a_lire")
    return cmd

#execute the command line correspondin to the information put in the argument.
def execute_command(path, function, estimator, true_ctc, true_str, list_ctc, list_str, modules):
    if (modules == 'desc'):
        cmd = create_command_desc(path, name, function, estimator)
    elif (modules == 'bgsu'):
        cmd = create_command_bgsu(path, name, function, estimator)
    elif (modules == 'rin'):
        cmd = create_command_rin(path, name, function, estimator)
    elif (modules == 'isaure'):
        cmd = create_command_isaure(path, name, function, estimator)

    os.system(cmd)

    """file_path = "results/test_" + name + ".json_pm" + function + "_" + estimator
    if os.path.isfile(file_path):
        tab = write_mcc_in_file(name, true_ctc, true_str, estimator, function)
        list_ctc.append(tab[0])
        list_str.append(tab[1])"""

#Retrieves the list of structures predicted by Biorseo for each sequence in the benchmark.txt file
def get_list_str_by_seq(name, estimator, function, list_str, true_str, modules):
    if modules == 'bgsu':
        extension = ".jar3d"
    elif modules == 'desc':
        extension = ".desc_pm"
    elif modules == 'rin':
        extension = ".rin_pm"
    elif modules == 'json':
        extension = ".json_pm"
    file_path = "results/test_" + name + extension + function + "_" + estimator
    if os.path.isfile(file_path):
        path_benchmark = "data/modules/ISAURE/benchmark.txt"
        max_mcc = get_mcc_structs_max(path_benchmark, name, estimator, function, extension, modules, true_str)
        list_str.append(max_mcc)

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
        #print("warning: division by zero! no contact in the prediction")
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

#Get the best MCC value for all prediction of the results file of the sequence in argument
def get_mcc_structs_max(path_benchmark, sequence_id, estimator, function, extension, modules, true_structure):
    read_prd = open("results/test_" + sequence_id + extension + function + "_" + estimator, "r")
    write = open("results/test_" + sequence_id + ".mcc_" + function + "_" + estimator + "_" + modules, "w")

    max_mcc_str = -1;

    title_exp = ">test_" + sequence_id + ": "
    write.write(title_exp)
    structure_exp = true_structure
    write.write("structure 2d attendue:\n" + structure_exp + "\n")

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
    write.write("max mcc 2D:" + str(max_mcc_str))
    read_prd.close()
    write.close()
    return max_mcc_str

#Create a file that store the information concerning the MCC value obtains for each prediction of the
#sequence in input
def write_mcc_in_file(sequence_id, true_contacts, true_structure, estimator, function):
    read_prd = open("results/test_" + sequence_id + ".json_pm" + function + "_" + estimator, "r")
    write = open("results/test_" + sequence_id + ".mcc_" + function + "_" + estimator, "w")

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

#Create one violin plot to see the distribution of the best MCC (of the base pairing)
def visualization_best_mcc_str(list_struct2d, estimator, function, modules, color, lines_color):

    print(estimator + " + " + function + ": ")
    np_struct2d = np.array(list_struct2d)

    data_to_plot = np_struct2d
    median_2d = np.median(np_struct2d)
    print("mediane 2D: " + str(median_2d) + "\n")

    fig = plt.figure()

    ax = fig.add_axes([0, 0, 1, 1])

    labels = ['structure 2D']
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
    plt.savefig('visualisation_28_09' + estimator + '_' + function + '_' + modules + '.png', bbox_inches='tight')

#Create 4 violin plot to see the distribution of the best MCC (of the base pairing)
def visualization_best_mcc_str_4_figures(list_struct2d, color, lines_color):

    np_struct2d_1 = np.array(list_struct2d[0])
    np_struct2d_2 = np.array(list_struct2d[1])
    np_struct2d_3 = np.array(list_struct2d[2])
    np_struct2d_4 = np.array(list_struct2d[3])

    data_to_plot = [np_struct2d_1, np_struct2d_2, np_struct2d_3, np_struct2d_4]
    fig = plt.figure()
    fig.set_size_inches(6, 3)
    ax = fig.add_axes([0, 0, 1, 1])

    labels = ['MEA + E', 'MEA + F', 'MFE + E', 'MFE + F']
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    #ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel("motifs d'Isaure")
    ax.set_ylabel('MCC (en fonction des appariements)')

    violins = ax.violinplot(data_to_plot, showmedians=True)

    for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        vp = violins[partname]
        vp.set_edgecolor(lines_color)
        vp.set_linewidth(1)

    for v in violins['bodies']:
        v.set_facecolor(color)
    plt.savefig('visualisation_28_09_Isaure_E_F.png', dpi=200, bbox_inches='tight')

#Create 2 violin plot to see the distribution of the best MCC (one for the MCC of
# the base pairing and one for the contacts)
def visualization_best_mcc(list_struct2d, list_contacts, estimator, function, modules, color, lines_color):

    print(estimator + " + " + function + ": ")
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
    plt.savefig('visualisation_16_06_' + estimator + '_' + function + '_' + modules + '.png', bbox_inches='tight')

#Return the list of names of all the sequence in the benchmark.txt, the list of structures and the list
# of contacts predicted in the result file of each sequence.
# This function is only use for the result in .json_pm format files
def get_list_structs_contacts_all(path_benchmark, estimator, function):
    myfile = open(path_benchmark, "r")
    list_name = []

    complete_list_struct2d = []
    complete_list_contacts = []

    name = myfile.readline()
    contacts = myfile.readline()
    seq = myfile.readline()
    structure2d = myfile.readline()
    count = 0
    while seq:
        name = name[6:].strip()
        count = count + 1
        file_path = "results/test_" + name + ".json_pm" + function + "_" + estimator
        if os.path.isfile(file_path):
            file_result = open(file_path, "r")
            list_struct2d = []
            list_contacts = []
            list_name.append(name)
            title_prd = file_result.readline()
            structure_prd = file_result.readline()
            sequence = structure_prd
            while structure_prd:
                structure_prd = file_result.readline()
                if (len(structure_prd) != 0):
                    mcc_tab = compare_two_structures(structure2d, structure_prd[:len(sequence)])
                    mcc_str = mattews_corr_coeff(mcc_tab[0], mcc_tab[1], mcc_tab[2], mcc_tab[3])
                    list_struct2d.append(mcc_str)

                    contacts_prd = file_result.readline()
                    if (len(contacts_prd) == len(contacts)):
                        mcc_tab = compare_two_contacts(contacts, contacts_prd)
                        mcc_ctc = mattews_corr_coeff(mcc_tab[0], mcc_tab[1], mcc_tab[2], mcc_tab[3])
                        list_contacts.append(mcc_ctc)
            complete_list_struct2d.append(list_struct2d)
            complete_list_contacts.append(list_contacts)
        name = myfile.readline()
        contacts = myfile.readline()
        seq = myfile.readline()
        structure2d = myfile.readline()
    return [list_name, complete_list_struct2d, complete_list_contacts]
    myfile.close()

#Return the list of names of all the sequence in the benchmark.txt, the list of structures
# predicted in the result file of each sequence.
def get_list_structs_all(path_benchmark, estimator, function, modules):
    if modules == 'bgsu':
        extension = ".jar3d"
    elif modules == 'desc':
        extension = ".desc_pm"
    elif modules == 'rin':
        extension = ".rin_pm"
    elif modules == 'json':
        extension = ".json_pm"

    myfile = open(path_benchmark, "r")
    list_name = []

    complete_list_struct2d = []
    complete_list_contacts = []

    name = myfile.readline()
    contacts = myfile.readline()
    seq = myfile.readline()
    structure2d = myfile.readline()
    count = 0
    while seq:
        name = name[6:].strip()
        count = count + 1
        file_path = "results/test_" + name + extension + function + "_" + estimator
        if os.path.isfile(file_path):
            file_result = open(file_path, "r")
            list_struct2d = []
            list_name.append(name)
            title_prd = file_result.readline()
            structure_prd = file_result.readline()
            sequence = structure_prd
            while structure_prd:
                structure_prd = file_result.readline()
                if (len(structure_prd) != 0):
                    mcc_tab = compare_two_structures(structure2d, structure_prd[:len(sequence)])
                    mcc_str = mattews_corr_coeff(mcc_tab[0], mcc_tab[1], mcc_tab[2], mcc_tab[3])
                    list_struct2d.append(mcc_str)

                    contacts_prd = file_result.readline()
            complete_list_struct2d.append(list_struct2d)
        name = myfile.readline()
        contacts = myfile.readline()
        seq = myfile.readline()
        structure2d = myfile.readline()
    return [list_name, complete_list_struct2d]
    myfile.close()

#Return the list in argument in two lists, each list containing half of the list in argument
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

#Create a boxplot with all the MCC (for the base pairing) obtains with all the prediction for each sequence,
# divide in two figures.
def visualization_all_mcc_str(path_benchmark, estimator, function, modules):

    list_name = get_list_structs_all(path_benchmark, estimator, function, modules)[0]
    tab_struct2d = get_list_structs_all(path_benchmark, estimator, function, modules)[1]

    min = 20
    max = 0
    max_i = 0
    min_i = 0
    for i in range(len(tab_struct2d)):
        if (len(tab_struct2d[i]) > max):
            max = len(tab_struct2d[i])
            max_i = i
        if (len(tab_struct2d[i]) < min):
            min = len(tab_struct2d[i])
            min_i = i

    print("max: " + list_name[max_i] + " " + str(max) + "     min: " + list_name[min_i] + " " + str(min) + "\n")

    np_struct2d = np.array(tab_struct2d)
    size = len(tab_struct2d)
    list_median_str = []
    for i in range(size):
        list_median_str.append(np.median(np_struct2d[i]))

    all_str = []
    for i in range(size):
        for j in range(len(np_struct2d[i])):
            all_str.append(np_struct2d[i][j])

    """print("mediane struct" + estimator + " + " + function + " : " + str(np.median(all_str)))
    print("ecart struct" + estimator + " + " + function + " : " + str(np.std(all_str)) + "\n")"""

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
    plt.savefig('visualisation_all_mcc_' + estimator + "_" + function + "_" + modules + '.png', bbox_inches='tight')

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
    plt.savefig('visualisation_all_mcc_' + estimator + "_" + function + "_" + modules + '_2.png', bbox_inches='tight')

# Create a boxplot with all the MCC (for the base pairing and the contacts) obtains with all
# the prediction for each sequence, divide in two figures.
def visualization_all_mcc(path_benchmark, estimator, function, modules):

    list_name = get_list_structs_contacts_all(path_benchmark, estimator, function)[0]
    tab_struct2d = get_list_structs_contacts_all(path_benchmark, estimator, function)[1]
    tab_contacts = get_list_structs_contacts_all(path_benchmark, estimator, function)[2]

    min = 20
    max = 0
    max_i = 0
    min_i = 0
    for i in range(len(tab_struct2d)):
        if (len(tab_struct2d[i]) > max):
            max = len(tab_struct2d[i])
            max_i = i
        if (len(tab_struct2d[i]) < min):
            min = len(tab_struct2d[i])
            min_i = i

    print("max: " + list_name[max_i] + " " + str(max) + "     min: " + list_name[min_i] + " " + str(min) + "\n")

    np_struct2d = np.array(tab_struct2d)
    size = len(tab_struct2d)
    list_median_str = []
    for i in range(size):
        list_median_str.append(np.median(np_struct2d[i]))

    all_str = []
    for i in range(size):
        for j in range(len(np_struct2d[i])):
            all_str.append(np_struct2d[i][j])

    """print("mediane struct" + estimator + " + " + function + " : " + str(np.median(all_str)))
    print("ecart struct" + estimator + " + " + function + " : " + str(np.std(all_str)) + "\n")"""

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
    plt.savefig('visualisation_128arn_structure2d_' + estimator + "_" + function + "_" + modules + '.png', bbox_inches='tight')

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
    plt.savefig('visualisation_128arn_structure2d_' + estimator + "_" + function + "_" + modules + '_2.png', bbox_inches='tight')

    np_contacts = np.array(tab_contacts)
    size = len(tab_contacts)
    list_median_ctc = []
    for i in range(size):
        list_median_ctc.append(np.median(np_contacts[i]))

    all_ctc = []
    for i in range(size):
        for j in range(len(np_contacts[i])):
            all_ctc.append(np_contacts[i][j])

    """print("mediane ctc" + estimator + " + " + function + " : " + str(np.median(all_ctc)))
    print("ecart ctc" + estimator + " + " + function + " : " + str(np.std(all_ctc)) + "\n")"""

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

# Most of those commands lines and this script work only for the benchmark_16-07-2021.json
# file provides by Isaure Chauvot de Beauchene. This script is means to automize the execution
# of BiORSEO for each sequence and the creation of figures.

# after compiling create_file.cpp in 'scripts' directory (compiling once is enough) :
# clang++ create_files.cpp -o create
#cmd = ("scripts/create")

# after compiling create_file.cpp in 'scripts' directory (compiling once is enough) :
# clang++ add_delimiter.cpp -o addDelimiter
#cmd0 = ("cppsrc/Scripts/addDelimiter")

# after compiling count_pattern.cpp in 'scripts' directory (compiling once is enough) :
# clang++ count_pattern.cpp -o countPattern
#cmd1 = ("cppsrc/Scripts/countPattern")

myfile = open("data/modules/ISAURE/benchmark.txt", "r")
name = myfile.readline()
contacts = myfile.readline()
seq = myfile.readline()
structure2d = myfile.readline()

list_struct2d = [[],[],[],[]]

# source path to the directory (for my computer)
path = "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo"
path2 = "/local/local/BiorseoNath"
while seq:
    name = name[6:].strip()
    print(name)

    # after compiling delete_same_pdb.cpp in 'scripts' directory (compiling once is enough) :
    # clang++ delete_same_pdb.cpp -o deletePdb
    #cmd2 = ("cppsrc/Scripts/deletePdb " + name)
    #os.system(cmd2)

    """get_list_str_by_seq(name, 'MEA', 'E', list_struct2d[0], structure2d, 'json')
    get_list_str_by_seq(name, 'MEA', 'F', list_struct2d[1], structure2d, 'json')
    get_list_str_by_seq(name, 'MFE', 'E', list_struct2d[2], structure2d, 'json')
    get_list_str_by_seq(name, 'MFE', 'F', list_struct2d[3], structure2d, 'json')"""

    name = myfile.readline()
    contacts = myfile.readline()
    seq = myfile.readline()
    structure2d = myfile.readline()

#visualization_best_mcc_str_4_figures(list_struct2d, 'red', '#900C3F')

"""visualization_best_mcc(list_struct2d_A_MFE, list_contacts_A_MFE, 'MFE', 'A', 'rin', 'red', '#900C3F')
visualization_best_mcc(list_struct2d_B_MFE, list_contacts_B_MFE, 'MFE', 'B', 'rin', 'blue', '#0900FF')
visualization_best_mcc(list_struct2d_A_MEA, list_contacts_A_MEA, 'MEA', 'A', 'rin', 'red', '#900C3F')
visualization_best_mcc(list_struct2d_B_MEA, list_contacts_B_MEA, 'MEA', 'B', 'rin', 'blue', '#0900FF')"""


myfile.close()

path_benchmark = "data/modules/ISAURE/benchmark.txt"
visualization_all_mcc_str(path_benchmark, 'MEA', 'C', 'bgsu')
visualization_all_mcc_str(path_benchmark, 'MEA', 'D', 'bgsu')
visualization_all_mcc_str(path_benchmark, 'MFE', 'C', 'bgsu')
visualization_all_mcc_str(path_benchmark, 'MFE', 'D', 'bgsu')