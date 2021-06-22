import time
import subprocess
import os
from math import sqrt, ceil


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

def create_command(name):
    cmd = ("python3 /mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/biorseo.py -i " +
      "/mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/fasta/" +
      name + ".fa  " +
      "-O results/ " +
      "--contacts " +
      "--patternmatch " +
      "--func E  -v " +
      "--biorseo-dir /mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo " +
      "--modules-path /mnt/c/Users/natha/Documents/IBISC/biorseo2/biorseo/data/modules/ISAURE/Motifs_derniere_version ")
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
    #print(str(tp) + " " + str(tn) + " " + str(fp) + " " + str(fn) + "\n")
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
        print("warning: division by zero!")
        return 0
    elif (tp + fp == 0):
        print("We have an issue : no positives detected ! (linear structure)")
    return (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

def f1_score(tp, tn, fp, fn):
    return 2 * tp / (2 * tp + fp + fn)

def specificity(tp, tn, fp, fn):
    return tn / (tn + fp)

# ================== Code from Louis Beckey Benchark.py ==============================

def write_mcc_in_file(sequence_id, true_contacts, true_structure):
    read_prd = open("results/test_" + sequence_id + ".json_pmE", "r")
    write = open("results/test_" + sequence_id + ".mcc", "w")

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
            write.write("mcc: " + str(mcc_str) + "\n")
            contacts_prd = read_prd.readline()
            write.write("\ncontacts predits:\n" + contacts_prd)
            if (len(contacts_prd) == len(contacts_exp)):
                mcc_tab = compare_two_contacts(contacts_exp, contacts_prd)
                mcc_ctc = mattews_corr_coeff(mcc_tab[0], mcc_tab[1], mcc_tab[2], mcc_tab[3])
                write.write("mcc: " + str(mcc_ctc) + "\n\n")
            else:
                write.write("mcc: no expected contacts sequence or not same length between expected and predicted\n\n")

    read_prd.close()
    write.close()


cmd = ("cppsrc/Scripts/addDelimiter")
cmd1 = ("cppsrc/Scripts/countPattern")
cmd2 = ("cppsrc/Scripts/deletePdb")

myfile = open("data/modules/ISAURE/Motifs_version_initiale/benchmark.txt", "r")
name = myfile.readline()
contacts = myfile.readline()
seq = myfile.readline()
structure2d = myfile.readline()

"""run_test(cmd2 + " 1JJ2" + ".fa", log)
print(cmd2 + " 1JJ2" + ".fa")
cmd3 = create_command("1JJ2")
print(cmd3)
os.system(cmd3)
read_exp = open("data/modules/ISAURE/Motifs_version_initiale/benchmark.txt", "r")
read_prd = open("results/test_1JJ2.json_pmE", "r")
write = open("results/test_1JJ2.mcc", "w")

title_exp = read_exp.readline()
write.write(title_exp)
contacts_exp = read_exp.readline()
structure_exp = read_exp.readline()
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
        write.write("mcc: " + str(mcc_str) + "\n")
        contacts_prd = read_prd.readline()
        write.write("\ncontacts predits:\n" + contacts_prd)
        mcc_tab = compare_two_contacts(contacts_exp, contacts_prd)
        mcc_ctc = mattews_corr_coeff(mcc_tab[0], mcc_tab[1], mcc_tab[2], mcc_tab[3])
        write.write("mcc: " + str(mcc_ctc) + "\n\n")

read_prd.close()
write.close()"""

while seq:
    name = name[6:].strip()
    print(name)
    """run_test(cmd, log)
    run_test(cmd1, log)
    run_test(cmd2 + " " + name + ".fa", log)
    print(cmd2 + " " + name + ".fa")
    cmd3 = create_command(name)
    os.system(cmd3)"""
    write_mcc_in_file(name, contacts, structure2d)
    name = myfile.readline()
    contacts = myfile.readline()
    seq = myfile.readline()
    structure2d = myfile.readline()
myfile.close()
