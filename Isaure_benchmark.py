import time
import subprocess
import os

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
    if (tp + fp == 0):
        print("We have an issue : no positives detected ! (linear structure)")
    return (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

def f1_score(tp, tn, fp, fn):
    return 2 * tp / (2 * tp + fp + fn)

def specificity(tp, tn, fp, fn):
    return tn / (tn + fp)

cmd = ("cppsrc/Scripts/addDelimiter")
cmd1 = ("cppsrc/Scripts/countPattern")
cmd2 = ("cppsrc/Scripts/deletePdb")

myfile = open("data/fasta/benchmark.fa", "r")
name = myfile.readline()
seq = myfile.readline()

run_test(cmd2 + " 1JJ2" + ".fa", log)
print(cmd2 + " 1JJ2" + ".fa")
cmd2 = create_command("1JJ2")
print(cmd3)
os.system(cmd3)
filetest = open("test_1L9A.json_pmE", "r")
line = filetest.readline()
while line:
    print(line)
    line = filetest.readline()
filetest.close()

"""while seq:
    name = name[6:].strip()
    print(name)
    run_test(cmd, log)
    run_test(cmd1, log)
    run_test(cmd2 + " " + name + ".fa", log)
    print(cmd2 + " " + name + ".fa")
    cmd3 = create_command(name)
    os.system(cmd3)
    name = myfile.readline()
    seq = myfile.readline()"""
myfile.close()
