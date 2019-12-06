#!/usr/bin/python3
# Created by Louis Becquey, louis.becquey@univ-evry.fr, Oct 2019
# This script processes files containing RNA structures obtained from bi-objective
# optimization programs, and a dot-bracket database of reference structures, to plot
# where are the best solutions in the Pareto set.
#
# The result files should follow this kind of format:
# for Biokop: (option --biokop)
# Structure        Free energy score       Expected accuracy score
# (((...(((...)))))) <tab> obj1_value <tab> obj2_value
# (((............))) <tab> obj1_value <tab> obj2_value
# ((((((...)))...))) <tab> obj1_value <tab> obj2_value
# ...
#
# for BiORSEO: (options --biorseo_**stuff**)
# >Header of the sequence
# GGCACAGAGUUAUGUGCC
# (((...(((...)))))) + Motif1 + Motif2 <tab> obj1_value <tab> obj2_value
# (((............))) <tab> obj1_value <tab> obj2_value
# ((((((...)))...))) + Motif1 <tab> obj1_value <tab> obj2_value
#
# typical Biokop usage:
# python3 pareto_visualizer.py --biokop --folder path/to/your/results/folder --database path/to/the/database_file.dbn
# typical Biorseo usage:
# python3 pareto_visualizer.py --folder path/to/your/results/folder --database path/to/the/database_file.dbn
#

from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm 
import scipy.stats as st
import sys
import os
import subprocess
import getopt

class SecStruct:
    def __init__(self, dot_bracket, obj1_value, obj2_value):
        self.dbn = dot_bracket
        self.objectives = [ obj1_value, obj2_value ]
        self.basepair_list = self.get_basepairs()
        self.length = len(dot_bracket)

    def get_basepairs(self):
        parenthesis = []
        brackets = []
        braces = []
        rafters = []
        basepairs = []
        As = []
        Bs = []
        for i, c in enumerate(self.dbn):
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

    def get_MCC_with(self, reference_structure):
        # Get true and false positives and negatives
        tp = 0
        fp = 0
        tn = 0
        fn = 0
        for bp in reference_structure.basepair_list:
            if bp in self.basepair_list:
                tp += 1
            else:
                fn += 1
        for bp in self.basepair_list:
            if bp not in reference_structure.basepair_list:
                fp += 1
        tn = reference_structure.length * (reference_structure.length - 1) * 0.5 - fp - fn - tp

        # Compute MCC
        if (tp+fp == 0):
            print("We have an issue : no positives detected ! (linear structure)")
        return (tp*tn-fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))


class Pareto:
    def __init__(self, list_of_structs, reference):
        self.predictions = list_of_structs
        self.true_structure = reference
        self.n_pred = len(list_of_structs)
        self.max_obj1 = max([ s.objectives[0] for s in self.predictions ])
        self.max_obj2 = max([ s.objectives[1] for s in self.predictions ])
        self.index_of_best = self.find_best_solution()
        
    def find_best_solution(self):
        # returns the index of the solution of the Pareto set which is the closest
        # to the real 2D structure (the one with the max MCC)
        max_i = -1
        max_mcc = -1
        for i,s in enumerate(self.predictions):
            mcc = s.get_MCC_with(self.true_structure)
            if mcc > max_mcc:
                max_mcc = mcc
                max_i = i
        return max_i

    def get_normalized_coords(self):
        # retrieves the objective values of the best solution and normlizes them
        coords = self.predictions[self.index_of_best].objectives
        if self.max_obj1: # avoid divide by zero if all solutions are 0
            x = coords[0]/self.max_obj1
        else:
            x = 0.5
        if self.max_obj2: # avoid divide by zero if all solutions are 0
            y = coords[1]/self.max_obj2
        else:
            y = 0.5
        return ( x,y )


class RNA:
    def __init__(self, filename, header, seq, struct):
        self.seq_ = seq
        self.header_ = header
        self.struct_ = struct
        self.basename_ = filename


ignored_nt_dict = {}
def is_canonical_nts(seq):
    for c in seq[:-1]:
        if c not in "ACGU":
            if c in ignored_nt_dict.keys():
                ignored_nt_dict[c] += 1
            else:
                ignored_nt_dict[c] = 1
            return False
    return True

def is_canonical_bps(struct):
    if "()" in struct:
        return False
    if "(.)" in struct:
        return False
    if "(..)" in struct:
        return False
    if "[]" in struct:
        return False
    if "[.]" in struct:
        return False
    if "[..]" in struct:
        return False
    return True

def load_from_dbn(file, header_style=1):
    container = []
    counter = 0

    db = open(file, "r")
    c = 0
    header = ""
    seq = ""
    struct = ""
    while True:
        l = db.readline()
        if l == "":
            break
        c += 1
        c = c % 3
        if c == 1:
            header = l[:-1]
        if c == 2:
            seq = l[:-1].upper()
        if c == 0:
            struct = l[:-1]
            n = len(seq)

            if n < 10 or n > 100:
                continue  # ignore too short and too long RNAs
            if is_canonical_nts(seq) and is_canonical_bps(struct):
                if header_style == 1: container.append(RNA(header.replace('/', '_').split('(')[-1][:-1], header, seq, struct))
                if header_style == 2: container.append(RNA(header.replace('/', '_').split('[')[-1][:-41], header, seq, struct))
                if '[' in struct: counter += 1
    db.close()
    return container, counter

def parse_biokop(folder, basename, ext=".biok"):
    solutions = []
    if os.path.isfile(os.path.join(folder, basename + ext)):
        rna = open(os.path.join(folder, basename + ext), "r")
        lines = rna.readlines()
        rna.close()
        different_2ds = []
        for s in lines[1:]:
            if s == '\n':
                continue
            splitted = s.split('\t')
            db2d = splitted[0]
            if db2d not in different_2ds:
                different_2ds.append(db2d)
            # here is a negative sign because Biokop actually minimizes -MEA instead
            # of maximizing MEA : we switch back to MEA
            solutions.append(SecStruct(db2d, -float(splitted[1]), -float(splitted[2][:-1])))

        # check the range of MEA in this pareto set
        min_mea = solutions[0].objectives[1]
        max_mea = min_mea
        for s in solutions:
            mea = s.objectives[1]
            if mea < min_mea:
                min_mea = mea
            if mea > max_mea:
                max_mea = mea

        # normalize so the minimum MEA of the set is 0
        for i in range(len(solutions)):
            solutions[i].objectives[1] -= min_mea

        if len(different_2ds) > 1:
            return solutions
        else:
            print("[%s] \033[36mWARNING: ignoring this RNA, only one 2D solution is found.\033[0m" % (basename))
    else:
        print("[%s] \033[36mWARNING: file not found !\033[0m" % (basename))

def parse_biorseo(folder, basename, ext):
    solutions = []
    if os.path.isfile(os.path.join(folder, basename + ext)):
        rna = open(os.path.join(folder, basename + ext), "r")
        lines = rna.readlines()
        rna.close()
        different_2ds = []
        for s in lines[2:]:
            if s == '\n':
                continue
            splitted = s.split('\t')
            db2d = splitted[0].split(' ')[0]
            if db2d not in different_2ds:
                different_2ds.append(db2d)
            solutions.append(SecStruct(db2d, float(splitted[1]), float(splitted[2][:-1])))
        if len(different_2ds) > 1:
            return solutions
        else:
            print("[%s] \033[36mWARNING: ignoring this RNA, only one 2D solution is found.\033[0m" % (basename))
    else:
        print("[%s] \033[36mWARNING: file not found !\033[0m" % (basename))
    return None

def prettify_biorseo(code):
    name = ""
    if "bgsu" in code:
        name += "RNA 3D Motif Atlas + "
    else:
        name += "Rna3Dmotifs + "
    if "raw" in code:
        name += "Direct P.M."
    if "byp" in code:
        name += "BPairing"
    if "jar3d" in code:
        name += "Jar3d"
    # name += " + $f_{1" + code[-1] + "}$"
    return name

# Parse options
try:
    opts, args = getopt.getopt( sys.argv[1:], "", 
                             [  "biokop", 
                                "biorseo_desc_byp_A",
                                "biorseo_desc_byp_B",
                                "biorseo_desc_byp_C",
                                "biorseo_desc_byp_D",
                                "biorseo_bgsu_byp_A",
                                "biorseo_bgsu_byp_B",
                                "biorseo_bgsu_byp_C",
                                "biorseo_bgsu_byp_D",
                                "biorseo_desc_raw_A",
                                "biorseo_desc_raw_B",
                                "biorseo_bgsu_jar3d_A",
                                "biorseo_bgsu_jar3d_B",
                                "biorseo_bgsu_jar3d_C",
                                "biorseo_bgsu_jar3d_D",
                                "folder=",
                                "database=",
                                "output="
                             ])
except getopt.GetoptError as err:
    print(err)
    sys.exit(2)

results_folder = "."
extension = "all"
outputf = ""
for opt, arg in opts:
    if opt == "--biokop":
        extension = ".biok"
        parse = parse_biokop
    elif opt == "--folder":
        results_folder = arg
    elif opt == "--database":
        database = arg
    elif opt == "--output":
        outputf = arg
    else:
        extension = '.' + opt[2:]
        parse = parse_biorseo

RNAcontainer, _ = load_from_dbn(database)

if results_folder[-1] != '/':
    results_folder = results_folder + '/'
if outputf == "":
    outputf = results_folder
if outputf[-1] != '/':
    outputf = outputf + '/'

def process_extension(ax, pos, ext, nsolutions=False, xlabel="Best solution performs\nwell on obj1", ylabel="Best solution performs\n well on obj2"):
    points = []
    sizes = []
    for rna in RNAcontainer:
        # Extracting the predictions from the results file
        solutions = parse(results_folder, rna.basename_, ext)
        reference = SecStruct(rna.struct_, float("inf"), float("inf"))
        if solutions is None:
            continue
        pset = Pareto(solutions, reference)
        points.append(pset.get_normalized_coords())
        sizes.append(pset.n_pred)
        print("[%s] Loaded %d solutions in a Pareto set, max(obj1)=%f, max(obj2)=%f" % (rna.basename_, pset.n_pred, pset.max_obj1, pset.max_obj2))
    print("Loaded %d points on %d." % (len(points), len(RNAcontainer)))

    x = np.array([ p[0] for p in points ])
    y = np.array([ p[1] for p in points ])
    xmin, xmax = 0, 1
    ymin, ymax = 0, 1
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([x, y])
    kernel = st.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    ax[pos].axhline(y=0, alpha=0.2, color='black')
    ax[pos].axhline(y=1, alpha=0.2, color='black')
    ax[pos].axvline(x=0, alpha=0.2, color='black')
    ax[pos].axvline(x=1, alpha=0.2, color='black')
    ax[pos].contourf(xx, yy, f, cmap=cm.Blues, alpha=0.5)
    ax[pos].scatter(x, y, s=25, alpha=0.1)
    ax[pos].set_xlim((-0.1,1.1))
    ax[pos].set_ylim((-0.1,1.1))
    ax[pos].set_title(prettify_biorseo(ext[1:]), fontsize=10)
    ax[pos].annotate("("+str(len(points))+'/'+str(len(RNAcontainer))+" RNAs)", (0.08,0.15))
    ax[pos].set_xlabel(xlabel)
    ax[pos].set_ylabel(ylabel)

    if nsolutions:
        ax[pos+1].hist(sizes, bins=range(0, max(sizes)+1, 2), histtype='bar')
        ax[pos+1].set_xlim((0,max(sizes)+2))
        ax[pos+1].set_xticks(range(0, max(sizes), 10))
        ax[pos+1].set_xticklabels(range(0, max(sizes), 10), rotation=90)
        ax[pos+1].set_xlabel("# solutions")
        ax[pos+1].set_ylabel("# RNAs")

if extension == "all":
    parse = parse_biorseo
    fig, ax = plt.subplots(1,4,figsize=(10,3),sharey=True)
    ax = ax.flatten()
    process_extension(ax, 0, ".biorseo_desc_raw_A", xlabel="Normalized $f_{1A}$", ylabel="Normalized MEA")
    process_extension(ax, 1, ".biorseo_desc_byp_A", xlabel="Normalized $f_{1A}$", ylabel="Normalized MEA")
    process_extension(ax, 2, ".biorseo_bgsu_byp_A", xlabel="Normalized $f_{1A}$", ylabel="Normalized MEA")
    process_extension(ax, 3, ".biorseo_bgsu_jar3d_A", xlabel="Normalized $f_{1A}$", ylabel="Normalized MEA")
    for a in ax:
        a.label_outer()
    plt.subplots_adjust(bottom=0.2, top=0.9, left=0.07, right=0.98, hspace=0.05, wspace = 0.05)
    plt.show()

    fig, ax = plt.subplots(1,4,figsize=(10,3),sharey=True)
    ax = ax.flatten()
    process_extension(ax, 0, ".biorseo_desc_raw_B", xlabel="Normalized $f_{1B}$", ylabel="Normalized MEA")
    process_extension(ax, 1, ".biorseo_desc_byp_B", xlabel="Normalized $f_{1B}$", ylabel="Normalized MEA")
    process_extension(ax, 2, ".biorseo_bgsu_byp_B", xlabel="Normalized $f_{1B}$", ylabel="Normalized MEA")
    process_extension(ax, 3, ".biorseo_bgsu_jar3d_B", xlabel="Normalized $f_{1B}$", ylabel="Normalized MEA")
    for a in ax:
        a.label_outer()
    plt.subplots_adjust(bottom=0.2, top=0.9, left=0.07, right=0.98, hspace=0.05, wspace = 0.05)
    plt.show()

    fig, ax = plt.subplots(1,4,figsize=(10,3),sharey=True)
    ax = ax.flatten()
    process_extension(ax, 1, ".biorseo_desc_byp_C", xlabel="Normalized $f_{1C}$", ylabel="Normalized MEA")
    process_extension(ax, 2, ".biorseo_bgsu_byp_C", xlabel="Normalized $f_{1C}$", ylabel="Normalized MEA")
    process_extension(ax, 3, ".biorseo_bgsu_jar3d_C", xlabel="Normalized $f_{1C}$", ylabel="Normalized MEA")
    ax[0].axis("off")
    for a in ax:
        a.label_outer()
    plt.subplots_adjust(bottom=0.2, top=0.9, left=0.07, right=0.98, hspace=0.05, wspace = 0.05)
    plt.show()

    fig, ax = plt.subplots(1,4,figsize=(10,3),sharey=True)
    ax = ax.flatten()
    process_extension(ax, 1, ".biorseo_desc_byp_D", xlabel="Normalized $f_{1D}$", ylabel="Normalized MEA")
    process_extension(ax, 2, ".biorseo_bgsu_byp_D", xlabel="Normalized $f_{1D}$", ylabel="Normalized MEA")
    process_extension(ax, 3, ".biorseo_bgsu_jar3d_D", xlabel="Normalized $f_{1D}$", ylabel="Normalized MEA")
    ax[0].axis("off")    
    for a in ax:
        a.label_outer()
    plt.subplots_adjust(bottom=0.2, top=0.9, left=0.07, right=0.98, hspace=0.05, wspace = 0.05)
    plt.show()
    

    
else:
    fig, ax = plt.subplots(2,1, figsize=(6,5))
    plt.subplots_adjust(bottom=0.12, top=0.9, left=0.15, right=0.9, hspace=0.4)
    if extension == ".biok":
        process_extension(ax, 0, extension, nsolutions=True, xlabel="Normalized MFE", ylabel="Normalized MEA")
    else:
        process_extension(ax, 0, extension, nsolutions=True)
    plt.show()