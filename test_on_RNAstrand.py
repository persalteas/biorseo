#!/usr/bin/python3
#coding=utf-8
from sys import argv
from scipy import stats
import subprocess
from os import path, makedirs, getcwd, chdir, devnull
import matplotlib.pyplot as plt
from matplotlib import colors
from math import sqrt
from multiprocessing import Pool, cpu_count, Manager
import multiprocessing
import ast

# ================== DEFINITION OF THE PATHS ==============================

# Retrieve Paths from file EditMe
jar3dexec = ""
HLmotifDir = ""
ILmotifDir = ""
descfolder = ""
bypdir = ""
biorseoDir = "."
exec(compile(open(biorseoDir+"/EditMe").read(), '', 'exec'))
runDir = path.dirname(path.realpath(__file__))
dataFile = argv[1]
outputDir = biorseoDir + "/results/"

# Create some folders to store the results
subprocess.call(["mkdir", "-p", outputDir])
subprocess.call(["mkdir", "-p", outputDir + "PK/"])
subprocess.call(["mkdir", "-p", outputDir + "noPK/"])

m = Manager()
running_stats = m.list()
running_stats.append(0) # n_launched
running_stats.append(0) # n_finished
running_stats.append(0) # n_skipped
fails = m.list()

# ================== CLASSES AND FUNCTIONS ================================

class Job:
    def __init__(self, command=[], function=None, args=[], how_many_in_parallel=0, priority=1, timeout=None, checkFunc=None, checkArgs=[]):
        self.cmd_ = command
        self.func_ = function
        self.args_ = args
        self.checkFunc_ = checkFunc
        self.checkArgs_ = checkArgs
        self.priority_ = priority
        self.timeout_ = timeout
        if not how_many_in_parallel:
            self.nthreads = cpu_count()
        elif how_many_in_parallel == -1:
            self.nthreads = cpu_count() - 1
        else:
            self.nthreads = how_many_in_parallel


class NoDaemonProcess(multiprocessing.Process):
    @property
    def daemon(self):
        return False

    @daemon.setter
    def daemon(self, value):
        pass


class NoDaemonContext(type(multiprocessing.get_context())):
    Process = NoDaemonProcess

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(multiprocessing.pool.Pool):
    def __init__(self, *args, **kwargs):
        kwargs['context'] = NoDaemonContext()
        super(MyPool, self).__init__(*args, **kwargs)


def execute_job(j):

    if j.checkFunc_ is not None:
        if j.checkFunc_(*j.checkArgs_):
            running_stats[2] += 1
            print("["+str(running_stats[0]+running_stats[2])+'/'+str(jobcount)+"]\tSkipping a finished job")
            return 0
    running_stats[0] += 1
    if len(j.cmd_):
        logfile = open("log_of_the_run.sh", 'a')
        logfile.write(" ".join(j.cmd_))
        logfile.write("\n")
        logfile.close()
        print("["+str(running_stats[0]+running_stats[2])+'/'+str(jobcount)+"]\t"+" ".join(j.cmd_))
        r = subprocess.call(j.cmd_, timeout=j.timeout_)
    elif j.func_ is not None:
        print("["+str(running_stats[0]+running_stats[2])+'/'+str(jobcount)+"]\t"+j.func_.__name__+'('+", ".join([a for a in j.args_])+')')
        # try:
        r = j.func_(*j.args_)
        # except:
        #     r = 1
        #     pass
    #if r:
    #    fails.append(j)
    running_stats[1] += 1
    return r

def check_RNAsubopt(basename):
    return path.isfile(outputDir + basename + ".subopt")

def check_biorseoBGSUJAR3DA(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".jar3dA") 

def check_biorseoBGSUJAR3DC(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".jar3dC")

def check_biorseoBGSUJAR3DD(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".jar3dD")

def check_biorseoBGSUJAR3DB(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".jar3dB")

def check_biorseoBGSUBayesPairA(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".bgsubypA")

def check_biorseoBGSUBayesPairC(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".bgsubypC")

def check_biorseoBGSUBayesPairD(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".bgsubypD")

def check_biorseoBGSUBayesPairB(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".bgsubypB")

def check_biorseoBayesPairA(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".bypA")

def check_biorseoBayesPairC(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".bypC")

def check_biorseoBayesPairD(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".bypD")

def check_biorseoBayesPairB(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".bypB")

def check_biorseoRawA(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".rawA")

def check_biorseoRawB(basename, with_PK):
    folder = outputDir+"PK/" if with_PK else outputDir+"noPK/"
    return path.isfile(folder + basename + ".rawB")

def check_JAR3D(basename):
    return path.isfile(outputDir + basename + ".sites.csv")

def check_BayesPairing(basename):
    return path.isfile(outputDir + basename + ".byp.csv")

def check_BGSUBayesPairing(basename):
    return path.isfile(outputDir + basename + ".bgsubyp.csv")

def check_biokop(basename):
    return path.isfile(outputDir + basename + ".biok")

def check_RNAMoIP(basename):
    return path.isfile(outputDir + basename + ".moip")

def launch_JAR3D_worker(loop):
    # write motif to a file
    newpath = getcwd()+'/'+loop.header[1:]
    if not path.exists(newpath):
        makedirs(newpath)
    chdir(newpath)
    filename = loop.header[1:]+".fasta"
    fasta = open(filename, 'w')
    fasta.write('>'+loop.get_header()+'\n'+loop.subsequence()+'\n')
    fasta.close()

    # Launch Jar3D on it
    if loop.type == 'h':
        cmd = ["java", "-jar", jar3dexec, filename, HLmotifDir+"/all.txt", loop.header[1:]+".HLloop.csv", loop.header[1:]+".HLseq.csv"]
    else:
        cmd = ["java", "-jar", jar3dexec, filename, ILmotifDir+"/all.txt", loop.header[1:]+".ILloop.csv", loop.header[1:]+".ILseq.csv"]
    nowhere = open(devnull, 'w')
    logfile = open("log_of_the_run.sh", 'a')
    logfile.write(' '.join(cmd))
    logfile.write("\n")
    logfile.close()
    subprocess.call(cmd, stdout=nowhere)
    nowhere.close()

    # Retrieve results
    insertion_sites = []
    if loop.type == 'h':
        capstype = "HL"
    else:
        capstype = "IL"
    csv = open(loop.header[1:]+".%sseq.csv" % capstype, 'r')
    l = csv.readline()
    while l:
        if "true" in l:
            insertion_sites.append(InsertionSite(loop, l))
        l = csv.readline()
    csv.close()

    # Cleaning
    chdir("..")
    subprocess.call(["rm", "-r", loop.header[1:]])
    return insertion_sites

def launch_JAR3D(seq_, basename):
    rnasubopt_preds = []
    # Extracting probable loops from RNA-subopt structures
    rna = open(outputDir + basename + ".subopt", "r")
    lines = rna.readlines()
    rna.close()
    for i in range(2, len(lines)):
        ss = lines[i].split(' ')[0]
        if ss not in rnasubopt_preds:
            rnasubopt_preds.append(ss)
    HLs = []
    ILs = []
    for ss in rnasubopt_preds:
        loop_candidates = enumerate_loops(ss)
        for loop_candidate in loop_candidates:
            if len(loop_candidate) == 1 and loop_candidate not in HLs:
                HLs.append(loop_candidate)
            if len(loop_candidate) == 2 and loop_candidate not in ILs:
                ILs.append(loop_candidate)
    # Retrieve subsequences corresponding to the possible loops
    loops = []
    for i, l in enumerate(HLs):
        loops.append(
            Loop(">HL%d" % (i+1), seq_[l[0][0]-1:l[0][1]], "h", l))
    for i, l in enumerate(ILs):
        loops.append(
            Loop(">IL%d" % (i+1), seq_[l[0][0]-1:l[0][1]]+'*'+seq_[l[1][0]-1:l[1][1]], "i", l))
    # Scanning loop subsequences against motif database
    pool = MyPool(processes=cpu_count())
    insertion_sites = [x for y in pool.map(launch_JAR3D_worker, loops) for x in y]
    insertion_sites.sort(reverse=True)
    # Writing results to CSV file
    c = 0
    resultsfile = open(outputDir+basename+".sites.csv", "w")
    resultsfile.write("Motif,Rotation,Score,Start1,End1,Start2,End2\n")
    for site in insertion_sites:
        if site.score > 10:
            c += 1
            string = "FOUND with score %d:\t\t possible insertion of motif " % site.score + site.atlas_id
            if site.rotation:
                string += " (reversed)"
            string += (" on " + site.loop.get_header() + " at positions")
        resultsfile.write(site.atlas_id+',' + str(bool(site.rotation))+",%d" % site.score+',')
        positions = [','.join([str(y) for y in x]) for x in site.position]
        if len(positions) == 1:
            positions.append("-,-")
        resultsfile.write(','.join(positions)+'\n')
    resultsfile.close()

def launch_BayesPairing(module_type, seq_, header_, basename):
    chdir(bypdir)

    cmd = ["python3","parse_sequences.py","-seq",outputDir + basename + ".fa", "-d", module_type, "-interm","1"]

    logfile = open("log_of_the_run.sh", 'a')
    logfile.write(" ".join(cmd))
    logfile.write("\n")
    logfile.close()

    out = subprocess.check_output(cmd).decode('utf-8')
    BypLog = out.split('\n')
    idx = 0
    l = BypLog[idx]
    while l[:3] != "PUR":
        idx += 1
        l = BypLog[idx]
    insertion_sites = [ x for x in ast.literal_eval(l.split(":")[1][1:])]
    if module_type=="rna3dmotif":
        rna = open(outputDir + basename + ".byp.csv", "w")
    else:
        rna = open(outputDir + basename + ".bgsubyp.csv", "w")
    rna.write("Motif,Score,Start1,End1,Start2,End2...\n")
    for i,module in enumerate(insertion_sites):
        if len(module):
            for (score, positions, sequence) in zip(*[iter(module)]*3):
                pos = []
                q = -2
                for p in positions:
                    if p-q > 1:
                        pos.append(q)
                        pos.append(p)
                    q = p
                pos.append(q)
                rna.write(module_type+str(i)+','+str(int(score)))
                for (p,q) in zip(*[iter(pos[1:])]*2):
                    if q>p:
                        rna.write(','+str(p)+','+str(q))
                rna.write('\n')
    rna.close()

def launch_RNAMoIP_worker(x):
    RNAMoIP = "../RNAMoIP/RNAMoIP.py"
    logfile = open("log_of_the_run.sh", 'a')
    logfile.write(" ".join(["gurobi.sh", RNAMoIP, "-s", '"' +x[1]+'"', "-ss", '"'+x[0].strip()+'"', "-d", descfolder]))
    logfile.write("\n")
    logfile.close()


    out = subprocess.check_output(["gurobi.sh", RNAMoIP, "-s", x[1], "-ss", x[0].strip(), "-d", descfolder]).decode("utf-8")
    gurobiLog = out.split('\n')
    idx = 0
    l = gurobiLog[idx]
    solution = ""
    while l != "Corrected secondary structure:" and l != " NO SOLUTIONS!":
        idx += 1
        l = gurobiLog[idx]
    if l == "Corrected secondary structure:":
        idx+=1
        solution =  gurobiLog[idx][1:]
        idx += 1
        motifs = []
        while gurobiLog[idx].count('.'):
            motif = gurobiLog[idx].split('-')[1]
            if motif not in motifs:
                motifs.append(motif)
            idx += 1
        nmotifs = len(motifs)
        score = float(gurobiLog[-2][1:-1])
    else:
        solution = ""
        nmotifs = 0
        score = 0
    return (solution, nmotifs, score)

def launch_RNAMoIP(seq_, header_, basename):
    rnasubopt_preds = []
    rna = open(outputDir + basename + ".subopt", "r")
    lines = rna.readlines()
    rna.close()
    for i in range(2, len(lines)):
        ss = lines[i].split(' ')[0]
        if ss not in rnasubopt_preds:
            rnasubopt_preds.append(ss)
    pool = MyPool(processes=cpu_count())
    results = [x for x in pool.map(launch_RNAMoIP_worker, zip([p for p in rnasubopt_preds], [seq_[:-1] for p in rnasubopt_preds]))]
    predictions = [ t[0] for t in results if t[0] != ""] 
    ninsertions = [ t[1] for t in results if t[0] != ""]
    scores = [ t[2] for t in results if t[0] != ""]
    rna = open(outputDir + basename + ".moip", "w")
    rna.write(header_)
    rna.write(seq_)
    for p,n,s in zip(predictions, ninsertions, scores):
        rna.write(p+'\t'+str(n)+'\t'+str(s)+'\n')
    rna.close()

def launch_pKiss(seq_, header_, basename):
    json = "{\"pkiss_input_rna_sequences\":\">%s\r\n%s\",\"paramset\":{\"pkiss_parameter_absoluteDeviation\":\"0.5\",\"pkiss_parameter_maxKnotSize\":\"3.0\",\"pkiss_parameter_windowSize\":\"1.0\",\"pkiss_parameter_param\":\"rna_andronescu2007\"}}" %(header_, seq_)
    cmd = "curl -X POST -d @[[%s]] http://bibiserv2.cebitec.uni-bielefeld.de:80/rest/pkiss/pkiss_function_subopt/request -H \"Content-Type: application/json\"" % json
    logfile = open("log_of_the_run.sh", 'a')
    logfile.write(cmd+"\n")
    logfile.close()
    print(cmd)

def mattews_corr_coeff(tp, tn, fp, fn):
    if (tp+fp == 0):
        print("We have an issue : no positives detected ! (linear structure)")
    return (tp*tn-fp*fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))

def accuracy(tp, tn, fp, fn):
    return (tp+tn)/(tp+fp+tn+fn)

def recall_sensitivity(tp, tn, fp, fn):
    return tp/(tp+fn)

def specificity(tp, tn, fp, fn):
    return tn/(tn+fp)

def precision_ppv(tp, tn, fp, fn):
    return tp/(tp+fp)

def npv(tp, tn, fp, fn):
    return tn/(tn+fn)

def f1_score(tp, tn, fp, fn):
    return 2*precision_ppv(tp, tn, fp, fn)*recall_sensitivity(tp, tn, fp, fn)/(precision_ppv(tp, tn, fp, fn)+recall_sensitivity(tp, tn, fp, fn))

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

def enumerate_loops(s):
    def resort(unclosedLoops):
        loops.insert(len(loops)-1-unclosedLoops, loops[-1])
        loops.pop(-1)

    opened = []
    openingStart = []
    closingStart = []
    loops = []
    loopsUnclosed = 0
    consecutiveOpenings = []
    if s[0] == '(':
        consecutiveOpenings.append(1)
    consecutiveClosings = 0

    lastclosed = -1
    previous = ''
    for i in range(len(s)):

        # If we arrive on an unpaired segment
        if s[i] == '.':
            if previous == '(':
                openingStart.append(i-1)
            if previous == ')':
                closingStart.append(i-1)

        # Opening basepair
        if s[i] == '(':
            if previous == '(':
                consecutiveOpenings[-1] += 1
            else:
                consecutiveOpenings.append(1)
            if previous == ')':
                closingStart.append(i-1)

            # We have something like (...(
            if len(openingStart) and openingStart[-1] == opened[-1]:
                # Create a new loop starting with this component.
                loops.append([(openingStart[-1], i)])
                openingStart.pop(-1)
                loopsUnclosed += 1
            # We have something like )...( or even )(
            if len(closingStart) and closingStart[-1] == lastclosed:
                # Append a component to existing multiloop
                loops[-1].append((closingStart[-1], i))
                closingStart.pop(-1)

            opened.append(i)

        # Closing basepair
        if s[i] == ')':
            if previous == ')':
                consecutiveClosings += 1
            else:
                consecutiveClosings = 1
            # This is not supposed to happen in real data, but whatever.
            if previous == '(':
                openingStart.append(i-1)

            # We have something like (...) or ()
            if len(openingStart) and openingStart[-1] == opened[-1]:
                # Create a new loop, and save it as already closed (HL)
                loops.append([(openingStart[-1], i)])
                openingStart.pop(-1)
                resort(loopsUnclosed)
            # We have something like )...)
            if len(closingStart) and closingStart[-1] == lastclosed:
                # Append a component to existing multiloop and close it.
                loops[-1].append((closingStart[-1], i))
                closingStart.pop(-1)
                loopsUnclosed -= 1
                resort(loopsUnclosed)

            if i+1 < len(s):
                if s[i+1] != ')':  # We are on something like: ).
                    # an openingStart has not been correctly detected, like in ...((((((...)))...)))
                    if consecutiveClosings < consecutiveOpenings[-1]:
                        # Create a new loop (uncompleted)
                        loops.append([(opened[-2], opened[-1])])
                        loopsUnclosed += 1

                    # We just completed an HL+stem, like ...(((...))).., we can forget its info
                    if consecutiveClosings == consecutiveOpenings[-1]:
                        consecutiveClosings = 0
                        consecutiveOpenings.pop(-1)
                    else:  # There are still several basepairs to remember, forget only the processed ones, keep the others
                        consecutiveOpenings[-1] -= consecutiveClosings
                        consecutiveClosings = 0

                else:  # We are on something like: ))
                    # we are on an closingStart that cannot be correctly detected, like in ...(((...(((...))))))
                    if consecutiveClosings == consecutiveOpenings[-1]:
                        # Append a component to the uncomplete loop and close it.
                        loops[-1].append((i, i+1))
                        loopsUnclosed -= 1
                        resort(loopsUnclosed)
                        # Forget the info about the processed stem.
                        consecutiveClosings = 0
                        consecutiveOpenings.pop(-1)

            opened.pop(-1)
            lastclosed = i

        previous = s[i]
        # print(i,"=",s[i],"\t", "consec. Op=", consecutiveOpenings,"Cl=",consecutiveClosings)

    return(loops)

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

def is_all(n, tot):
    if n == tot:
        return "\033[32m%d\033[0m/%d" % (n, tot)
    else:
        return "\033[91m%d\033[0m/%d" % (n, tot)

class Loop:
    def __init__(self, header, subsequence, looptype, position):
        self.header = header
        self.seq = subsequence
        self.type = looptype
        self.position = position

    def get_header(self):
        return self.header

    def subsequence(self):
        return self.seq


class InsertionSite:
    def __init__(self, loop, csv_line):
        # BEWARE : jar3d csv output is crap because of java's locale settings.
        # On french OSes, it uses commas to delimit the fields AND as floating point delimiters !!
        # Parse with caution, and check what the csv output files look like on your system...
        info = csv_line.split(',')
        self.loop = loop  # the Loop object that has been searched with jar3d
        # position of the loop's components, so the motif's ones, in the query sequence.
        self.position = loop.position
        # Motif model identifier of the RNA 3D Motif Atlas
        self.atlas_id = info[2]
        # alignment score of the subsequence to the motif model
        self.score = int(float(info[4]))
        # should the motif model be inverted to fit the sequence ?
        self.rotation = int(info[-2])

    def __lt__(self, other):
        return self.score < other.score

    def __gt__(self, other):
        return self.score > other.score


class Method:
    def __init__(self):
        self.predictions = []
        self.scores = []
        self.ninsertions = []
        self.max_mcc = 0
        self.min_mcc = 0
        self.avg_mcc = 0
        self.best_pred = ""
        self.n_pred = 0
        self.ratio = 0  # ratio of the number of inserted motifs in the best solution on the max number of inserted motifs for this RNA


class RNA:
    def __init__(self, filename, header, seq, struct):
        self.seq_ = seq
        self.header_ = header
        self.true2d = struct
        self.basename = filename

        self.rnasubopt = Method()
        self.biokop = Method()
        self.rnamoip = Method()
        self.biorseoRawA = Method()
        self.biorseoRawB = Method()
        self.biorseoBGSUJAR3DA = Method()
        self.biorseoBGSUJAR3DC = Method()
        self.biorseoBGSUJAR3DD = Method()
        self.biorseoBGSUJAR3DB = Method()
        self.biorseoBayesPairA = Method()
        self.biorseoBayesPairC = Method()
        self.biorseoBayesPairD = Method()
        self.biorseoBayesPairB = Method()
        self.biorseoBGSUBayesPairA = Method()
        self.biorseoBGSUBayesPairC = Method()
        self.biorseoBGSUBayesPairD = Method()
        self.biorseoBGSUBayesPairB = Method()

        if not path.isfile(outputDir + self.basename + ".fa"):
            rna = open(outputDir + self.basename + ".fa", "w")
            rna.write(">"+self.header_+'\n')
            rna.write(self.seq_+'\n')
            rna.close()
        rna = open(outputDir + "allsequences.fa", "a")
        rna.write(">"+self.header_+'\n')
        rna.write(self.seq_+'\n')
        rna.close()
        rna = open(outputDir + "allsequences.dbn", "a")
        rna.write(">"+self.header_+'\n')
        rna.write(self.seq_+'\n')
        rna.write(self.true2d + '\n')
        rna.close()
        if "TRNA" in self.header_:
            rna = open(outputDir + "tRNAs.dbn", "a")
            rna.write(">"+self.header_+'\n')
            rna.write(self.seq_+'\n')
            rna.write(self.true2d + '\n')
            rna.close()

    def evaluate(self):

        methods = [self.rnasubopt, self.biokop, self.rnamoip,
                   self.biorseoBayesPairA, 
                   self.biorseoBayesPairB, 
                   self.biorseoBayesPairC, 
                   self.biorseoBayesPairD,
                   self.biorseoRawA, 
                   self.biorseoRawB,
                   self.biorseoBGSUJAR3DA, 
                   self.biorseoBGSUJAR3DB, 
                   self.biorseoBGSUJAR3DC, 
                   self.biorseoBGSUJAR3DD,
                   self.biorseoBGSUBayesPairA, 
                   self.biorseoBGSUBayesPairB, 
                   self.biorseoBGSUBayesPairC, 
                   self.biorseoBGSUBayesPairD 
        ]

        for m in methods:
            if len(m.predictions):
                mccs = []
                m.n_pred = len(m.predictions)
                sec_structs = [] # store the dot-brackets to check for redundancy
                for p in m.predictions:
                    if not ')' in p:
                        m.n_pred -= 1
                        continue
                    ss = p.split('\t')[0].split(' ')[0]
                    if ss not in sec_structs:
                        sec_structs.append(p.split('\t')[0])
                    else:
                        m.n_pred -= 1
                        continue                    
                    mccs.append(mattews_corr_coeff(*compare_two_structures(self.true2d, p)))

                if len(mccs):
                    m.max_mcc = max(mccs)
                    m.min_mcc = min(mccs)
                    m.avg_mcc = sum(mccs)/float(len(mccs))
                    m.best_pred = sec_structs[mccs.index(m.max_mcc)]
                for p,n in zip(m.predictions, m.ninsertions):
                    if not ')' in p:
                        continue
                    if m.max_mcc == mattews_corr_coeff(*compare_two_structures(self.true2d, p)):
                        m.best_pred = p
                        if max(m.ninsertions) > 0 and float(n)/max(m.ninsertions) > m.ratio:
                            m.ratio = float(n)/max(m.ninsertions)
            
    def get_biokop_results(self):
        if path.isfile(outputDir + self.basename + ".biok"):
            rna = open(outputDir + self.basename + ".biok", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(1, len(lines)-1):
                ss = lines[i].split(' ')[0]
                if ss not in self.biokop.predictions:
                    self.biokop.predictions.append(ss)

    def get_RNAsubopt_results(self):
        rna = open(outputDir + self.basename + ".subopt", "r")
        lines = rna.readlines()
        rna.close()
        for i in range(2, len(lines)):
            ss = lines[i].split(' ')[0]
            if ss not in self.rnasubopt.predictions:
                self.rnasubopt.predictions.append(ss)

    def get_RNAMoIP_results(self):
        rna = open(outputDir + self.basename + ".moip", "r")
        lines = rna.readlines()
        rna.close()
        for i in range(2, len(lines)):
            self.rnamoip.predictions.append(lines[i].split('\t')[0])
            self.rnamoip.ninsertions.append(int(lines[i].split('\t')[1]))
            self.rnamoip.scores.append(float(lines[i].split('\t')[2][:-1]))

    def get_biorseoBayesPairA_results(self, targetdir):  
        if path.isfile(targetdir+ self.basename + ".bypA"):
            rna = open(targetdir+ self.basename + ".bypA", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBayesPairA.predictions:
                    self.biorseoBayesPairA.predictions.append(ss)
                self.biorseoBayesPairA.ninsertions.append(lines[i].count('+'))
    
    def get_biorseoBayesPairB_results(self, targetdir):
        if path.isfile(targetdir+ self.basename + ".bypB"):
            rna = open(targetdir+ self.basename + ".bypB", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBayesPairB.predictions:
                    self.biorseoBayesPairB.predictions.append(ss)
                self.biorseoBayesPairB.ninsertions.append(lines[i].count('+'))

    def get_biorseoBayesPairC_results(self, targetdir):  
        if path.isfile(targetdir+ self.basename + ".bypC"):
            rna = open(targetdir+ self.basename + ".bypC", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBayesPairC.predictions:
                    self.biorseoBayesPairC.predictions.append(ss)
                self.biorseoBayesPairC.ninsertions.append(lines[i].count('+'))

    def get_biorseoBayesPairD_results(self, targetdir):  
        if path.isfile(targetdir+ self.basename + ".bypD"):
            rna = open(targetdir+ self.basename + ".bypD", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBayesPairD.predictions:
                    self.biorseoBayesPairD.predictions.append(ss)
                self.biorseoBayesPairD.ninsertions.append(lines[i].count('+'))

    def get_biorseoRawA_results(self, targetdir):
        if path.isfile(targetdir+ self.basename + ".rawA"):
            rna = open(targetdir+ self.basename + ".rawA", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoRawA.predictions:
                    self.biorseoRawA.predictions.append(ss)
                self.biorseoRawA.ninsertions.append(lines[i].count('+'))

    def get_biorseoRawB_results(self, targetdir):
        if path.isfile(targetdir+ self.basename + ".rawB"):
            rna = open(targetdir+ self.basename + ".rawB", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoRawB.predictions:
                    self.biorseoRawB.predictions.append(ss)
                self.biorseoRawB.ninsertions.append(lines[i].count('+'))

    def get_biorseoBGSUJAR3DA_results(self, targetdir):  
        if path.isfile(targetdir+ self.basename + ".jar3dA"):
            rna = open(targetdir+ self.basename + ".jar3dA", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBGSUJAR3DA.predictions:
                    self.biorseoBGSUJAR3DA.predictions.append(ss)
                self.biorseoBGSUJAR3DA.ninsertions.append(lines[i].count('+'))
    
    def get_biorseoBGSUJAR3DB_results(self, targetdir):
        if path.isfile(targetdir+ self.basename + ".jar3dB"):
            rna = open(targetdir+ self.basename + ".jar3dB", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBGSUJAR3DB.predictions:
                    self.biorseoBGSUJAR3DB.predictions.append(ss)
                self.biorseoBGSUJAR3DB.ninsertions.append(lines[i].count('+'))

    def get_biorseoBGSUJAR3DC_results(self, targetdir):  
        if path.isfile(targetdir+ self.basename + ".jar3dC"):
            rna = open(targetdir+ self.basename + ".jar3dC", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBGSUJAR3DC.predictions:
                    self.biorseoBGSUJAR3DC.predictions.append(ss)
                self.biorseoBGSUJAR3DC.ninsertions.append(lines[i].count('+'))

    def get_biorseoBGSUJAR3DD_results(self, targetdir):  
        if path.isfile(targetdir+ self.basename + ".jar3dD"):
            rna = open(targetdir+ self.basename + ".jar3dD", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBGSUJAR3DD.predictions:
                    self.biorseoBGSUJAR3DD.predictions.append(ss)
                self.biorseoBGSUJAR3DD.ninsertions.append(lines[i].count('+'))

    def get_biorseoBGSUBayesPairA_results(self, targetdir):  
        if path.isfile(targetdir+ self.basename + ".bgsubypA"):
            rna = open(targetdir+ self.basename + ".bgsubypA", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBGSUBayesPairA.predictions:
                    self.biorseoBGSUBayesPairA.predictions.append(ss)
                self.biorseoBGSUBayesPairA.ninsertions.append(lines[i].count('+'))
        else:
            print(targetdir+ self.basename + ".bgsubypA not found !")
    
    def get_biorseoBGSUBayesPairB_results(self, targetdir):
        if path.isfile(targetdir+ self.basename + ".bgsubypB"):
            rna = open(targetdir+ self.basename + ".bgsubypB", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBGSUBayesPairB.predictions:
                    self.biorseoBGSUBayesPairB.predictions.append(ss)
                self.biorseoBGSUBayesPairB.ninsertions.append(lines[i].count('+'))
        else:
            print(targetdir+ self.basename + ".bgsubypB not found !")

    def get_biorseoBGSUBayesPairC_results(self, targetdir):  
        if path.isfile(targetdir+ self.basename + ".bgsubypC"):
            rna = open(targetdir+ self.basename + ".bgsubypC", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBGSUBayesPairC.predictions:
                    self.biorseoBGSUBayesPairC.predictions.append(ss)
                self.biorseoBGSUBayesPairC.ninsertions.append(lines[i].count('+'))
        else:
            print(targetdir+ self.basename + ".bgsubypC not found !")

    def get_biorseoBGSUBayesPairD_results(self, targetdir):  
        if path.isfile(targetdir+ self.basename + ".bgsubypD"):
            rna = open(targetdir+ self.basename + ".bgsubypD", "r")
            lines = rna.readlines()
            rna.close()
            for i in range(2, len(lines)):
                ss = lines[i].split(' ')[0].split('\t')[0]
                if ss not in self.biorseoBGSUBayesPairD.predictions:
                    self.biorseoBGSUBayesPairD.predictions.append(ss)
                self.biorseoBGSUBayesPairD.ninsertions.append(lines[i].count('+'))
        else:
            print(targetdir+ self.basename + ".bgsubypD not found !")

    def load_results_from(self, targetDir):
        self.get_biokop_results()
        self.get_RNAsubopt_results()
        self.get_RNAMoIP_results()
        self.get_biorseoBayesPairA_results(targetDir)
        self.get_biorseoBayesPairB_results(targetDir)
        self.get_biorseoBayesPairC_results(targetDir)
        self.get_biorseoBayesPairD_results(targetDir)
        self.get_biorseoRawA_results(targetDir)
        self.get_biorseoRawB_results(targetDir)
        self.get_biorseoBGSUJAR3DA_results(targetDir)
        self.get_biorseoBGSUJAR3DB_results(targetDir)
        self.get_biorseoBGSUJAR3DC_results(targetDir)
        self.get_biorseoBGSUJAR3DD_results(targetDir)
        self.get_biorseoBGSUBayesPairA_results(targetDir)
        self.get_biorseoBGSUBayesPairB_results(targetDir)
        self.get_biorseoBGSUBayesPairC_results(targetDir)
        self.get_biorseoBGSUBayesPairD_results(targetDir)

    def has_complete_results(self, with_PK):
        if not with_PK and not check_RNAsubopt(self.basename): return False
        if not with_PK and not check_RNAMoIP(self.basename): return False
        if with_PK and not check_biokop(self.basename): return False
        if not check_biorseoBayesPairA(self.basename, with_PK): return False
        if not check_biorseoBayesPairB(self.basename, with_PK): return False
        if not check_biorseoBayesPairC(self.basename, with_PK): return False
        if not check_biorseoBayesPairD(self.basename, with_PK): return False
        if not check_biorseoRawA(self.basename, with_PK): return False
        if not check_biorseoRawB(self.basename, with_PK): return False
        if not check_biorseoBGSUJAR3DA(self.basename, with_PK): return False
        if not check_biorseoBGSUJAR3DB(self.basename, with_PK): return False
        if not check_biorseoBGSUJAR3DC(self.basename, with_PK): return False
        if not check_biorseoBGSUJAR3DD(self.basename, with_PK): return False
        if not check_biorseoBGSUBayesPairA(self.basename, with_PK): return False
        if not check_biorseoBGSUBayesPairB(self.basename, with_PK): return False
        if not check_biorseoBGSUBayesPairC(self.basename, with_PK): return False
        if not check_biorseoBGSUBayesPairD(self.basename, with_PK): return False
        return True

# ================= EXTRACTION OF STRUCTURES FROM DATABASE ===============================
RNAcontainer = []
pk_counter = 0
print("loading files...")

db = open(dataFile, "r")
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
            RNAcontainer.append(RNA(header.replace('/', '_').split('(')[-1][:-1], header, seq, struct))
            #RNAcontainer.append(RNA(header.replace('/', '_').split('[')[-1][:-41], header, seq, struct))
            if '[' in struct: pk_counter += 1
db.close()

for nt, number in ignored_nt_dict.items():
    print("ignored %d sequences because of char %c" % (number, nt))
tot = len(RNAcontainer)
print("Loaded %d RNAs of length between 10 and 100. %d of them contain pseudoknots." % (tot, pk_counter))

# # ================= PREDICTION OF STRUCTURES ===============================

# # define job list
# joblist = []
# for instance in RNAcontainer:
#     basename = instance.basename
#     # RNAsubopt
#     joblist.append(Job(command=["RNAsubopt", "-i", outputDir + basename + ".fa", "--outfile="+ basename + ".subopt"], priority=1, checkFunc=check_RNAsubopt, checkArgs=[basename]))
#     joblist.append(Job(command=["mv", basename + ".subopt", outputDir], priority=2, checkFunc=check_RNAsubopt, checkArgs=[basename]))
#     # JAR3D
#     joblist.append(Job(function=launch_JAR3D, args=[instance.seq_, basename], priority=3, how_many_in_parallel=1, checkFunc=check_JAR3D, checkArgs=[basename]))
#     # BayesPairing and BGSUBayesPairing
#     joblist.append(Job(function=launch_BayesPairing, args=["rna3dmotif", instance.seq_, instance.header_, basename], how_many_in_parallel=-1, priority=3, checkFunc=check_BayesPairing, checkArgs=[basename]))
#     joblist.append(Job(function=launch_BayesPairing, args=["3dmotifatlas", instance.seq_, instance.header_, basename], how_many_in_parallel=-1, priority=3, checkFunc=check_BGSUBayesPairing, checkArgs=[basename]))
#     # bmotinBGSUJAR3DA-D
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--jar3dcsv", outputDir+basename+".sites.csv", "-o", outputDir+"noPK/"+basename+".jar3dA", "--type", "A", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUJAR3DA, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--jar3dcsv", outputDir+basename+".sites.csv", "-o", outputDir+"noPK/"+basename+".jar3dB", "--type", "B", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUJAR3DB, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--jar3dcsv", outputDir+basename+".sites.csv", "-o", outputDir+"noPK/"+basename+".jar3dC", "--type", "C", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUJAR3DC, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--jar3dcsv", outputDir+basename+".sites.csv", "-o", outputDir+"noPK/"+basename+".jar3dD", "--type", "D", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUJAR3DD, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--jar3dcsv", outputDir+basename+".sites.csv", "-o", outputDir+"PK/"+basename+".jar3dA", "--type", "A"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUJAR3DA, checkArgs=[basename, True]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--jar3dcsv", outputDir+basename+".sites.csv", "-o", outputDir+"PK/"+basename+".jar3dB", "--type", "B"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUJAR3DB, checkArgs=[basename, True]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--jar3dcsv", outputDir+basename+".sites.csv", "-o", outputDir+"PK/"+basename+".jar3dC", "--type", "C"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUJAR3DC, checkArgs=[basename, True]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--jar3dcsv", outputDir+basename+".sites.csv", "-o", outputDir+"PK/"+basename+".jar3dD", "--type", "D"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUJAR3DD, checkArgs=[basename, True]))
#     # bmotinBGSUBayesPairA-D
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".bgsubyp.csv", "-o", outputDir+"noPK/"+basename+".bgsubypA", "--type", "A", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUBayesPairA, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".bgsubyp.csv", "-o", outputDir+"noPK/"+basename+".bgsubypB", "--type", "B", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUBayesPairB, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".bgsubyp.csv", "-o", outputDir+"noPK/"+basename+".bgsubypC", "--type", "C", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUBayesPairC, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".bgsubyp.csv", "-o", outputDir+"noPK/"+basename+".bgsubypD", "--type", "D", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUBayesPairD, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".bgsubyp.csv", "-o", outputDir+"PK/"+basename+".bgsubypA", "--type", "A"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUBayesPairA, checkArgs=[basename, True]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".bgsubyp.csv", "-o", outputDir+"PK/"+basename+".bgsubypB", "--type", "B"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUBayesPairB, checkArgs=[basename, True]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".bgsubyp.csv", "-o", outputDir+"PK/"+basename+".bgsubypC", "--type", "C"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUBayesPairC, checkArgs=[basename, True]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".bgsubyp.csv", "-o", outputDir+"PK/"+basename+".bgsubypD", "--type", "D"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBGSUBayesPairD, checkArgs=[basename, True]))
#     # bmotinBayesPairA-D
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".byp.csv", "-o", outputDir+"noPK/"+basename+".bypA", "--type", "A", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBayesPairA, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".byp.csv", "-o", outputDir+"noPK/"+basename+".bypB", "--type", "B", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBayesPairB, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".byp.csv", "-o", outputDir+"noPK/"+basename+".bypC", "--type", "C", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBayesPairC, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".byp.csv", "-o", outputDir+"noPK/"+basename+".bypD", "--type", "D", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBayesPairD, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".byp.csv", "-o", outputDir+"PK/"+basename+".bypA", "--type", "A"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBayesPairA, checkArgs=[basename, True]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".byp.csv", "-o", outputDir+"PK/"+basename+".bypB", "--type", "B"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBayesPairB, checkArgs=[basename, True]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".byp.csv", "-o", outputDir+"PK/"+basename+".bypC", "--type", "C"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBayesPairC, checkArgs=[basename, True]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir+basename+".fa", "--bayespaircsv", outputDir+basename+".byp.csv", "-o", outputDir+"PK/"+basename+".bypD", "--type", "D"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoBayesPairD, checkArgs=[basename, True]))
#     # biorseoRawA,B
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir + basename + ".fa", "-d", descfolder, "-o", outputDir+"noPK/" + basename + ".rawA", "--type", "A", "-n"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoRawA, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir + basename + ".fa", "-d", descfolder, "-o", outputDir+"noPK/" + basename + ".rawB", "--type", "B", "-n"], priority=4, timeout=3600,  how_many_in_parallel=3, checkFunc=check_biorseoRawB, checkArgs=[basename, False]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir + basename + ".fa", "-d", descfolder, "-o", outputDir+"PK/" + basename + ".rawA", "--type", "A"], priority=4, timeout=3600, how_many_in_parallel=3, checkFunc=check_biorseoRawA, checkArgs=[basename, True]))
#     joblist.append(Job(command=[biorseoDir+"/bin/biorseo", "-s", outputDir + basename + ".fa", "-d", descfolder, "-o", outputDir+"PK/" + basename + ".rawB", "--type", "B"], priority=4, timeout=3600,  how_many_in_parallel=3, checkFunc=check_biorseoRawB, checkArgs=[basename, True]))
#     # RNA MoIP
#     joblist.append(Job(function=launch_RNAMoIP, args=[instance.seq_, instance.header_, basename], priority=3, timeout=3600, checkFunc=check_RNAMoIP, checkArgs=[basename]))
#     # Biokop
#     joblist.append(Job(command=[biorseoDir + "/../biokop/biokop", "-n1", "-i", outputDir + basename + ".fa", "-o", outputDir + basename + ".biok"], priority=5, timeout=15000, how_many_in_parallel=3, checkFunc=check_biokop, checkArgs=[basename]))


# # execute jobs
# jobs = {}
# jobcount = len(joblist)
# for job in joblist:
#     if job.priority_ not in jobs.keys():
#         jobs[job.priority_] = {}
#     if job.nthreads not in jobs[job.priority_].keys():
#         jobs[job.priority_][job.nthreads] = []
#     jobs[job.priority_][job.nthreads].append(job)
# nprio = max(jobs.keys())


# for i in range(1,nprio+1):
#     if not len(jobs[i].keys()): continue

#     # check the thread numbers
#     different_thread_numbers = [n for n in jobs[i].keys()]
#     different_thread_numbers.sort()

#     for n in different_thread_numbers:
#         bunch = jobs[i][n]
#         if not len(bunch): continue
#         pool = MyPool(processes=n)
#         results = pool.map(execute_job, bunch)
#         pool.close()
#         pool.join()

# if len(fails):
#     print()
#     print("Some jobs failed! :")
#     print()
#     for j in fails:
#         print(j.cmd_)

# exit()



# ================= Statistics (without pseudoknots) ========================

# load results in objects 
for instance in RNAcontainer:
    instance.load_results_from(outputDir + "noPK/")
    instance.evaluate()

RNAs_fully_predicted = [ x for x in RNAcontainer if x.has_complete_results(False)]

x_noPK = [
    [ rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.rnasubopt.predictions)],
    [ rna.rnamoip.max_mcc for rna in RNAcontainer if len(rna.rnamoip.predictions)],
    [ rna.biorseoBGSUJAR3DA.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DA.predictions)],
    [ rna.biorseoBGSUJAR3DB.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DB.predictions)],
    [ rna.biorseoBGSUJAR3DC.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DC.predictions)],
    [ rna.biorseoBGSUJAR3DD.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DD.predictions)],
    [ rna.biorseoBGSUBayesPairA.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairA.predictions)],
    [ rna.biorseoBGSUBayesPairB.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairB.predictions)],
    [ rna.biorseoBGSUBayesPairC.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairC.predictions)],
    [ rna.biorseoBGSUBayesPairD.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairD.predictions)],
    [ rna.biorseoRawA.max_mcc for rna in RNAcontainer if len(rna.biorseoRawA.predictions)],
    [ rna.biorseoRawB.max_mcc for rna in RNAcontainer if len(rna.biorseoRawB.predictions)],
    [ rna.biorseoBayesPairA.max_mcc  for rna in RNAcontainer if len(rna.biorseoBayesPairA.predictions)],
    [ rna.biorseoBayesPairB.max_mcc  for rna in RNAcontainer if len(rna.biorseoBayesPairB.predictions)],
    [ rna.biorseoBayesPairC.max_mcc  for rna in RNAcontainer if len(rna.biorseoBayesPairC.predictions)],
    [ rna.biorseoBayesPairD.max_mcc  for rna in RNAcontainer if len(rna.biorseoBayesPairD.predictions)],
]

x_noPK_fully = [
    [ rna.rnasubopt.max_mcc for rna in RNAs_fully_predicted],
    [ rna.rnamoip.max_mcc for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUJAR3DA.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUJAR3DB.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUJAR3DC.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUJAR3DD.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUBayesPairA.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUBayesPairB.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUBayesPairC.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUBayesPairD.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoRawA.max_mcc for rna in RNAs_fully_predicted],
    [ rna.biorseoRawB.max_mcc for rna in RNAs_fully_predicted],
    [ rna.biorseoBayesPairA.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBayesPairB.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBayesPairC.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBayesPairD.max_mcc  for rna in RNAs_fully_predicted],
]  # We ensure having the same number of RNAs in every sample by discarding the one for which computations did not ended/succeeded.



print()
print("Without PK:")
print("%s RNAsubopt predictions" % is_all(len(x_noPK[0]), tot))
print("%s RNA MoIP predictions" % is_all(len(x_noPK[1]), tot))
print("%s biorseo + BGSU + JAR3D + f1A predictions" % is_all(len(x_noPK[2]), tot))
print("%s biorseo + BGSU + JAR3D + f1B predictions" % is_all(len(x_noPK[3]), tot))
print("%s biorseo + BGSU + JAR3D + f1C predictions" % is_all(len(x_noPK[4]), tot))
print("%s biorseo + BGSU + JAR3D + f1D predictions" % is_all(len(x_noPK[5]), tot))
print("%s biorseo + BGSU + BayesPairing + f1A predictions" % is_all(len(x_noPK[6]), tot))
print("%s biorseo + BGSU + BayesPairing + f1B predictions predictions" % is_all(len(x_noPK[7]), tot))
print("%s biorseo + BGSU + BayesPairing + f1C predictions predictions" % is_all(len(x_noPK[8]), tot))
print("%s biorseo + BGSU + BayesPairing + f1D predictions predictions" % is_all(len(x_noPK[9]), tot))
print("%s biorseo + Patternmatch + f1A predictions predictions" % is_all(len(x_noPK[10]), tot))
print("%s biorseo + Patternmatch + f1B predictions predictions" % is_all(len(x_noPK[11]), tot))
print("%s biorseo + BayesPairing + f1A predictions predictions" % is_all(len(x_noPK[12]), tot))
print("%s biorseo + BayesPairing + f1B predictions predictions" % is_all(len(x_noPK[13]), tot))
print("%s biorseo + BayesPairing + f1C predictions predictions" % is_all(len(x_noPK[14]), tot))
print("%s biorseo + BayesPairing + f1D predictions predictions" % is_all(len(x_noPK[15]), tot))
print("==> %s ARN were predicted with all methods successful." % is_all(len(x_noPK_fully[0]), tot) )

# stat tests
# Search if all methods are equal in positions with Friedman test:
test = stats.friedmanchisquare(*x_noPK_fully)
print("Friedman test without PK: H0 = 'The position parameter of all distributions is equal', p-value = ", test.pvalue)
# ==> No they are not, but none does better, no need to test one further.

# ================= Statistics (with pseudoknots) ========================

# load results in objects
for instance in RNAcontainer:
    instance.load_results_from(outputDir + "PK/")
    instance.evaluate()

x_PK = [
    [ rna.biokop.max_mcc for rna in RNAcontainer if len(rna.biokop.predictions)],
    [ rna.biokop.max_mcc for rna in RNAcontainer if len(rna.biokop.predictions)],
    [ rna.biorseoRawA.max_mcc for rna in RNAcontainer if len(rna.biorseoRawA.predictions)],
    [ rna.biorseoRawB.max_mcc for rna in RNAcontainer if len(rna.biorseoRawB.predictions)],
    [ rna.biorseoBayesPairA.max_mcc  for rna in RNAcontainer if len(rna.biorseoBayesPairA.predictions)],
    [ rna.biorseoBayesPairB.max_mcc  for rna in RNAcontainer if len(rna.biorseoBayesPairB.predictions)],
    [ rna.biorseoBayesPairC.max_mcc  for rna in RNAcontainer if len(rna.biorseoBayesPairC.predictions)],
    [ rna.biorseoBayesPairD.max_mcc  for rna in RNAcontainer if len(rna.biorseoBayesPairD.predictions)],
    [ rna.biorseoBGSUJAR3DA.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DA.predictions)],
    [ rna.biorseoBGSUJAR3DB.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DB.predictions)],
    [ rna.biorseoBGSUJAR3DC.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DC.predictions)],
    [ rna.biorseoBGSUJAR3DD.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DD.predictions)],
    [ rna.biorseoBGSUBayesPairA.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairA.predictions)],
    [ rna.biorseoBGSUBayesPairB.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairB.predictions)],
    [ rna.biorseoBGSUBayesPairC.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairC.predictions)],
    [ rna.biorseoBGSUBayesPairD.max_mcc  for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairD.predictions)]
]

RNAs_fully_predicted = [ x for x in RNAcontainer if x.has_complete_results(True)]

x_PK_fully = [
    [ rna.biokop.max_mcc for rna in RNAs_fully_predicted],
    [ rna.biokop.max_mcc for rna in RNAs_fully_predicted],
    [ rna.biorseoRawA.max_mcc for rna in RNAs_fully_predicted],
    [ rna.biorseoRawB.max_mcc for rna in RNAs_fully_predicted],
    [ rna.biorseoBayesPairA.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBayesPairB.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBayesPairC.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBayesPairD.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUJAR3DA.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUJAR3DB.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUJAR3DC.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUJAR3DD.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUBayesPairA.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUBayesPairB.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUBayesPairC.max_mcc  for rna in RNAs_fully_predicted],
    [ rna.biorseoBGSUBayesPairD.max_mcc  for rna in RNAs_fully_predicted],
]  # We ensure having the same number of RNAs in every sample by discarding the one for which computations did not ended/succeeded.

print()
print("With PK:")
print("%s Biokop predictions" % is_all(len(x_PK[1]), tot))
print("%s biorseo + Patternmatch + f1A predictions predictions" % is_all(len(x_PK[2]), tot))
print("%s biorseo + Patternmatch + f1B predictions predictions" % is_all(len(x_PK[3]), tot))
print("%s biorseo + BayesPairing + f1A predictions predictions" % is_all(len(x_PK[4]), tot))
print("%s biorseo + BayesPairing + f1B predictions predictions" % is_all(len(x_PK[5]), tot))
print("%s biorseo + BayesPairing + f1C predictions predictions" % is_all(len(x_PK[6]), tot))
print("%s biorseo + BayesPairing + f1D predictions predictions" % is_all(len(x_PK[7]), tot))
print("%s biorseo + BGSU + JAR3D + f1A predictions" % is_all(len(x_PK[8]), tot))
print("%s biorseo + BGSU + JAR3D + f1B predictions" % is_all(len(x_PK[9]), tot))
print("%s biorseo + BGSU + JAR3D + f1C predictions" % is_all(len(x_PK[10]), tot))
print("%s biorseo + BGSU + JAR3D + f1D predictions" % is_all(len(x_PK[11]), tot))
print("%s biorseo + BGSU + BayesPairing + f1A predictions" % is_all(len(x_PK[12]), tot))
print("%s biorseo + BGSU + BayesPairing + f1B predictions predictions" % is_all(len(x_PK[13]), tot))
print("%s biorseo + BGSU + BayesPairing + f1C predictions predictions" % is_all(len(x_PK[14]), tot))
print("%s biorseo + BGSU + BayesPairing + f1D predictions predictions" % is_all(len(x_PK[15]), tot))
print("==> %s ARN were predicted with all methods successful." % is_all(len(x_PK_fully[0]), tot) )


# stat tests
# First, search if all methods are equal in positions with Friedman test:
test = stats.friedmanchisquare(*x_PK_fully)
print("Friedman test with PK: H0 = 'The position parameter of all distributions is equal', p-value = ", test.pvalue)
# it looks like some methods do better. Let's test the difference:

test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[2])
print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and RawA are equal', p-value = ", test.pvalue)

test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[3])
print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and RawB are equal', p-value = ", test.pvalue)

test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[8])
print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dA are equal', p-value = ", test.pvalue)

test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[9])
print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dB are equal', p-value = ", test.pvalue)

test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[10])
print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dC are equal', p-value = ", test.pvalue)

test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[11])
print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dD are equal', p-value = ", test.pvalue)


# ================== Print results for application cases =====================

labels = ["Biokop","Biokop","RawA","RawB","BayesPairingA","BayesPairingB","BayesPairingC","BayesPairingD","JAR3DA","JAR3DB","JAR3DC","JAR3DD","BGSUBayesPairingA","BGSUBayesPairingB","BGSUBayesPairingC","BGSUBayesPairingD"]
print("RNAsubopt",":",x_noPK[0])
print("RNA-MOIP",":",x_noPK[1])
for data, name in zip(x_PK, labels):
    print(name,":",data)
labels = ["RNAsubopt","Biokop\t", "RNA MoIP\t","RawA\t","RawB\t","BayesPairingA","BayesPairingB","BayesPairingC","BayesPairingD","JAR3DA\t","JAR3DB\t","JAR3DC\t","JAR3DD\t","BGSUBPairingA","BGSUBPairingB","BGSUBPairingC","BGSUBPairingD"]
for r in RNAcontainer:
    print("\n",r.header_,"\nTrue structure:\t", r.true2d)
    for m, name in zip([r.rnasubopt, r.biokop, r.rnamoip,
                   r.biorseoRawA, 
                   r.biorseoRawB,
                   r.biorseoBayesPairA, 
                   r.biorseoBayesPairB, 
                   r.biorseoBayesPairC, 
                   r.biorseoBayesPairD,
                   r.biorseoBGSUJAR3DA, 
                   r.biorseoBGSUJAR3DB, 
                   r.biorseoBGSUJAR3DC, 
                   r.biorseoBGSUJAR3DD,
                   r.biorseoBGSUBayesPairA, 
                   r.biorseoBGSUBayesPairB, 
                   r.biorseoBGSUBayesPairC, 
                   r.biorseoBGSUBayesPairD ], labels):
        print(name+":\t",m.best_pred)

# # ================= PLOTS OF RESULTS =======================================

# merge = [   x_PK_fully[0], # Biokop
#             x_noPK_fully[0], # RNA subopt
#             x_noPK_fully[1], # RNA MoIP
#             x_noPK_fully[2], x_PK_fully[2], #biorseoRawA
#             x_noPK_fully[3], x_PK_fully[3], #biorseoRawB
#             x_noPK_fully[4], x_PK_fully[4], #biorseoBayesPairA
#             x_noPK_fully[5], x_PK_fully[5], #biorseoBayesPairB
#             x_noPK_fully[6], x_PK_fully[6], #biorseoBayesPairC
#             x_noPK_fully[7], x_PK_fully[7], #biorseoBayesPairD
#             x_noPK_fully[8], x_PK_fully[8], #biorseoBGSUJAR3DA
#             x_noPK_fully[9], x_PK_fully[9], #biorseoBGSUJAR3DB
#             x_noPK_fully[10], x_PK_fully[10], #biorseoBGSUJAR3DC
#             x_noPK_fully[11], x_PK_fully[11], #biorseoBGSUJAR3DD
#             x_noPK_fully[12], x_PK_fully[12], #biorseoBGSUBayesPairA
#             x_noPK_fully[13], x_PK_fully[13], #biorseoBGSUBayesPairB
#             x_noPK_fully[14], x_PK_fully[14], #biorseoBGSUBayesPairC
#             x_noPK_fully[15], x_PK_fully[15], #biorseoBGSUBayesPairD
# ]

# colors = [  'green', 'blue', 'goldenrod',
#             'darkturquoise', 'darkturquoise', 
#             'red', 'red', 
#             'firebrick', 'firebrick',
#             'limegreen', 'limegreen',
#             'olive', 'olive',
#             'forestgreen', 'forestgreen', 
#             'lime', 'lime',
#             'darkcyan', 'darkcyan', 
#             'royalblue', 'royalblue', 
#             'navy', 'navy', 
#             'limegreen', 'limegreen',
#             'olive', 'olive', 
#             'forestgreen', 'forestgreen', 
#             'lime', 'lime'
# ]
# labels = [  "Biokop", "RNAsubopt",
#             "RNA MoIP",
#             "$f_{1A}$",
#             "$f_{1B}$",
#             "$f_{1A}$",
#             "$f_{1B}$",
#             "$f_{1C}$",
#             "$f_{1D}$",
#             "$f_{1A}$",
#             "$f_{1B}$",
#             "$f_{1C}$",
#             "$f_{1D}$",
#             "$f_{1A}$",
#             "$f_{1B}$",
#             "$f_{1C}$",
#             "$f_{1D}$"        
# ]

# ax = plt.subplot(211)
# ax.tick_params(labelsize=12)
# for y in [ i/10 for i in range(11) ]:
#     plt.axhline(y=y, color="grey", linestyle="--", linewidth=1)
# colors = [  'blue','goldenrod',
#             'red', 'firebrick','limegreen','olive', 'forestgreen', 'lime',
#             'darkturquoise', 'darkcyan', 'royalblue', 'navy', 'limegreen','olive', 'forestgreen', 'lime'
#          ]
# bplot = plt.boxplot(x_noPK_fully, vert=True, patch_artist=True, notch=False, whis=[3,97])
# for patch, color in zip(bplot['boxes'], colors):
#     patch.set_facecolor(color)
# # plt.axhline(y=0, color="black", linewidth=1)
# # plt.axhline(y=1, color="black", linewidth=1)
# plt.xticks([1.0+i for i in range(16)], labels[1:])
# plt.ylim((0, 1.01))
# plt.ylabel("MCC", fontsize=12)
# plt.subplots_adjust(left=0.05, right=0.95)
# # plt.title("Performance without pseudoknots (%d RNAs included)" % len(x_noPK_fully[0]))


# ax = plt.subplot(212)
# ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, labelsize=12)
# ax.xaxis.set_label_position('top')
# for y in [ i/10 for i in range(11) ]:
#     plt.axhline(y=y, color="grey", linestyle="--", linewidth=1)
# colors = [  'green','green', 
#             'red', 'firebrick','limegreen','olive', 'forestgreen', 'lime',
#             'darkturquoise', 'darkcyan', 'royalblue', 'navy', 'limegreen','olive', 'forestgreen', 'lime'
#          ]
# labels = [  "Biokop"]
# bplot = plt.boxplot(x_PK_fully, vert=True, patch_artist=True, notch=False, whis=[3,97])
# for patch, color in zip(bplot['boxes'], colors):
#     patch.set_facecolor(color)
# # plt.axhline(y=0, color="black", linewidth=1)
# # plt.axhline(y=1, color="black", linewidth=1)
# plt.xticks([1.0+i for i in range(16)], labels)
# plt.ylim((0, 1.01))
# plt.ylabel("MCC", fontsize=12)
# plt.subplots_adjust(left=0.05, right=0.95)
# # plt.text(6.2,-0.3,"Performance with pseudoknots (%d RNAs included)" % len(x_PK_fully[0]), fontsize=12)


# plt.show()


# # ================== MCC performance ====================================
# # plt.subplot(141)
# x = [
#     [ rna.rnasubopt.max_mcc for rna in RNAs_fully_predicted],
#     # [ rna.rnamoip.max_mcc for rna in RNAs_fully_predicted],
#     [ rna.biorseoRawA.max_mcc for rna in RNAs_fully_predicted],
#     # [ rna.biorseoRawB.max_mcc for rna in RNAs_fully_predicted],
#     [ rna.biokop.max_mcc for rna in RNAs_fully_predicted],
#     [ rna.biorseoBGSUJAR3DA.max_mcc for rna in RNAs_fully_predicted]
# ]
# colors = ['xkcd:blue','xkcd:red', 'green', 'cyan']
# labels = ["Best RNAsubopt prediction",  "Best RawA prediction", "Best Biokop prediction", "Best JAR3DA prediction"]
# for y, col, lab in zip(x, colors, labels):
#     x_data = [ i for i in range(len(y)) if y[i]]
#     y_data = [ i for i in y if i]
#     plt.scatter(x_data, y_data, color=col, label=lab, marker='o', s=2.5)
# plt.axhline(y=0, color='black', linewidth=1)
# plt.axvline(x=0, color='black', linewidth=1)
# plt.xlabel("RNA Strand verified tRNA structures (10 < |nt| < 100)")
# plt.ylabel("Mattews Correlation Coefficient")
# plt.title("Performance of the prediction method")
# plt.legend(loc="lower right")
# plt.show()



# # ========================= Number of solutions ===========================

# plt.subplot(231)
# x = [
#     [ rna.biorseoBayesPairA.n_pred for rna in RNAs_fully_predicted],
#     [ rna.biorseoBayesPairB.n_pred for rna in RNAs_fully_predicted],
#     [ rna.biorseoBayesPairC.n_pred for rna in RNAs_fully_predicted],
#     [ rna.biorseoBayesPairD.n_pred for rna in RNAs_fully_predicted]
# ]
# colors = ['red', 'black', 'blue', 'limegreen']
# labels = ["$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$"]
# # plt.hist(x, max([ max(x[i]) for i in range(len(x))])0, color=colors, align="mid", density=False, fill=False, histtype="step", stacked=False, label=labels)
# plt.hist(x, max([ max(x[i]) for i in range(len(x))]),  align="mid", density=False, stacked=False, label=labels)
# for i in range(0, max([ max(x[i]) for i in range(len(x))]), 10):
#     plt.axvline(x=i, linestyle='--', color='gray')
# plt.xlabel("Size of Pareto set")
# plt.xticks([i for i in range(max([ max(x[i]) for i in range(len(x))]))])
# plt.ylabel("Number of RNAs")
# plt.ylim((0,265))
# plt.legend(loc="upper right")
# plt.title("(A) Rna3Dmotifs + BayesPairing")

# plt.subplot(232)
# x = [
#     [ rna.biorseoBGSUBayesPairA.n_pred for rna in RNAs_fully_predicted],
#     [ rna.biorseoBGSUBayesPairB.n_pred for rna in RNAs_fully_predicted],
#     [ rna.biorseoBGSUBayesPairC.n_pred for rna in RNAs_fully_predicted],
#     [ rna.biorseoBGSUBayesPairD.n_pred for rna in RNAs_fully_predicted]
# ]
# plt.hist(x, max([ max(x[i]) for i in range(len(x))]),  align="mid", density=False, stacked=False, label=labels)
# for i in range(0, max([ max(x[i]) for i in range(len(x))]), 10):
#     plt.axvline(x=i, linestyle='--', color='gray')
# plt.xticks([i for i in range(max([ max(x[i]) for i in range(len(x))]))])
# plt.xlabel("Size of Pareto set")
# # plt.ylabel("Number of RNAs")
# plt.ylim((0,265))
# plt.legend(loc="upper right")
# plt.title("(B) The RNA Motif Atlas 3.2 + BayesPairing")


# plt.subplot(233)
# x = [
#     [ rna.biorseoRawA.n_pred for rna in RNAs_fully_predicted],
#     [ rna.biorseoRawB.n_pred for rna in RNAs_fully_predicted if rna.biorseoRawB.n_pred < 55],
# ]
# # colors = ['red', 'firebrick']
# colors = ['red', 'black']
# labels = ["$f_{1A}$", "$f_{1B}$"]
# plt.hist(x, 55,  align="mid", density=False, stacked=False, label=labels)
# for i in range(0, 55, 10):
#     plt.axvline(x=i, linestyle='--', color='gray')
# plt.xticks([i for i in range(55)])
# plt.ylim((0,265))
# plt.xlabel("Size of Pareto set")
# # plt.ylabel("Number of RNAs")
# plt.legend(loc="upper right")
# plt.title("(C) Rna3Dmotifs + Simple pattern matching")


# plt.subplot(234)
# x = [
#     [ rna.biorseoBGSUJAR3DA.n_pred for rna in RNAs_fully_predicted],
#     [ rna.biorseoBGSUJAR3DB.n_pred for rna in RNAs_fully_predicted],
#     [ rna.biorseoBGSUJAR3DC.n_pred for rna in RNAs_fully_predicted],
#     [ rna.biorseoBGSUJAR3DD.n_pred for rna in RNAs_fully_predicted]
# ]
# # colors = ['darkturquoise', 'darkcyan', 'royalblue', 'navy']
# colors = ['red', 'black', 'blue', 'limegreen']
# labels = ["$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$"]
# plt.hist(x, max([ max(x[i]) for i in range(len(x))]),  align="mid", density=False, stacked=False, label=labels)
# for i in range(0, max([ max(x[i]) for i in range(len(x))]), 10):
#     plt.axvline(x=i, linestyle='--', color='gray')
# plt.xticks([i for i in range(max([ max(x[i]) for i in range(len(x))]))])
# plt.xlabel("Size of Pareto set")
# # plt.ylabel("Number of RNAs")
# plt.legend(loc="upper right")
# plt.title("(D) The RNA Motif Atlas 3.2 + JAR3D")


# plt.subplot(235)
# x = [
#     [ rna.rnasubopt.n_pred for rna in RNAs_fully_predicted],
#     [ rna.rnamoip.n_pred for rna in RNAs_fully_predicted],
# ]
# colors = ['blue', 'goldenrod']
# labels = ["RNAsubopt", "RNA-MoIP"]
# plt.hist(x, max([ max(x[i]) for i in range(len(x))]), color=colors, align="mid", density=False, stacked=False, label=labels)
# for i in range(0, max([ max(x[i]) for i in range(len(x))]), 10):
#     plt.axvline(x=i, linestyle='--', color='gray')
# plt.xticks([i for i in range(max([ max(x[i]) for i in range(len(x))]))])
# plt.xlabel("Size of results set")
# plt.ylim((0,265))
# # plt.ylabel("Number of RNAs")
# plt.legend(loc="upper right")
# plt.title("(E) Other methods")

# plt.subplot(236)
# x = [
#     [ rna.biokop.n_pred for rna in RNAs_fully_predicted],
# ]
# colors = ['green']
# labels = [ "Biokop"]
# plt.hist(x, max([ max(x[i]) for i in range(len(x))]), color=colors, align="mid", density=False, stacked=False, label=labels)
# for i in range(0, max([ max(x[i]) for i in range(len(x))]), 10):
#     plt.axvline(x=i, linestyle='--', color='gray')
# plt.xticks([i for i in range(max([ max(x[i]) for i in range(len(x))]))])
# plt.xlabel("Size of Pareto set")
# plt.ylim((0,265))
# # plt.ylabel("Number of RNAs")
# plt.legend(loc="upper right")
# plt.title("(F) Biokop")

# plt.show()


# # MCC boost compared to RNA subopt
# plt.subplot(143)
# x = [
#     [ rna.rnamoip.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.rnamoip.predictions)],
#     [ rna.biorseoRawA.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoRawA.predictions)],
#     [ rna.biorseoRawB.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoRawB.predictions)],
#     [ rna.biokop.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biokop.predictions)],
#]
# colors = ['xkcd:goldenrod', 'xkcd:red', 'firebrick', 'limegreen']
# labels = ["$\Delta$MCC(RNAsubopt,RNA MoIP)","$\Delta$MCC(RNAsubopt,RNA MoBOIP)",
#     "$\Delta$MCC(RNAsubopt,RNA MoBOIP++)","$\Delta$MCC(RNAsubopt,Biokop)"]
# bplot = plt.boxplot(x, vert=False, patch_artist=True, notch=False, whis=[3,97])
# for patch, color in zip(bplot['boxes'], colors):
#     patch.set_facecolor(color)
# plt.axvline(x=0, color="black", linewidth=1)
# plt.yticks([1.0+i for i in range(4)], labels)
# plt.xlim((-1.1, 1.1))
# plt.xlabel("Improvement in MCC")
# plt.title("MCC performance relatively to RNAsubopt")
# plt.show()


# plt.subplot(222)
# x = [
#     [ rna.biorseoBGSUBayesPairA.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairA.predictions)],
#     [ rna.biorseoBGSUBayesPairB.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairB.predictions)],
#     [ rna.biorseoBGSUBayesPairC.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairC.predictions)],
#     [ rna.biorseoBGSUBayesPairD.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairD.predictions)],
#]
# bplot = plt.boxplot(x, vert=False, patch_artist=True, notch=False, whis=[3,97])
# for patch, color in zip(bplot['boxes'], colors):
#     patch.set_facecolor(color)
# plt.axvline(x=0, color="black", linewidth=1)
# plt.yticks([1.0+i for i in range(4)], labels)
# plt.xlim((-1.1, 1.1))
# # plt.xlabel("Improvement in MCC")
# plt.title("(B) The RNA Motif Atlas 3.2 + BayesPairing")


# plt.subplot(223)
# x = [
#     [ rna.biorseoRawA.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoRawA.predictions)],
#     [ rna.biorseoRawB.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoRawB.predictions)],
#]
# colors = ['red', 'firebrick']
# labels = ["$f_{1A}$", "$f_{1B}$"]
# bplot = plt.boxplot(x, vert=False, patch_artist=True, notch=False, whis=[3,97])
# for patch, color in zip(bplot['boxes'], colors):
#     patch.set_facecolor(color)
# plt.axvline(x=0, color="black", linewidth=1)
# plt.yticks([1.0+i for i in range(2)], labels)
# plt.xlabel("Improvement in MCC")
# plt.xlim((-1.1, 1.1))
# plt.title("(C) Rna3Dmotifs + Simple pattern matching")


# plt.subplot(224)
# x = [
#     [ rna.biorseoBGSUJAR3DA.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DA.predictions)],
#     [ rna.biorseoBGSUJAR3DB.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DB.predictions)],
#     [ rna.biorseoBGSUJAR3DC.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DC.predictions)],
#     [ rna.biorseoBGSUJAR3DD.max_mcc - rna.rnasubopt.max_mcc for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DD.predictions)],
#]
# colors = ['darkturquoise', 'darkcyan', 'royalblue', 'navy']
# labels = ["$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$"]
# bplot = plt.boxplot(x, vert=False, patch_artist=True, notch=False, whis=[3,97])
# for patch, color in zip(bplot['boxes'], colors):
#     patch.set_facecolor(color)
# plt.axvline(x=0, color="black", linewidth=1)
# plt.yticks([1.0+i for i in range(4)], labels)
# plt.xlabel("Improvement in MCC")
# plt.xlim((-1.1, 1.1))
# plt.title("(D) The RNA Motif Atlas 3.2 + JAR3D")
# plt.show()


# # insertion ratio of the best structure
# plt.subplot(221)
# x = [
#     [ rna.biorseoBayesPairA.ratio for rna in RNAcontainer if len(rna.biorseoBayesPairA.predictions)],
#     [ rna.biorseoBayesPairC.ratio for rna in RNAcontainer if len(rna.biorseoBayesPairC.predictions)],
#     [ rna.biorseoBayesPairD.ratio for rna in RNAcontainer if len(rna.biorseoBayesPairD.predictions)],
#     [ rna.biorseoBayesPairB.ratio for rna in RNAcontainer if len(rna.biorseoBayesPairB.predictions)]
#]
# colors = ['olive', 'forestgreen', 'lime', 'limegreen']
# labels = ["$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$"]
# plt.hist(x, 30, color=colors, align="mid", density=True, fill=False, histtype="step", stacked=False, label=labels)
# plt.xlim(0, 1)
# # plt.xlabel("Ratio $Ninserted_{best} / Ninserted_{max}$")
# plt.ylabel("Percentage of RNAs")
# plt.yticks([])
# plt.title("(A) Rna3Dmotifs + BayesPairing")
# plt.legend(loc="upper left")

# plt.subplot(222)
# x = [
#     [ rna.biorseoBGSUBayesPairA.ratio for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairA.predictions)],
#     [ rna.biorseoBGSUBayesPairC.ratio for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairC.predictions)],
#     [ rna.biorseoBGSUBayesPairD.ratio for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairD.predictions)],
#     [ rna.biorseoBGSUBayesPairB.ratio for rna in RNAcontainer if len(rna.biorseoBGSUBayesPairB.predictions)]
#]
# plt.hist(x, 30, color=colors, align="mid", density=True, fill=False, histtype="step", stacked=False, label=labels)
# plt.xlim(0, 1)
# # plt.xlabel("Ratio $Ninserted_{best} / Ninserted_{max}$")
# plt.ylabel("Percentage of RNAs")
# plt.yticks([])
# plt.title("(B) The RNA Motif Atlas 3.2 + BayesPairing")
# plt.legend(loc="upper left")

# plt.subplot(223)
# x = [
#     [ rna.biorseoRawA.ratio for rna in RNAcontainer if len(rna.biorseoRawA.predictions)],
#     [ rna.biorseoRawB.ratio for rna in RNAcontainer if len(rna.biorseoRawB.predictions)],
#]
# colors = ['red', 'firebrick']
# labels = ["$f_{1A}$", "$f_{1B}$"]
# plt.hist(x, 30, color=colors, align="mid", density=True, fill=False, histtype="step", stacked=False, label=labels)
# plt.xlim(0, 1)
# plt.xlabel("Ratio $Ninserted_{best} / Ninserted_{max}$")
# plt.ylabel("Percentage of RNAs")
# plt.yticks([])
# plt.title("(C) Rna3Dmotifs + Simple pattern matching")
# plt.legend(loc="upper left")

# plt.subplot(224)
# x = [
#     [ rna.biorseoBGSUJAR3DA.ratio for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DA.predictions)],
#     [ rna.biorseoBGSUJAR3DC.ratio for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DC.predictions)],
#     [ rna.biorseoBGSUJAR3DD.ratio for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DD.predictions)],
#     [ rna.biorseoBGSUJAR3DB.ratio for rna in RNAcontainer if len(rna.biorseoBGSUJAR3DB.predictions)]
#]
# colors = ['darkturquoise', 'darkcyan', 'royalblue', 'navy']
# labels = ["$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$"]
# plt.hist(x, 30, color=colors, align="mid", density=True, fill=False, histtype="step", stacked=False, label=labels)
# plt.xlim(0, 1)
# plt.xlabel("Ratio $Ninserted_{best} / Ninserted_{max}$")
# plt.ylabel("Percentage of RNAs")
# plt.yticks([])
# plt.title("(D) The RNA Motif Atlas 3.2 + JAR3D")
# plt.legend(loc="upper left")

# plt.show()