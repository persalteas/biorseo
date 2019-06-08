#!/usr/bin/python3
# coding=utf-8
import sys
import getopt
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
bypdir = ""
biorseoDir = "."
exec(compile(open(biorseoDir+"/EditMe").read(), '', 'exec'))
runDir = path.dirname(path.realpath(__file__))
self.outputf = biorseoDir + "/results/"
tempDir = biorseoDir + "/temp/"
HLmotifDir = biorseoDir + "/data/modules/BGSU/HL/3.2/lib"
ILmotifDir = biorseoDir + "/data/modules/BGSU/IL/3.2/lib"
descfolder = biorseoDir + "/data/modules/DESC"


# ================== CLASSES AND FUNCTIONS ================================

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


class RNA:
    def __init__(self, header, seq):
        self.seq_ = seq
        self.header_ = header
        self.length = len(seq)

        self.rnasubopt = []
        self.biorseoRawA = []
        self.biorseoRawB = []
        self.biorseoBGSUJAR3DA = []
        self.biorseoBGSUJAR3DC = []
        self.biorseoBGSUJAR3DD = []
        self.biorseoBGSUJAR3DB = []
        self.biorseoBayesPairA = []
        self.biorseoBayesPairC = []
        self.biorseoBayesPairD = []
        self.biorseoBayesPairB = []
        self.biorseoBGSUBayesPairA = []
        self.biorseoBGSUBayesPairC = []
        self.biorseoBGSUBayesPairD = []
        self.biorseoBGSUBayesPairB = []

    def get_RNAsubopt_results(self):
        rna = open(self.outputf + self.basename + ".subopt", "r")
        lines = rna.readlines()
        rna.close()
        for i in range(2, len(lines)):
            ss = lines[i].split(' ')[0]
            if ss not in self.rnasubopt.predictions:
                self.rnasubopt.predictions.append(ss)

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
        # else:
        #     print(targetdir+ self.basename + ".bgsubypA not found !")
    
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
        # else:
        #     print(targetdir+ self.basename + ".bgsubypB not found !")

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
        # else:
        #     print(targetdir+ self.basename + ".bgsubypC not found !")

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
        # else:
        #     print(targetdir+ self.basename + ".bgsubypD not found !")


class BiorseoInstance:
    def __init__(self, argv):
        # set default options
        self.type = "dpm"
        self.modules = "desc"
        self.func = 'B'
        self.outputf = self.outputf
        self.jobcount = 0

        # Parse options
        try:
            opts, args = getopt.getopt(
                argv, "hi:o:", ["type=", "func=", "modules=", "rna3dmotifs", "3dmotifatlas", "jar3d", "bayespairing", "patternmatching"])
        except getopt.GetoptError:
            print("Please provide arguments !")
            sys.exit(2)
        for opt, arg in opts:
            if opt == "-h":
                print("biorseo.py -i myRNA.fa -o myRNA.jar3dB --type jar3d --func B")
                sys.exit()
            elif opt == "-i":
                self.inputfile = arg
                self.mode = 0  # single sequence mode
            elif opt == "-o":
                self.outputf = arg  # output file or folder...
            elif opt == "--func":
                if arg in ['A', 'B', 'C', 'D']:
                    self.func = arg
                else:
                    raise "Unknown scoring function " + arg
            elif opt == "--type":
                if arg in ['dpm', 'jar3d', 'byp']:
                    self.type = arg
                else:
                    raise "Unknown pattern matching method " + arg
            elif opt == "--modules":
                if arg in ['desc', 'bgsu']:
                    self.modules = arg
                else:
                    raise "Unsupported module model type " + arg
            else:
                raise "Unknown option " + opt

        # create jobs
        self.list_jobs()

        if self.mode:
            # Create a job manager
            self.manager = Manager()
            self.running_stats = self.manager.list()
            self.running_stats.append(0)  # n_launched
            self.running_stats.append(0)  # n_finished
            self.running_stats.append(0)  # n_skipped
            self.fails = self.manager.list()

            # Create the output folder
            subprocess.call(["mkdir", "-p", self.outputf])

    def enumerate_loops(self, s):
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

    def launch_JAR3D_worker(self, loop):
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
            cmd = ["java", "-jar", jar3dexec, filename, HLmotifDir+"/all.txt",
                   loop.header[1:]+".HLloop.csv", loop.header[1:]+".HLseq.csv"]
        else:
            cmd = ["java", "-jar", jar3dexec, filename, ILmotifDir+"/all.txt",
                   loop.header[1:]+".ILloop.csv", loop.header[1:]+".ILseq.csv"]
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

    def launch_JAR3D(self, seq_, basename):
        rnasubopt_preds = []
        # Extracting probable loops from RNA-subopt structures
        rna = open(self.outputf + basename + ".subopt", "r")
        lines = rna.readlines()
        rna.close()
        for i in range(2, len(lines)):
            ss = lines[i].split(' ')[0]
            if ss not in rnasubopt_preds:
                rnasubopt_preds.append(ss)
        HLs = []
        ILs = []
        for ss in rnasubopt_preds:
            loop_candidates = self.enumerate_loops(ss)
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
        insertion_sites = [x for y in pool.map(
            self.launch_JAR3D_worker, loops) for x in y]
        insertion_sites.sort(reverse=True)
        # Writing results to CSV file
        c = 0
        resultsfile = open(self.outputf+basename+".sites.csv", "w")
        resultsfile.write("Motif,Rotation,Score,Start1,End1,Start2,End2\n")
        for site in insertion_sites:
            if site.score > 10:
                c += 1
                string = "FOUND with score %d:\t\t possible insertion of motif " % site.score + site.atlas_id
                if site.rotation:
                    string += " (reversed)"
                string += (" on " + site.loop.get_header() + " at positions")
            resultsfile.write(site.atlas_id+',' +
                              str(bool(site.rotation))+",%d" % site.score+',')
            positions = [','.join([str(y) for y in x]) for x in site.position]
            if len(positions) == 1:
                positions.append("-,-")
            resultsfile.write(','.join(positions)+'\n')
        resultsfile.close()

    def launch_BayesPairing(self, module_type, seq_, header_, basename):
        chdir(bypdir)

        cmd = ["python3", "parse_sequences.py", "-seq", self.outputf +
               basename + ".fa", "-d", module_type, "-interm", "1"]

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
        insertion_sites = [x for x in ast.literal_eval(l.split(":")[1][1:])]
        if module_type == "rna3dmotif":
            rna = open(self.outputf + basename + ".byp.csv", "w")
        else:
            rna = open(self.outputf + basename + ".bgsubyp.csv", "w")
        rna.write("Motif,Score,Start1,End1,Start2,End2...\n")
        for i, module in enumerate(insertion_sites):
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
                    for (p, q) in zip(*[iter(pos[1:])]*2):
                        if q > p:
                            rna.write(','+str(p)+','+str(q))
                    rna.write('\n')
        rna.close()

    def execute_job(self, j):
        if j.checkFunc_ is not None:
            if j.checkFunc_(*j.checkArgs_):
                self.running_stats[2] += 1
                print("["+str(self.running_stats[0]+self.running_stats[2]) +
                      '/'+str(self.jobcount)+"]\tSkipping a finished job")
                return 0
        self.running_stats[0] += 1
        if len(j.cmd_):
            logfile = open("log_of_the_run.sh", 'a')
            logfile.write(" ".join(j.cmd_))
            logfile.write("\n")
            logfile.close()
            print("["+str(self.running_stats[0]+self.running_stats[2]) +
                  '/'+str(self.jobcount)+"]\t"+" ".join(j.cmd_))
            r = subprocess.call(j.cmd_, timeout=j.timeout_)
        elif j.func_ is not None:
            print("["+str(self.running_stats[0]+self.running_stats[2])+'/'+str(self.jobcount) +
                  "]\t"+j.func_.__name__+'('+", ".join([a for a in j.args_])+')')
            try:
                r = j.func_(*j.args_)
            except:
                r = 1
                pass
        if r:
            self.fails.append(j)
        self.running_stats[1] += 1
        return r

    def check_result_existence(self, datatype, method, function, with_PK, basename):
        folder = self.outputf+"PK/" if with_PK else self.outputf+"noPK/"
        if datatype == "bgsu":
            if method == "jar3d":
                extension = ".jar3d"
            elif method == "byp":
                extension = ".bgsubyp"
            else:
                raise "Unknown method !"
        elif datatype == "desc":
            if method == "dpm":
                extension = ".raw"
            elif method == "byp":
                extension = ".byp"
            else:
                raise "Unknown method !"
        else:
            raise "Unknown data type !"
        return path.isfile(folder + basename + extension + function)

    def check_csv_existence(self, datatype, method, basename):
        if datatype == "bgsu":
            if method == "jar3d":
                extension = ".sites.csv"
            elif method == "byp":
                extension = ".bgsubyp.csv"
            else:
                raise "Unknown method !"
        elif datatype == "desc":
            if method == "byp":
                extension = ".byp.csv"
            else:
                raise "You cannot use " + method + " with " + datatype + " data !"
        else:
            raise "Unknown data type !"
        return path.isfile(self.outputf + basename + extension)

    def list_jobs(self):

        # Read fasta file, which can contain one or several RNAs
        RNAcontainer = []
        print("loading file(s)...")

        db = open(self.inputfile, "r")
        c = 0
        header = ""
        seq = ""
        while True:
            l = db.readline()
            if l == "":
                break
            c += 1
            c = c % 2
            if c == 1:
                if header != "": # This is our second RNA in the fasta file
                    self.mode = 1
                header = l[:-1]
            if c == 0:
                seq = l[:-1].upper()
                if is_canonical_nts(seq):
                    header = header.replace('/', '_')
                    RNAcontainer.append(RNA(header, seq))
                    if not path.isfile(self.outputf + header + ".fa"):
                        rna = open(self.outputf + header + ".fa", "w")
                        rna.write(">" + header +'\n')
                        rna.write(seq +'\n')
                        rna.close()
                    db.close()

        for nt, number in ignored_nt_dict.items():
            print("ignored %d sequences because of char %c" % (number, nt))
        tot = len(RNAcontainer)
        print("Loaded %d RNAs." % (tot))

        #define job list
        joblist = []
        for instance in RNAcontainer:
            
            executable = biorseoDir + "/bin/biorseo"
            fastafile = self.outputf+instance.header+".fa"
            method_type = ""
            ext = ".raw"

            if self.type == "jar3d":
                ext = ".jar3d"
                method_type = "--jar3dcsv"
                csv = self.outputf + instance.header + ".sites.csv"

                # RNAsubopt
                joblist.append(Job(command=["RNAsubopt", "-i", fastafile, "--outfile="+ instance.header + ".subopt"], priority=1, checkFunc=check_RNAsubopt, checkArgs=[instance.header]))
                joblist.append(Job(command=["mv", instance.header + ".subopt", self.outputf], priority=2, checkFunc=check_RNAsubopt, checkArgs=[instance.header]))
                # JAR3D
                joblist.append(Job(function=self.launch_JAR3D, args=[instance.seq_, instance.header], priority=3, how_many_in_parallel=1, checkFunc=check_JAR3D, checkArgs=[instance.header]))

            if self.type == "byp":
                method_type = "--bayespaircsv"
                if self.modules == "desc":
                    ext = ".byp"
                    csv = self.outputf + instance.header + ".byp.csv"
                    joblist.append(Job(function=self.launch_BayesPairing, args=["rna3dmotif", instance.seq_, instance.header_, instance.header], how_many_in_parallel=-1, priority=1, checkFunc=check_BayesPairing, checkArgs=[instance.header]))
                elif self.modules == "bgsu":
                    ext = ".bgsubyp"
                    csv = self.outputf + instance.header + ".bgsubyp.csv"
                    joblist.append(Job(function=self.launch_BayesPairing, args=["3dmotifatlas", instance.seq_, instance.header_, instance.header], how_many_in_parallel=-1, priority=1, checkFunc=check_BGSUBayesPairing, checkArgs=[instance.header]))

            command = [executable, "-s", fastafile ]
            if method_type:
                command += [ method_type, csv ]
            command += [ "-o", self.outputf + instance.header + ext + self.func, "--type", self.func ]
            joblist.append(Job(command=command, priority=4, timeout=3600, how_many_in_parallel=3))
            

if __name__ == "__main__":
    BiorseoInstance(sys.argv)
