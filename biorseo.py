#!/usr/bin/python3
# coding=utf-8
import sys
import getopt
from scipy import stats
import subprocess
from os import path, makedirs, getcwd, chdir, devnull, remove, walk
import matplotlib.pyplot as plt
from matplotlib import colors
from math import sqrt
from multiprocessing import cpu_count, Manager
import multiprocessing
import ast
from shutil import move


# ================== DEFINITION OF THE PATHS ==============================

# Parse options
try:
    opts, args = getopt.getopt( sys.argv[1:], 
                                "bc:f:hi:jl:no:O:pt:v", 
                             [  "verbose","rna3dmotifs","3dmotifatlas","jar3d","bayespairing","patternmatch","func=",
                                "help","version","seq=","modules-path=", "jar3dexec=", "bypdir=", "biorseodir=", "first-objective=","output=","theta=",
                                "interrupt-limit=", "outputf="])
except getopt.GetoptError as err:
    print(err)
    sys.exit(2)

m = Manager()
running_stats = m.list()
running_stats.append(0) # n_launched
running_stats.append(0) # n_finished
running_stats.append(0) # n_skipped

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


class MyPool(multiprocessing.pool.Pool):
    # We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
    # because the latter is only a wrapper function, not a proper class.
    def __init__(self, *args, **kwargs):
        kwargs['context'] = NoDaemonContext()
        super(MyPool, self).__init__(*args, **kwargs)


class Loop:
    def __init__(self, header, subsequence, looptype, position):
        self.header = header
        self.seq = subsequence
        self.type = looptype
        self.position = position


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
    def __init__(self, command=[], function=None, args=[], how_many_in_parallel=0, priority=1, timeout=None):
        self.cmd_ = command
        self.func_ = function
        self.args_ = args
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
        self.header = header
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


class BiorseoInstance:
    def __init__(self, opts):
        # set default options
        self.type = "dpm"
        self.modules = "desc"
        self.func = 'B'
        self.inputfile = ""
        self.finalname = ""
        self.outputf = ""
        if path.exists("/biorseo/results"): # docker image default
            self.outputf = "/biorseo/results"
        self.output = ""
        self.jobcount = 0
        self.joblist = []
        self.mode = 0 # default is single sequence mode
        self.forward_options = []
        self.jar3dexec = "/jar3d_2014-12-11.jar"
        self.bypdir = "/byp/src"
        self.biorseoDir = "/biorseo"
        self.HLmotifDir = "/modules/BGSU/HL/3.2/lib"
        self.ILmotifDir = "/modules/BGSU/IL/3.2/lib"
        self.descfolder = "/modules/DESC"
        self.runDir = path.dirname(path.realpath(__file__))
        self.tempDir = "temp/"

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(  "Biorseo, Bi-Objective RNA Structure Efficient Optimizer\n"
                        "Bio-objective integer linear programming framework to predict RNA secondary structures by including known RNA modules.\n"
                        "Developped by Louis Becquey (louis.becquey@univ-evry.fr), 2019\n\n")
                print("Usage:\tYou must provide:\n\t1) a FASTA input file with -i,\n\t2) a module type with --rna3dmotifs or --3dmotifatlas"
                      "\n\t3) one module placement method in { --patternmatch, --jar3d, --bayespairing }\n\t4) one scoring function with --func A, B, C or D"
                      "\n\n\tIf you are not using the Docker image: \n\t5) --modules-path, --biorseodir and (--jar3dexec or --bypdir)")
                print()
                print("Options:")
                print("-h [ --help ]\t\t\tPrint this help message")
                print("--version\t\t\tPrint the program version")
                print("-i [ --seq=… ]\t\t\tFASTA file with the query RNA sequence")
                print("-p [ --patternmatch ]\t\tUse regular expressions to place modules in the sequence")
                print("-j [ --jar3d ]\t\t\tUse JAR3D to place modules in the sequence (requires --3dmotifatlas)")
                print("-b [ --bayespairing ]\t\tUse BayesPairing to place modules in the sequence")
                print("-o [ --output=… ]\t\tFile to summarize the results")
                print("-O [ --outputf=… ]\t\tFolder where to output result and temp files")
                print("-f [ --func=… ]\t\t\t(A, B, C or D, default is B)"
                      " Objective function to score module insertions:\n\t\t\t\t  (A) insert big modules (B) insert light, high-order modules"
                      "\n\t\t\t\t  (c) insert modules which score well with the sequence\n\t\t\t\t  (D) insert light, high-order modules which score well with the sequence."
                      "\n\t\t\t\t  C and D require cannot be used with --patternmatch.")
                print("-c [ --first-objective=… ]\t(default 1) Objective to solve in the mono-objective portions of the algorithm."
                      "\n\t\t\t\t  (1) is the module objective given by --func, (2) is the expected accuracy of the structure.")
                print("-l [ --limit=… ]\t\t(default 500) Intermediate number of solutions in the Pareto set from which"
                      " we give up the computation.")
                print("-t [ --theta=… ]\t\t(default 0.001) Pairing-probability threshold to consider or not the possibility of pairing")
                print("-n [ --disable-pseudoknots ]\tAdd constraints to explicitly forbid the formation of pseudoknots")
                print("-v [ --verbose ]\t\tPrint what is happening to stdout")
                print("--modules-path=…\t\tPath to the modules data.\n\t\t\t\t  The folder should contain modules in the DESC format as output by Djelloul & Denise's"
                      "\n\t\t\t\t  'catalog' program for use with --rna3dmotifs, or should contain the IL/ and HL/ folders from a release of\n\t\t\t\t  the RNA 3D Motif Atlas"
                      "for use with --3dmotifatlas.\n\t\t\t\t  Consider placing these files on a fast I/O device (NVMe SSD, ...)")
                print("--jar3dexec=…\t\t\tPath to the jar3d executable.\n\t\t\t\t  Default is /jar3d_2014-12-11.jar, you should use this option if you run"
                      "\n\t\t\t\t  BiORSEO from outside the docker image.")
                print("--bypdir=…\t\t\tPath to the BayesParing src folder.\n\t\t\t\t  Default is /byp/src, you should use this option if you run"
                      "\n\t\t\t\t  BiORSEO from outside the docker image.")
                print("--biorseodir=…\t\t\tPath to the BiORSEO root directory.\n\t\t\t\t  Default is /biorseo, you should use this option if you run"
                      "\n\t\t\t\t  BiORSEO from outside the docker image. Use the FULL path.")
                print("\nExamples:")
                print("biorseo.py -i myRNA.fa -O myResultsFolder/ --rna3dmotifs --patternmatch --func B")
                print("biorseo.py -i myRNA.fa -O myResultsFolder/ --3dmotifatlas --jar3d --func B -l 800")
                print("biorseo.py -i myRNA.fa -v --3dmotifatlas --bayespairing --func D")
                sys.exit()
            elif opt == "-i" or opt == "--seq":
                self.inputfile = arg
            elif opt == "-O" or opt == "--outputf":
                self.outputf = arg # output folder
                if self.outputf[0] != '/':
                    self.outputf = getcwd() + '/' + self.outputf
                if self.outputf[-1] != '/':
                    self.outputf = self.outputf + '/'
            elif opt == "-o" or opt == "--output":
                self.output = arg # output file 
                if self.output[0] != '/':
                    self.output = getcwd() + self.output
            elif opt == "-f" or opt == "--func":
                if arg in ['A', 'B', 'C', 'D']:
                    self.func = arg
                else:
                    raise "Unknown scoring function " + arg
            elif opt == "-p" or opt == "--patternmatch":
                self.type = "dpm"
            elif opt == "-j" or opt == "--jar3d":
                self.type = "jar3d"
            elif opt == "-b" or opt == "--bayespairing":
                self.type = "byp"
            elif opt == "--rna3dmotifs":
                self.modules = "desc"
            elif opt == "--3dmotifatlas":
                self.modules = "bgsu"
            elif opt == "--modules-path":
                self.HLmotifDir = arg + "/HL/3.2/lib"
                self.ILmotifDir = arg + "/IL/3.2/lib"
                self.descfolder = arg
                print("Looking for modules in", arg)
            elif opt == "--jar3dexec":
                self.jar3dexec = arg
                print("Using ", arg)
            elif opt == "--bypdir":
                self.bypdir = arg
                print("Using trained BayesPairing in", arg)
            elif opt == "--biorseodir":
                self.biorseoDir = arg
            elif opt == "--version":
                subprocess.call([self.biorseoDir+"/bin/biorseo", "--version"])
                exit(0)
            elif opt == "-l" or opt == "--interrupt-limit":
                self.forward_options.append("-l")
                self.forward_options.append(arg)
            elif opt == "-v" or opt == "--verbose":
                self.forward_options.append("-v")
            elif opt == "-n" or opt == "--disable-pseudoknots":
                self.forward_options.append("-n")
            elif opt == "-t" or opt == "--theta":
                self.forward_options.append("-t")
                self.forward_options.append(arg)
            elif opt == "-c" or opt == "--first-objective":
                self.forward_options.append("-c")
                self.forward_options.append(arg)

        if self.outputf != "":
            print("saving files to", self.outputf)
    
        # create jobs
        self.list_jobs()

        # run them
        self.execute_jobs()         

        # locate the results at the right place
        if self.output != "":
            subprocess.call(["mv", self.tempDir+self.finalname.split('/')[-1], self.output])
        if self.outputf != "":
            for src_dir, dirs, files in walk(self.tempDir):
                dst_dir = src_dir.replace(self.tempDir, self.outputf, 1)
                if not path.exists(dst_dir):
                    makedirs(dst_dir)
                for file_ in files:
                    src_file = path.join(src_dir, file_)
                    dst_file = path.join(dst_dir, file_)
                    if path.exists(dst_file):
                        # in case of the src and dst are the same file
                        if path.samefile(src_file, dst_file):
                            continue
                        remove(dst_file)
                    move(src_file, dst_dir)
        subprocess.call(["rm", "-rf", self.tempDir])  # remove the temp folder  

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
        modulefolder = self.tempDir + loop.header[1:] + '/'
        if not path.exists(modulefolder):
            makedirs(modulefolder)
        filename = modulefolder + loop.header[1:]+".fasta"
        fasta = open(filename, 'w')
        fasta.write('>'+loop.header+'\n'+loop.seq+'\n')
        fasta.close()

        # Launch Jar3D on it
        if loop.type == 'h':
            cmd = ["java", "-jar", self.jar3dexec, loop.header[1:]+".fasta", self.HLmotifDir+"/all.txt",
                   loop.header[1:]+".HLloop.csv", loop.header[1:]+".HLseq.csv"]
        else:
            cmd = ["java", "-jar", self.jar3dexec, loop.header[1:]+".fasta", self.ILmotifDir+"/all.txt",
                   loop.header[1:]+".ILloop.csv", loop.header[1:]+".ILseq.csv"]
        nowhere = open(devnull, 'w')
        logfile = open(self.tempDir + "log_of_the_run.sh", 'a')
        logfile.write(' '.join(cmd))
        logfile.write("\n")
        logfile.close()
        chdir(modulefolder)
        subprocess.call(cmd, stdout=nowhere)
        chdir(self.biorseoDir)
        nowhere.close()

        # Retrieve results
        insertion_sites = []
        if loop.type == 'h':
            capstype = "HL"
        else:
            capstype = "IL"
        csv = open(modulefolder + loop.header[1:] +".%sseq.csv" % capstype, 'r')
        l = csv.readline()
        while l:
            if "true" in l:
                insertion_sites.append(InsertionSite(loop, l))
            l = csv.readline()
        csv.close()

        return insertion_sites

    def launch_JAR3D(self, seq_, basename):
        rnasubopt_preds = []
        # Extracting probable loops from RNA-subopt structures
        rna = open(self.tempDir + basename + ".subopt", "r")
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
            loops.append(Loop(">HL%d" % (i+1), seq_[l[0][0]-1:l[0][1]], "h", l))
        for i, l in enumerate(ILs):
            loops.append(Loop(">IL%d" % (i+1), seq_[l[0][0]-1:l[0][1]]+'*'+seq_[l[1][0]-1:l[1][1]], "i", l))
        # Scanning loop subsequences against motif database
        pool = MyPool(processes=cpu_count())
        insertion_sites = [x for y in pool.map(self.launch_JAR3D_worker, loops) for x in y]
        insertion_sites.sort(reverse=True)
        # Writing results to CSV file
        c = 0
        resultsfile = open(self.biorseoDir + "/" + self.tempDir+basename+".sites.csv", "w")
        resultsfile.write("Motif,Rotation,Score,Start1,End1,Start2,End2\n")
        for site in insertion_sites:
            if site.score > 10:
                c += 1
                string = "FOUND with score %d:\t\t possible insertion of motif " % site.score + site.atlas_id
                if site.rotation:
                    string += " (reversed)"
                string += (" on " + site.loop.header + " at positions")
            resultsfile.write(site.atlas_id+',' +
                              str(bool(site.rotation))+",%d" % site.score+',')
            positions = [','.join([str(y) for y in x]) for x in site.position]
            if len(positions) == 1:
                positions.append("-,-")
            resultsfile.write(','.join(positions)+'\n')
        resultsfile.close()

    def launch_BayesPairing(self, module_type, seq_, header_):
        cmd = ["python3", "parse_sequences.py", "-seq", self.biorseoDir + '/' + self.tempDir +
               header_ + ".fa", "-d", module_type, "-interm", "1"]

        logfile = open(self.tempDir + "log_of_the_run.sh", 'a')
        logfile.write(" ".join(cmd))
        logfile.write("\n")
        logfile.close()

        chdir(self.bypdir)
        out = subprocess.check_output(cmd).decode('utf-8')
        BypLog = out.split('\n')
        idx = 0
        l = BypLog[idx]
        while l[:3] != "PUR":
            idx += 1
            l = BypLog[idx]
        insertion_sites = [x for x in ast.literal_eval(l.split(":")[1][1:])]
        if module_type == "rna3dmotif":
            rna = open(self.biorseoDir + "/" + self.tempDir + header_ + ".byp.csv", "w")
        else:
            rna = open(self.biorseoDir + "/" + self.tempDir + header_ + ".bgsubyp.csv", "w")
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


    def launch_BayesPairing2(module_type, seq_, header_, basename):

        cmd = ["python3.7", "parse_sequences.py", "-seq", self.biorseoDir + '/' + self.tempDir + header_ + ".fa", "-samplesize", 1000, "-d", module_type]

        logfile = open(self.tempDir + "log_of_the_run.sh", 'a')
        logfile.write(" ".join(cmd))
        logfile.write("\n")
        logfile.close()

        chdir(self.bypdir)
        out = subprocess.check_output(cmd).decode('utf-8')
        Byp2Log = out.split('\n')

        #remove the 2 first lines of output, and the last one
        Byp2Log.pop(0)
        Byp2Log.pop(0)
        Byp2Log.pop()

        lines = []
        for i in range(len(Byp2Log)):
            line = Byp2Log[i].split()

            #remove "|",  ",", "-", and the sequence
            while "|" in line:
                line.remove("|")
            while "," in line:
                line.remove(",")
            while "-" in line:
                line.remove("-")
            line.pop()

            lines.append(line)


        if module_type=="rna3dmotif":
            rna = open(self.biorseoDir + "/" + self.tempDir + header_ + ".byp2.csv", "w")
        else:
            rna = open(self.biorseoDir + "/" + self.tempDir + header_ + ".bgsubyp2.csv", "w")
        rna.write("Motif,Score,Start1,End1,Start2,End2...\n")

        for line in lines:
            rna.write(module_type)
            for i in range(len(line)-1):
                rna.write(line[i] + ",")
            rna.write(line[-1] + "\n")

        rna.close()



    def execute_job(self, j):
        
        running_stats[0] += 1
        if len(j.cmd_):
            logfile = open(self.tempDir + "log_of_the_run.sh", 'a')
            logfile.write(" ".join(j.cmd_))
            logfile.write("\n")
            logfile.close()
            print("["+str(running_stats[0]+running_stats[2]) +
                  '/'+str(self.jobcount)+"]\t"+" ".join(j.cmd_))
            r = subprocess.call(j.cmd_, timeout=j.timeout_)
        elif j.func_ is not None:
            print("["+str(running_stats[0]+running_stats[2])+'/'+str(self.jobcount) +
                  "]\t"+j.func_.__name__+'('+", ".join([a for a in j.args_])+')')
            try:
                r = j.func_(*j.args_)
            except:
                r = 1
                pass
        running_stats[1] += 1
        return r

    def execute_jobs(self):
        jobs = {}
        self.jobcount = len(self.joblist)
        for job in self.joblist:
            if job.priority_ not in jobs.keys():
                jobs[job.priority_] = {}
            if job.nthreads not in jobs[job.priority_].keys():
                jobs[job.priority_][job.nthreads] = []
            jobs[job.priority_][job.nthreads].append(job)
        nprio = max(jobs.keys())

        for i in range(1,nprio+1):
            if not len(jobs[i].keys()): continue

            # check the thread numbers
            different_thread_numbers = [n for n in jobs[i].keys()]
            different_thread_numbers.sort()

            for n in different_thread_numbers:
                bunch = jobs[i][n]
                if not len(bunch): continue
                pool = MyPool(processes=n)
                pool.map(self.execute_job, bunch)
                pool.close()
                pool.join()

    def list_jobs(self):

        # Read fasta file, which can contain one or several RNAs
        RNAcontainer = []
        if self.outputf != "":
            subprocess.call(["mkdir", "-p", self.outputf])  # Create the output folder
        subprocess.call(["mkdir", "-p", self.tempDir])  # Create the temp folder
        print("loading file %s..." % self.inputfile)
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
                    self.mode = 1 # we switch to batch mode
                header = l[1:-1]
            if c == 0:
                seq = l[:-1].upper()
                if is_canonical_nts(seq):
                    header = header.replace('/', '_').replace('\'','').replace('(','').replace(')','').replace(' ','_').replace('>','')
                    RNAcontainer.append(RNA(header, seq))
                    if not path.isfile(self.tempDir + header + ".fa"):
                        rna = open(self.tempDir + header + ".fa", "w")
                        rna.write(">" + header +'\n')
                        rna.write(seq +'\n')
                        rna.close()
        db.close()

        for nt, number in ignored_nt_dict.items():
            print("ignored %d sequences because of char %c" % (number, nt))
        tot = len(RNAcontainer)
        print("Loaded %d RNAs." % (tot))

        # define job list
        for instance in RNAcontainer:
            
            executable = self.biorseoDir + "/bin/biorseo"
            fastafile = self.tempDir+instance.header+".fa"
            method_type = ""
            ext = ".raw"
            priority = 1

            if self.type == "jar3d":
                ext = ".jar3d"
                method_type = "--jar3dcsv"
                csv = self.tempDir + instance.header + ".sites.csv"

                # RNAsubopt
                self.joblist.append(Job(command=["RNAsubopt", "-i", fastafile, "--outfile="+ instance.header + ".subopt"], priority=1))
                self.joblist.append(Job(command=["mv", instance.header + ".subopt", self.tempDir], priority=2))
                # JAR3D
                self.joblist.append(Job(function=self.launch_JAR3D, args=[instance.seq_, instance.header], priority=3, how_many_in_parallel=1))
                priority = 4
            if self.type == "byp":

                method_type = "--bayespaircsv"

                if self.modules == "desc":
                    ext = ".byp"
                    csv = self.tempDir + instance.header + ".byp.csv"
                    self.joblist.append(Job(function=self.launch_BayesPairing, args=["rna3dmotif", instance.seq_, instance.header], how_many_in_parallel=-1, priority=1))

                    ext = ".byp2"
                    csv = self.tempDir + instance.header + ".byp2.csv"
                    self.joblist.append(Job(function=self.launch_BayesPairing2, args=["rna3dmotif", instance.seq_, instance.header], how_many_in_parallel=-1, priority=1))

                elif self.modules == "bgsu":
                    ext = ".bgsubyp"
                    csv = self.tempDir + instance.header + ".bgsubyp.csv"
                    self.joblist.append(Job(function=self.launch_BayesPairing, args=["3dmotifatlas", instance.seq_, instance.header], how_many_in_parallel=-1, priority=1))

                    ext = ".bgsubyp2"
                    csv = self.tempDir + instance.header + ".bgsubyp2.csv"
                    self.joblist.append(Job(function=self.launch_BayesPairing2, args=["3dmotifatlas", instance.seq_, instance.header], how_many_in_parallel=-1, priority=1))

                priority = 2

            if self.type == "dpm":
                method_type = "--descfolder"
                csv = self.descfolder
            command = [executable, "-s", fastafile ]
            if method_type:
                command += [ method_type, csv ]
            self.finalname =  self.tempDir + instance.header + ext + self.func
            command += [ "-o", self.finalname, "--function", self.func ]
            command += self.forward_options
            self.joblist.append(Job(command=command, priority=priority, timeout=3600, how_many_in_parallel=3))
            

BiorseoInstance(opts)
