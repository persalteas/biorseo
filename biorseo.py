#!/usr/bin/python3
# coding=utf-8
import getopt, multiprocessing, subprocess, sys
import matplotlib.pyplot as plt
from scipy import stats
from os import path, makedirs, getcwd, chdir, devnull, remove, walk
from matplotlib import colors
from math import sqrt
from multiprocessing import cpu_count, Manager
from shutil import move
from ast import literal_eval

# Parse options
try:
    cmd_opts, cmd_args = getopt.getopt( sys.argv[1:],
                                "bc:f:hi:jl:no:O:pt:v",
                             [  "verbose","rna3dmotifs","3dmotifatlas","carnaval","contacts","jar3d","bayespairing","patternmatch","func=",
                                "help","version","seq=","modules-path=", "jar3d-exec=", "bypdir=", "biorseo-dir=", "first-objective=","output=","theta=",
                                "interrupt-limit=", "outputf="])
except getopt.GetoptError as err:
    print(err)
    sys.exit(2)

m = Manager()
running_stats = m.list()
running_stats.append(0) # n_launched
running_stats.append(0) # n_finished
running_stats.append(0) # n_skipped


ignored_nt_dict = {}
def is_canonical_nts(seq):
    """Checks if a nucleotide is ACGU, and stores it in a dictionnary if it's not,
    with an occurrence count."""

    for c in seq[:-1]:
        if c not in "ACGU":
            if c in ignored_nt_dict.keys():
                ignored_nt_dict[c] += 1
            else:
                ignored_nt_dict[c] = 1
            return False
    return True


def absolutize_path(p, directory=False):
    if p[0] != '/':
        p = getcwd() + '/' + p
    if directory and p[-1] != '/':
        p = p + '/'
    return p



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
    """Just a data structure to store module informations."""

    def __init__(self, header, subsequence, looptype, position):
        self.header = header
        self.seq = subsequence
        self.type = looptype
        self.position = position


class InsertionSite:
    """An interval of an RNA sequence where a particular module could be inserted.
    The philosophy of Biorseo is : if this portion can fold like this module, 
    then it may be a loop of the corresponding type."""

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
        """compare two insertion sites scores"""

        return self.score < other.score

    def __gt__(self, other):
        """compare two insertion sites scores"""

        return self.score > other.score


class Job:
    """A class to store the properties of a tool execution, in order to run similar jobs in parallel."""

    def __init__(self, command=None, function=None, args=None, how_many_in_parallel=0, priority=1, timeout=None):
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
    """Just a data structure gathering header, sequence and length."""

    def __init__(self, header, seq):
        self.seq_ = seq
        self.header = header
        self.length = len(seq)


class BiorseoInstance:
    """A run of the biorseo tool, to predict one or several RNA sequences' folding(s),
    including all the necessary previous run of other tools."""

    def __init__(self, opts):
        # set default run type options
        self.type = "dpm"       # direct pattern mathcing
        self.modules = "desc"   # ...with Rna3dMotifs "DESC" modules
        self.func = 'B'         # ...and function B
        self.forward_options = []   # options to pass to the C++ biorseo
        self.jobcount = 0
        self.joblist = []

        # set default file output locations
        self.finalname = ""     
        self.outputf = ""       # A folder where to output the computation files
        if path.exists("/biorseo/results"): # docker image default
            self.outputf = "/biorseo/results"
        self.output = ""        # A file to store the solutions

        # set default data input locations
        self.mode = 0           # default is single sequence mode
        self.inputfile = ""
        self.jar3d_exec = "/jar3d_2014-12-11.jar"
        self.bypdir = "/byp/src" # Docker image locations
        self.HL_motif_dir = "/modules/BGSU/HL/3.2/lib"
        self.IL_motif_dir = "/modules/BGSU/IL/3.2/lib"
        self.desc_folder = "/modules/DESC"
        self.rin_folder = "/modules/RIN/Subfiles"
        self.json_folder = "/modules/ISAURE"
        self.biorseo_dir = "/biorseo"
        self.run_dir = path.dirname(path.realpath(__file__))
        self.temp_dir = "temp/"

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(  "Biorseo, Bi-Objective RNA Structure Efficient Optimizer\n"
                        "Bio-objective integer linear programming framework to predict RNA secondary structures by including known RNA modules.\n"
                        "Developped by Louis Becquey (louis.becquey@univ-evry.fr), 2018-2020\n\n")
                print("Usage:\tYou must provide:\n\t1) a FASTA input file with -i,\n\t2) a module type with --rna3dmotifs, --carnaval, --contacts or --3dmotifatlas"
                      "\n\t3) one module placement method in { --patternmatch, --jar3d, --bayespairing }\n\t4) one scoring function with --func A, B, C or D"
                      "\n\n\tIf you are not using the Docker image: \n\t5) --modules-path, --biorseo-dir and (--jar3d-exec or --bypdir)")
                print()
                print("Options:")
                print("-h [ --help ]\t\t\tPrint this help message")
                print("--version\t\t\tPrint the program version")
                print("-i [ --seq=… ]\t\t\tFASTA file with the query RNA sequence")
                print("--rna3dmotifs\t\t\tUse DESC modules from Djelloul & Denise, 2008")
                print("--carnaval\t\t\tUse RIN modules from Reinharz & al, 2018")
                print("--3dmotifatlas\t\t\tUse the HL and IL loops from BGSU's 3D Motif Atlas (updated)")
                print("--contacts\t\t\tUse .json motif from Isaure")

                print("-p [ --patternmatch ]\t\tUse regular expressions to place modules in the sequence (requires --rna3dmotifs or --carnaval)")
                print("-j [ --jar3d ]\t\t\tUse JAR3D to place modules in the sequence (requires --3dmotifatlas)")
                print("-b [ --bayespairing ]\t\tUse BayesPairing2 to place modules in the sequence (requires --rna3dmotifs or --3dmotifatlas)")
                print("-o [ --output=… ]\t\tFile to summarize the results")
                print("-O [ --outputf=… ]\t\tFolder where to output result and temp files")
                print("-f [ --func=… ]\t\t\t(A, B, C or D, default is B)"
                      " Objective function to score module insertions:\n\t\t\t\t  (A) insert big modules (B) insert light, high-order modules"
                      "\n\t\t\t\t  (c) insert modules which score well with the sequence\n\t\t\t\t  (D) insert light, high-order modules which score well with the sequence."
                      "\n\t\t\t\t  C and D require cannot be used with --patternmatch.")
                print("-c [ --first-objective=… ]\t(default 1) Objective to solve in the mono-objective portions of the algorithm."
                      "\n\t\t\t\t  (1) is the module objective given by --func, (2) is the expected accuracy of the structure.")
                print("-l [ --limit=… ]\t\t(default 500) Number of solutions in the Pareto set from which"
                      "\n\t\t\t\t  we give up the computation, before eliminating secondary structure doublons.")
                print("-t [ --theta=… ]\t\t(default 0.001) Pairing-probability threshold to consider or not the possibility of pairing")
                print("-n [ --disable-pseudoknots ]\tAdd constraints to explicitly forbid the formation of pseudoknots")
                print("-v [ --verbose ]\t\tPrint what is happening to stdout")
                print("--modules-path=…\t\tPath to the modules data.\n\t\t\t\t  The folder should contain modules in the DESC format as output by Djelloul & Denise's"
                      "\n\t\t\t\t  'catalog' program for use with --rna3dmotifs, or the IL/ and HL/ folders\n\t\t\t\t  from a release of the RNA 3D Motif Atlas "
                      "for use with --3dmotifatlas, or the\n\t\t\t\t  data/modules/RIN/Subfiles/ folder for use with --carnaval.\n\t\t\t\t  Consider placing these files on a fast I/O device (NVMe SSD, ...)")
                print("--jar3d-exec=…\t\t\tPath to the jar3d executable.\n\t\t\t\t  Default is /jar3d_2014-12-11.jar, you should use this option if you run"
                      "\n\t\t\t\t  BiORSEO from outside the docker image.")
                print("--bypdir=…\t\t\tPath to the BayesParing src folder.\n\t\t\t\t  Default is /byp/src, you should use this option if you run"
                      "\n\t\t\t\t  BiORSEO from outside the docker image.")
                print("--biorseo-dir=…\t\t\tPath to the BiORSEO root directory.\n\t\t\t\t  Default is /biorseo, you should use this option if you run"
                      "\n\t\t\t\t  BiORSEO from outside the docker image. Use the FULL path.")
                print("\nExamples:")
                print("biorseo.py -i myRNA.fa -O myResultsFolder/ --rna3dmotifs --patternmatch --func B")
                print("biorseo.py -i myRNA.fa -O myResultsFolder/ --3dmotifatlas --jar3d --func B -l 800")
                print("biorseo.py -i myRNA.fa -v --3dmotifatlas --bayespairing --func D")
                print("\nThe allowed module/placement-method/function combinations are:\n")
                print("                --patternmatch  --bayespairing    --jar3d")
                print("--rna3dmotifs     A. B.           A. B. C. D.")
                print("--3dmotifatlas                    A. B. C. D.     A. B. C. D.")
                print("--carnaval        A. B.\n")
                sys.exit()
            elif opt == "-i" or opt == "--seq":
                self.inputfile = arg
            elif opt == "-O" or opt == "--outputf":
                self.outputf = absolutize_path(arg, directory=True) # output folder
            elif opt == "-o" or opt == "--output":
                self.output = absolutize_path(arg) # output file
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
            elif opt == "--carnaval":
                self.modules = "rin"
            elif opt == "--contacts":
                self.modules = "json"
            elif opt == "--rna3dmotifs":
                self.modules = "desc"
            elif opt == "--3dmotifatlas":
                self.modules = "bgsu"
            elif opt == "--modules-path":
                self.HL_motif_dir = absolutize_path(arg, directory=True) + "HL/3.2/lib"
                self.IL_motif_dir = absolutize_path(arg, directory=True) + "IL/3.2/lib"
                self.desc_folder = absolutize_path(arg, directory=True)
                self.rin_folder = absolutize_path(arg, directory=True)
                self.json_folder = absolutize_path(arg, directory=True)
                print("Looking for modules in", arg)
            elif opt == "--jar3d-exec":
                self.jar3d_exec = absolutize_path(arg)
                print("Using ", arg)
            elif opt == "--bypdir":
                self.bypdir = absolutize_path(arg, directory=True)
                print("Using trained BayesPairing in", arg)
            elif opt == "--biorseo-dir":
                self.biorseo_dir = absolutize_path(arg, directory=True)
            elif opt == "--version":
                subprocess.run([self.biorseo_dir+"bin/biorseo", "--version"])
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

        # Check the argument combination is OK
        self.check_args()

        if self.outputf != "":
            print("saving files to", self.outputf)

        # create jobs
        self.list_jobs()

        # run them
        self.execute_jobs()

        # locate the results at the right place
        if self.output != "":
            subprocess.run(["mv", self.temp_dir+self.finalname.split('/')[-1], self.output], check=True)
        if self.outputf != "":
            for src_dir, _, files in walk(self.temp_dir):
                dst_dir = src_dir.replace(self.temp_dir, self.outputf, 1)
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
        subprocess.run(["rm", "-rf", self.temp_dir], check=True)  # remove the temp folder

    def check_args(self):
        """Checks that the argument combination passed by user is a realistic project"""

        warning = "ERROR: The argument list you passed contains errors:"
        issues = False
        if self.modules == "desc" and self.type == "jar3d":
            issues = True
            print(warning)
            print("/!\\ Using jar3d requires the 3D Motif Atlas modules. Use --3dmotifatlas instead of --rna3dmotifs or --carnaval.")
        if self.modules == "rin" and self.type != "dpm":
            issues = True
            print(warning)
            print("/!\\ CaRNAval does not support placement tools (yet). Please use it with --patternmatch, not --jar3d nor --bayespairing.")
        if self.modules == "json" and self.type != "dpm":
            issues = True
            print(warning)
            print("/!\\ Contacts does not support placement tools (yet). Please use it with --patternmatch, not --jar3d nor --bayespairing.")
        if self.modules == "bgsu" and self.type == "dpm":
            issues = True
            print(warning)
            print("/!\\ Cannot place the Atlas loops by direct pattern matching. Please use a dedicated tool --jar3d or --bayespairing to do so.")

        if issues:
            print("\nUsage:\tYou must provide:\n\t1) a FASTA input file with -i,\n\t2) one module type in { --rna3dmotifs, --carnaval, --3dmotifatlas, --contacts }"
                  "\n\t3) one module placement method in { --patternmatch, --jar3d, --bayespairing }\n\t4) one scoring function with --func A, B, C or D"
                  "\n\n\tIf you are not using the Docker image: \n\t5) --modules-path, --biorseo-dir and (--jar3d-exec or --bypdir)")
            print("\nThe allowed module/placement-method/function combinations are:\n")
            print("                --patternmatch  --bayespairing    --jar3d")
            print("--rna3dmotifs     A. B.           A. B. C. D.")
            print("--3dmotifatlas                    A. B. C. D.     A. B. C. D.")
            print("--carnaval        A. B.\n")
            exit(1)

    def enumerate_loops(self, s):
        """Lists all the loop positions in a dot-bracket notation"""

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

        return loops

    def launch_JAR3D_worker(self, loop):
        """Launches a jar3d search of specific loop types (IL or HL) on a RNA subsequence (fasta file)"""

        # write motif to a file
        modulefolder = self.temp_dir + loop.header[1:] + '/'
        if not path.exists(modulefolder):
            makedirs(modulefolder)
        filename = modulefolder + loop.header[1:]+".fasta"
        fasta = open(filename, 'w')
        fasta.write('>'+loop.header+'\n'+loop.seq+'\n')
        fasta.close()

        # Launch Jar3D on it
        if loop.type == 'h':
            cmd = ["java", "-jar", self.jar3d_exec, loop.header[1:]+".fasta", self.HL_motif_dir+"/all.txt",
                   loop.header[1:]+".HLloop.csv", loop.header[1:]+".HLseq.csv"]
        else:
            cmd = ["java", "-jar", self.jar3d_exec, loop.header[1:]+".fasta", self.IL_motif_dir+"/all.txt",
                   loop.header[1:]+".ILloop.csv", loop.header[1:]+".ILseq.csv"]

        with open(self.temp_dir + "log_of_the_run.sh", 'a') as logfile:
            logfile.write(' '.join(cmd) + '\n')

        chdir(modulefolder)

        with open(devnull, 'w') as nowhere:
            subprocess.run(cmd, stdout=nowhere)

        chdir(self.biorseo_dir)

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
        """Identify loops in a RNA sequence from RNAsubopt results and search
        for 3D motif Atlas modules to place on it using jar3d."""

        rnasubopt_preds = []
        
        # Extracting probable loops from RNA-subopt structures
        rna = open(self.temp_dir + basename + ".subopt", "r")
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
        insertion_sites = [ x for y in pool.map(self.launch_JAR3D_worker, loops) for x in y ]
        pool.close()
        pool.join()
        insertion_sites.sort(reverse=True)

        # Writing results to CSV file
        c = 0
        resultsfile = open(self.biorseo_dir + "/" + self.temp_dir+basename+".sites.csv", "w")
        resultsfile.write("Motif,Rotation,Score,Start1,End1,Start2,End2\n")
        for site in insertion_sites:
            if site.score > 10:
                c += 1
                string = "FOUND with score %d:\t\t possible insertion of motif " % site.score + site.atlas_id
                if site.rotation:
                    string += " (reversed)"
                string += (" on " + site.loop.header + " at positions")
            resultsfile.write(site.atlas_id+',' + str(bool(site.rotation))+",%d" % site.score+',')
            positions = [','.join([str(y) for y in x]) for x in site.position]
            if len(positions) == 1:
                positions.append("-,-")
            resultsfile.write(','.join(positions)+'\n')
        resultsfile.close()

    def launch_BayesPairing(self, module_type, seq_, header_):
        """DEPRECATED: now we are using Bayespairing 2.
        Searches for module occurrences in the provided RNA sequence using BayesPairing 1."""

        # Run BayePairing
        cmd = ["python3", "parse_sequences.py", "-seq", self.biorseo_dir + '/' + self.temp_dir + header_ + ".fa", "-d", module_type, "-interm", "1"]

        logfile = open(self.biorseo_dir + "/" + self.temp_dir + "log_of_the_run.sh", 'a')
        logfile.write(" ".join(cmd))
        logfile.write("\n")
        logfile.close()

        chdir(self.bypdir)
        out = subprocess.check_output(cmd).decode('utf-8')
        BypLog = out.split('\n')

        # Extract results from command output to file
        idx = 0
        l = BypLog[idx]
        while l[:3] != "PUR":
            idx += 1
            l = BypLog[idx]
        insertion_sites = [x for x in literal_eval(l.split(":")[1][1:])]
        if module_type == "rna3dmotif":
            rna = open(self.biorseo_dir + "/" + self.temp_dir + header_ + ".byp.csv", "w")
        else:
            rna = open(self.biorseo_dir + "/" + self.temp_dir + header_ + ".bgsubyp.csv", "w")
        rna.write("Motif,Score,Start1,End1,Start2,End2...\n")
        for i, module in enumerate(insertion_sites):
            if len(module):
                for (score, positions, _) in zip(*[iter(module)]*3):
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

    def launch_BayesPairing2(self, module_type, seq_, header_):
        """Search for module occurrences in the provided RNA sequence using BayesPairing 2"""

        # Run BayesPairing 2
        if module_type=="rna3dmotif":
            BP2_type = "rna3dmotif"
        else:
            BP2_type = "3dMotifAtlas_ALL"

        cmd = ["python3", "parse_sequences.py", "-seq", self.biorseo_dir + '/' + self.temp_dir + header_ + ".fa", "-samplesize", "1000", "-d", BP2_type ]
        logfile = open(self.temp_dir + "log_of_the_run.sh", 'a')
        logfile.write(" ".join(cmd))
        logfile.write("\n")
        logfile.close()

        chdir(self.bypdir)
        out = subprocess.check_output(cmd).decode('utf-8')
        Byp2Log = out.splitlines()

        # remove what is not in the original input
        Byp2Log.pop(0)
        Byp2Log.pop(0)
        Byp2Log.pop()
        Byp2Log.pop()

        # remove the 2 first lines of output
        Byp2Log.pop(0)
        Byp2Log.pop(0)

        # Extract results from command line output to file
        lines = []
        for i in range(len(Byp2Log)):
            line = Byp2Log[i].replace("|", ' ').replace(",", ' ').split()

            if line != []:
                if "=" in line[0]: #skip the "| MODULE  N HITS  PERCENTAGE  |" part
                    break

                module_sequence = line.pop().split("&") #remove the sequence

                if line != []:
                    new_line = [line[0], line[1]]

                    num_comp = 0
                    position_index = 2

                    while num_comp < len(module_sequence):

                        len_comp = 0
                        comp = module_sequence[num_comp]

                        element = line[position_index].split("-")
                        new_line.append(element[0])

                        if position_index >= len(line):
                            print("!!! Skipped one BP2 result : positions not matching sequence length\n")
                            new_line = []
                            break

                        while len_comp < len(comp):

                            element = line[position_index].split("-")

                            if len(element)==1:
                                len_comp += 1

                            else:
                                len_comp += int(element[1]) - int(element[0]) + 1

                            position_index += 1

                        new_line.append(element[-1])
                        num_comp += 1

                    if new_line != [] :
                        lines.append(new_line)

        if module_type=="rna3dmotif":
            rna = open(self.biorseo_dir + "/" + self.temp_dir + header_ + ".byp2.csv", "w")
        else:
            rna = open(self.biorseo_dir + "/" + self.temp_dir + header_ + ".bgsubyp2.csv", "w")

        rna.write("Motif,Score,Start1,End1,Start2,End2...\n")

        for line in lines:
            rna.write(module_type)
            for i in range(len(line)-1):
                rna.write(line[i] + ",")
            rna.write(line[-1] + "\n")

        rna.close()

    def execute_job(self, j):
        """Execute the command or function stored in a Job class object."""

        running_stats[0] += 1
        if j.cmd_ is not None:
            logfile = open(self.biorseo_dir + "/" + self.temp_dir + "log_of_the_run.sh", 'a')
            logfile.write(" ".join(j.cmd_))
            logfile.write("\n")
            logfile.close()
            print("["+str(running_stats[0]+running_stats[2]) +
                  '/'+str(self.jobcount)+"]\t"+" ".join(j.cmd_))
            r = subprocess.run(j.cmd_, timeout=j.timeout_)
        elif j.func_ is not None:
            print("["+str(running_stats[0]+running_stats[2])+'/'+str(self.jobcount) +
                  "]\t"+j.func_.__name__+'('+", ".join([a for a in j.args_])+')')
            try:
                if j.args_ is not None:
                    r = j.func_(*j.args_)
                else:
                    r = j.func_()
            except Exception as e:
                r = 1
                print("\033[31m", e, "\033[0m")
                pass
        running_stats[1] += 1
        return r

    def execute_jobs(self):
        """Groups job by priority and ability to be run in parallel,
        and runs them."""

        jobs = {}
        self.jobcount = len(self.joblist)
        for job in self.joblist:
            if job.priority_ not in jobs.keys():
                jobs[job.priority_] = {}
            if job.nthreads not in jobs[job.priority_].keys():
                jobs[job.priority_][job.nthreads] = []
            jobs[job.priority_][job.nthreads].append(job)
        nprio = max(jobs.keys())

        if len(jobs) > 1 :
            for i in range(1,nprio+1):
                if not len(jobs[i].keys()): continue

                # check the thread numbers
                different_thread_numbers = [ n for n in jobs[i].keys() ] 
                different_thread_numbers.sort()
                for n in different_thread_numbers:
                    bunch = jobs[i][n]
                    if not len(bunch): continue
                    pool = MyPool(processes=n)
                    pool.map(self.execute_job, bunch)
                    pool.close()
                    pool.join()
        else:
            for job in self.joblist:
                self.execute_job(job)

    def list_jobs(self):
        """Determines the required tool runs we need before running the C++ Biorseo,
        and creates a job list to run."""

        if self.outputf != "":
            subprocess.run(["mkdir", "-p", self.outputf])  # Create the output folder
        subprocess.run(["mkdir", "-p", self.temp_dir])  # Create the temp folder

        # Read fasta file, which can contain one or several RNAs
        RNAcontainer = []
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
                    if not path.isfile(self.temp_dir + header + ".fa"):
                        rna = open(self.temp_dir + header + ".fa", "w")
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

            executable = self.biorseo_dir + "bin/biorseo"
            fastafile = self.temp_dir+instance.header+".fa"
            method_type = ""
            priority = 1

            if self.type == "jar3d":
                ext = ".jar3d"
                method_type = "--jar3dcsv"
                csv = self.temp_dir + instance.header + ".sites.csv"

                # RNAsubopt
                self.joblist.append(Job(command=["RNAsubopt", "-i", fastafile, "--outfile="+ instance.header + ".subopt"], priority=1))
                self.joblist.append(Job(command=["mv", instance.header + ".subopt", self.temp_dir], priority=2))
                # JAR3D
                self.joblist.append(Job(function=self.launch_JAR3D, args=[instance.seq_, instance.header], priority=3, how_many_in_parallel=1))
                priority = 4

            if self.type == "byp":

                method_type = "--bayespaircsv"

                if self.modules == "desc":
                    ext = ".byp2"
                    csv = self.temp_dir + instance.header + ".byp2.csv"
                    self.joblist.append(Job(function=self.launch_BayesPairing2, args=["rna3dmotif", instance.seq_, instance.header], how_many_in_parallel=1, priority=1))

                elif self.modules == "bgsu":
                    ext = ".bgsubyp2"
                    csv = self.temp_dir + instance.header + ".bgsubyp2.csv"
                    self.joblist.append(Job(function=self.launch_BayesPairing2, args=["3dmotifatlas", instance.seq_, instance.header], how_many_in_parallel=1, priority=1))

                priority = 2

            if self.type == "dpm":

                if self.modules == "desc":
                    method_type = "--descfolder"
                    csv = self.desc_folder
                    ext = ".desc_pm"
                
                elif self.modules == "rin":
                    method_type = "--rinfolder"
                    csv = self.rin_folder
                    ext = ".rin_pm"

                elif self.modules == "json":
                    method_type = "--jsonfolder"
                    csv = self.json_folder
                    ext = ".json_pm"

            command = [ executable, "-s", fastafile ]
            if method_type:
                command += [ method_type, csv ]
            self.finalname =  self.temp_dir + instance.header + ext + self.func
            command += [ "-o", self.finalname, "--function", self.func ]
            command += self.forward_options
            self.joblist.append(Job(command=command, priority=priority, timeout=3600, how_many_in_parallel=3))


BiorseoInstance(cmd_opts)
