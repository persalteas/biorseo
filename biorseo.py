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
outputDir = biorseoDir + "/results/"
HLmotifDir = biorseoDir + "/data/modules/BGSU/HL/3.2/lib"
ILmotifDir = biorseoDir + "/data/modules/BGSU/IL/3.2/lib"
descfolder = biorseoDir + "/data/modules/DESC"


# ================== CLASSES AND FUNCTIONS ================================


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


class BiorseoInstance:
    def __init__(self, argv):
        self.type = "dpm"
        self.modules = "desc"
        self.func = 'B'
        self.outputf = outputDir
        try:
            opts, args = getopt.getopt(
                argv, "hil::o:", ["type=", "func=", "modules="])
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
            elif opt == "-l":
                self.inputfile = arg
                self.mode = 1  # batch mode
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
        insertion_sites = [x for y in pool.map(
            launch_JAR3D_worker, loops) for x in y]
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
            resultsfile.write(site.atlas_id+',' +
                              str(bool(site.rotation))+",%d" % site.score+',')
            positions = [','.join([str(y) for y in x]) for x in site.position]
            if len(positions) == 1:
                positions.append("-,-")
            resultsfile.write(','.join(positions)+'\n')
        resultsfile.close()

    def launch_BayesPairing(self, module_type, seq_, header_, basename):
        chdir(bypdir)

        cmd = ["python3", "parse_sequences.py", "-seq", outputDir +
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
            rna = open(outputDir + basename + ".byp.csv", "w")
        else:
            rna = open(outputDir + basename + ".bgsubyp.csv", "w")
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
                running_stats[2] += 1
                print("["+str(running_stats[0]+running_stats[2]) +
                      '/'+str(jobcount)+"]\tSkipping a finished job")
                return 0
        running_stats[0] += 1
        if len(j.cmd_):
            logfile = open("log_of_the_run.sh", 'a')
            logfile.write(" ".join(j.cmd_))
            logfile.write("\n")
            logfile.close()
            print("["+str(running_stats[0]+running_stats[2]) +
                  '/'+str(jobcount)+"]\t"+" ".join(j.cmd_))
            r = subprocess.call(j.cmd_, timeout=j.timeout_)
        elif j.func_ is not None:
            print("["+str(running_stats[0]+running_stats[2])+'/'+str(jobcount) +
                  "]\t"+j.func_.__name__+'('+", ".join([a for a in j.args_])+')')
            try:
                r = j.func_(*j.args_)
            except:
                r = 1
                pass
        if r:
            fails.append(j)
        running_stats[1] += 1
        return r

    def check_existence(self, datatype, method, function, with_PK, basename):
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

    def check_existence(self, datatype, method, basename):
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


if __name__ == "__main__":
    BiorseoInstance(sys.argv)
