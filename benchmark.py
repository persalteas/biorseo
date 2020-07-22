#!/usr/bin/python3.8
#coding=utf-8

# typical usage : ./benchmark.py data/sec_structs/verified_secondary_structures_database.dbn data/sec_structs/pseudoknots.dbn data/sec_structs/applications.dbn

# the .dbn files should be formatted the following way:
# > header of the sequence (somecode)
# ACGUACGUACGUACGUACGU
# ...(((...((...))))).
# > header of the next sequence (somecode2)
# ACGUACGUACGGGCGUACGU
# ...(((..........))).

from sys import argv
from scipy import stats
import subprocess
from os import path, makedirs, getcwd, chdir, devnull
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from math import sqrt, ceil
from multiprocessing import Pool, cpu_count, Manager
import multiprocessing
import ast, time

# ================== DEFINITION OF THE PATHS ==============================

biorseoDir = path.realpath(".")
jar3dexec = "/opt/jar3d_2014-12-11.jar"
bypdir = "/opt/BayesPairing/bayespairing/src"
byp2dir = "/opt/rnabayespairing2.git/bayespairing/src"
moipdir = "/opt/RNAMoIP/Src/RNAMoIP.py"
biokopdir = "/opt/biokop/biokop"
runDir = path.dirname(path.realpath(__file__))
RNAStrandFile = argv[1]
PseudobaseFile = argv[2]
StudyCaseFile = argv[3]
outputDir = biorseoDir + "/benchmark_results/"
HLmotifDir = biorseoDir + "/data/modules/BGSU/HL/3.2/lib"
ILmotifDir = biorseoDir + "/data/modules/BGSU/IL/3.2/lib"
descfolder = biorseoDir + "/data/modules/DESC"
rinfolder = biorseoDir + "/data/modules/CaRNAval/Subfiles/"

# Create some folders to store the results
subprocess.call(["mkdir", "-p", outputDir])
subprocess.call(["mkdir", "-p", outputDir + "PK/"])
subprocess.call(["mkdir", "-p", outputDir + "noPK/"])

m = Manager()
running_stats = m.list()
running_stats.append(0) # n_launched
running_stats.append(0) # n_finished
running_stats.append(0) # n_skipped

# ================== CLASSES AND FUNCTIONS ================================

"""
class NoDaemonProcess(multiprocessing.Process):
	# make 'daemon' attribute always return False
	def _get_daemon(self):
		return False
	def _set_daemon(self, value):
		pass
	daemon = property(_get_daemon, _set_daemon)

# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(multiprocessing.pool.Pool):
	Process = NoDaemonProcess
"""


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
	def __init__(self, results, command=[], function=None, args=[], how_many_in_parallel=0, priority=1, timeout=None, checkFunc=None, checkArgs=[], label=""):
		self.cmd_ = command
		self.func_ = function
		self.args_ = args
		self.checkFunc_ = checkFunc
		self.checkArgs_ = checkArgs
		self.results_file = results
		self.priority_ = priority
		self.timeout_ = timeout
		self.comp_time = -1 # -1 is not executed yet
		self.label = label
		if not how_many_in_parallel:
			self.nthreads = cpu_count()
		elif how_many_in_parallel == -1:
			self.nthreads = cpu_count() - 1
		else:
			self.nthreads = how_many_in_parallel
		self.nthreads = ceil(self.nthreads/2)
		self.useless_bool = False

	def __str__(self):
		if self.func_ is None:
			s = f"{self.priority_}({self.nthreads}) [{self.comp_time}]\t{j.label:25}" + " ".join(self.cmd_)
		else:
			s = f"{self.priority_}({self.nthreads}) [{self.comp_time}]\t{j.label:25}{self.func_.__name__}(" + " ".join([str(a) for a in self.args_]) + ")"
		return s

def execute_job(j):
	# Check if you really need to execute it
	if path.isfile(j.results_file) or ((j.checkFunc_ is not None) and j.checkFunc_(*j.checkArgs_)):
		running_stats[2] += 1
		print(f"[{running_stats[0]+running_stats[2]}/{jobcount}]\tSkipping {j.label} (already finished)")
		return (0, 0)

	running_stats[0] += 1

	# Add the job to log file and run
	if len(j.cmd_):
		logfile = open(runDir + "/log_of_the_run.sh", 'a')
		logfile.write(" ".join(j.cmd_))
		logfile.write("\n")
		logfile.close()
		print(f"[{running_stats[0]+running_stats[2]}/{jobcount}]\t{j.label}")
		start_time = time.time()
		r = subprocess.call(j.cmd_, timeout=j.timeout_)
		end_time = time.time()
	elif j.func_ is not None:
		print(f"[{running_stats[0]+running_stats[2]}/{jobcount}]\t{j.func_.__name__}({', '.join([str(a) for a in j.args_])})")
		start_time = time.time()
		r = j.func_(*j.args_)
		end_time = time.time()

	# Job is finished
	running_stats[1] += 1
	t = end_time - start_time
	return (t,r)

def launch_JAR3D_worker(loop):
	# write motif to a file
	newpath = getcwd()+'/' + loop.rna_name + '/'+ loop.header[1:]
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

	logfile = open(runDir + "/log_of_the_run.sh", 'a')
	logfile.write(' '.join(cmd)+"\n")
	logfile.close()
	subprocess.run(cmd, stdout=subprocess.DEVNULL)

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

	return insertion_sites

def launch_JAR3D(seq_, basename):
	rnasubopt_preds = []

	# Extracting probable loops from RNA-subopt structures
	rna = open(outputDir + f"noPK/{basename}.subopt", "r")
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
		loops.append(Loop(basename, ">HL%d" % (i+1), seq_[l[0][0]-1:l[0][1]], "h", l))
	for i, l in enumerate(ILs):
		loops.append(Loop(basename, ">IL%d" % (i+1), seq_[l[0][0]-1:l[0][1]]+'*'+seq_[l[1][0]-1:l[1][1]], "i", l))

	# Scanning loop subsequences against motif database
	if not path.exists(basename):
		makedirs(basename)
	p = MyPool(processes=cpu_count())
	insertion_sites = [x for y in p.map(launch_JAR3D_worker, loops) for x in y]
	p.close()
	p.join()
	insertion_sites.sort(reverse=True)
	subprocess.call(["rm", "-r", basename])

	# Writing results to CSV file
	c = 0
	resultsfile = open(outputDir+basename+".bgsu_jar3d.csv", "w")
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

	cmd = ["python3.8","parse_sequences.py","-seq",outputDir + basename + ".fa", "-d", module_type, "-interm","1"]

	logfile = open(runDir + "/log_of_the_run.sh", 'a')
	logfile.write(" ".join(cmd))
	logfile.write("\n")
	logfile.close()

	chdir(bypdir)
	out = subprocess.check_output(cmd).decode('utf-8')
	BypLog = out.split('\n')
	idx = 0
	l = BypLog[idx]
	while l[:3] != "PUR":
		idx += 1
		l = BypLog[idx]
	insertion_sites = [ x for x in ast.literal_eval(l.split(":")[1][1:])]
	if module_type=="carnaval":
		rna = open(outputDir + basename + ".rin_byp.csv", "w")
	elif module_type=="rna3dmotif":
		rna = open(outputDir + basename + ".desc_byp.csv", "w")
	else:
		rna = open(outputDir + basename + ".bgsu_byp.csv", "w")
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



def launch_BayesPairing2(module_type, seq_, header_, basename):

	if module_type=="rna3dmotif":
		BP2_type = "rna3dmotif"
	else:
		BP2_type = "3DmotifAtlas_ALL"

	cmd = ["python3.7", "parse_sequences.py", "-seq", outputDir+basename+".fa", "-samplesize", "1000", "-d", BP2_type]

	logfile = open(runDir + "/log_of_the_run.sh", 'a')
	logfile.write(" ".join(cmd))
	logfile.write("\n")
	logfile.close()

	chdir(byp2dir)
	out = subprocess.check_output(cmd).decode('utf-8')
	Byp2Log = out.splitlines()

    #remove what is not in the original input
    Byp2Log.pop(0)
    Byp2Log.pop(0)
    Byp2Log.pop()
    Byp2Log.pop()

    #remove the 2 first lines of output
    Byp2Log.pop(0)
    Byp2Log.pop(0)
    
	lines = []
	for i in range(len(Byp2Log)):
		line = Byp2Log[i].replace("|", ' ').replace(",", ' ').replace("-", ' ').split()

		if line != []:
			if "=" in line[0]: #skip the "| MODULE  N HITS  PERCENTAGE  |" part
				break
			line.pop() #remove the sequence

			if line != []:
				lines.append(line)
				#print(line)


	if module_type=="rna3dmotif":
		rna = open(outputDir + basename + ".desc_byp2.csv", "w")
	else:
		rna = open(outputDir + basename + ".bgsu_byp2.csv", "w")
	rna.write("Motif,Score,Start1,End1,Start2,End2...\n")

	for line in lines:
		rna.write(module_type)
		for i in range(len(line)-1):
			rna.write(line[i] + ",")
		rna.write(line[-1] + "\n")

	rna.close()



def launch_RNAMoIP_worker(c):
	# launch gurobi
	try:
		out = subprocess.check_output(c).decode("utf-8")
	except subprocess.CalledProcessError as e:
		print(e.output)
		exit()
	gurobiLog = out.split('\n')

	# parse output
	idx = 0
	l = gurobiLog[idx]
	solution = ""
	nsolutions = 0
	while l != "Corrected secondary structure:" and l != " NO SOLUTIONS!":
		if l[:19] == "Optimal solution nb:":
			nsolutions = int(l.split(' ')[-1])
		idx += 1
		l = gurobiLog[idx]
	if nsolutions > 1:
		print("WARNING: RNA-MoIP found several solutions !")
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

	return solution, nmotifs, score

def launch_RNAMoIP(seq_, header_, basename, one_by_one):
	RNAMoIP = moipdir

	# read RNAsubopt predictions
	rnasubopt_preds = []
	rna = open(outputDir + f"noPK/{basename}.subopt", "r")
	lines = rna.readlines()
	rna.close()
	for i in range(2, len(lines)):
		ss = lines[i].split(' ')[0]
		if ss not in rnasubopt_preds:
			rnasubopt_preds.append(ss)
	if one_by_one:
		logfile = open(runDir + "/log_of_the_run.sh", 'a')
		rna = open(outputDir + f"noPK/{basename}.moip", "w")
		rna.write(header_+'\n')
		rna.write(seq_+'\n')
		for ss in rnasubopt_preds:
			c = ["gurobi.sh", RNAMoIP, "-s", f"{seq_}", "-ss", f"{ss}", "-d", descfolder]
			logfile.write(" ".join(c)+'\n')
			solution, nmotifs, score = launch_RNAMoIP_worker(c)
			rna.write("{}\t{}\t{}\n".format(solution, nmotifs, score))
		rna.close()
		logfile.close()
	else:
		subopts = open(f"{runDir}/{basename}.temp_subopt", "w")
		for ss in rnasubopt_preds:
			subopts.write(ss + '\n')
		subopts.close()
		c = ["gurobi.sh", RNAMoIP, "-s", f"{seq_}", "-ss", f"{runDir}/{basename}.temp_subopt", "-d", descfolder]

		logfile = open(runDir + "/log_of_the_run.sh", 'a')
		logfile.write(" ".join(c))
		logfile.write("\n")
		logfile.close()

		solution, nmotifs, score = launch_RNAMoIP_worker(c)
		rna = open(outputDir + f"noPK/{basename}.moipc", "w")
		rna.write(header_+'\n')
		rna.write(seq_+'\n')
		rna.write(f"{solution}\t{nmotifs}\t{score}\n")
		rna.close()
		subprocess.call(["rm", f"{runDir}/{basename}.temp_subopt"])

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
	return 2*tp/(2*tp+fp+fn)

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


class Loop:
	def __init__(self, rna_name, header, subsequence, looptype, position):
		self.rna_name = rna_name
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
	def __init__(self, parent_rna, tool, data_source=None, placement_method=None, obj_func=None, PK=True, flat=False):
		self.parent_rna = parent_rna

		# defintion of the method (theoretical)
		self.tool = tool
		self.data_source = data_source
		self.placement_method = placement_method
		self.func = obj_func
		self.allow_pk = PK
		self.label = self.get_label()

		# things related to execution
		self.joblist = []
		self.flat = flat
		self.build_job_list()

		# descriptors of the results set:
		self.predictions = []
		self.scores = []
		self.ninsertions = []
		self.max_mcc = 0
		self.min_mcc = 0
		self.avg_mcc = 0
		self.max_f1 = 0
		self.min_f1 = 0
		self.avg_f1 = 0
		self.best_pred = ""
		self.n_pred = 0
		self.ratio = 0  # ratio of the number of inserted motifs in the best solution on the max number of inserted motifs for this RNA

	def get_label(self):
		if self.tool == "biorseo":
			if self.allow_pk:
				return f"{self.data_source}-{self.placement_method}-{self.func}"
			else:
				return f"{self.data_source}-{self.placement_method}-{self.func}-noPK"
		else:
			return self.tool

	def build_job_list(self):
		basename = self.parent_rna.basename
		fasta = outputDir+basename+".fa"

		# Things that require RNAsubopt calculations
		if self.tool in ["RNAsubopt", "RNA-MoIP (1by1)", "RNA-MoIP (chunk)"] or self.placement_method == "Jar3d":
			self.joblist.append(Job(command=["RNAsubopt", "-i", fasta, "--outfile="+ basename + ".subopt"], priority=1, how_many_in_parallel=1 if self.flat else 0, results = outputDir + "noPK/" +  basename + ".subopt", label=f"{basename} RNAsubopt"))
			self.joblist.append(Job(command=["mv", basename + ".subopt", outputDir + "noPK/"], priority=2, how_many_in_parallel=1 if self.flat else 0, results = outputDir + "noPK/" +  basename + ".subopt", label=f"{basename} mv"))

		# Find modules using Jar3d or BayesPairing:
		if self.placement_method == "Jar3d":
			self.joblist.append(Job(function=launch_JAR3D, args=[instance.seq_, basename],
								priority=3, how_many_in_parallel=1, results = outputDir + basename + ".bgsu_jar3d.csv", label=f"{basename} BGSU-Jar3d"))

		if self.placement_method == "ByP":
			if self.data_source == "DESC" :
				module_type_arg = "rna3dmotif"
			elif self.data_source == "RIN" :
				module_type_arg = "carnaval"
			else:
				module_type_arg = "3dmotifatlas"

			self.joblist.append(Job(function=launch_BayesPairing, args=[module_type_arg, instance.seq_, instance.header_, basename],
								how_many_in_parallel=1 if self.flat else -1, priority=3, results = outputDir + basename + f".{self.data_source.lower()}_byp.csv", label=f"{basename} {self.data_source}-ByP"))

			if module_type_arg != "carnaval":
				self.joblist.append(Job(function=launch_BayesPairing2, args=[module_type_arg, instance.seq_, instance.header_, basename],
								how_many_in_parallel=1 if self.flat else -1, priority=3, results = outputDir + basename + f".{self.data_source.lower()}_byp2.csv", label=f"{basename} {self.data_source}-ByP"))

		if self.tool == "biorseo":
			c = [ biorseoDir+"/bin/biorseo", "-s", fasta ]
			if self.placement_method == "D.P.":
				if self.data_source == "RIN" :
					results_file = outputDir+f"{'' if self.allow_pk else 'no'}PK/"+basename+f".biorseo_rin_raw_{self.func}"
					c += [ "-x", rinfolder]
				else:
					results_file = outputDir+f"{'' if self.allow_pk else 'no'}PK/"+basename+f".biorseo_desc_raw_{self.func}"
					c += [ "-d", descfolder]
			else:
				results_file = outputDir+f"{'' if self.allow_pk else 'no'}PK/"+basename+f".biorseo_{self.data_source.lower()}_{self.placement_method.lower()}_{self.func}"
				c += ["--bayespaircsv", outputDir+basename+f".{self.data_source.lower()}_{self.placement_method.lower()}.csv"]
			c += ["-o", results_file, "--func", self.func]
			if not self.allow_pk:
				c += ["-n"]
			self.joblist.append(Job(command=c, priority=4, timeout=3600, how_many_in_parallel=1 if self.flat else 3, results = results_file, label=f"{basename} {self.label}"))
		
		if self.tool == "RNA-MoIP (chunk)":
			self.joblist.append(Job(function=launch_RNAMoIP, args=[instance.seq_, instance.header_, basename, False],
								priority=3, how_many_in_parallel=1 if self.flat else 0, timeout=3600, results = outputDir + f"noPK/{basename}.moipc", label=f"{basename} {self.label}"))
		
		if self.tool == "RNA-MoIP (1by1)":
			self.joblist.append(Job(function=launch_RNAMoIP, args=[instance.seq_, instance.header_, basename, True],
								priority=3, how_many_in_parallel=1 if self.flat else 0, timeout=3600, results = outputDir + f"noPK/{basename}.moip", label=f"{basename} {self.label}"))
		
		if self.tool == "Biokop":
			self.joblist.append(Job(command=[biokopdir, "-n1", "-i", fasta, "-o", outputDir + f"PK/{basename}.biok"],
								priority=5, timeout=15000, how_many_in_parallel=1 if self.flat else 3, results = outputDir + f"PK/{basename}.biok", label=f"{basename} {self.label}"))

	def get_comp_times(self):
		s = ""
		for j in self.joblist:
			s += f"{j.comp_time:.1f} + "
		return s[:-3]


class RNA:
	def __init__(self, filename, header, seq, struct):
		self.seq_ = seq
		self.header_ = header
		self.true2d = struct
		self.basename = filename
		self.methods = []
		self.meth_idx = {}
		# subprocess.call(["rm", "-f", outputDir + self.basename + "*"])
		# subprocess.call(["rm", "-f", outputDir + "PK/" + self.basename + "*"])
		# subprocess.call(["rm", "-f", outputDir + "noPK/" + self.basename + "*"])

		if not path.isfile(outputDir + self.basename + ".fa"):
			rna = open(outputDir + self.basename + ".fa", "w")
			rna.write(">"+self.header_+'\n')
			rna.write(self.seq_+'\n')
			rna.close()

	def add_method_evaluation(self, *args, **kwargs):
		new_m = Method(*args, **kwargs)
		self.meth_idx[new_m.label] = len(self.methods)
		self.methods.append(new_m)

	def evaluate(self, verbose=False):
		if verbose:
			print("{:<24}{}".format("True 2D", self.true2d))

		for m in self.methods:
			if len(m.predictions):
				mccs = []
				f1s = []
				m.n_pred = len(m.predictions)

				# List unique solutions
				sec_structs = []
				for p in m.predictions:
					if not ')' in p: # ignore flat solutions
						m.n_pred -= 1
						continue
					ss = p.split('\t')[0].split(' ')[0]
					if ss not in sec_structs:
						sec_structs.append(p.split('\t')[0])
					else:
						m.n_pred -= 1
						continue
					f1s.append(f1_score(*compare_two_structures(self.true2d, p)))
					mccs.append(mattews_corr_coeff(*compare_two_structures(self.true2d, p)))

				if len(mccs):
					m.max_mcc = max(mccs)
					m.min_mcc = min(mccs)
					m.avg_mcc = sum(mccs)/float(len(mccs))
					m.best_pred = sec_structs[mccs.index(m.max_mcc)]
				if len(f1s):
					m.max_f1 = max(f1s)
					m.min_f1 = min(f1s)
					m.avg_f1 = sum(f1s)/float(len(f1s))

				for p,n in zip(m.predictions, m.ninsertions):
					if not ')' in p: # ignore linear structures
						continue
					# if several structures have the max_MCC
					if m.max_mcc == mattews_corr_coeff(*compare_two_structures(self.true2d, p)):
						# if one of them has a higher ratio, update the ratio
						if max(m.ninsertions) > 0 and float(n)/max(m.ninsertions) > m.ratio:
							m.ratio = float(n)/max(m.ninsertions)
							m.best_pred = p

				if verbose:
					print(f"{m.label:<21}\t{m.best_pred}\t{m.max_mcc:.2f}\t{m.n_pred}\t{m.get_comp_times()}")

	def get_results(self, label):
		return self.methods[self.meth_idx[label]]

	def load_biokop_results(self):
		if path.isfile(outputDir + "PK/" + self.basename + ".biok"):
			rna = open(outputDir + "PK/" + self.basename + ".biok", "r")
			lines = rna.readlines()
			rna.close()
			for i in range(1, len(lines)-1):
				ss = lines[i].split(' ')[0]
				if ss not in self.get_results("Biokop").predictions:
					self.get_results("Biokop").predictions.append(ss)

	def load_RNAsubopt_results(self):
		rna = open(outputDir + "noPK/" + self.basename + ".subopt", "r")
		lines = rna.readlines()
		rna.close()
		for i in range(2, len(lines)):
			ss = lines[i].split(' ')[0]
			if ss not in self.get_results("RNAsubopt").predictions:
				self.get_results("RNAsubopt").predictions.append(ss)

	def load_RNAMoIP_results(self):
		rna = open(outputDir + "noPK/" + self.basename + ".moipc", "r")
		lines = rna.readlines()
		rna.close()
		method = self.get_results("RNA-MoIP (chunk)")
		for i in range(2, len(lines)):
			method.predictions.append(lines[i].split('\t')[0])
			method.ninsertions.append(int(lines[i].split('\t')[1]))
			method.scores.append(float(lines[i].split('\t')[2][:-1]))
		rna = open(outputDir + "noPK/" + self.basename + ".moip", "r")
		lines = rna.readlines()
		rna.close()
		method = self.get_results("RNA-MoIP (1by1)")
		for i in range(2, len(lines)):
			method.predictions.append(lines[i].split('\t')[0])
			method.ninsertions.append(int(lines[i].split('\t')[1]))
			method.scores.append(float(lines[i].split('\t')[2][:-1]))

	def load_biorseo_results(self, filename, method):
		if path.isfile(filename):
			rna = open(filename, "r")
			lines = rna.readlines()
			rna.close()
			for i in range(2, len(lines)):
				ss = lines[i].split(' ')[0].split('\t')[0]
				# if ss not in method.predictions:
				method.predictions.append(ss)
				method.ninsertions.append(lines[i].count('+'))
		# else:
		#     print(filename, "not found !")

	def load_results(self, include_noPK=False):
		self.load_RNAsubopt_results()
		self.load_RNAMoIP_results()
		#self.load_biokop_results()
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_desc_raw_A", self.get_results("DESC-D.P.-A"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_desc_raw_B", self.get_results("DESC-D.P.-B"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_desc_byp_A", self.get_results("DESC-ByP-A"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_desc_byp_B", self.get_results("DESC-ByP-B"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_desc_byp_C", self.get_results("DESC-ByP-C"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_desc_byp_D", self.get_results("DESC-ByP-D"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_bgsu_byp_A", self.get_results("BGSU-ByP-A"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_bgsu_byp_B", self.get_results("BGSU-ByP-B"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_bgsu_byp_C", self.get_results("BGSU-ByP-C"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_bgsu_byp_D", self.get_results("BGSU-ByP-D"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_bgsu_jar3d_A", self.get_results("BGSU-Jar3d-A"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_bgsu_jar3d_B", self.get_results("BGSU-Jar3d-B"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_bgsu_jar3d_C", self.get_results("BGSU-Jar3d-C"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_bgsu_jar3d_D", self.get_results("BGSU-Jar3d-D"))

		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_rin_raw_A", self.get_results("RIN-D.P.-A"))
		self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_rin_raw_B", self.get_results("RIN-D.P.-B"))

		if include_noPK:
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_desc_raw_A", self.get_results("DESC-D.P.-A-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_desc_raw_B", self.get_results("DESC-D.P.-B-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_desc_byp_A", self.get_results("DESC-ByP-A-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_desc_byp_B", self.get_results("DESC-ByP-B-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_desc_byp_C", self.get_results("DESC-ByP-C-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_desc_byp_D", self.get_results("DESC-ByP-D-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_byp_A", self.get_results("BGSU-ByP-A-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_byp_B", self.get_results("BGSU-ByP-B-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_byp_C", self.get_results("BGSU-ByP-C-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_byp_D", self.get_results("BGSU-ByP-D-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_jar3d_A", self.get_results("BGSU-Jar3d-A-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_jar3d_B", self.get_results("BGSU-Jar3d-B-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_jar3d_C", self.get_results("BGSU-Jar3d-C-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_jar3d_D", self.get_results("BGSU-Jar3d-D-noPK"))

			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_rin_raw_A", self.get_results("RIN-D.P.-A-noPK"))
			self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_rin_raw_B", self.get_results("RIN-D.P.-B-noPK"))


	def has_complete_results(self, with_PK):
		if not with_PK:
			if not path.isfile(outputDir + "noPK/" + self.basename + ".subopt"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".moip"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".moipc"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_desc_raw_A"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_desc_raw_B"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_desc_byp_A"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_desc_byp_B"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_desc_byp_C"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_desc_byp_D"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_byp_A"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_byp_B"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_byp_C"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_byp_D"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_jar3d_A"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_jar3d_B"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_jar3d_C"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_bgsu_jar3d_D"): return False

			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_rin_raw_A"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_rin_raw_B"): return False

			return True
		else:
			#if not path.isfile(outputDir + "PK/" + self.basename + ".biok"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_desc_raw_A"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_desc_raw_B"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_desc_byp_A"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_desc_byp_B"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_desc_byp_C"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_desc_byp_D"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_bgsu_byp_A"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_bgsu_byp_B"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_bgsu_byp_C"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_bgsu_byp_D"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_bgsu_jar3d_A"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_bgsu_jar3d_B"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_bgsu_jar3d_C"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_bgsu_jar3d_D"): return False

			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_rin_raw_A"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_rin_raw_B"): return False

			return True

# ================= EXTRACTION OF STRUCTURES FROM FILES ===============================

if __name__ == '__main__':

	print("loading files...")
	RNAStrandContainer, RNAStrand_pk_counter = load_from_dbn(RNAStrandFile)
	PseudobaseContainer, Pseudobase_pk_counter = load_from_dbn(PseudobaseFile)
	StudycaseContainer, StudyCase_pk_counter = load_from_dbn(StudyCaseFile)

	for nt, number in ignored_nt_dict.items():
		print("ignored %d sequences because of char %c" % (number, nt))

	RNAStrand_tot = len(RNAStrandContainer)
	Pseudobase_tot = len(PseudobaseContainer)
	StudyCase_tot = len(StudycaseContainer)
	print("Loaded %d RNAs of length between 10 and 100 from RNA Strand. %d of them contain pseudoknots." % (RNAStrand_tot, RNAStrand_pk_counter))
	print("Loaded %d RNAs of length between 10 and 100 from Pseudobase. %d of them contain pseudoknots." % (Pseudobase_tot, Pseudobase_pk_counter))
	print("Loaded %d RNAs of length between 10 and 100 from study case. %d of them contain pseudoknots." % (StudyCase_tot, StudyCase_pk_counter))

	#================= PREDICTION OF STRUCTURES ===============================

	#define job list
	fulljoblist = []
	joblabel_list = []

	for instance in PseudobaseContainer + RNAStrandContainer:
		instance.add_method_evaluation(instance, "RNAsubopt")
		#instance.add_method_evaluation(instance, "Biokop")
		instance.add_method_evaluation(instance, "RNA-MoIP (1by1)")
		instance.add_method_evaluation(instance, "RNA-MoIP (chunk)")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="D.P.", obj_func="A")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="D.P.", obj_func="B")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="A")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="B")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="C")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="D")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="A")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="B")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="C")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="D")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="A")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="B")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="C")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="D")

		instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", placement_method="D.P.", obj_func="A")
		instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", placement_method="D.P.", obj_func="B")

		for method in instance.methods:
			for i in range(len(method.joblist)):
				j = method.joblist[i]
				if j.label not in joblabel_list: # look for a duplicate job (Jar3d, BayesPairing, RNAsubopt...)
					fulljoblist.append(j)
					joblabel_list.append(j.label)
				else:
					index = joblabel_list.index(j.label)
					method.joblist[i] = fulljoblist[index] # point to the previous occurrence

	for instance in RNAStrandContainer:
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="D.P.", obj_func="A", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="D.P.", obj_func="B", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="A", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="B", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="C", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="D", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="A", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="B", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="C", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="D", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="A", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="B", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="C", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="D", PK=False)

		instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", placement_method="D.P.", obj_func="A", PK=False)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", placement_method="D.P.", obj_func="B", PK=False)

		for method in instance.methods:
			for i in range(len(method.joblist)):
				j = method.joblist[i]
				if j.label not in joblabel_list: # look for a duplicate job (Jar3d, BayesPairing, RNAsubopt...)
					fulljoblist.append(j)
					joblabel_list.append(j.label)
				else:
					index = joblabel_list.index(j.label)
					method.joblist[i] = fulljoblist[index] # point to the previous occurrence

	
	for instance in StudycaseContainer: # We need to define these separately because we do not want concurrency, to measure proper run times.
		instance.add_method_evaluation(instance, "RNAsubopt", flat=True)
		#instance.add_method_evaluation(instance, "Biokop", flat=True)
		instance.add_method_evaluation(instance, "RNA-MoIP (1by1)", flat=True)
		instance.add_method_evaluation(instance, "RNA-MoIP (chunk)", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="D.P.", obj_func="A", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="D.P.", obj_func="B", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="A", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="B", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="C", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", placement_method="ByP", obj_func="D", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="A", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="B", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="C", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="ByP", obj_func="D", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="A", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="B", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="C", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="BGSU", placement_method="Jar3d", obj_func="D", flat=True)

		instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", placement_method="D.P.", obj_func="A", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", placement_method="D.P.", obj_func="B", flat=True)

		for method in instance.methods:
			for i in range(len(method.joblist)):
				j = method.joblist[i]
				if j.label not in joblabel_list: # look for a duplicate job (Jar3d, BayesPairing, RNAsubopt...)
					fulljoblist.append(j)
					joblabel_list.append(j.label)
				else:
					index = joblabel_list.index(j.label)
					method.joblist[i] = fulljoblist[index] # point to the previous occurrence
	

	# sort jobs in a tree structure
	"""
	jobs = {}
	jobcount = len(fulljoblist)
	for job in fulljoblist:
		if job.priority_ not in jobs.keys():
			jobs[job.priority_] = {}
		if job.nthreads not in jobs[job.priority_].keys():
			jobs[job.priority_][job.nthreads] = []
		jobs[job.priority_][job.nthreads].append(job)
	nprio = max(jobs.keys())
	# for each priority level
	for i in range(1,nprio+1):
		if i not in jobs.keys(): continue # ignore this priority level if no job available
		different_thread_numbers = [n for n in jobs[i].keys()]
		different_thread_numbers.sort()
		print("processing jobs of priority", i)
		# jobs should be processed 1 by 1, 2 by 2, or n by n depending on their definition
		for n in different_thread_numbers:
			bunch = jobs[i][n]
			if not len(bunch): continue # ignore if no jobs should be processed n by n
			print("using", n, "processes:")

			try :
				# execute jobs of priority i that should be processed n by n:
				p = MyPool(processes=n, maxtasksperchild=10)
				raw_results = p.map(execute_job, bunch)
				p.close()
				p.join()

				# extract computation times
				times = [ r[0] for r in raw_results ]
				for j, t in zip(bunch, times):
					j.comp_time = t

			except (subprocess.TimeoutExpired) :
				print("Skipping, took more than 3600s")
				pass
	"""


	# ================= Statistics ========================

	def get_RNAStrand_statistics(include_noPK=True):
		print("\nLoading RNA Strand results from files...")
		# load results in objects
		for instance in RNAStrandContainer:
			instance.load_results(include_noPK=True)
			instance.evaluate()

		RNAs_fully_predicted_noPK = [ x for x in RNAStrandContainer if x.has_complete_results(False)]

		x_noPK = [
			[ rna.get_results("RNAsubopt").max_mcc for rna in RNAStrandContainer if rna.get_results("RNAsubopt").n_pred],
			[ rna.get_results("RNA-MoIP (1by1)").max_mcc for rna in RNAStrandContainer if rna.get_results("RNA-MoIP (1by1)").n_pred],
			[ rna.get_results("RNA-MoIP (chunk)").max_mcc for rna in RNAStrandContainer if rna.get_results("RNA-MoIP (chunk)").n_pred],
			[ rna.get_results("DESC-D.P.-A-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-A-noPK").n_pred],
			[ rna.get_results("DESC-D.P.-B-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-B-noPK").n_pred],
			[ rna.get_results("DESC-ByP-A-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("DESC-ByP-A-noPK").n_pred],
			[ rna.get_results("DESC-ByP-B-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("DESC-ByP-B-noPK").n_pred],
			[ rna.get_results("DESC-ByP-C-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("DESC-ByP-C-noPK").n_pred],
			[ rna.get_results("DESC-ByP-D-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("DESC-ByP-D-noPK").n_pred],
			[ rna.get_results("BGSU-Jar3d-A-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-A-noPK").n_pred],
			[ rna.get_results("BGSU-Jar3d-B-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-B-noPK").n_pred],
			[ rna.get_results("BGSU-Jar3d-C-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-C-noPK").n_pred],
			[ rna.get_results("BGSU-Jar3d-D-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-D-noPK").n_pred],
			[ rna.get_results("BGSU-ByP-A-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-A-noPK").n_pred],
			[ rna.get_results("BGSU-ByP-B-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-B-noPK").n_pred],
			[ rna.get_results("BGSU-ByP-C-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-C-noPK").n_pred],
			[ rna.get_results("BGSU-ByP-D-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-D-noPK").n_pred],

			[ rna.get_results("RIN-D.P.-A-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-A-noPK").n_pred],
			[ rna.get_results("RIN-D.P.-B-noPK").max_mcc  for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-B-noPK").n_pred],
		]

		x_noPK_fully = [
			[ rna.get_results("RNAsubopt").max_mcc for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("RNA-MoIP (1by1)").max_mcc for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("RNA-MoIP (chunk)").max_mcc for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("DESC-D.P.-A-noPK").max_mcc for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("DESC-D.P.-B-noPK").max_mcc for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("DESC-ByP-A-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("DESC-ByP-B-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("DESC-ByP-C-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("DESC-ByP-D-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("BGSU-Jar3d-A-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("BGSU-Jar3d-B-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("BGSU-Jar3d-C-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("BGSU-Jar3d-D-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("BGSU-ByP-A-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("BGSU-ByP-B-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("BGSU-ByP-C-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("BGSU-ByP-D-noPK").max_mcc  for rna in RNAs_fully_predicted_noPK],

			[ rna.get_results("RIN-D.P.-A-noPK").max_mcc for rna in RNAs_fully_predicted_noPK],
			[ rna.get_results("RIN-D.P.-B-noPK").max_mcc for rna in RNAs_fully_predicted_noPK],
		]  # We ensure having the same number of RNAs in every sample by discarding the ones for which computations did not ended/succeeded.


		print()
		print("Without PK:")
		print("%s RNAsubopt predictions" % is_all(len(x_noPK[0]), RNAStrand_tot))
		print("%s RNA-MoIP 1 by 1 predictions" % is_all(len(x_noPK[1]), RNAStrand_tot))
		print("%s RNA-MoIP chunk predictions" % is_all(len(x_noPK[2]), RNAStrand_tot))
		print("%s biorseo + DESC + Patternmatch + f1A predictions" % is_all(len(x_noPK[3]), RNAStrand_tot))
		print("%s biorseo + DESC + Patternmatch + f1B predictions" % is_all(len(x_noPK[4]), RNAStrand_tot))
		print("%s biorseo + DESC + BayesPairing + f1A predictions" % is_all(len(x_noPK[5]), RNAStrand_tot))
		print("%s biorseo + DESC + BayesPairing + f1B predictions" % is_all(len(x_noPK[6]), RNAStrand_tot))
		print("%s biorseo + DESC + BayesPairing + f1C predictions" % is_all(len(x_noPK[7]), RNAStrand_tot))
		print("%s biorseo + DESC + BayesPairing + f1D predictions" % is_all(len(x_noPK[8]), RNAStrand_tot))
		print("%s biorseo + BGSU + JAR3D + f1A predictions" % is_all(len(x_noPK[9]), RNAStrand_tot))
		print("%s biorseo + BGSU + JAR3D + f1B predictions" % is_all(len(x_noPK[10]), RNAStrand_tot))
		print("%s biorseo + BGSU + JAR3D + f1C predictions" % is_all(len(x_noPK[11]), RNAStrand_tot))
		print("%s biorseo + BGSU + JAR3D + f1D predictions" % is_all(len(x_noPK[12]), RNAStrand_tot))
		print("%s biorseo + BGSU + BayesPairing + f1A predictions" % is_all(len(x_noPK[13]), RNAStrand_tot))
		print("%s biorseo + BGSU + BayesPairing + f1B predictions" % is_all(len(x_noPK[14]), RNAStrand_tot))
		print("%s biorseo + BGSU + BayesPairing + f1C predictions" % is_all(len(x_noPK[15]), RNAStrand_tot))
		print("%s biorseo + BGSU + BayesPairing + f1D predictions" % is_all(len(x_noPK[16]), RNAStrand_tot))

		print("%s biorseo + RIN + Patternmatch + f1A predictions" % is_all(len(x_noPK[17]), RNAStrand_tot))
		print("%s biorseo + RIN + Patternmatch + f1B predictions" % is_all(len(x_noPK[18]), RNAStrand_tot))

		print("==> %s ARN were predicted with all methods successful." % is_all(len(x_noPK_fully[0]), RNAStrand_tot) )

		# stat tests
		# Search if all methods are equal in positions with Friedman test:
		test = stats.friedmanchisquare(*x_noPK_fully)
		print("Friedman test without PK: H0 = 'The position parameter of all distributions is equal', p-value = ", test.pvalue)
		# ==> No they are not, but none does better, no need to test one further.
		test = stats.wilcoxon(x_noPK_fully[1], x_noPK_fully[3])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of RNA-MoIP and RawA are equal', p-value = ", test.pvalue)
		test = stats.wilcoxon(x_noPK_fully[1], x_noPK_fully[4])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of RNA-MoIP and RawB are equal', p-value = ", test.pvalue)


		RNAs_fully_predicted_PK = [ x for x in RNAStrandContainer if x.has_complete_results(True)]

		x_PK = [
			#[ rna.get_results("Biokop").max_mcc for rna in RNAStrandContainer if rna.get_results("Biokop").n_pred],
			[],
			[ rna.get_results("DESC-D.P.-A").max_mcc for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-A").n_pred],
			[ rna.get_results("DESC-D.P.-B").max_mcc for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-B").n_pred],
			[ rna.get_results("DESC-ByP-A").max_mcc  for rna in RNAStrandContainer if rna.get_results("DESC-ByP-A").n_pred],
			[ rna.get_results("DESC-ByP-B").max_mcc  for rna in RNAStrandContainer if rna.get_results("DESC-ByP-B").n_pred],
			[ rna.get_results("DESC-ByP-C").max_mcc  for rna in RNAStrandContainer if rna.get_results("DESC-ByP-C").n_pred],
			[ rna.get_results("DESC-ByP-D").max_mcc  for rna in RNAStrandContainer if rna.get_results("DESC-ByP-D").n_pred],
			[ rna.get_results("BGSU-Jar3d-A").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-A").n_pred],
			[ rna.get_results("BGSU-Jar3d-B").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-B").n_pred],
			[ rna.get_results("BGSU-Jar3d-C").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-C").n_pred],
			[ rna.get_results("BGSU-Jar3d-D").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-D").n_pred],
			[ rna.get_results("BGSU-ByP-A").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-A").n_pred],
			[ rna.get_results("BGSU-ByP-B").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-B").n_pred],
			[ rna.get_results("BGSU-ByP-C").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-C").n_pred],
			[ rna.get_results("BGSU-ByP-D").max_mcc  for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-D").n_pred],

			[ rna.get_results("RIN-D.P.-A").max_mcc for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-A").n_pred],
			[ rna.get_results("RIN-D.P.-B").max_mcc for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-B").n_pred],
		]

		# We ensure having the same number of RNAs in every sample by discarding the one for which computations did not ended/succeeded.
		x_PK_fully = [
			#[ rna.get_results("Biokop").max_mcc for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("DESC-D.P.-A").max_mcc for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("DESC-D.P.-B").max_mcc for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("DESC-ByP-A").max_mcc  for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("DESC-ByP-B").max_mcc  for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("DESC-ByP-C").max_mcc  for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("DESC-ByP-D").max_mcc  for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("BGSU-Jar3d-A").max_mcc  for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("BGSU-Jar3d-B").max_mcc  for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("BGSU-Jar3d-C").max_mcc  for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("BGSU-Jar3d-D").max_mcc  for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("BGSU-ByP-A").max_mcc  for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("BGSU-ByP-B").max_mcc  for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("BGSU-ByP-C").max_mcc  for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("BGSU-ByP-D").max_mcc  for rna in RNAs_fully_predicted_PK],

			[ rna.get_results("RIN-D.P.-A").max_mcc for rna in RNAs_fully_predicted_PK],
			[ rna.get_results("RIN-D.P.-B").max_mcc for rna in RNAs_fully_predicted_PK],
		]



		print()
		print("With PK:")
		#print("%s Biokop predictions" % is_all(len(x_PK[0]), RNAStrand_tot))
		print("%s biorseo + DESC + Patternmatch + f1A predictions" % is_all(len(x_PK[1]), RNAStrand_tot))
		print("%s biorseo + DESC + Patternmatch + f1B predictions" % is_all(len(x_PK[2]), RNAStrand_tot))
		print("%s biorseo + DESC + BayesPairing + f1A predictions" % is_all(len(x_PK[3]), RNAStrand_tot))
		print("%s biorseo + DESC + BayesPairing + f1B predictions" % is_all(len(x_PK[4]), RNAStrand_tot))
		print("%s biorseo + DESC + BayesPairing + f1C predictions" % is_all(len(x_PK[5]), RNAStrand_tot))
		print("%s biorseo + DESC + BayesPairing + f1D predictions" % is_all(len(x_PK[6]), RNAStrand_tot))
		print("%s biorseo + BGSU + JAR3D + f1A predictions" % is_all(len(x_PK[7]), RNAStrand_tot))
		print("%s biorseo + BGSU + JAR3D + f1B predictions" % is_all(len(x_PK[8]), RNAStrand_tot))
		print("%s biorseo + BGSU + JAR3D + f1C predictions" % is_all(len(x_PK[9]), RNAStrand_tot))
		print("%s biorseo + BGSU + JAR3D + f1D predictions" % is_all(len(x_PK[10]), RNAStrand_tot))
		print("%s biorseo + BGSU + BayesPairing + f1A predictions" % is_all(len(x_PK[11]), RNAStrand_tot))
		print("%s biorseo + BGSU + BayesPairing + f1B predictions" % is_all(len(x_PK[12]), RNAStrand_tot))
		print("%s biorseo + BGSU + BayesPairing + f1C predictions" % is_all(len(x_PK[13]), RNAStrand_tot))
		print("%s biorseo + BGSU + BayesPairing + f1D predictions" % is_all(len(x_PK[14]), RNAStrand_tot))

		print("%s biorseo + RIN + Patternmatch + f1A predictions" % is_all(len(x_PK[15]), RNAStrand_tot))
		print("%s biorseo + RIN + Patternmatch + f1B predictions" % is_all(len(x_PK[16]), RNAStrand_tot))

		print("==> %s ARN were predicted with all methods successful." % is_all(len(x_PK_fully[0]), RNAStrand_tot) )

		# stat tests
		# First, search if all methods are equal in positions with Friedman test:
		test = stats.friedmanchisquare(*x_PK_fully)
		print("Friedman test with PK: H0 = 'The position parameter of all distributions is equal', p-value = ", test.pvalue)
		# it looks like some methods do better. Let's test the difference:
		test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[1])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and RawA are equal', p-value = ", test.pvalue)
		test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[2])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and RawB are equal', p-value = ", test.pvalue)
		test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[7])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dA are equal', p-value = ", test.pvalue)
		test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[8])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dB are equal', p-value = ", test.pvalue)
		test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[9])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dC are equal', p-value = ", test.pvalue)
		test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[10])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dD are equal', p-value = ", test.pvalue)
		return x_noPK_fully, x_PK_fully

	def get_Pseudobase_statistics():
		# load results in objects
		print("\nLoading Pseudobase results from files...")
		for instance in PseudobaseContainer:
			instance.load_results()
			instance.evaluate()

		RNAs_fully_predicted_Pseudobase = [ x for x in PseudobaseContainer if x.has_complete_results(True)]

		x_pseudobase = [
			#[ rna.get_results("Biokop").max_mcc for rna in PseudobaseContainer if rna.get_results("Biokop").n_pred],
			[],
			[ rna.get_results("RNAsubopt").max_mcc for rna in PseudobaseContainer if rna.get_results("RNAsubopt").n_pred],
			[ rna.get_results("RNA-MoIP (1by1)").max_mcc for rna in PseudobaseContainer if rna.get_results("RNA-MoIP (1by1)").n_pred],
			[ rna.get_results("RNA-MoIP (chunk)").max_mcc for rna in PseudobaseContainer if rna.get_results("RNA-MoIP (chunk)").n_pred],
			[ rna.get_results("DESC-D.P.-A").max_mcc for rna in PseudobaseContainer if rna.get_results("DESC-D.P.-A").n_pred],
			[ rna.get_results("DESC-D.P.-B").max_mcc for rna in PseudobaseContainer if rna.get_results("DESC-D.P.-B").n_pred],
			[ rna.get_results("DESC-ByP-A").max_mcc  for rna in PseudobaseContainer if rna.get_results("DESC-ByP-A").n_pred],
			[ rna.get_results("DESC-ByP-B").max_mcc  for rna in PseudobaseContainer if rna.get_results("DESC-ByP-B").n_pred],
			[ rna.get_results("DESC-ByP-C").max_mcc  for rna in PseudobaseContainer if rna.get_results("DESC-ByP-C").n_pred],
			[ rna.get_results("DESC-ByP-D").max_mcc  for rna in PseudobaseContainer if rna.get_results("DESC-ByP-D").n_pred],
			[ rna.get_results("BGSU-Jar3d-A").max_mcc  for rna in PseudobaseContainer if rna.get_results("BGSU-Jar3d-A").n_pred],
			[ rna.get_results("BGSU-Jar3d-B").max_mcc  for rna in PseudobaseContainer if rna.get_results("BGSU-Jar3d-B").n_pred],
			[ rna.get_results("BGSU-Jar3d-C").max_mcc  for rna in PseudobaseContainer if rna.get_results("BGSU-Jar3d-C").n_pred],
			[ rna.get_results("BGSU-Jar3d-D").max_mcc  for rna in PseudobaseContainer if rna.get_results("BGSU-Jar3d-D").n_pred],
			[ rna.get_results("BGSU-ByP-A").max_mcc  for rna in PseudobaseContainer if rna.get_results("BGSU-ByP-A").n_pred],
			[ rna.get_results("BGSU-ByP-B").max_mcc  for rna in PseudobaseContainer if rna.get_results("BGSU-ByP-B").n_pred],
			[ rna.get_results("BGSU-ByP-C").max_mcc  for rna in PseudobaseContainer if rna.get_results("BGSU-ByP-C").n_pred],
			[ rna.get_results("BGSU-ByP-D").max_mcc  for rna in PseudobaseContainer if rna.get_results("BGSU-ByP-D").n_pred],

			[ rna.get_results("RIN-D.P.-A").max_mcc for rna in PseudobaseContainer if rna.get_results("RIN-D.P.-A").n_pred],
			[ rna.get_results("RIN-D.P.-B").max_mcc for rna in PseudobaseContainer if rna.get_results("RIN-D.P.-B").n_pred],
		]

		# We ensure having the same number of RNAs in every sample by discarding the one for which computations did not ended/succeeded.
		x_pseudobase_fully = [
			#[ rna.get_results("Biokop").max_mcc for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("RNAsubopt").max_mcc for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("RNA-MoIP (1by1)").max_mcc for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("RNA-MoIP (chunk)").max_mcc for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("DESC-D.P.-A").max_mcc for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("DESC-D.P.-B").max_mcc for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("DESC-ByP-A").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("DESC-ByP-B").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("DESC-ByP-C").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("DESC-ByP-D").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("BGSU-Jar3d-A").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("BGSU-Jar3d-B").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("BGSU-Jar3d-C").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("BGSU-Jar3d-D").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("BGSU-ByP-A").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("BGSU-ByP-B").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("BGSU-ByP-C").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("BGSU-ByP-D").max_mcc  for rna in RNAs_fully_predicted_Pseudobase],

			[ rna.get_results("RIN-D.P.-A").max_mcc for rna in RNAs_fully_predicted_Pseudobase],
			[ rna.get_results("RIN-D.P.-B").max_mcc for rna in RNAs_fully_predicted_Pseudobase],
		]


		print()
		print("With PK:")
		#print("%s Biokop predictions" % is_all(len(x_pseudobase[0]), Pseudobase_tot))
		print("%s RNAsubopt predictions" % is_all(len(x_pseudobase[1]), Pseudobase_tot))
		print("%s RNA-MoIP 1 by 1 predictions" % is_all(len(x_pseudobase[2]), Pseudobase_tot))
		print("%s RNA-MoIP chunk predictions" % is_all(len(x_pseudobase[3]), Pseudobase_tot))
		print("%s biorseo + DESC + Patternmatch + f1A predictions" % is_all(len(x_pseudobase[4]), Pseudobase_tot))
		print("%s biorseo + DESC + Patternmatch + f1B predictions" % is_all(len(x_pseudobase[5]), Pseudobase_tot))
		print("%s biorseo + DESC + BayesPairing + f1A predictions" % is_all(len(x_pseudobase[6]), Pseudobase_tot))
		print("%s biorseo + DESC + BayesPairing + f1B predictions" % is_all(len(x_pseudobase[7]), Pseudobase_tot))
		print("%s biorseo + DESC + BayesPairing + f1C predictions" % is_all(len(x_pseudobase[8]), Pseudobase_tot))
		print("%s biorseo + DESC + BayesPairing + f1D predictions" % is_all(len(x_pseudobase[9]), Pseudobase_tot))
		print("%s biorseo + BGSU + JAR3D + f1A predictions" % is_all(len(x_pseudobase[10]), Pseudobase_tot))
		print("%s biorseo + BGSU + JAR3D + f1B predictions" % is_all(len(x_pseudobase[11]), Pseudobase_tot))
		print("%s biorseo + BGSU + JAR3D + f1C predictions" % is_all(len(x_pseudobase[12]), Pseudobase_tot))
		print("%s biorseo + BGSU + JAR3D + f1D predictions" % is_all(len(x_pseudobase[13]), Pseudobase_tot))
		print("%s biorseo + BGSU + BayesPairing + f1A predictions" % is_all(len(x_pseudobase[14]), Pseudobase_tot))
		print("%s biorseo + BGSU + BayesPairing + f1B predictions" % is_all(len(x_pseudobase[15]), Pseudobase_tot))
		print("%s biorseo + BGSU + BayesPairing + f1C predictions" % is_all(len(x_pseudobase[16]), Pseudobase_tot))
		print("%s biorseo + BGSU + BayesPairing + f1D predictions" % is_all(len(x_pseudobase[17]), Pseudobase_tot))

		print("%s biorseo + RIN + Patternmatch + f1A predictions" % is_all(len(x_pseudobase[18]), Pseudobase_tot))
		print("%s biorseo + RIN + Patternmatch + f1B predictions" % is_all(len(x_pseudobase[19]), Pseudobase_tot))

		print("==> %s ARN were predicted with all methods successful." % is_all(len(x_pseudobase_fully[0]), Pseudobase_tot) )

		# stat tests
		# First, search if all methods are equal in positions with Friedman test:
		test = stats.friedmanchisquare(*x_pseudobase_fully)
		print("Friedman test with PK: H0 = 'The position parameter of all distributions is equal', p-value = ", test.pvalue)
		# it looks like some methods do better. Let's test the difference:
		test = stats.wilcoxon(x_pseudobase_fully[0], x_pseudobase_fully[4])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and RawA are equal', p-value = ", test.pvalue)
		test = stats.wilcoxon(x_pseudobase_fully[0], x_pseudobase_fully[5])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and RawB are equal', p-value = ", test.pvalue)
		test = stats.wilcoxon(x_pseudobase_fully[0], x_pseudobase_fully[10])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dA are equal', p-value = ", test.pvalue)
		test = stats.wilcoxon(x_pseudobase_fully[0], x_pseudobase_fully[11])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dB are equal', p-value = ", test.pvalue)
		test = stats.wilcoxon(x_pseudobase_fully[0], x_pseudobase_fully[12])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dC are equal', p-value = ", test.pvalue)
		test = stats.wilcoxon(x_pseudobase_fully[0], x_pseudobase_fully[13])
		print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and Jar3dD are equal', p-value = ", test.pvalue)
		return x_pseudobase_fully

	
	def print_StudyCase_results():
		print("\nLoading study case results from files...")

		# load results in objects
		for instance in StudycaseContainer:
			instance.load_results()
			instance.evaluate(verbose=True)
	

	x_noPK_fully, x_PK_fully = get_RNAStrand_statistics()
	x_pseudobase_fully = get_Pseudobase_statistics()
	print_StudyCase_results()
	"""
	nbetter = [
			len([ rna for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("DESC-ByP-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("DESC-ByP-B").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("DESC-ByP-C").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("DESC-ByP-D").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-B").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-C").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-D").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-B").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-C").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-D").max_mcc > rna.get_results("Biokop").max_mcc ]),

			len([ rna for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
		]
	print("Number of RNAs that are better predicted by BiORSEO than Biokop on RNAStrand:")
	print(nbetter, "/", len(RNAStrandContainer))

	nbetter = [
			len([ rna for rna in PseudobaseContainer if rna.get_results("DESC-D.P.-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("DESC-D.P.-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("DESC-ByP-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("DESC-ByP-B").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("DESC-ByP-C").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("DESC-ByP-D").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("BGSU-ByP-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("BGSU-ByP-B").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("BGSU-ByP-C").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("BGSU-ByP-D").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("BGSU-Jar3d-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("BGSU-Jar3d-B").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("BGSU-Jar3d-C").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("BGSU-Jar3d-D").max_mcc > rna.get_results("Biokop").max_mcc ]),

			len([ rna for rna in PseudobaseContainer if rna.get_results("RIN-D.P.-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
			len([ rna for rna in PseudobaseContainer if rna.get_results("RIN-D.P.-A").max_mcc > rna.get_results("Biokop").max_mcc ]),
		]
	print("Number of RNAs that are better predicted by BiORSEO than Biokop on Pseudobase:")
	print(nbetter, "/", len(PseudobaseContainer))
	"""

	# ================= PLOTS OF RESULTS =======================================
	colors = [
		'#911eb4', #purple
		'#000075', #navy
		'#ffe119', '#ffe119', # yellow
		'#e6194B', '#e6194B', #red
		'#3cb44b', '#3cb44b', '#3cb44b', '#3cb44b', #green
		'#4363d8', '#4363d8', '#4363d8', '#4363d8', #blue
		'#3cb44b', '#3cb44b', '#3cb44b', '#3cb44b' # green
	]

	def plot_best_MCCs(x_noPK_fully, x_PK_fully, x_pseudobase_fully):

		#"Biokop\n", 
		labels = [
			"RNA\nsubopt", "RNA-\nMoIP\n1by1", "RNA-\nMoIP\nchunk",
			"$f_{1A}$", "$f_{1B}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$",
			"$f_{1A}$", "$f_{1B}$",
		]

		fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,3.5), dpi=80)
		fig.suptitle(" \n ")
		fig.subplots_adjust(left=0.1, right=0.99, top=0.83, bottom=0.02)

		for ax in axes:
			ax.set_ylim((0.5, 1.01))
			ax.set_xlim((-1,19))
			yticks = [ i/10 for i in range(5,11) ]
			ax.set_yticks(yticks)
			for y in yticks:
				ax.axhline(y=y, color="grey", linestyle="--", linewidth=1)
			ax.tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)
			ax.set_xticks([ 1.0+i for i in range(18)])
		axes[0].tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
		axes[0].set_xticklabels(labels)
		for i, tick in enumerate(axes[0].xaxis.get_major_ticks()):
			if i<4: # Reduce size of Biokop, RNAsubopt and RNA-MoIP labels to stay readable
				tick.label2.set_fontsize(10)
			else:
				tick.label2.set_fontsize(12)

		# Line 1 : no Pseudoknots
		# bplot = axes[0].boxplot([[]] + x_noPK_fully, vert=True, patch_artist=True, notch=False, whis=[3,97], medianprops=dict(color="white"))
		# for patch, color in zip(bplot['boxes'], colors):
		#     patch.set_facecolor(color)
		#xpos = [ x for x in range(2,19) ]
		xpos = [ x for x in range(len(x_noPK_fully)) ]
		#print(len(x_noPK_fully), len(xpos))
		vplot = axes[0].violinplot(x_noPK_fully, showmeans=False, showmedians=False, showextrema=False, points=len(x_noPK_fully[0]), positions=xpos)
		axes[0].set_xticks(xpos)
		for patch, color in zip(vplot['bodies'], colors[1:]):
			patch.set_facecolor(color)
			patch.set_edgecolor(color)
			patch.set_alpha(0.5)
		quartile1, medians, quartile3 = np.percentile(x_noPK_fully, [25, 50, 75], axis=1)
		axes[0].scatter(xpos, medians, marker='o', color='k', s=30, zorder=3)
		axes[0].vlines(xpos, quartile1, quartile3, color='k', linestyle='-', lw=1)
		for x, y1, y2 in zip(xpos, quartile1, quartile3):
			bar1 = Line2D([x-0.1, x+0.1], [y1, y1], color="k", lw=1)
			bar2 = Line2D([x-0.1, x+0.1], [y2, y2], color="k", lw=1)
			axes[0].add_line(bar1)
			axes[0].add_line(bar2)
		axes[0].set_ylabel("(A)\nmax MCC\n(%d RNAs)" % (len(x_noPK_fully[0])), fontsize=12)

		# Line 2 : Pseudoknots supported
		# bplot = axes[1].boxplot([x_PK_fully[0],[],[]] + x_PK_fully[1:], vert=True, patch_artist=True, notch=False, whis=[3,97], medianprops=dict(color="white"))
		# for patch, color in zip(bplot['boxes'], colors):
		#     patch.set_facecolor(color)
		#xpos = [1] + [ i for i in range(5,19) ]
		xpos = [ i for i in range(3,len(x_PK_fully)+3) ]
		#print(len(x_PK_fully), len(xpos))
		vplot = axes[1].violinplot(x_PK_fully, showmeans=False, showmedians=False, showextrema=False, points=len(x_PK_fully[0]), positions=xpos)
		#print(len(vplot['bodies']), len(xpos))
		for patch, color in zip(vplot['bodies'], colors[4:]):
			patch.set_facecolor(color)
			patch.set_edgecolor(color)
			patch.set_alpha(0.5)
		quartile1, medians, quartile3 = np.percentile(x_PK_fully, [25, 50, 75], axis=1)
		print(len(medians), len(xpos))
		axes[1].scatter(xpos, medians, marker='o', color='k', s=30, zorder=3)
		axes[1].vlines(xpos, quartile1, quartile3, color='k', linestyle='-', lw=1)
		for x, y1, y2 in zip(xpos, quartile1, quartile3):
			bar1 = Line2D([x-0.1, x+0.1], [y1, y1], color="k", lw=1)
			bar2 = Line2D([x-0.1, x+0.1], [y2, y2], color="k", lw=1)
			axes[1].add_line(bar1)
			axes[1].add_line(bar2)
		axes[1].set_ylabel("(B)\nmax MCC\n(%d RNAs)" % (len(x_PK_fully[0])), fontsize=12)

		# Line 3 : all methods on pseudoknotted dataset
		# bplot = axes[2].boxplot(x_pseudobase_fully, vert=True, patch_artist=True, notch=False, whis=[3,97], medianprops=dict(color="white"))
		# for patch, color in zip(bplot['boxes'], colors):
		#     patch.set_facecolor(color)
		#xpos = [ x for x in range(1,19) ]
		xpos = [ x for x in range(len(x_pseudobase_fully)) ]
		#print(len(x_pseudobase_fully), len(xpos))
		vplot = axes[2].violinplot(x_pseudobase_fully, showmeans=False, showmedians=False, showextrema=False, points=len(x_pseudobase_fully[0]), positions=xpos)
		for patch, color in zip(vplot['bodies'], colors[1:]):
			patch.set_facecolor(color)
			patch.set_edgecolor(color)
			patch.set_alpha(0.5)
		quartile1, medians, quartile3 = np.percentile(x_pseudobase_fully, [25, 50, 75], axis=1)
		axes[2].scatter(xpos, medians, marker='o', color='k', s=30, zorder=3)
		axes[2].vlines(xpos, quartile1, quartile3, color='k', linestyle='-', lw=1)
		for x, y1, y2 in zip(xpos, quartile1, quartile3):
			bar1 = Line2D([x-0.1, x+0.1], [y1, y1], color="k", lw=1)
			bar2 = Line2D([x-0.1, x+0.1], [y2, y2], color="k", lw=1)
			axes[2].add_line(bar1)
			axes[2].add_line(bar2)
		axes[2].set_ylabel("(C)\nmax MCC\n(%d RNAs)" % (len(x_pseudobase_fully[0])), fontsize=12)

	def plot_more_info():
		# ======= number of solutions, insertion ratio, etc ========================

		n = [
			#[ rna.get_results("Biokop").n_pred for rna in RNAStrandContainer if rna.get_results("Biokop").n_pred ],
			[ rna.get_results("RNAsubopt").n_pred for rna in RNAStrandContainer if rna.get_results("RNAsubopt").n_pred ],
			[ rna.get_results("RNA-MoIP (1by1)").n_pred for rna in RNAStrandContainer if rna.get_results("RNA-MoIP (1by1)").n_pred ],
			[ rna.get_results("RNA-MoIP (chunk)").n_pred for rna in RNAStrandContainer if rna.get_results("RNA-MoIP (chunk)").n_pred ],
			[ rna.get_results("DESC-D.P.-A").n_pred for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-A").n_pred ],
			[ rna.get_results("DESC-D.P.-B").n_pred for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-B").n_pred ],
			[ rna.get_results("DESC-ByP-A").n_pred for rna in RNAStrandContainer if rna.get_results("DESC-ByP-A").n_pred ],
			[ rna.get_results("DESC-ByP-B").n_pred for rna in RNAStrandContainer if rna.get_results("DESC-ByP-B").n_pred ],
			[ rna.get_results("DESC-ByP-C").n_pred for rna in RNAStrandContainer if rna.get_results("DESC-ByP-C").n_pred ],
			[ rna.get_results("DESC-ByP-D").n_pred for rna in RNAStrandContainer if rna.get_results("DESC-ByP-D").n_pred ],
			[ rna.get_results("BGSU-Jar3d-A").n_pred for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-A").n_pred ],
			[ rna.get_results("BGSU-Jar3d-B").n_pred for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-B").n_pred ],
			[ rna.get_results("BGSU-Jar3d-C").n_pred for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-C").n_pred ],
			[ rna.get_results("BGSU-Jar3d-D").n_pred for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-D").n_pred ],
			[ rna.get_results("BGSU-ByP-A").n_pred for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-A").n_pred ],
			[ rna.get_results("BGSU-ByP-B").n_pred for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-B").n_pred ],
			[ rna.get_results("BGSU-ByP-C").n_pred for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-C").n_pred ],
			[ rna.get_results("BGSU-ByP-D").n_pred for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-D").n_pred ],

			[ rna.get_results("RIN-D.P.-A").n_pred for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-A").n_pred ],
			[ rna.get_results("RIN-D.P.-B").n_pred for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-B").n_pred ],
		]

		r = [
			[ rna.get_results("RNA-MoIP (1by1)").ratio for rna in RNAStrandContainer if rna.get_results("RNA-MoIP (1by1)").n_pred > 1 ],
			[ rna.get_results("DESC-D.P.-A").ratio for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-A").n_pred > 1 ],
			[ rna.get_results("DESC-D.P.-B").ratio for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-B").n_pred > 1 ],
			[ rna.get_results("DESC-ByP-A").ratio for rna in RNAStrandContainer if rna.get_results("DESC-ByP-A").n_pred > 1 ],
			[ rna.get_results("DESC-ByP-B").ratio for rna in RNAStrandContainer if rna.get_results("DESC-ByP-B").n_pred > 1 ],
			[ rna.get_results("DESC-ByP-C").ratio for rna in RNAStrandContainer if rna.get_results("DESC-ByP-C").n_pred > 1 ],
			[ rna.get_results("DESC-ByP-D").ratio for rna in RNAStrandContainer if rna.get_results("DESC-ByP-D").n_pred > 1 ],
			[ rna.get_results("BGSU-Jar3d-A").ratio for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-A").n_pred > 1 ],
			[ rna.get_results("BGSU-Jar3d-B").ratio for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-B").n_pred > 1 ],
			[ rna.get_results("BGSU-Jar3d-C").ratio for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-C").n_pred > 1 ],
			[ rna.get_results("BGSU-Jar3d-D").ratio for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-D").n_pred > 1 ],
			[ rna.get_results("BGSU-ByP-A").ratio for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-A").n_pred > 1 ],
			[ rna.get_results("BGSU-ByP-B").ratio for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-B").n_pred > 1 ],
			[ rna.get_results("BGSU-ByP-C").ratio for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-C").n_pred > 1 ],
			[ rna.get_results("BGSU-ByP-D").ratio for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-D").n_pred > 1 ],

			[ rna.get_results("RIN-D.P.-A").ratio for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-A").n_pred > 1 ],
			[ rna.get_results("RIN-D.P.-B").ratio for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-B").n_pred > 1 ],
		]

		max_i = [
			[ max(rna.get_results("RNA-MoIP (1by1)").ninsertions) for rna in RNAStrandContainer if rna.get_results("RNA-MoIP (1by1)").n_pred ],
			[ max(rna.get_results("RNA-MoIP (chunk)").ninsertions) for rna in RNAStrandContainer if rna.get_results("RNA-MoIP (chunk)").n_pred ],
			[ max(rna.get_results("DESC-D.P.-A").ninsertions) for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-A").n_pred ],
			[ max(rna.get_results("DESC-D.P.-B").ninsertions) for rna in RNAStrandContainer if rna.get_results("DESC-D.P.-B").n_pred ],
			[ max(rna.get_results("DESC-ByP-A").ninsertions) for rna in RNAStrandContainer if rna.get_results("DESC-ByP-A").n_pred ],
			[ max(rna.get_results("DESC-ByP-B").ninsertions) for rna in RNAStrandContainer if rna.get_results("DESC-ByP-B").n_pred ],
			[ max(rna.get_results("DESC-ByP-C").ninsertions) for rna in RNAStrandContainer if rna.get_results("DESC-ByP-C").n_pred ],
			[ max(rna.get_results("DESC-ByP-D").ninsertions) for rna in RNAStrandContainer if rna.get_results("DESC-ByP-D").n_pred ],
			[ max(rna.get_results("BGSU-Jar3d-A").ninsertions) for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-A").n_pred ],
			[ max(rna.get_results("BGSU-Jar3d-B").ninsertions) for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-B").n_pred ],
			[ max(rna.get_results("BGSU-Jar3d-C").ninsertions) for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-C").n_pred ],
			[ max(rna.get_results("BGSU-Jar3d-D").ninsertions )for rna in RNAStrandContainer if rna.get_results("BGSU-Jar3d-D").n_pred ],
			[ max(rna.get_results("BGSU-ByP-A").ninsertions) for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-A").n_pred ],
			[ max(rna.get_results("BGSU-ByP-B").ninsertions) for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-B").n_pred ],
			[ max(rna.get_results("BGSU-ByP-C").ninsertions) for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-C").n_pred ],
			[ max(rna.get_results("BGSU-ByP-D").ninsertions) for rna in RNAStrandContainer if rna.get_results("BGSU-ByP-D").n_pred ],

			[ max(rna.get_results("RIN-D.P.-A").ninsertions) for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-A").n_pred ],
			[ max(rna.get_results("RIN-D.P.-B").ninsertions) for rna in RNAStrandContainer if rna.get_results("RIN-D.P.-B").n_pred ],
		]

		# Figure : number of solutions
		plt.figure(figsize=(6,1.5), dpi=80)
		plt.suptitle(" \n ")
		plt.subplots_adjust(left=0.05, right=0.99, top=0.5, bottom=0.05)
		labels = [
			#"Biokop",
			"RNAsubopt","RNA-MoIP\n1by1", "RNA-MoIP\nchunk",
			"$f_{1A}$", "$f_{1B}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$",
			"$f_{1A}$", "$f_{1B}$"
		]
		plt.xlim((0,19))
		plt.tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)
		plt.xticks([ 1.0+i for i in range(18)], labels)
		plt.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
		for i, tick in enumerate(plt.gca().xaxis.get_major_ticks()):
			if i<4: # Reduce size of RNA-MoIP labels to stay readable
				tick.label2.set_fontsize(8)
				tick.label2.set_rotation(90)
			else:
				tick.label2.set_fontsize(12)
		#xpos = [ x for x in range(1,19) ]
		xpos = [ x for x in range(len(n)) ]
		#print(len(n), len(xpos))
		plt.yticks([ 20*x for x in range(3) ])
		plt.ylim((0,40))
		for y in [ 10*x for x in range(8) ]:
			plt.axhline(y=y, color="grey", linestyle="-", linewidth=0.5)
		plt.axhline(y=1, color="grey", linestyle="-", linewidth=0.5)
		vplot = plt.violinplot(n, showmeans=False, showmedians=False, showextrema=False, points=len(n[0]), positions=xpos)
		for patch, color in zip(vplot['bodies'], colors):
			patch.set_facecolor(color)
			patch.set_edgecolor(color)
			patch.set_alpha(0.5)




		fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6,2.5), dpi=80)
		fig.suptitle(" \n ")
		fig.subplots_adjust(left=0.09, right=0.99, top=0.7, bottom=0.05)
		labels = [
			"RNA-MoIP\n1by1", "RNA-MoIP\nchunk",
			"$f_{1A}$", "$f_{1B}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$"
		]
		for ax in axes:
			ax.set_xlim((0,17))
			ax.tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)
			ax.set_xticks([ 1.0+i for i in range(16)])
		axes[0].tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
		axes[0].set_xticklabels(labels)
		for i, tick in enumerate(axes[0].xaxis.get_major_ticks()):
			if i<2: # Reduce size of RNA-MoIP labels to stay readable
				tick.label2.set_fontsize(9)
				tick.label2.set_rotation(90)
			else:
				tick.label2.set_fontsize(12)

		# Figure : max inserted
		#xpos = [ x for x in range(1,17) ]
		xpos = [ x for x in range(18) ]
		#print(len(max_i), len(xpos))
		axes[0].set_yticks([ 5*x for x in range(3) ])
		for y in [ 2*x for x in range(7) ]:
			axes[0].axhline(y=y, color="grey", linestyle="-", linewidth=0.5)
			# axes[0].axhline(y=y-1, color="lightgray", linestyle="-", linewidth=0.5)
		vplot = axes[0].violinplot(max_i, showmeans=False, showmedians=False, showextrema=False, points=len(max_i[0]), positions=xpos)
		for patch, color in zip(vplot['bodies'], colors[2:]):
			patch.set_facecolor(color)
			patch.set_edgecolor(color)
			patch.set_alpha(0.5)
		axes[0].set_ylabel("(A)", fontsize=12)

		# Figure : insertion ratio
		#xpos = [1] + [ x for x in range(3,17) ]
		xpos = [ x for x in range(len(r)) ]
		#print(len(r), len(xpos))
		axes[1].set_ylim((-0.01, 1.01))
		yticks = [ 0, 0.5, 1.0 ]
		axes[1].set_yticks(yticks)
		for y in yticks:
			axes[1].axhline(y=y, color="grey", linestyle="-", linewidth=0.5)
		vplot = axes[1].violinplot(r, showmeans=False, showmedians=False, showextrema=False, points=len(r[0]), positions=xpos)
		for patch, color in zip(vplot['bodies'], [colors[2]] + colors[4:]):
			patch.set_facecolor(color)
			patch.set_edgecolor(color)
			patch.set_alpha(0.5)
		for i,x in enumerate(xpos):
			axes[1].annotate(str(len(r[i])), (x-0.25, 0.05), fontsize=8)
		axes[1].set_ylabel("(B)", fontsize=12)

	def compare_subopt_MoIP():
		# ================== MCC performance ====================================

		plt.figure(figsize=(10,4), dpi=80)
		RNAStrandContainer.sort(key=lambda x: x.get_results("RNA-MoIP (chunk)").max_mcc)

		x = [
			[ rna.get_results("RNA-MoIP (chunk)").max_mcc for rna in RNAStrandContainer ],
			[ rna.get_results("RNA-MoIP (1by1)").max_mcc for rna in RNAStrandContainer ],
			[ rna.get_results("RNAsubopt").max_mcc for rna in RNAStrandContainer ]
		]
		diffs = [
			[ x[1][i] - x[0][i] for i in range(len(x[0])) ], # 1by1 - chunk
			[ x[1][i] - x[2][i] for i in range(len(x[0])) ], # 1by1 - subopt
			[ x[0][i] - x[2][i] for i in range(len(x[0])) ]  # chunk - subopt
		]

		plt.subplot(121)
		colors = [ 'firebrick','goldenrod', 'xkcd:blue']
		labels = ["RNA-MoIP 'chunk' MCC (1 solution)", "Best RNA-MoIP '1by1' MCC", "Best RNAsubopt MCC" ]
		for y, col, lab in zip(x, colors, labels):
			plt.scatter(range(len(y)), y, color=col, label=lab, marker='o', s=2)
		plt.axvline(x=0, color='black', linewidth=1)
		plt.xlabel("RNA Strand verified structures (10 <  nt  < 100)")
		plt.ylabel("Mattews Correlation Coefficient")
		plt.ylim((-0.05,1.05))
		plt.title("(a) Performance of the prediction method")
		plt.legend(loc="lower right")

		plt.subplot(122)
		plt.axhline(y=0, color="black")
		plt.boxplot(diffs)
		plt.ylabel("Difference in max MCC")
		plt.title("(b) Difference between prediction methods")
		plt.xticks([1,2,3], ["MoIP '1by1'\n-\nMoIP 'chunk'", "MoIP '1by1'\n-\nRNAsubopt", "MoIP 'chunk'\n-\nRNAsubopt"])
		plt.subplots_adjust(wspace=0.25, bottom=0.2, left=0.1, right=0.99)

	# savefig() rajoutés par Louis Samedi soir
	plot_best_MCCs(x_noPK_fully, x_PK_fully, x_pseudobase_fully)
	plt.savefig("best_MCCs.png")
	plot_more_info()
	plt.savefig("detailed_stats.png")
	compare_subopt_MoIP()
	plt.savefig("compare_subopt_MOIP.png")
