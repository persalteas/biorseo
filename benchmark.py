#!/usr/bin/python3
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
from tqdm import tqdm
import subprocess
from os import path, makedirs, getcwd, chdir, devnull
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from math import sqrt, ceil
from multiprocessing import Pool, cpu_count, Manager, Value
import multiprocessing
import multiprocessing.pool
import ast, time
import pickle

# ================== DEFINITION OF THE PATHS ==============================

biorseoDir = path.realpath(".")
runDir = path.dirname(path.realpath(__file__))
bpRNAFile = argv[1]
PseudobaseFile = argv[2]
StudyCaseFile = argv[3]
outputDir = biorseoDir + "/benchmark_results/"
descfolder = biorseoDir + "/data/modules/DESC"
rinfolder = biorseoDir + "/data/modules/RIN/Subfiles/"
jsonfile = biorseoDir + "/data/modules/ISAURE/bibliotheque_a_lire/motifs_final.json"

# Create some folders to store the results
subprocess.call(["mkdir", "-p", outputDir])
subprocess.call(["mkdir", "-p", outputDir + "PK/"])
subprocess.call(["mkdir", "-p", outputDir + "noPK/"])

n_launched = Value('i', 0)
n_finished = Value('i', 0)
n_skipped = Value('i', 0)
ncores = cpu_count()

# ================== CLASSES AND FUNCTIONS ================================

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
			self.nthreads = ncores
		elif how_many_in_parallel == -1:
			self.nthreads = ncores - 1
		else:
			self.nthreads = how_many_in_parallel
		self.useless_bool = False

	def __str__(self):
		if self.func_ is None:
			s = f"{self.priority_}({self.nthreads}) [{self.comp_time}]\t{j.label:25}" + " ".join(self.cmd_)
		else:
			s = f"{self.priority_}({self.nthreads}) [{self.comp_time}]\t{j.label:25}{self.func_.__name__}(" + " ".join([str(a) for a in self.args_]) + ")"
		return s


class Method:
	def __init__(self, parent_rna, tool, data_source=None, obj_func=None, PK=True, flat=False):
		self.parent_rna = parent_rna

		# defintion of the method (theoretical)
		self.tool = tool
		self.data_source = data_source
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
				return f"{self.data_source}-{self.func}"
			else:
				return f"{self.data_source}-{self.func}-noPK"
		else:
			return self.tool

	def build_job_list(self):
		basename = self.parent_rna.basename
		fasta = outputDir+basename+".fa"

		# Things that require RNAsubopt calculations
		if self.tool == "RNAsubopt":
			if f"{basename} RNAsubopt" in issues:
				return
			self.joblist.append(Job(command=["RNAsubopt", "-i", fasta, "--outfile="+ basename + ".subopt"],
									priority=1, how_many_in_parallel=1 if self.flat else 0,
									results = outputDir + "noPK/" +  basename + ".subopt",
									label=f"{basename} RNAsubopt"))
			self.joblist.append(Job(command=["mv", basename + ".subopt", outputDir + "noPK/"],
									priority=2, how_many_in_parallel=1 if self.flat else 0,
									results = outputDir + "noPK/" +  basename + ".subopt",
									label=f"{basename} mv"))

		if f"{basename} {self.label}" in issues:
			return

		if self.tool == "biorseo":
			c = [ biorseoDir+"/bin/biorseo", "-s", fasta ]
			if self.data_source == "RIN":
				results_file = outputDir+f"{'' if self.allow_pk else 'no'}PK/"+basename+f".biorseo_rin_{self.func}"
				c += [ "-r", rinfolder]
			elif self.data_source == "DESC":
				results_file = outputDir+f"{'' if self.allow_pk else 'no'}PK/"+basename+f".biorseo_desc_{self.func}"
				c += [ "-d", descfolder]
			elif self.data_source == "JSON":
				results_file = outputDir+f"{'' if self.allow_pk else 'no'}PK/"+basename+f".biorseo_json_{self.func}"
				c += [ "-j", jsonfile]
			
			c += ["-o", results_file, "--func", self.func]
			if not self.allow_pk:
				c += ["-n"]
			self.joblist.append(Job(command=c, priority=3, timeout=3600,
									how_many_in_parallel=1 if self.flat else 10,
									results = results_file,
									label=f"{basename} {self.label}"))

		if self.tool == "Biokop-mode":
			results_file = outputDir+"PK/"+basename+".biok"
			c = [ biorseoDir+"/bin/biorseo", "-s", fasta, "-o", results_file, "--mfe", "--mea"]
			if not self.allow_pk:
				c += ["-n"]
			self.joblist.append(Job(command=c, priority=4, timeout=3600,
									how_many_in_parallel=1 if self.flat else 10,
									results = results_file,
									label=f"{basename} {self.label}"))
	
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

	def get_method(self, label):
		return self.methods[self.meth_idx[label]]

	def load_RNAsubopt_results(self):
		if not path.isfile(outputDir + "noPK/" + self.basename + ".subopt"):
			return
		rna = open(outputDir + "noPK/" + self.basename + ".subopt", "r")
		lines = rna.readlines()
		rna.close()
		for i in range(2, len(lines)):
			ss = lines[i].split(' ')[0]
			if ss not in self.get_method("RNAsubopt").predictions:
				self.get_method("RNAsubopt").predictions.append(ss)

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

	def load_biokop_results(self):
		filename = outputDir+"PK/"+basename+".biok"
		if path.isfile(filename):
			rna = open(filename, "r")
			lines = rna.readlines()
			rna.close()
			for i in range(2, len(lines)):
				ss = lines[i].split(' ')[0].split('\t')[0]
				method.predictions.append(ss)
				method.ninsertions.append(lines[i].count('+'))

	def load_results(self, include_noPK=False):
		if "Biokop-mode" in self.meth_idx.keys():
			self.load_biokop_results()
		if "RNAsubopt" in self.meth_idx.keys():
			self.load_RNAsubopt_results()
		if "DESC-A" in self.meth_idx.keys():
			self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_desc_A", self.get_method("DESC-A"))
		if "DESC-B" in self.meth_idx.keys():
			self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_desc_B", self.get_method("DESC-B"))
		if "RIN-A" in self.meth_idx.keys():
			self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_rin_A", self.get_method("RIN-A"))
		if "RIN-B" in self.meth_idx.keys():
			self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_rin_B", self.get_method("RIN-B"))
		if "JSON-A" in self.meth_idx.keys():
			self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_json_A", self.get_method("JSON-A"))
		if "JSON-B" in self.meth_idx.keys():
			self.load_biorseo_results(outputDir + "PK/" + self.basename + ".biorseo_json_B", self.get_method("JSON-B"))

		if include_noPK:
			if "DESC-A-noPK" in self.meth_idx.keys():
				self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_desc_A", self.get_method("DESC-A-noPK"))
			if "DESC-B-noPK" in self.meth_idx.keys():
				self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_desc_B", self.get_method("DESC-B-noPK"))
			if "RIN-A-noPK" in self.meth_idx.keys():
				self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_rin_A", self.get_method("RIN-A-noPK"))
			if "RIN-B-noPK" in self.meth_idx.keys():
				self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_rin_B", self.get_method("RIN-B-noPK"))
			if "JSON-A-noPK" in self.meth_idx.keys():
				self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_json_A", self.get_method("JSON-A-noPK"))
			if "JSON-B-noPK" in self.meth_idx.keys():
				self.load_biorseo_results(outputDir + "noPK/" + self.basename + ".biorseo_json_B", self.get_method("JSON-B-noPK"))

	def has_complete_results(self, with_PK):
		if not with_PK:
			if not path.isfile(outputDir + "noPK/" + self.basename + ".subopt"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_desc_A"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_desc_B"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_rin_A"): return False
			if not path.isfile(outputDir + "noPK/" + self.basename + ".biorseo_rin_B"): return False

			return True
		else:
			if not path.isfile(outputDir + "PK/" + self.basename + ".biok"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_desc_A"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_desc_B"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_rin_A"): return False
			if not path.isfile(outputDir + "PK/" + self.basename + ".biorseo_rin_B"): return False

			return True


def init(arg1, arg2, arg3):
	global n_launched, n_finished, n_skipped
	n_launched = arg1
	n_finished = arg2
	n_skipped = arg3

def execute_job(j):
	global n_launched, n_skipped, n_finished

	# Check if you really need to execute it
	if path.isfile(j.results_file) or ((j.checkFunc_ is not None) and j.checkFunc_(*j.checkArgs_)):
		# skip it
		with n_skipped.get_lock():
			n_skipped.value += 1
			n_finished.value += 1
		print(f"[{n_launched.value+n_skipped.value}/{jobcount}]\tSkipping {j.label} (already finished)")
		return (0, 0)


	# Add the job to log file and run
	with n_launched.get_lock():
		n_launched.value += 1
	if len(j.cmd_):
		logfile = open(runDir + "/log_of_the_run.sh", 'a')
		logfile.write(" ".join(j.cmd_))
		logfile.write("\n")
		logfile.close()
		print(f"[{n_launched.value+n_skipped.value}/{jobcount}]\t{j.label}")
		start_time = time.time()
		r = subprocess.run(j.cmd_, timeout=j.timeout_)
		end_time = time.time()
		if r.returncode != 0:
			if r.stderr is not None:
				print(r.stderr, flush=True)
			print(f"[{n_launched.value+n_skipped.value}/{jobcount}]\tIssue faced with {j.label}, skipping it and adding it to known issues (if not known).")
			with n_launched.get_lock():
				n_launched.value -= 1
			with n_skipped.get_lock():
				n_skipped.value += 1
			if j.label not in issues:
				issues.add(j.label)
				with open("benchmark_results/known_issues.txt", "a") as iss:
					iss.write(j.label+"\n")
	elif j.func_ is not None:
		print(f"[{n_launched.value+n_skipped.value}/{jobcount}]\t{j.func_.__name__}({', '.join([str(a) for a in j.args_])})")
		start_time = time.time()
		r = j.func_(*j.args_)
		end_time = time.time()

	# Job is finished
	with n_finished.get_lock():
		n_finished.value += 1
	t = end_time - start_time
	return (t,r)

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
	try:
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
	except IndexError: # pop from empty list
		print("Error in structure :", structure)
		exit(0)
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

def load_from_dbn(file, header_style=3):
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
			if not '(' in struct:
				continue # ignore linear structures
			if is_canonical_nts(seq) and is_canonical_bps(struct):
				# keeps what's inside brackets at the end as the filename
				if header_style == 1: container.append(RNA(header.replace('/', '_').split('(')[-1][:-1], header, seq, struct))
				# keeps what's inside square brackets at the end as the filename
				if header_style == 2: container.append(RNA(header.replace('/', '_').split('[')[-1][:-41], header, seq, struct))
				# keeps all the header as filename
				if header_style == 3: container.append(RNA(header[1:], header, seq, struct))
				if '[' in struct: counter += 1
	db.close()
	return container, counter

def get_bpRNA_statistics(include_noPK=True):

	print("\nLoading bpRNA results from files...")

	# load results in objects
	for instance in tqdm(bpRNAContainer, desc="bpRNA instances"):
		instance.load_results(include_noPK=True)
		instance.evaluate()

	RNAs_fully_predicted_noPK = [ x for x in bpRNAContainer if x.has_complete_results(with_PK=False) ]

	# Get max MCCs for each method without PK, and see who is complete
	x_noPK = [
		[ rna.get_method("RNAsubopt").max_mcc if rna.get_method("RNAsubopt").n_pred else print(rna.basename, "has no RNAsubopt structure (linear)") for rna in bpRNAContainer if rna.get_method("RNAsubopt").n_pred ],
		[ rna.get_method("DESC-A-noPK").max_mcc  if rna.get_method("DESC-A-noPK").n_pred else print(rna.basename, "has no DESC-A-noPK") for rna in bpRNAContainer if rna.get_method("DESC-A-noPK").n_pred ],
		[ rna.get_method("DESC-B-noPK").max_mcc  if rna.get_method("DESC-B-noPK").n_pred else print(rna.basename, "has no DESC-B-noPK") for rna in bpRNAContainer if rna.get_method("DESC-B-noPK").n_pred ],
		[ rna.get_method("RIN-A-noPK").max_mcc  if rna.get_method("RIN-A-noPK").n_pred else print(rna.basename, "has no RIN-A-noPK") for rna in bpRNAContainer if rna.get_method("RIN-A-noPK").n_pred ],
		[ rna.get_method("RIN-B-noPK").max_mcc  if rna.get_method("RIN-B-noPK").n_pred else print(rna.basename, "has no RIN-B-noPK") for rna in bpRNAContainer if rna.get_method("RIN-B-noPK").n_pred ],
		[ rna.get_method("JSON-A-noPK").max_mcc  if rna.get_method("JSON-A-noPK").n_pred else print(rna.basename, "has no JSON-A-noPK") for rna in bpRNAContainer if rna.get_method("JSON-A-noPK").n_pred ],
		[ rna.get_method("JSON-B-noPK").max_mcc  if rna.get_method("JSON-B-noPK").n_pred else print(rna.basename, "has no JSON-B-noPK") for rna in bpRNAContainer if rna.get_method("JSON-B-noPK").n_pred ],
	]

	x_noPK_fully = [
		[ rna.get_method("RNAsubopt").max_mcc for rna in RNAs_fully_predicted_noPK ],
		[ rna.get_method("DESC-A-noPK").max_mcc for rna in RNAs_fully_predicted_noPK ],
		[ rna.get_method("DESC-B-noPK").max_mcc for rna in RNAs_fully_predicted_noPK ],
		[ rna.get_method("RIN-A-noPK").max_mcc for rna in RNAs_fully_predicted_noPK ],
		[ rna.get_method("RIN-B-noPK").max_mcc for rna in RNAs_fully_predicted_noPK ],
		[ rna.get_method("JSON-A-noPK").max_mcc for rna in RNAs_fully_predicted_noPK ],
		[ rna.get_method("JSON-B-noPK").max_mcc for rna in RNAs_fully_predicted_noPK ],
	]  # We ensure having the same number of RNAs in every sample by discarding the ones for which computations did not ended/succeeded.


	print()
	print("Without PK:")
	print("%s RNAsubopt predictions" % is_all(len(x_noPK[0]), bpRNA_tot))
	print("%s biorseo + DESC + f1A predictions" % is_all(len(x_noPK[1]), bpRNA_tot))
	print("%s biorseo + DESC + f1B predictions" % is_all(len(x_noPK[2]), bpRNA_tot))
	print("%s biorseo + RIN + f1A predictions" % is_all(len(x_noPK[3]), bpRNA_tot))
	print("%s biorseo + RIN + f1B predictions" % is_all(len(x_noPK[4]), bpRNA_tot))
	print("%s biorseo + JSON + f1A predictions" % is_all(len(x_noPK[5]), bpRNA_tot))
	print("%s biorseo + JSON + f1B predictions" % is_all(len(x_noPK[6]), bpRNA_tot))

	print("==> %s ARN were predicted with all methods successful." % is_all(len(x_noPK_fully[0]), bpRNA_tot) )

	# Stat tests
	# Search if all methods are equal in positions with Friedman test:
	test = stats.friedmanchisquare(*x_noPK_fully)
	print("Friedman test without PK: H0 = 'The position parameter of all distributions is equal', p-value = ", test.pvalue)

	RNAs_fully_predicted_PK = [ x for x in bpRNAContainer if x.has_complete_results(with_PK=True) ]

	# Get max MCCs for each method with PK, and see who is complete
	x_PK = [
		[ rna.get_method("Biokop-mode").max_mcc if rna.get_method("Biokop-mode").n_pred else print(rna.basename, "has no Biokop-mode") for rna in bpRNAContainer if rna.get_method("Biokop-mode").n_pred ],
		[ rna.get_method("DESC-A").max_mcc if rna.get_method("DESC-A").n_pred else print(rna.basename, "has no DESC-A") for rna in bpRNAContainer if rna.get_method("DESC-A").n_pred ],
		[ rna.get_method("DESC-B").max_mcc if rna.get_method("DESC-B").n_pred else print(rna.basename, "has no DESC-B") for rna in bpRNAContainer if rna.get_method("DESC-B").n_pred ],
		[ rna.get_method("RIN-A").max_mcc if rna.get_method("RIN-A").n_pred else print(rna.basename, "has no RIN-A") for rna in bpRNAContainer if rna.get_method("RIN-A").n_pred ],
		[ rna.get_method("RIN-B").max_mcc if rna.get_method("RIN-B").n_pred else print(rna.basename, "has no RIN-B") for rna in bpRNAContainer if rna.get_method("RIN-B").n_pred ],
		[ rna.get_method("JSON-A").max_mcc if rna.get_method("JSON-A").n_pred else print(rna.basename, "has no JSON-A") for rna in bpRNAContainer if rna.get_method("JSON-A").n_pred ],
		[ rna.get_method("JSON-B").max_mcc if rna.get_method("JSON-B").n_pred else print(rna.basename, "has no JSON-B") for rna in bpRNAContainer if rna.get_method("JSON-B").n_pred ],
	]

	# We ensure having the same number of RNAs in every sample by discarding the one for which computations did not ended/succeeded.
	x_PK_fully = [
		[ rna.get_method("Biokop-mode").max_mcc for rna in RNAs_fully_predicted_PK ],
		[ rna.get_method("DESC-A").max_mcc for rna in RNAs_fully_predicted_PK ],
		[ rna.get_method("DESC-B").max_mcc for rna in RNAs_fully_predicted_PK ],
		[ rna.get_method("RIN-A").max_mcc for rna in RNAs_fully_predicted_PK ],
		[ rna.get_method("RIN-B").max_mcc for rna in RNAs_fully_predicted_PK ],
		[ rna.get_method("JSON-A").max_mcc for rna in RNAs_fully_predicted_PK ],
		[ rna.get_method("JSON-B").max_mcc for rna in RNAs_fully_predicted_PK ],
	]

	print()
	print("With PK:")
	print("%s biokop-mode predictions" % is_all(len(x_PK[0]), bpRNA_tot))
	print("%s biorseo + DESC + f1A predictions" % is_all(len(x_PK[1]), bpRNA_tot))
	print("%s biorseo + DESC + f1B predictions" % is_all(len(x_PK[2]), bpRNA_tot))
	print("%s biorseo + RIN + f1A predictions" % is_all(len(x_PK[3]), bpRNA_tot))
	print("%s biorseo + RIN + f1B predictions" % is_all(len(x_PK[4]), bpRNA_tot))
	print("%s biorseo + JSON + f1A predictions" % is_all(len(x_PK[5]), bpRNA_tot))
	print("%s biorseo + JSON + f1B predictions" % is_all(len(x_PK[6]), bpRNA_tot))

	print("==> %s ARN were predicted with all methods successful." % is_all(len(x_PK_fully[0]), bpRNA_tot) )

	# stat tests
	# First, search if all methods are equal in positions with Friedman test:
	test = stats.friedmanchisquare(*x_PK_fully)
	print("Friedman test with PK: H0 = 'The position parameter of all distributions is equal', p-value = ", test.pvalue)
	# it looks like some methods do better. Let's test the difference:
	test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[1])
	print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and DESC-A are equal', p-value = ", test.pvalue)
	test = stats.wilcoxon(x_PK_fully[0], x_PK_fully[2])
	print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and DESC-B are equal', p-value = ", test.pvalue)
	
	n = [
		[ rna.get_method("Biokop-mode").n_pred  for rna in bpRNAContainer if rna.get_method("Biokop-mode").n_pred ],
		[ rna.get_method("RNAsubopt").n_pred  for rna in bpRNAContainer if rna.get_method("RNAsubopt").n_pred ],
		[ rna.get_method("DESC-A").n_pred for rna in bpRNAContainer if rna.get_method("DESC-A").n_pred ],
		[ rna.get_method("DESC-B").n_pred for rna in bpRNAContainer if rna.get_method("DESC-B").n_pred ],
		[ rna.get_method("RIN-A").n_pred for rna in bpRNAContainer if rna.get_method("RIN-A").n_pred ],
		[ rna.get_method("RIN-B").n_pred for rna in bpRNAContainer if rna.get_method("RIN-B").n_pred ],
		[ rna.get_method("JSON-A").n_pred for rna in bpRNAContainer if rna.get_method("JSON-A").n_pred ],
		[ rna.get_method("JSON-B").n_pred for rna in bpRNAContainer if rna.get_method("JSON-B").n_pred ],
	]

	r = [
		[ rna.get_method("DESC-A").ratio for rna in bpRNAContainer if rna.get_method("DESC-A").n_pred > 1 ],
		[ rna.get_method("DESC-B").ratio for rna in bpRNAContainer if rna.get_method("DESC-B").n_pred > 1 ],
		[ rna.get_method("RIN-A").ratio for rna in bpRNAContainer if rna.get_method("RIN-A").n_pred > 1 ],
		[ rna.get_method("RIN-B").ratio for rna in bpRNAContainer if rna.get_method("RIN-B").n_pred > 1 ],
		[ rna.get_method("JSON-A").ratio for rna in bpRNAContainer if rna.get_method("JSON-A").n_pred > 1 ],
		[ rna.get_method("JSON-B").ratio for rna in bpRNAContainer if rna.get_method("JSON-B").n_pred > 1 ],
	]

	max_i = [
		[ max(rna.get_method("DESC-A").ninsertions) for rna in bpRNAContainer if rna.get_method("DESC-A").n_pred ],
		[ max(rna.get_method("DESC-B").ninsertions) for rna in bpRNAContainer if rna.get_method("DESC-B").n_pred ],
		[ max(rna.get_method("RIN-A").ninsertions) for rna in bpRNAContainer if rna.get_method("RIN-A").n_pred ],
		[ max(rna.get_method("RIN-B").ninsertions) for rna in bpRNAContainer if rna.get_method("RIN-B").n_pred ],
		[ max(rna.get_method("JSON-A").ninsertions) for rna in bpRNAContainer if rna.get_method("JSON-A").n_pred ],
		[ max(rna.get_method("JSON-B").ninsertions) for rna in bpRNAContainer if rna.get_method("JSON-B").n_pred ],
	]
	
	return x_noPK_fully, x_PK_fully, n, r, max_i

def get_Pseudobase_statistics():

	# load results in objects
	print("\nLoading Pseudobase results from files...")
	for instance in tqdm(PseudobaseContainer, desc="Pseudobase instances"):
		instance.load_results()
		instance.evaluate()

	RNAs_fully_predicted_Pseudobase = [ x for x in PseudobaseContainer if x.has_complete_results(with_PK=True)]

	x_pseudobase = [
		[ rna.get_method("Biokop-mode").max_mcc if rna.get_method("Biokop-mode").n_pred else print(rna.basename, "has no") for rna in PseudobaseContainer if rna.get_method("Biokop-mode").n_pred ],
		[ rna.get_method("RNAsubopt").max_mcc if rna.get_method("RNAsubopt").n_pred else print(rna.basename, "has no") for rna in PseudobaseContainer if rna.get_method("RNAsubopt").n_pred ],
		[ rna.get_method("DESC-A").max_mcc if rna.get_method("DESC-A").n_pred else print(rna.basename, "has no") for rna in PseudobaseContainer if rna.get_method("DESC-A").n_pred ],
		[ rna.get_method("DESC-B").max_mcc if rna.get_method("DESC-B").n_pred else print(rna.basename, "has no") for rna in PseudobaseContainer if rna.get_method("DESC-B").n_pred ],
		[ rna.get_method("RIN-A").max_mcc if rna.get_method("RIN-A").n_pred else print(rna.basename, "has no") for rna in PseudobaseContainer if rna.get_method("RIN-A").n_pred ],
		[ rna.get_method("RIN-B").max_mcc if rna.get_method("RIN-B").n_pred else print(rna.basename, "has no") for rna in PseudobaseContainer if rna.get_method("RIN-B").n_pred ],
		[ rna.get_method("JSON-A").max_mcc if rna.get_method("JSON-A").n_pred else print(rna.basename, "has no") for rna in PseudobaseContainer if rna.get_method("JSON-A").n_pred ],
		[ rna.get_method("JSON-B").max_mcc if rna.get_method("JSON-B").n_pred else print(rna.basename, "has no") for rna in PseudobaseContainer if rna.get_method("JSON-B").n_pred ],
	]

	# We ensure having the same number of RNAs in every sample by discarding the one for which computations did not ended/succeeded.
	x_pseudobase_fully = [
		[ rna.get_method("Biokop-mode").max_mcc for rna in RNAs_fully_predicted_Pseudobase ],
		[ rna.get_method("RNAsubopt").max_mcc for rna in RNAs_fully_predicted_Pseudobase ],
		[ rna.get_method("DESC-A").max_mcc for rna in RNAs_fully_predicted_Pseudobase ],
		[ rna.get_method("DESC-B").max_mcc for rna in RNAs_fully_predicted_Pseudobase ],
		[ rna.get_method("RIN-A").max_mcc for rna in RNAs_fully_predicted_Pseudobase ],
		[ rna.get_method("RIN-B").max_mcc for rna in RNAs_fully_predicted_Pseudobase ],
		[ rna.get_method("JSON-A").max_mcc for rna in RNAs_fully_predicted_Pseudobase ],
		[ rna.get_method("JSON-B").max_mcc for rna in RNAs_fully_predicted_Pseudobase ],
	]


	print()
	print("With PK:")
	print("%s Biokop-mode predictions" % is_all(len(x_pseudobase[0]), Pseudobase_tot))
	print("%s RNAsubopt predictions" % is_all(len(x_pseudobase[1]), Pseudobase_tot))
	print("%s biorseo + DESC + f1A predictions" % is_all(len(x_pseudobase[2]), Pseudobase_tot))
	print("%s biorseo + DESC + f1B predictions" % is_all(len(x_pseudobase[3]), Pseudobase_tot))
	print("%s biorseo + RIN + f1A predictions" % is_all(len(x_pseudobase[4]), Pseudobase_tot))
	print("%s biorseo + RIN + f1B predictions" % is_all(len(x_pseudobase[5]), Pseudobase_tot))
	print("%s biorseo + JSON + f1A predictions" % is_all(len(x_pseudobase[6]), Pseudobase_tot))
	print("%s biorseo + JSON + f1B predictions" % is_all(len(x_pseudobase[7]), Pseudobase_tot))

	print("==> %s ARN were predicted with all methods successful." % is_all(len(x_pseudobase_fully[0]), Pseudobase_tot) )

	# stat tests
	# First, search if all methods are equal in positions with Friedman test:
	test = stats.friedmanchisquare(*x_pseudobase_fully)
	print("Friedman test with PK: H0 = 'The position parameter of all distributions is equal', p-value = ", test.pvalue)
	# it looks like some methods do better. Let's test the difference:
	test = stats.wilcoxon(x_pseudobase_fully[0], x_pseudobase_fully[2])
	print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and DESC-A are equal', p-value = ", test.pvalue)
	test = stats.wilcoxon(x_pseudobase_fully[0], x_pseudobase_fully[3])
	print("Wilcoxon signed rank test with PK: H0 = 'The position parameter of Biokop and DESC-B are equal', p-value = ", test.pvalue)
	return x_pseudobase_fully

def print_StudyCase_results():
	print("\nLoading study case results from files...")

	# load results in objects
	for instance in StudycaseContainer:
		instance.load_results()
		instance.evaluate(verbose=True)

# ================= EXTRACTION OF STRUCTURES FROM FILES ===============================

if __name__ == '__main__':

	print("> Loading files...", flush=True)
	bpRNAContainer, bpRNA_pk_counter = load_from_dbn(bpRNAFile, header_style=1)
	PseudobaseContainer, Pseudobase_pk_counter = load_from_dbn(PseudobaseFile, header_style=3)
	StudycaseContainer, StudyCase_pk_counter = load_from_dbn(StudyCaseFile, header_style=1)

	for nt, number in ignored_nt_dict.items():
		print("\t> ignored %d sequences because of char %c" % (number, nt))

	bpRNA_tot = len(bpRNAContainer)
	Pseudobase_tot = len(PseudobaseContainer)
	StudyCase_tot = len(StudycaseContainer)
	print("\t> Loaded %d RNAs of length between 10 and 100 from RNA Strand. %d of them contain pseudoknots." % (bpRNA_tot, bpRNA_pk_counter))
	print("\t> Loaded %d RNAs of length between 10 and 100 from Pseudobase. %d of them contain pseudoknots." % (Pseudobase_tot, Pseudobase_pk_counter))
	print("\t> Loaded %d RNAs of length between 10 and 100 from study case. %d of them contain pseudoknots." % (StudyCase_tot, StudyCase_pk_counter))

	issues = set()
	if path.isfile("benchmark_results/known_issues.txt"):
		with open("benchmark_results/known_issues.txt") as f:
			issues = set([ j[:-1] for j in f.readlines() ])
		print(f"\t> Ignoring {len(issues)} known failing jobs.")

	#================= PREDICTION OF STRUCTURES ===============================

	#define job list
	print("> Defining jobs...")
	fulljoblist = []
	joblabel_list = set()

	if path.isfile("containers.pickle"):
		with open("containers.pickle", "rb") as cont:
			bpRNAContainer, PseudobaseContainer = pickle.load(cont)
	else:
		for instance in tqdm(bpRNAContainer, desc="bpRNA jobs"):
			instance.add_method_evaluation(instance, "Biokop-mode", flat=False)
			instance.add_method_evaluation(instance, "RNAsubopt", flat=False)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", obj_func="A", PK=False)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", obj_func="A", PK=True)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", obj_func="B", PK=False)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", obj_func="B", PK=True)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", obj_func="A", PK=False)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", obj_func="A", PK=True)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", obj_func="B", PK=False)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", obj_func="B", PK=True)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="JSON", obj_func="A", PK=False)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="JSON", obj_func="A", PK=True)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="JSON", obj_func="B", PK=False)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="JSON", obj_func="B", PK=True)

			for method in instance.methods:
				for i in range(len(method.joblist)):
					j = method.joblist[i]
					if j.label in joblabel_list: # look for a duplicate job
						# for index, job in enumerate(fulljoblist):
						# 	if job.label == j.label:
						# 		method.joblist[i] = fulljoblist[index] # point to the previous occurrence
						# 		break
						continue
					else:
						fulljoblist.append(j)
						joblabel_list.add(j.label)

		for instance in tqdm(PseudobaseContainer, desc="Pseudobase jobs"):
			instance.add_method_evaluation(instance, "Biokop-mode", flat=False)
			instance.add_method_evaluation(instance, "RNAsubopt", flat=False)
			instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", obj_func="A")
			instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", obj_func="B")
			instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", obj_func="A")
			instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN", obj_func="B")
			instance.add_method_evaluation(instance, tool="biorseo", data_source="JSON", obj_func="A")
			instance.add_method_evaluation(instance, tool="biorseo", data_source="JSON", obj_func="B")

			for method in instance.methods:
				for i in range(len(method.joblist)):
					j = method.joblist[i]
					if j.label in joblabel_list: # look for a duplicate job
						# for index, job in enumerate(fulljoblist):
						# 	if job.label == j.label:
						# 		method.joblist[i] = fulljoblist[index] # point to the previous occurrence
						# 		break
						continue
					else:
						fulljoblist.append(j)
						joblabel_list.add(j.label)

		with open("containers.pickle", "wb") as cont:
			pickle.dump((bpRNAContainer, PseudobaseContainer), cont)

	for instance in StudycaseContainer: # We need to define these separately because we do not want concurrency, to measure proper run times.
		instance.add_method_evaluation(instance, "Biokop-mode", flat=True)
		instance.add_method_evaluation(instance, "RNAsubopt", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", obj_func="A", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="DESC", obj_func="B", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN",  obj_func="A", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="RIN",  obj_func="B", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="JSON",  obj_func="A", flat=True)
		instance.add_method_evaluation(instance, tool="biorseo", data_source="JSON",  obj_func="B", flat=True)

		for method in instance.methods:
			for i in range(len(method.joblist)):
				j = method.joblist[i]
				if j.label in joblabel_list: # look for a duplicate job
					# for index, job in enumerate(fulljoblist):
					# 	if job.label == j.label:
					# 		method.joblist[i] = fulljoblist[index] # point to the previous occurrence
					# 		break
					continue
				else:
					fulljoblist.append(j)
					joblabel_list.add(j.label)

	# sort jobs in a tree structure
	jobs = {}
	jobcount = len(fulljoblist)
	for job in fulljoblist:
		if job.priority_ not in jobs.keys():
			jobs[job.priority_] = {}
		if job.nthreads not in jobs[job.priority_].keys():
			print(f"New job priority/concurrency: {job.priority_} {job.nthreads}")
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
				p = Pool(initializer = init, initargs = (n_launched, n_finished, n_skipped), processes=n, maxtasksperchild=10)
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


	# ================= Statistics ========================

	if path.isfile("pickleresults.pickle"):
		with open("pickleresults.pickle", "rb") as rf:
			t = pickle.load(rf)
		x_noPK_fully, x_PK_fully, n, r, max_i, x_pseudobase_fully = t
	else:
		x_noPK_fully, x_PK_fully, n, r, max_i = get_bpRNA_statistics()
		x_pseudobase_fully = get_Pseudobase_statistics()
		with open("pickleresults.pickle", "wb") as rf:
			pickle.dump((x_noPK_fully, x_PK_fully, n, r, max_i, x_pseudobase_fully), rf)
	print_StudyCase_results()

	# ================= PLOTS OF RESULTS =======================================

	colors = [
		'#911eb4', #purple
		'#000075', #navy
		'#ffe119', '#ffe119', # yellow
		'#e6194B', '#e6194B', #red
		'#3cb44b', '#3cb44b', #green
		'#4363d8', '#4363d8', #blue
	]

	def plot_best_MCCs(x_noPK_fully, x_PK_fully, x_pseudobase_fully):

		print("Best MCCs...")
		labels = [
			"Biokop-mode\n", "RNAsubopt",
			"$f_{1A}$", "$f_{1B}$",
			"$f_{1A}$", "$f_{1B}$",
			"$f_{1A}$", "$f_{1B}$",
			"$f_{1A}$", "$f_{1B}$",
			"$f_{1A}$", "$f_{1B}$",
		]

		fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,5), dpi=150)
		fig.suptitle(" \n ")
		fig.subplots_adjust(left=0.1, right=0.97, top=0.83, bottom=0.05)



		# Line 1 : no Pseudoknots
		xpos = [ 1+x for x in range(len(x_noPK_fully)) ] # skip Biokop's column
		vplot = axes[0].violinplot(x_noPK_fully, showmeans=False, showmedians=False, showextrema=False,
								   points=len(x_noPK_fully[0]), positions=xpos)
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
		xpos = [ 0 ] + [ i for i in range(4,20) ]
		vplot = axes[1].violinplot(x_PK_fully, showmeans=False, showmedians=False, showextrema=False,
								   points=len(x_PK_fully[0]), positions=xpos)
		for patch, color in zip(vplot['bodies'], colors[:1] + colors[4:]):
			patch.set_facecolor(color)
			patch.set_edgecolor(color)
			patch.set_alpha(0.5)
		quartile1, medians, quartile3 = np.percentile(x_PK_fully, [25, 50, 75], axis=1)
		axes[1].scatter(xpos, medians, marker='o', color='k', s=30, zorder=3)
		axes[1].vlines(xpos, quartile1, quartile3, color='k', linestyle='-', lw=1)
		for x, y1, y2 in zip(xpos, quartile1, quartile3):
			bar1 = Line2D([x-0.1, x+0.1], [y1, y1], color="k", lw=1)
			bar2 = Line2D([x-0.1, x+0.1], [y2, y2], color="k", lw=1)
			axes[1].add_line(bar1)
			axes[1].add_line(bar2)
		axes[1].set_ylabel("(B)\nmax MCC\n(%d RNAs)" % (len(x_PK_fully[0])), fontsize=12)

		# Line 3 : all methods on pseudoknotted dataset
		xpos = [ x for x in range(len(x_pseudobase_fully)) ]
		vplot = axes[2].violinplot(x_pseudobase_fully, showmeans=False, showmedians=False, showextrema=False,
								   points=len(x_pseudobase_fully[0]), positions=xpos)
		for patch, color in zip(vplot['bodies'], colors):
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

		for ax in axes:
			ax.set_ylim((0.0, 1.01))
			ax.set_xlim((-1, 20))
			yticks = [ i/10 for i in range(0, 11, 2) ]
			ax.set_yticks(yticks)
			for y in yticks:
				ax.axhline(y=y, color="grey", linestyle="--", linewidth=1)
			ax.tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)
			ax.set_xticks([i for i in range(20)])
		axes[0].tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
		axes[0].set_xticklabels(labels)
		for i, tick in enumerate(axes[0].xaxis.get_major_ticks()):
			if i<4: # Reduce size of Biokop, RNAsubopt and RNA-MoIP labels to stay readable
				tick.label2.set_fontsize(10)
			else:
				tick.label2.set_fontsize(12)

	def plot_more_info():
		# ======= number of solutions, insertion ratio, etc ========================

		# Figure : number of solutions
		print("Number of solutions...")
		plt.figure(figsize=(9,2.5), dpi=80)
		plt.suptitle(" \n ")
		plt.subplots_adjust(left=0.05, right=0.97, top=0.6, bottom=0.05)
		xpos = [ x for x in range(len(n)) ]
		for y in [ 10*x for x in range(8) ]:
			plt.axhline(y=y, color="grey", linestyle="-", linewidth=0.5)
		plt.axhline(y=1, color="grey", linestyle="-", linewidth=0.5)
		vplot = plt.violinplot(n, showmeans=False, showmedians=False, showextrema=False, points=len(n[0]), positions=xpos)
		for patch, color in zip(vplot['bodies'], colors):
			patch.set_facecolor(color)
			patch.set_edgecolor(color)
			patch.set_alpha(0.5)
		labels = [
			"Biokop",
			"RNAsubopt","RNA-MoIP\n1by1", "RNA-MoIP\nchunk",
			"$f_{1A}$", "$f_{1B}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$",
			"$f_{1A}$", "$f_{1B}$", "$f_{1C}$", "$f_{1D}$",
			"$f_{1A}$", "$f_{1B}$"
		]
		plt.xlim((-1,20))
		plt.tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)
		plt.xticks([ i for i in range(len(labels))], labels)
		plt.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
		for i, tick in enumerate(plt.gca().xaxis.get_major_ticks()):
			if i<4: # Reduce size of RNA-MoIP labels to stay readable
				# tick.label2.set_fontsize(8)
				tick.label2.set_rotation(90)
			else:
				tick.label2.set_fontsize(12)
		plt.yticks([ 20*x for x in range(3) ])
		plt.ylim((0,40))
		plt.savefig("number_of_solutions.png")

		# Figure : max number of insertions and ratio
		fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10,4), dpi=80)
		fig.suptitle(" \n ")
		fig.subplots_adjust(left=0.09, right=0.99, top=0.7, bottom=0.05)
		
		# Figure : max inserted
		print("Max inserted...")
		xpos = [ x for x in range(18) ]
		axes[0].set_yticks([ 5*x for x in range(3) ])
		for y in [ 2*x for x in range(7) ]:
			axes[0].axhline(y=y, color="grey", linestyle="-", linewidth=0.5)
		vplot = axes[0].violinplot(max_i, showmeans=False, showmedians=False, showextrema=False, points=len(max_i[0]), positions=xpos)
		for patch, color in zip(vplot['bodies'], colors[2:]):
			patch.set_facecolor(color)
			patch.set_edgecolor(color)
			patch.set_alpha(0.5)
		axes[0].set_ylabel("(A)", fontsize=12)

		# Figure : insertion ratio
		print("Ratio of insertions...")
		xpos = [ 0 ] + [ x for x in range(2, 1+len(r)) ]
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

		labels = labels[2:]
		for ax in axes:
			ax.set_xlim((-1,18))
			ax.tick_params(top=False, bottom=False, labeltop=False, labelbottom=False)
			ax.set_xticks([ i for i in range(18)])
		axes[0].tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
		axes[0].set_xticklabels(labels)
		for i, tick in enumerate(axes[0].xaxis.get_major_ticks()):
			if i<2: # Reduce size of RNA-MoIP labels to stay readable
				# tick.label2.set_fontsize(9)
				tick.label2.set_rotation(90)
			else:
				tick.label2.set_fontsize(12)

	plot_best_MCCs(x_noPK_fully, x_PK_fully, x_pseudobase_fully)
	plt.savefig("best_MCCs.png")
	plot_more_info()
	plt.savefig("detailed_stats.png")
	compare_subopt_MoIP()
	plt.savefig("compare_subopt_MOIP.png")
