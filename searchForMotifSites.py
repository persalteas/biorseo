#!/usr/bin/python3

from sys import argv
import subprocess
import inspect
from multiprocessing import Pool, TimeoutError, cpu_count
from os import path, makedirs, getcwd, chdir, devnull

# Retrieve Jar3D Paths from file EditMe
jar3dexec = ""
HLmotifDir = ""
ILmotifDir = ""
bminDir = path.dirname(path.abspath(inspect.getfile(inspect.currentframe())))
exec(compile(open(bminDir+"/EditMe").read(), '', 'exec'))


class Loop:
    def __init__(self, header, subsequence, looptype, position):
        self.header = header
        self.seq = subsequence
        self.type = looptype
        self.position = position

    def get_header(self):
        return self.header[1:]

    def subsequence(self):
        return self.seq


class InsertionSite:
    def __init__(self, loop, csv_line):
        # BEWARE : jar3d csv output is crap because of java's locale settings.
        # On french OSes, it uses commas to delimit the fields AND as floating point delimiters !!
        # Parse with caution, and check what the csv output files look like on your system...
        info = csv_line.split(',')
        self.loop = loop                    # the Loop object that has been searched with jar3d
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


def enumerate_loops(s):
    def printLoops(s, loops):
        print(s)
        for l in loops:
            i = 0
            m = ''
            for c in l:
                while i < c[0]:
                    m += ' '
                    i += 1
                while i < c[1]+1:
                    m += '-'
                    i += 1
            while i < len(s):
                m += ' '
                i += 1
            print(m, "\tfound in position", l)

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

    if (verbose):
        printLoops(s, loops)
    return(loops)


def launchJar3d(loop):
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
    if (verbose):
        print(' '.join(cmd))
    nowhere = open(devnull, 'w')
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


filename = argv[1]
verbose = int(argv[3])
basename = filename[0:filename.index('.')]

# Retrieving possible 2D structrures from RNAsubopt
if (verbose):
    print("Retrieving possible 2D structures from RNAsubopt...")
dbn = open(basename+"_temp.dbn", "w")
subprocess.call(["RNAsubopt", "-i", filename], stdout=dbn)
dbn.close()
dbn = open(basename+"_temp.dbn", "r")
dbn.readline()
s = dbn.readline().split(' ')[0]
structures = []
l = dbn.readline()
while l:
    structures.append(l.split(' ')[0])
    l = dbn.readline()
dbn.close()
subprocess.call(["rm", basename+"_temp.dbn"])
for ss in structures:
    if (verbose):
        print(ss)
if (verbose):
    print()


# Extracting probable loops from these structures
if (verbose):
    print("Extracting probable loops from these structures...")
HLs = []
ILs = []
for ss in structures:
    loop_candidates = enumerate_loops(ss)
    for loop_candidate in loop_candidates:
        if len(loop_candidate) == 1 and loop_candidate not in HLs:
            HLs.append(loop_candidate)
        if len(loop_candidate) == 2 and loop_candidate not in ILs:
            ILs.append(loop_candidate)
if (verbose):
    print("TOTAL:")
    print(len(HLs), "probable hairpin loops found")
    print(len(ILs), "probable internal loops")
    print()

# Retrieve subsequences corresponding to the possible loops
if (verbose):
    print("Retrieving subsequences corresponding to the possible loops...")
loops = []
for i, l in enumerate(HLs):
    loops.append(Loop(">HL%d" % (i+1), s[l[0][0]-1:l[0][1]], "h", l))
    if (verbose):
        print()
        print(loops[-1].get_header(), "\t\t", l)
        print(loops[-1].subsequence())
for i, l in enumerate(ILs):
    loops.append(
        Loop(">IL%d" % (i+1), s[l[0][0]-1:l[0][1]]+'*'+s[l[1][0]-1:l[1][1]], "i", l))
    if (verbose):
        print()
        print(loops[-1].get_header(), "\t\t", l)
        print(loops[-1].subsequence())
if (verbose):
    print()

# Scanning loop subsequences against motif database
if (verbose):
    print("Scanning loop subsequences against motif database (using %d threads)..." %
          (cpu_count()-1))
pool = Pool(processes=cpu_count()-1)
insertion_sites = [x for y in pool.map(launchJar3d, loops) for x in y]
insertion_sites.sort(reverse=True)
if (verbose):
    print(len(insertion_sites), "insertions found:")

# Writing results to file
c = 0
resultsfile = open(basename+".sites.csv", "w")
resultsfile.write("Motif,Rotation,Score,Start1,End1,Start2,End2\n")
for site in insertion_sites:
    if site.score > 10:
        c += 1
        string = "FOUND with score %d:\t\t possible insertion of motif " % site.score + site.atlas_id
        if site.rotation:
            string += " (reversed)"
        string += (" on " + site.loop.get_header() + " at positions")
        if (verbose):
            print(string, site.loop.subsequence(),
                  '*'.join([str(x) for x in site.position]))
    resultsfile.write(site.atlas_id+',' +
                      str(bool(site.rotation))+",%d" % site.score+',')
    positions = [','.join([str(y) for y in x]) for x in site.position]
    if len(positions) == 1:
        positions.append("-,-")
    resultsfile.write(','.join(positions)+'\n')
if c < len(insertion_sites):
    if (verbose):
        print("... and %d more with score(s) below 10." %
              (len(insertion_sites)-c))
resultsfile.close()

# Lauching biominserter to get 2D predictions
subprocess.call([bminDir+"/bin/biominserter",
                 filename, basename+".sites.csv", argv[2], str(verbose)])
