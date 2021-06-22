import time
import subprocess
import os

log_path = "test.log"
log = open(log_path, 'a')

def run_test(cmd,log):
    log.write(time.asctime(time.localtime(time.time())) + " : Run process \"" + cmd +"\"\n")
    log.flush()
    process = subprocess.Popen(cmd.split(' '),shell=False,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
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

cmd = ("cppsrc/Scripts/addDelimiter")
cmd1 = ("cppsrc/Scripts/countPattern")
cmd2 = ("cppsrc/Scripts/deletePdb")


myfile = open("data/fasta/benchmark.fa", "r")
name = myfile.readline()
seq = myfile.readline()

"""run_test(cmd1 + " 1JJ2" + ".fa", log)
print(cmd1 + " 1JJ2" + ".fa")
cmd2 = create_command("1JJ2")
print(cmd2)
os.system(cmd2)"""
while seq:
    name = name[6:].strip()
    print(name)
    run_test(cmd, log)
    run_test(cmd1, log)
    run_test(cmd2 + " " + name + ".fa", log)
    print(cmd2 + " " + name + ".fa")
    cmd3 = create_command(name)
    os.system(cmd3)
    name = myfile.readline()
    seq = myfile.readline()
myfile.close()
