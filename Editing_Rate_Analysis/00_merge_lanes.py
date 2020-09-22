import sys
import subprocess as sp
from collections import defaultdict

fq_files = str.strip(sp.check_output("ls *.fastq.gz", shell = True)).split("\n")

fq_keys = defaultdict(list)
for fq in fq_files:
    guide = fq.split("_")[0]
    read = fq.split("_")[3]
    fq_keys[guide + "_" + read].append(fq)

for k,v in fq_keys.items():
    #print k + ":\t" + str(v)
    if len(v) > 1:
        cmd = "cat " +  " ".join(v) + " > " + k + ".fastq.gz"
        sp.call(cmd, shell = True)
