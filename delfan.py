#!/usr/bin/python3
import random, sys
import numpy as np
infile=sys.argv[1]
outfile=sys.argv[2]

d=dict()
f=open(infile)
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))
with f as fp:
    for name, seq in read_fasta(fp):
        d[name]=seq
f.close()
myseq=[]
mylen=[]
fw=open(outfile,'w')
for key in d.keys():
    fw.write(key+'\n')
    fw.write(d[key]+'\n')
fw.close()
