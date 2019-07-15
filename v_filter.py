#!/usr/bin/python3
import random, sys, os
infile=sys.argv[1]
outfile=sys.argv[2]
f=open(infile)
fw=open(outfile,'w')
while True:
    l=f.readline()
    if l.startswith('#') or ':PS' in l:
       fw.write(l)
    if not l: break
f.close()
fw.close()
