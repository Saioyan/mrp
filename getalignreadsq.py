import sys, os, subprocess
import pandas as pd
import numpy as np
import io,gzip

argv=sys.argv
if '-m' in argv:
    mappingfile=argv[argv.index('-m')+1]
if '-i' in argv:
    fqfile=argv[argv.index('-i')+1]


mappingfile=os.path.abspath(mappingfile)
fqile=os.path.abspath(fqfile)

d=dict()
f=open(fqfile)

while True:
    h=f.readline()
    if not h: break
    h=h.replace('\n','')
    #print (h)
    rID=h.split()[0]
    rID=rID.replace('@','')
    #print (rID)
    seq=f.readline().replace('\n','')
    qh=f.readline().replace('\n','')
    qual=f.readline().replace('\n','')
    d[rID]=[]
    d[rID].append(h)
    d[rID].append(seq)
    d[rID].append(qual)
f.close()

df=pd.read_csv(mappingfile,names=['Qname','Qlen','Qstart','Qend','strand','Tname','Tlen','Tstart','Tend','Nmatch','Alen','MQ1','MQ2','MQ3','MQ4','MQ5','MQ6','MQ7'])
#print(df)
df=df.drop_duplicates('Qname')
#print(df)
getID = np.unique(df[['Qname']].values)
fw=open('reads_f.fastq','w')
for header in d.keys():
    if header in getID:
        fw.write('@'+header+'\n')
        fw.write(d[header][1]+'\n')
        fw.write('+'+'\n')
        fw.write(d[header][2]+'\n')

fw.close()
