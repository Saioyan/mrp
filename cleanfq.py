#!/usr/bin/python3
import sys,subprocess
import io, gzip
import pandas as pd

infq=sys.argv[1]
indexfile=sys.argv[2]

threads=100
argv=sys.argv
if '-t' in argv:
    threads=int(argv[agrv.index('-t')+1])

if 'gz' in infq:
    f=io.TextIOWrapper(gzip.open(infq,'r'))
else:
    f=open(infq)

outfile='half.fa'
d=dict()
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

fw=open(outfile,'w')
for ID in d.keys():
    fw.write('>'+ID+'\n')
    seq=d[ID][1]
    start=int(len(seq)*0.05)
    end=int(len(seq)*0.45)
    fw.write(seq[start:end]+'\n')
fw.close()


comm='minimap2 -x ava-ont -t {0} {1} {2} > mapreads.paf'.format(threads,infq,outfile)
print (comm)
subprocess.getoutput(comm)

comm="cat mapreads.paf | awk '$1==$6' > mapreadssub.paf"
print (comm)
subprocess.getoutput(comm)

comm="cat mapreadssub.paf | awk '$5=="+'"-"'+"&&$11>$2*0.9 {print $1}'"
print (comm)
stdout=subprocess.getoutput(comm)
removeID=stdout.splitlines()
if (len(removeID)>0):
    print ('Removing reads......')
    print (removeID)

outfq=infq.replace('.f','_clean.f')
outfq=outfq.replace('.gz','')
fw=open(outfq,'w')
for ID in d.keys():
    if ID not in removeID:
        fw.write(d[ID][0]+'\n')
        fw.write(d[ID][1]+'\n')
        fw.write('+'+'\n')
        fw.write(d[ID][2]+'\n')
    else:
        print (ID)
fw.close()

comm='rm -f map*.paf half.fa'
print (comm)
subprocess.getoutput(comm)

outindex=indexfile.replace('.txt','_clean.txt')
df=pd.read_table(indexfile)
df=df[['filename','read_id','run_id','sequence_length_template', 'mean_qscore_template','barcode_arrangement']]
print (df)
for id in removeID:
    df=df[df.read_id!=id]
print (df)
df.to_csv(outindex,sep='\t')
