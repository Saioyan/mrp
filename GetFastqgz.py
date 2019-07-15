#!/usr/bin/python3
import sys, os, subprocess
import pandas as pd
import numpy as np
import io,gzip

minseqlen=5000
outdir='output'
Qqcut=0.25
lowerlen=5000

argv=sys.argv
if '-i' in argv:
    indexfile=argv[argv.index('-i')+1]
    ifile=indexfile.split('/')[-1]
    print (ifile)
    bc=ifile.replace('.txt','')
    bc=bc.replace('_clean','')
if '-q' in argv:
    fqfile=argv[argv.index('-q')+1]
    fqfile=os.path.abspath(fqfile)
if '-o' in argv:
    outdir=argv[argv.index('-o')+1]
    outdir=os.path.abspath(outdir)
if '-l' in argv:
    minseqlen=int(argv[argv.index('-l')+1])
if '-t' in argv:
    mintotal=int(argv[argv.index('-t')+1])

if not os.path.exists(outdir):
    os.mkdir(outdir)

cwd=os.getcwd()
#print(fqapath)
os.chdir(cwd)
#df=pd.read_table(indexfile)
df=pd.read_csv(indexfile,sep='\t',index_col=False)
#print (df)
df=df.drop_duplicates('read_id')
df=df[['filename','read_id','run_id','sequence_length_template', 'mean_qscore_template','barcode_arrangement']]
#df=df[['filename','read_id','run_id','sequence_length_template', 'mean_qscore_template']]
df=df[df['barcode_arrangement']==bc]
df=df[df['mean_qscore_template']>7]
print (df)
Lcutvalue=df['sequence_length_template'].quantile(q=Qqcut)

if Lcutvalue>=minseqlen:
    Lcutvalue=minseqlen
if Lcutvalue<lowerlen:
    Lcutvalue=lowerlen

Qcutvalue=df['mean_qscore_template'].quantile(q=Qqcut)
print (bc)
print ('Desired bases: {0}'.format(mintotal))
print ('    The minimum quality: {0}'.format(Qcutvalue))
#print (df)

os.chdir(outdir)
#comm="sed '/^\s*$/d' {0}/joinedreads.fastq > fastq_runid.fastq".format(fqapath)
#print (comm)
#stdout=subprocess.getoutput(comm)
d=dict()

if '.gz' in fqfile:
    f=io.TextIOWrapper(gzip.open(fqfile,'r'))
else:
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

fdfL=df[df['mean_qscore_template']>=Qcutvalue]  #Only consider reads with mininmun quality at Qqcut quantile
#To sort by length
fdfL=fdfL.sort_values(['sequence_length_template'], ascending=False)

sumbasesL=fdfL['sequence_length_template'].cumsum()

#print (sumbasesL)

lineidx1=0
for i in sumbasesL.index:
    lineidx1+=1
    total1=int(sumbasesL.ix[i])
    #print (total)
    if total1 > mintotal*1/2:
        tempi=i
        break

set1=fdfL.iloc[0:lineidx1,:]
set1.to_csv(ifile+'A',sep='\t')
set1left=pd.concat([fdfL.iloc[lineidx1:,:],df[df['mean_qscore_template']<Qcutvalue]])
print ('    The minimum length: {:,} bp'.format(set1left.iloc[0,3]))
print ('  Get long read: {:,} bases'.format(total1))

fw=open('readsA.fastq','w')
for ID in set1['read_id']:
    if ID in d:
        fw.write(d[ID][0]+'\n')
        fw.write(d[ID][1]+'\n')
        fw.write('+'+'\n')
        fw.write(d[ID][2]+'\n')
    else:
        print (ID+'not in A')
fw.close()


#To sort by quality
#fdfS=set1left[set1left['sequence_length_template']>=minseqlen]
fdfS=set1left[set1left['sequence_length_template']>=Lcutvalue]
fdfS=fdfS.sort_values(['mean_qscore_template'], ascending=False)

sumbasesS=fdfS['sequence_length_template'].cumsum()
#print (fdfS)

lineidx2=0
for i in sumbasesS.index:
    lineidx2+=1
    total2=int(sumbasesS.ix[i])
    #print (total)
    if total2 > mintotal-total1:
        break
set2=fdfS.iloc[0:lineidx2,:]
set2.to_csv(ifile+'B',sep='\t')
print ('    The minimum quality: {:05.2f}'.format(set2.iloc[-1,4]))
print ('    The minimun length: {:,} bp'.format(int(Lcutvalue)))
print ('  Get total: {:,} bases'.format(total1+total2))

fw=open('readsB.fastq','w')
for ID in set2['read_id']:
    if ID in d:
        fw.write(d[ID][0]+'\n')
        fw.write(d[ID][1]+'\n')
        fw.write('+'+'\n')
        fw.write(d[ID][2]+'\n')
    else:
        print (ID+'not in B')
fw.close()


setcomb=pd.concat([set1,set2])
#print (setcomb)
setcomb.to_csv(ifile+'AB',sep='\t')
fw=open('reads.fastq','w')
for ID in setcomb['read_id']:
    if ID in d:
        fw.write(d[ID][0]+'\n')
        fw.write(d[ID][1]+'\n')
        fw.write('+'+'\n')
        fw.write(d[ID][2]+'\n')
    else:
        print (ID)
fw.close()
os.chdir (cwd)
