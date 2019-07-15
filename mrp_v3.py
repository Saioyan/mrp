#!/usr/bin/python3
import sys, os, subprocess
import threading
Illumina_R1 = sys.argv[1]
Illumina_R2 = sys.argv[2]
aa=Illumina_R1[-8:]
bb=Illumina_R2[-8:]
infq=sys.argv[3]
comm='pwd'
pd=subprocess.getoutput(comm)
print(pd)
if aa != "fastq.gz" or bb!='fastq.gz' :
    print("Illumina:error filename, 'fastq.gz' available")
    sys.exit()
listt=['fastq','gz']
strt=""
for i in listt:
    if i in infq and i=='fastq':
        strt+="fastq."
    if i in infq and i=='gz':
        strt+="gz"
if strt == "fastq.gz":
    numgz=infq.index(".gz")
    comm='gunzip '+infq+''
    subprocess.getoutput(comm)
    infq=infq[0:numgz]

comm='mv '+infq+' reads.fastq'
subprocess.getoutput(comm)
comm='mv wtasm.dmo.cns wtasm.dmo.cns.fa'
subprocess.getoutput(comm)
comm='grep ">" wtasm.dmo.cns.fa'
print(subprocess.getoutput(comm))
comm='trimOverlapseq.py wtasm.dmo.cns.fa'
subprocess.getoutput(comm)
print(comm)
comm='delfan.py wtasm.dmo.cns.fa wtasm.dmo.cns.fasta'
subprocess.getoutput(comm)
comm='mv wtasm.dmo.cns.fa wtasm.dmo.cns'
subprocess.getoutput(comm)
comm='grep ">" trimmedH.fa'
a=subprocess.getoutput(comm)
a=a.replace("\n","")
li=a.split('>')
li.remove("")
print(li)
if len(li) > 0:
    fw=open('mic.fasta','w')
    for i in li :
        comm='grep -n '+i+' wtasm.dmo.cns.fasta'
        strr=subprocess.getoutput(comm)
        num=strr.index(':')
        strr=strr[0:num]
        comm='sed -i "'+strr+' d" wtasm.dmo.cns.fasta'
        subprocess.getoutput(comm)
        comm='sed -i "'+strr+' d" wtasm.dmo.cns.fasta'
        subprocess.getoutput(comm)
        comm='grep -n '+i+' trimmedseqs.fa'
        micnum=subprocess.getoutput(comm)
        micnumid=micnum.index(':')
        micnum=micnum[0:micnumid]
        comm='sed -n "'+micnum+','+str(int(micnum)+1)+'p" trimmedseqs.fa'
        outp=subprocess.getoutput(comm)
        fw.write(outp)
    fw.close()
    print('found circlur ',li)
else:
    print('not found circlur')
    sys.exit()
comm='grep ">" wtasm.dmo.cns.fasta'
pr=subprocess.getoutput(comm)
print(pr)
comm='grep ">" wtasm.dmo.cns.fasta|wc -l'
tignum=subprocess.getoutput(comm)
print('number of contig ref: ',tignum)

comm='minimap2 -x map-ont -t 100 wtasm.dmo.cns.fasta reads.fastq > mapreads.paf'
subprocess.getoutput(comm)
print(comm)
comm="cat mapreads.paf | awk '$11>=$2*0.95' > mapreads_1.paf"
subprocess.getoutput(comm)
print(comm)
comm='racon -t 100 reads.fastq mapreads_1.paf wtasm.dmo.cns.fasta > polished_con1.fa'
subprocess.getoutput(comm)
print(comm)
comm='minimap2 -x map-ont -t 100 polished_con1.fa reads.fastq > mapreads_2.paf'
subprocess.getoutput(comm)
print(comm)
comm="cat mapreads_2.paf | awk '$11>=$2*0.95' > mapreads_3.paf"
subprocess.getoutput(comm)
print(comm)
comm='racon -t 100 reads.fastq mapreads_3.paf polished_con1.fa > polished_con2.fa'
subprocess.getoutput(comm)
print(comm)

#medaka times 1
comm='minimap2 -x map-ont -t 100 polished_con2.fa reads.fastq > mapreads_4.paf'
subprocess.getoutput(comm)
print(comm)
comm="cat mapreads_4.paf | awk '$11>=$2*0.95' > mapreads_5.paf"
subprocess.getoutput(comm)
print(comm)
comm='getalignreadsq.py -i reads.fastq -m mapreads_5.paf'
subprocess.getoutput(comm)
print(comm)
comm='medaka_consensus -i reads_f.fastq -d polished_con2.fa -o medaka -t 100'
subprocess.getoutput(comm)
print(comm)
comm='mv reads_f.fastq reads_f_NC.fastq'
subprocess.getoutput(comm)
comm='cp medaka/consensus.fasta '+pd+''
subprocess.getoutput(comm)
comm='grep ">" consensus.fasta'
pr=subprocess.getoutput(comm)
print(pr)
comm='grep ">" consensus.fasta|wc -l'
m1tig=subprocess.getoutput(comm)
print('number of contig run medaka times 1: ',m1tig)
comm='mv consensus.fasta consensus_o.fasta'
subprocess.getoutput(comm)

p1='polished_con2.fa'
p2='consensus_mix.fasta'
p3=''
if m1tig > tignum :
    comm = 'mv polished_con2.fa polished_con2'
    subprocess.getoutput(comm)
    print(comm)
    comm = 'cat polished_con2 mic.fasta > polished_con2.fa'
    subprocess.getoutput(comm)
    print(comm)
    p3=p1
else:
    ###########medaka times 2
    comm='minimap2 -x map-ont -t 100 consensus_o.fasta reads.fastq > mapreads_6.paf'
    subprocess.getoutput(comm)
    print(comm)
    comm="cat mapreads_6.paf | awk '$11>=$2*0.95' > mapreads_7.paf"
    subprocess.getoutput(comm)
    print(comm)
    comm='getalignreadsq.py -i reads.fastq -m mapreads_7.paf'
    subprocess.getoutput(comm)
    print(comm)
    comm='medaka_consensus -i reads_f.fastq -d consensus_o.fasta -o medaka_2 -t 100'
    subprocess.getoutput(comm)
    print(comm)
    comm='cp medaka_2/consensus.fasta '+pd+''
    subprocess.getoutput(comm)
    ##########
    comm='grep ">" consensus.fasta'
    pr=subprocess.getoutput(comm)
    print(pr)
    comm='grep ">" consensus.fasta|wc -l'
    m2tig=subprocess.getoutput(comm)
    print('number of contig run medaka times 2: ',m2tig)
    if m2tig > tignum :
        comm = 'mv consensus_o.fasta consensus_o'
        subprocess.getoutput(comm)
        print(comm)
        comm = 'cat consensus_o mic.fasta > consensus_mix.fasta'
        subprocess.getoutput(comm)
        print(comm)
    else:
        comm = 'cat consensus.fasta mic.fasta > consensus_mix.fasta'
        subprocess.getoutput(comm)
        print(comm)
    p3=p2
pnum=p3.index('.')
p4=p3[0:pnum]+'_id.fasta'
ref = p3
comm = 'cp '+ref+' pilon_0.fasta'
subprocess.getoutput(comm)
print(comm)

refname=''
numi=0
while numi < 1000 :
    adnum = numi+1
    comm = 'minimap2 -ax sr -t 100 pilon_'+str(numi)+'.fasta '+Illumina_R1+' '+Illumina_R2+'  > pilon_'+str(numi)+'.sam '
    pr=subprocess.getoutput(comm)
    print(comm)
    print(pr)
    comm = 'samtools view -T pilon_'+str(numi)+'.fasta -b pilon_'+str(numi)+'.sam | samtools sort -o pilon_'+str(numi)+'.bam -@ 100 -'
    pr=subprocess.getoutput(comm)
    print(comm)
    print(pr)
    comm = 'samtools index pilon_'+str(numi)+'.bam'
    subprocess.getoutput(comm)
    print(comm)
    comm = 'java -jar /saioyan/tools/pilon-1.22.jar --genome pilon_'+str(numi)+'.fasta --frags pilon_'+str(numi)+'.bam --outdir pilon_sdn'+str(adnum)+' --threads 100'
    pr = subprocess.getoutput(comm)
    print(comm)
    fa = open('pilon_times_'+str(adnum)+'_p.log', 'w')
    fa.write(pr)
    fa.write('\n')
    fa.close()
    comm = 'cp pilon_sdn'+str(adnum)+'/pilon.fasta ' + pd + ''
    subprocess.getoutput(comm)
    print(comm)
    comm = 'grep ">" pilon.fasta'
    pr = subprocess.getoutput(comm)
    print(pr)
    comm = 'renamefa.py pilon.fasta pilon_'+str(adnum)+'.fasta'
    subprocess.getoutput(comm)
    print(comm)
    refname = 'pilon_'+str(adnum)+'.fasta'
    print('ref filename is : ',refname)
    comm = 'rm pilon.fasta'
    subprocess.getoutput(comm)
    comm = 'find ./ -type f |grep _p.log|wc -l'
    numlog = subprocess.getoutput(comm)
    n=1
    li=[]
    while n < int(numlog+1):
        comm = 'grep "Corrected" pilon_times_' + str(n) + '_p.log'
        outp = subprocess.getoutput(comm)
        print(outp)
        num = outp.index("Corrected")
        num2 = outp.index("snps")
        outp2 = outp[num + 10:num2 - 1]
        li.append(outp2)
        if len(li) > 1:
            num2 = int(li[n - 2]) - int(li[n - 1])
            print(num2)
            if num2 < 5:
                numi = 1001
            else:
                print('pilon times ' + str(numi) + ' finished')
                print('continue')
        n += 1
    numi += 1


refname_num = refname.index(".")
refname_id = refname[0:refname_num]+'_id.fasta'
def job(ii):
  if  ii == 0:
      comm='graphmap align -r '+refname+' -d '+infq+' --max-error 0.15 -t 50 -o graphmap.sam'
      subprocess.getoutput(comm)
      print(comm)
      comm='samtools view -T '+refname+' -b graphmap.sam | samtools sort -o graphmap.bam -@ 50 -'
      subprocess.getoutput(comm)
      print(comm)
      comm='samtools index graphmap.bam'
      subprocess.getoutput(comm)
      print(comm)
  if ii == 1:
      comm = 'minimap2 -ax sr -t 50 '+refname+' '+Illumina_R1+' '+Illumina_R2+'  > pilon_out.sam '
      pr=subprocess.getoutput(comm)
      print(comm)
      print(pr)
      comm='samtools view -T '+refname+' -b pilon_out.sam | samtools sort -o pilon_out.bam -@ 50 -'
      subprocess.getoutput(comm)
      print(comm)
      comm='samtools index pilon_out.bam'
      subprocess.getoutput(comm)
      print(comm)
      comm='freebayes -f '+refname+' pilon_out.bam > pilon_out.vcf'
      subprocess.getoutput(comm)
      print(comm)
j=2
threads = []
for ii in range(j):
  threads.append(threading.Thread(target = job, args = (ii,)))
  threads[ii].start()
for ii in range(j):
  threads[ii].join()

comm='whatshap phase --reference '+refname+' -o pilon_new.vcf pilon_out.vcf graphmap.bam pilon_out.bam --ignore-read-groups --indel'
subprocess.getoutput(comm)
print(comm)
comm='v_filter.py pilon_new.vcf pilon_new_filter.vcf'
subprocess.getoutput(comm)
print(comm)

comm='bgzip pilon_new_filter.vcf'
subprocess.getoutput(comm)
print(comm)
comm='tabix pilon_new_filter.vcf.gz'
subprocess.getoutput(comm)
print(comm)


comm='bcftools consensus -H 1 -f '+refname+' pilon_new_filter.vcf.gz > pilon_new_filter1.fasta'
subprocess.getoutput(comm)
print(comm)
comm='bcftools consensus -H 2 -f '+refname+' pilon_new_filter.vcf.gz > pilon_new_filter2.fasta'
subprocess.getoutput(comm)
print(comm)

comm='renamefa.py pilon_new_filter1.fasta hap1.fasta'
subprocess.getoutput(comm)
print(comm)
comm='renamefa.py pilon_new_filter2.fasta hap2.fasta'
subprocess.getoutput(comm)
print(comm)

comm='gunzip pilon_new_filter.vcf.gz'
subprocess.getoutput(comm)
