# mrp_v3.py
雙套體分析程序mrp_v3.py
使用前設定以及安裝

設定export\
http_proxy=""\
https_proxy=""

下載更新必需套件\
https://github.com/Saioyan/mrp/blob/master/install.txt

下載其他必須程式\
https://github.com/Saioyan/mrp/blob/master/mrp_v3.py

https://github.com/Saioyan/mrp/blob/master/cleanfq.py

https://github.com/Saioyan/mrp/blob/master/GetFastqgz.py

https://github.com/Saioyan/mrp/blob/master/delfan.py

https://github.com/Saioyan/mrp/blob/master/getalignreadsq.py

https://github.com/Saioyan/mrp/blob/master/renamefa.py

https://github.com/Saioyan/mrp/blob/master/trimOverlapseq.py

https://github.com/Saioyan/mrp/blob/master/v_filter.py

過濾\
ex:\
cleanfq.py barcode01.fastq.gz barcode01.txt \
GetFastqgz.py -i barcode01_clean.txt -q barcode01_clean.fastq -t 1200000000 -o 自訂目錄


組裝\
smartdenovo.pl -c 1 reads.fastq > wtasm.mak\
make -f wtasm.mak

mrp_v3.py使用方法\
進入medaka環境\
. /saioyan/tools/medaka/venv/bin/activate\
執行mrp_v3.py\
mrp_v3.py Illumina_data_R1.fastq Illumina_data_R2.fastq reads.fastq

Output file\
雙套體序列\
hap1.fasta\
hap2.fasta



