# mrp
雙套體分析程序mrp_v3.py
使用前設定以及安裝

設定export
http_proxy=""
https_proxy=""

下載更新必需套件
pip3 install --upgrade pip
pip3 install pandas
apt-get update
apt-get install git
<<java>>
apt-get update
apt-get install default-jre
java -version
<<albacore-2.3.4>>
wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.3.4-cp35-cp35m-manylinux1_x86_64.whl
pip3 install ont_albacore-2.3.4-cp35-cp35m-manylinux1_x86_64.whl
read_fast5_basecaller.py -v
<<SPAdes>>
wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz
tar -xzf SPAdes-3.13.0-Linux.tar.gz
<<Smartdenovo>>
git clone https://github.com/ruanjue/smartdenovo.git && (cd smartdenovo; make)
<<MUMmer3.23>>
wget https://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz
tar -zxvf MUMmer3.23.tar.gz
cd MUMmer3.23/
make && make install
<<racon>>
git clone --recursive https://github.com/isovic/racon.git racon
cd racon/
mkdir build
cd build/
apt-get install cmake
cmake -DCMAKE_BUILD_TYPE=Release ..
make
<<pilon>>
wget https://github.com/broadinstitute/pilon/releases/download/v1.22/pilon-1.22.jar
<<samtools>>
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -xjf  samtools-1.9.tar.bz2
cd samtools-1.9
make && make install
<<bcftools>>
git clone git://github.com/samtools/htslib.git
git clone git://github.com/samtools/bcftools.git
cd htslib; git pull
cd ../bcftools; git pull
make clean
make
<<minimap2>>
git clone https://github.com/lh3/minimap2 && (cd minimap2 && make)
<<htslib-1.9>>
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar jxvf htslib-1.9.tar.bz2
cd htslib-1.9/
make
make prefix=/home/jxdong/biosoft/htslib-1.9 install
echo 'export PATH=/home/jxdong/biosoft/htslib-1.9/bin:$PATH' >>~/.bashrc
source ~/.bashrc

<<other>>
apt-get gcc zlib1g-dev 
apt-get libbz2-dev 
apt-get liblzma-dev 
apt-get libffi-dev 
apt-get libncurses5-dev 
apt-get libcurl4-gnutls-dev 
apt-get libssl-dev 
apt-get make wget python3-all-dev 
apt-get python-virtualenv

<<medaka>>
git clone https://github.com/nanoporetech/medaka.git
cd medaka
make install
. ./venv/bin/activate

<<freebyes>>
git config --global url.https://github.com/.insteadOf git://github.com/
git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes
make

<<whatshap>>
pip3 install --user whatshap
vi ~/.profile
PATH=$HOME/.local/bin:$PATH
source ~/.profile
whatshap --help

<<Graphmap>>
git clone https://github.com/isovic/graphmap.git
cd graphmap/
make modules
make
vi ~/.profile
PATH=$PATH:/XXX/Linux-x64

<<bowtie2>>
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip
apt-get install zip
unzip bowtie2-2.3.4.3-linux-x86_64.zip













