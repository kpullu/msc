#!/bin/sh
BWA_VERSION='bwa-0.7.15'
SAM_VERSION='samtools-1.4.1'

sudo apt-get update
sudo apt-get -y install python2.7 python-pip python-dev zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev

cd $HOME
mkdir biotools
cd biotools

# setup bwa
wget https://sourceforge.net/projects/bio-bwa/files/$BWA_VERSION.tar.bz2
bunzip2 $BWA_VERSION.tar.bz2
tar -xvf $BWA_VERSION.tar
cd $BWA_VERSION
make

# setup samtools
cd $HOME/biotools

wget https://sourceforge.net/projects/samtools/files/samtools/1.4.1/$SAM_VERSION.tar.bz2
bunzip2 $SAM_VERSION.tar.bz2
tar -xvf $SAM_VERSION.tar
cd $SAM_VERSION
./configure
make

# setup bwa and samtools in path
sudo cp $HOME/hd_biotools.sh /etc/profile.d/

# install python packages
sudo pip install biopython sh pysam
