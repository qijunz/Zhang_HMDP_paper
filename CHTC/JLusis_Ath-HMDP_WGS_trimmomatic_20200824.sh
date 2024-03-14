#!/bin/bash

# unpack trimmomatic and java installation
unzip Trimmomatic-0.39.zip
tar -xzvf jre-8u231-linux-x64.tar.gz

# make sure the script will use your trimmomatic and java installation
# adding Trimmomatic-0.39/ to PATH may not work
export PATH=$(pwd)/Trimmomatic-0.39:$PATH
export JAVA_HOME=$(pwd)/jre1.8.0_231
export PATH=$(pwd)/jre1.8.0_231/bin:$PATH

### $1 is the sample id, e.g., AA-0826
### $2 is the long sample id, e.g., AA-0826_S68

# both R1 and R2 fastq file
cp /staging/qzhang333/$1*.fastq.gz .          

# unzip raw reads
gzip -d $2_L001_R1_001.fastq.gz
gzip -d $2_L001_R2_001.fastq.gz

# clear .fastq.gz input
# rm *.fastq.gz

# trimming: quality trmming, for R1/R2 in lane1
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 1 \
$2_L001_R1_001.fastq \
$2_L001_R2_001.fastq \
$1_R1_trimmomaticTrimmed_paired.fastq \
$1_R1_trimmomaticTrimmed_unpaired.fastq \
$1_R2_trimmomaticTrimmed_paired.fastq \
$1_R2_trimmomaticTrimmed_unpaired.fastq \
TRAILING:20 MINLEN:30

# clean raw seq
rm $2_L001_R1_001.fastq
rm $2_L001_R2_001.fastq

# compress output
gzip -f $1_R1_trimmomaticTrimmed_paired.fastq
gzip -f $1_R1_trimmomaticTrimmed_unpaired.fastq
gzip -f $1_R2_trimmomaticTrimmed_paired.fastq
gzip -f $1_R2_trimmomaticTrimmed_unpaired.fastq

# transfer trimmed seq back to gluster
mv $1_R1_trimmomaticTrimmed_paired.fastq.gz /staging/qzhang333/
mv $1_R1_trimmomaticTrimmed_unpaired.fastq.gz /staging/qzhang333/
mv $1_R2_trimmomaticTrimmed_paired.fastq.gz /staging/qzhang333/
mv $1_R2_trimmomaticTrimmed_unpaired.fastq.gz /staging/qzhang333/

# clean individual seq data
# rm *.fastq
