#!/bin/bash

# transfer fastq reads files to working directory
cp /staging/qzhang333/$1_R1_trimmed_paired_mouseDNAremoved.fastq.gz .
cp /staging/qzhang333/$1_R2_trimmed_paired_mouseDNAremoved.fastq.gz .

# unpack compressed fastq files
gunzip $1_R1_trimmed_paired_mouseDNAremoved.fastq.gz
gunzip $1_R2_trimmed_paired_mouseDNAremoved.fastq.gz

# concat R1 and R2 reads
cat $1_R1_trimmed_paired_mouseDNAremoved.fastq $1_R2_trimmed_paired_mouseDNAremoved.fastq > $1.fastq

rm $1_R1_trimmed_paired_mouseDNAremoved.fastq
rm $1_R2_trimmed_paired_mouseDNAremoved.fastq

# transfer MetaPhlAn4 bowtie2db
cp /staging/qzhang333/metaphlan4_bowtie2db.tar.gz .

tar -xzf metaphlan4_bowtie2db.tar.gz
rm metaphlan4_bowtie2db.tar.gz

# make output dir
mkdir -p metaphlan4_out_$1

# modify this line to run desired Python script and any other work need to do
metaphlan $1.fastq --bowtie2db metaphlan4_bowtie2db \
                   --bowtie2out $1.bowtie2.bz2 \
                   --unclassified_estimation \
                   --input_type fastq \
                   --nproc 2 \
                   -o $1_metaphlan_output.txt


# transfer output filt to ourput dir
mv $1_metaphlan_output.txt metaphlan4_out_$1
mv $1.bowtie2.bz2 metaphlan4_out_$1

# compress output dir
tar -czvf metaphlan4_out_$1.tar.gz metaphlan4_out_$1

mv metaphlan4_out_$1.tar.gz /staging/qzhang333/

# clean space
rm $1.fastq