#!/bin/bash

#unpack RSEM, Bowtie
tar -xzvf RSEM-1.3.1.tar.gz
unzip bowtie2-2.3.4-linux-x86_64.zip

# set PATH
export PATH=$(pwd)/bowtie2-2.3.4-linux-x86_64:$PATH
export PATH=$(pwd)/RSEM-1.3.1/bin:$PATH
export PATH=$(pwd)/RSEM-1.3.1/bin/samtools-1.3:$PATH

#transfer and unpack reads
cp /staging/qzhang333/$1_R1_trimmed_paired_mouseDNAremoved.fastq.gz .
cp /staging/qzhang333/$1_R2_trimmed_paired_mouseDNAremoved.fastq.gz .

gzip -d *fastq.gz

# transfer RSEM bowtie2 indexed DO metagenes reference
cp /staging/qzhang333/296DOmice_1.9M_GeneSet_rsem-prepared-bowtie2-20181215.tar.gz .

tar -xzvf 296DOmice_1.9M_GeneSet_rsem-prepared-bowtie2-20181215.tar.gz
rm *.tar.gz

# make directory for expression calculating
mkdir rsem-DOmetaGenes-$1

# calculate expression
rsem-calculate-expression -p 1 \
                          --bowtie2 \
                          --paired-end \
                          --no-bam-output \
                          $1_R1_trimmed_paired_mouseDNAremoved.fastq \
                          $1_R2_trimmed_paired_mouseDNAremoved.fastq \
                          ref/296DOmice_1.9M_GeneSet \
                          rsem-DOmetaGenes-$1/$1
                    
# move gene result file out
tar -czvf rsem-DOmetaGenes-$1.tar.gz rsem-DOmetaGenes-$1
mv rsem-DOmetaGenes-$1.tar.gz /staging/qzhang333

rm -R ref
rm *.fastq
