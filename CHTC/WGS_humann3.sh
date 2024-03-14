#!/bin/bash

# transfer humann3 built database to working directory
cp /staging/qzhang333/humann3_db.tar.gz .
tar -xzvf humann3_db.tar.gz
rm humann3_db.tar.gz

# transfer MetaPhlAn database to working directory
cp /staging/qzhang333/metaphlan_bowtie2db.tar.gz .
tar -xzvf metaphlan_bowtie2db.tar.gz
rm metaphlan_bowtie2db.tar.gz

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

### run Humann3
# make output dir
mkdir -p humann3_out_$1

# run humann3 using db: chocophlan, uniref, metaphlan_bowtie2db
humann --input $1.fastq \
       --output humann3_out_$1 \
       --nucleotide-database ./humann3_db/chocophlan \
       --protein-database ./humann3_db/uniref \
       --metaphlan-options "--bowtie2db ./metaphlan_bowtie2db --index mpa_v30_CHOCOPhlAn_201901" \
       --threads 2

# make a dir to store all table outputs
mkdir humann3_out_tbl_$1
mv humann3_out_$1/*.tsv humann3_out_tbl_$1

# normalize gene families to counts per million
humann_renorm_table --input humann3_out_tbl_$1/*genefamilies.tsv --output humann3_out_tbl_$1/$1_genefamilies_cpm.tsv --units cpm --update-snames

# convert to EC, KO, Metacyc
humann_regroup_table --input humann3_out_tbl_$1/$1_genefamilies_cpm.tsv --output humann3_out_tbl_$1/$1_genefamilies_cpm_MetaCyc.tsv --g uniref90_rxn

# copy log file to table dir
cp humann3_out_$1/*temp/*.log humann3_out_tbl_$1

# move docker std error file to output dir
mv docker_stderror humann3_out_$1/docker_stderror

# compress output dir
tar -czvf humann3_out_$1.tar.gz humann3_out_$1
tar -czvf humann3_out_tbl_$1.tar.gz humann3_out_tbl_$1

mv humann3_out_$1.tar.gz /staging/qzhang333/
mv humann3_out_tbl_$1.tar.gz /staging/qzhang333/

# clean space
rm $1.fastq
#rm humann3_out_$1.tar.gz