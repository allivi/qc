#!/bin/bash
#USAGE: ./remove_contaminants.sh sample.fq sample_ref.fa.ind contaminant.fa.ind threads_num

sample=$1
ref=$2
cont_ref=$3
threads=$4

dir=$(dirname $0)

bowtie2 -q -p $threads --very-fast -x $cont_ref -U $sample > tmp.sam
grep -oE "^[^@]\S+" tmp.sam > contaminant_reads.txt
python $dir/filter_reads.py include $sample contaminant_reads.txt > shared.fq
bowtie2 -q -p $threads --very-fast -x $ref --un uniq_cont.fq -U shared.fq > /dev/null
sed -re "/^[^@]/ d; s/^@(.*)/\1/" uniq_cont.fq > contaminant_reads.txt
python $dir/filter_reads.py exclude $sample contaminant_reads.txt
#rm tmp.sam contaminant_reads.txt shared.fq uniq_cont.fq
